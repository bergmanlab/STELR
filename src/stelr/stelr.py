#!/usr/bin/env python3
__version__ = "2.0"

import sys
import os
import random
import argparse
import subprocess
import logging
import json
from shutil import rmtree
from time import perf_counter

def main():
    config = get_args()
    if config["resume"]:
        resume_params = config
        with open(f"{config['out']}/stelr_run_{config['resume']}/config.json","r") as config_from_file:
            config = json.load(config_from_file)
        resume_params.pop("out")
        config.update(resume_params)
    else:
        config["verbose"] = False #TODO: add as option later
        if_verbose = verbose(config["verbose"])
        config["sample_name"] = os.path.splitext(os.path.basename(config["reads"]))[0]
    
        config.update(handle_file_paths(config))
        config.update(setup_run(if_verbose, config))
    run_workflow(config)

def handle_file_paths(config):
    stelr_src = os.path.abspath(__file__)
    stelr_dir = os.path.split(stelr_src)[0]
    
    source_files = ["alignment","assembly","sv","te","liftover","utility","output"]
    for file in source_files:
        config[f"STELR_{file}"] = f"{stelr_dir}/STELR_{file}.py"
    config["smk"] = f"{stelr_dir}/STELR.smk"
    config["STELR_contig"] = f"{stelr_dir}/STELR_contig.smk"
    config["fix_ngmlr"] = f"{stelr_dir}/fix_ngmlr.py"
    env_dir = stelr_dir.split("src")[0] + "envs"
    config["envs"] = {}


    return config

def setup_run(if_verbose, config):
    run_id = random.randint(1000000,9999999) #generate a random run ID
    tmp_dir = f"{config['out']}/stelr_run_{run_id}"
    mkdir(if_verbose, tmp_dir)
    mkdir(if_verbose, os.path.join(tmp_dir, "input"))

    config.update(process_input_files(
        {"reads":config["reads"],
        "reference":config["reference"],
        "library":config["library"]},
        f"{tmp_dir}/input",
        config["sample_name"])
    )

    config["output"] = [
        f"reads.stelr.contig.fasta",
        f"reads.stelr.te.fasta",
        f"reads.stelr.bed",
        f"reads.stelr.json",
        f"reads.stelr.expanded.json",
        f"reads.stelr.vcf"
    ]

    config["run_id"] = run_id
    config["tmp_dir"] = tmp_dir
    config["config_path"] = f"{tmp_dir}/config.json" # path to config file
    with open(config["config_path"], "w") as conf:
        json.dump(config, conf, indent=4) #write config file as json    
    
    return config

def run_workflow(config):
    print(f"STELR version {__version__}")
    print(f"Run ID {config['run_id']}")
    command = [
        "snakemake", #"--use-conda",#"--conda-prefix",config["conda_yaml"],
        "-s",config["smk"],
        "--configfile", config["config_path"],
        "--cores", str(config["thread"]),# "--quiet"
        #"--unlock"
    ]
    try:
        if config["resume"] is None: 
            start = perf_counter()
            subprocess.run(command, cwd=config["tmp_dir"])
        else:
            subprocess.run(command + ["--unlock"], cwd=config["tmp_dir"])
            subprocess.run(command + ["--rerun-incomplete"] + ["--rerun-triggers","mtime"], cwd=config["tmp_dir"])
        for output_file in config["output"]:
            os.rename(f"{config['tmp_dir']}/{output_file}",f"{config['out']}/{output_file.replace('reads',config['sample_name'])}")
        if not config["keep_files"]:
            rmtree(config['tmp_dir'])
    except Exception as e:
        print(e)
        print("STELR failed!")
        sys.exit(1)
    if config["resume"] is None: print(f"STELR finished in {perf_counter()-start} seconds!")
    print(f"STELR run {config['run_id']} finished.")

def process_input_files(input_file_paths, input_dir, sample_name):
    def file_extension_of(file):
        return os.path.splitext(input_file_paths[file])[1]
    accepted_formats_for = { #valid file extensions for each type of input
        "reads":[".fasta",".fastq",".fa",".fq",".fna",".fa",".bam"],
        "library":[".fasta",".fastq",".fa",".fq",".fna",".fa"],
        "reference":[".fasta",".fastq",".fa",".fq",".fna",".fa"]
        }
    new_paths = {}
    for file in ["reads","reference","library"]:
        try: open(input_file_paths[file], "r")
        except Exception as e:
            print(e)
            logging.exception(f"Cannot open input file: {input_file_paths[file]}")
            sys.exit(1)
        if not file_extension_of(file) in accepted_formats_for[file]:
                print(f"Input {file} format not recognized, exiting...")
                logging.error("Input format not recognized")
                sys.exit(1)
        else: 
            new_paths[file] = f"{input_dir}/{file}{file_extension_of(file)}"
            symlink(input_file_paths[file], new_paths[file])
    if ".bam" in new_paths["reads"]: new_paths["fasta_reads"] = f"{input_dir}/reads.fasta"
    else: new_paths["fasta_reads"] = new_paths["reads"]
    return new_paths
        

def symlink(input_file, link): #create a symbolic link at the output location referencing the input path.
    input_file = os.path.abspath(input_file)
    #abspath needed because path for os.symlink is relative to output, not current directory.
    if os.path.islink(link):
        os.remove(link)
    try:
        os.symlink(input_file, link)
    except Exception:
        logging.exception(f"Create symbolic link for {input_file} failed")
        sys.exit(1)

def mkdir(if_verbose, dir):
    if os.path.isdir(dir):
        if_verbose.print(f"Directory {dir} exists") # verbose object print
        return
    try:
        os.mkdir(dir)
    except OSError:
        print(f"Creation of the directory {dir} failed")
    else:
        if_verbose.print(f"Successfully created the directory {dir}")

def install(args):
    threads = args["thread"]
    from STELR_utility import prefix, mkdir, abs_path
    from shutil import copyfile
    import glob
    stelr_dir = prefix(abs_path(__file__),"/src")
    envs_dir = f"{stelr_dir}/envs"
    snake_dir = f"{envs_dir}/snakemake"
    mkdir(snake_dir)
    stelr_env_github = glob.glob(f"{envs_dir}/stable_environment_v*")[0]
    stelr_env_installed = snake_dir + stelr_env_github.split(envs_dir)[1]
    conda_version=stelr_env_github.split("stable_environment_v")[-1].split(".yaml")[0]
    copyfile(stelr_env_github,stelr_env_installed)
    smk_config = f"{snake_dir}/installation_config"
    args.update({
        "output":".installation_complete",
        "conda":stelr_env_installed
        })
    for unnecessary_arg in ["fasta_reads","input","reference","library"]:
        args.update({unnecessary_arg:"not needed for install"})
    with open(smk_config,"w") as output:
        json.dump(args,output)
    command=[
        "snakemake","-s",f"{stelr_dir}/src/stelr/STELR.smk",#"--quiet",
        "--configfile",smk_config,
        "--cores",str(threads),
        "--use-conda",#stelr_env_installed,
        "--conda-create-envs-only","--use-singularity"
        ]
    subprocess.run(command, cwd=snake_dir)
    os.remove(smk_config)
    installed_version = {
        "version":conda_version,
        "installed_env":stelr_env_installed.split("envs/")[-1]
    }
    with open(f"{snake_dir}/installed_version","w") as output:
        json.dump(installed_version,output)
    print(f"STELR installation complete:\n\tSTELR version: {__version__}\n\tRunning environment version: {conda_version}")

def get_args():

    if "--resume" in sys.argv:
        config = {"out":"."}
        resume_params = {
            "--resume":"resume",
            "--out":"out",
            "-o":"out",
            "-t":"thread"
        }
        for param in resume_params:
            if param in sys.argv:
                config.update({resume_params[param]:sys.argv[sys.argv.index(param)+1]})
        return config
    parser = argparse.ArgumentParser(
        description="Program for detecting non-reference TEs in long read data"
    )
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")

    # required
    if not "--install" in sys.argv:
        required.add_argument(
            "-i",
            "--reads",
            type=str,
            help="reads in fasta/fastq format or read alignments in bam format",
            required=True,
        )
        required.add_argument(
            "-r",
            "--reference",
            type=str,
            help="reference genome in fasta format",
            required=True,
        )
        required.add_argument(
            "-l",
            "--library",
            type=str,
            help="TE consensus sequences in fasta format",
            required=True,
        )
        

    # optional
    optional.add_argument(
        "--aligner",
        type=str,
        help="choose method for read alignment, please provide 'ngmlr' or 'minimap2' (default = 'ngmlr')",
        required=False,
    )
    optional.add_argument(
        "--assembler",
        type=str,
        help="Choose the method to be used for local contig assembly step, please provide 'wtdbg2' or 'flye' (default = 'wtdbg2')",
        required=False,
    )
    optional.add_argument(
        "--polisher",
        type=str,
        help="Choose the method to be used for local contig polishing step, please provide 'wtdbg2' or 'flye' (default = 'wtdbg2')",
        required=False,
    )
    optional.add_argument(
        "-x",
        "--presets",
        type=str,
        help="parameter presets for different sequencing technologies, please provide 'pacbio' or 'ont' (default = 'pacbio')",
        required=False,
    )
    optional.add_argument(
        "-p",
        "--polish_iterations",
        type=int,
        help="iterations of contig polishing (default = 1)",
        required=False,
    )
    optional.add_argument(
        "-o",
        "--out",
        type=str,
        help="directory to output data (default = '.')",
        required=False,
    )
    optional.add_argument(
        "-t",
        "--thread",
        type=int,
        help="max cpu threads to use (default = '1')",
        required=False,
    )
    optional.add_argument(
        "-g",
        "--gap",
        type=int,
        help="max gap size for flanking sequence alignment (default = '20')",
        required=False,
    )
    optional.add_argument(
        "-v",
        "--overlap",
        type=int,
        help="max overlap size for flanking sequence alignment (default = '20')",
        required=False,
    )
    optional.add_argument(
        "--flank_len",
        type=int,
        help="flanking sequence length (default = '500')",
        required=False,
    )
    optional.add_argument(
        "--af_flank_interval",
        type=int,
        help="5' and 3'flanking sequence interval size used for allele frequency estimation (default = '100')",
        required=False,
    )
    optional.add_argument(
        "--af_flank_offset",
        type=int,
        help="5' and 3' flanking sequence offset size used for allele frequency estimation (default = '200')",
        required=False,
    )
    optional.add_argument(
        "--af_te_interval",
        type=int,
        help="5' and 3' te sequence interval size used for allele frequency estimation (default: '50')",
        required=False,
    )
    optional.add_argument(
        "--af_te_offset",
        type=int,
        help="5' and 3' te sequence offset size used for allele frequency estimation (default: '50')",
        required=False,
    )
    optional.add_argument(
        "--different_contig_name",
        action="store_true",
        help="If provided then STELR does not require the contig name to match before and after annotation liftover (default: require contig name to be the same before and after liftover)",
        required=False,
    )
    optional.add_argument(
        "--minimap2_family",
        action="store_true",
        help="If provided then minimap2 will be used to annotate TE families in the assembled contigs (default: use repeatmasker for contig TE annotation)",
        required=False,
    )
    optional.add_argument(
        "-k",
        "--keep_files",
        action="store_true",
        help="If provided then all intermediate files will be kept (default: remove intermediate files)",
        required=False,
    )
    optional.add_argument(
        "--resume",
        type=int,
        help="Resume a previous run by ID (default: begin new run)",
        required=False,
    )
    optional.add_argument(
        "--install",
        action="store_true",
        help="Use this option to install the stable conda environments required to run STELR.",
        required=False,
    )
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # check if optional arguments are valid
    if args.aligner is None:
        args.aligner = "ngmlr"
    elif args.aligner not in ["ngmlr", "minimap2"]:
        print("Please provide a valid alignment method (ngmlr/minimap2), exiting...")
        sys.exit(1)

    if args.assembler is None:
        args.assembler = "wtdbg2"
    elif args.assembler not in ["wtdbg2", "flye"]:
        print("Please provide a valid assembly method (wtdbg2/flye), exiting...")
        sys.exit(1)

    if args.polisher is None:
        args.polisher = "wtdbg2"
    elif args.polisher not in ["wtdbg2", "flye"]:
        print("Please provide a valid polish method (wtdbg2/flye), exiting...")
        sys.exit(1)

    if args.presets is None:
        args.presets = "pacbio"
    elif args.presets not in ["pacbio", "ont"]:
        print("Please provide a valid preset option (pacbio/ont), exiting...")
        sys.exit(1)

    if args.polish_iterations is None:
        args.polish_iterations = 1
    elif args.polish_iterations < 1:
        print("Please provide a valid number of iterations for polishing, exiting...")

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    mkdir(verbose(False), args.out)

    if args.thread is None:
        args.thread = 1

    if args.flank_len is None:
        args.flank_len = 500

    if args.af_flank_interval is None:
        args.af_flank_interval = 100
    else:
        if args.af_flank_interval <= 0:
            print(
                "Please provide a valid flanking sequence interval size (positive integer) for allele frequency estimation, exiting..."
            )
            sys.exit(1)

    if args.af_flank_offset is None:
        args.af_flank_offset = 200
    else:
        if args.af_flank_offset < 0:
            print(
                "Please provide a valid flanking sequence offset size (positive integer) for allele frequency estimation, exiting..."
            )

    if args.af_te_interval is None:
        args.af_te_interval = 50
    else:
        if args.af_te_interval <= 0:
            print(
                "Please provide a valid TE interval size (positive integer) for allele frequency estimation, exiting..."
            )

    if args.af_te_offset is None:
        args.af_te_offset = 50
    else:
        if args.af_te_offset < 0:
            print(
                "Please provide a valid TE offset size (positive integer) for allele frequency estimation, exiting..."
            )

    if args.gap is None:
        args.gap = 20

    if args.overlap is None:
        args.overlap = 20
    

    if "--install" in sys.argv:
        install(vars(args))
        quit()
    else: return vars(args)

class verbose:
    def __init__(self, verbose=False):
        self.verbose = verbose
    
    def print(self, string):
        if(self.verbose):
            print(string)

if __name__ == "__main__":
    main()