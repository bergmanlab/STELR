import os
import sys
import random
import json
import subprocess
import glob
import traceback
telr_dir = f"{__file__.split('/evaluation')[0]}/telr"
sys.path.insert(0,telr_dir)
from STELR_utility import getdict, setdict, memory_format, abs_path, mkdir, symlink, check_exist, string_to_bool
from interpret_config import config_from_file

def main():
    config = parse_args()
    config.update(get_file_paths())

    if not "resume" in config:
        #config.update(handle_file_paths(config))
        config.update(setup_run(config))
    else:
        run_id = config["resume"]
        tmp_dir = f"{config['out']}/stelr_eval_run_{run_id}"
        config_file = f"{tmp_dir}/config.json"
        with open(config_file,"r") as config_from_file:
            config = json.load(config_from_file)
        config["run_id"] = run_id
        config["resume"] = run_id
        config["tmp_dir"] = tmp_dir
        config["config_file"] = config_file
        config = update_config(config)

    setdict(config,("resources","estimated memory"),f'mem_mb={getdict(config,("resources","estimated memory"))}')
    snakefile = os.path.split(__file__)[0] + "/evaluation.smk"
    run_workflow(snakefile, config)


def get_file_paths():
    eval_src = abs_path(__file__)
    eval_src_dir = os.path.split(eval_src)[0]
    telr_dir = eval_src_dir.split("evaluation")[0] + "telr"
    env_dir = eval_src_dir.split("src")[0] + "envs"
    
    path_dict = {
        "telr":f"{telr_dir}/stelr.py",
        #"liftover":#TODO add a liftover script?
        "telr_dir":telr_dir,
        "eval_src_dir":eval_src_dir,
        "telr_conda":f"{env_dir}/stelr.yaml",
        "liftover_eval":f"{eval_src_dir}/liftover_evaluation.py",
        "af_eval":f"{eval_src_dir}/telr_af_eval.py",
        "seq_eval":f"{eval_src_dir}/telr_seq_eval.py"
    }

    return path_dict

def override_args():
    params = {("out",):abs_path()}
    args = {
        "-t":[("resources","threads"),int],
        "--threads":[("resources","threads"),int],
        "--mem":[("resources","estimated memory"),memory_format],
        "--memory":[("resources","estimated memory"),memory_format],
        "-o":[("out",),abs_path],
        "--resume":[("resume",),int],
        "--use-slurm":[("slurm options","use slurm"), string_to_bool]
    }
    unparsed = []
    for index in range(1,len(sys.argv)):
        prev = sys.argv[index-1]
        this = sys.argv[index]
        if prev in args:
            params[args[prev][0]] = args[prev][1](this)
        elif this in args:
            pass
        else:
            unparsed.append(this)
    
    return params, unparsed

def update_config(config, params=False):
    if not params:
        params = override_args()[0]
    for param in params:
        setdict(config,param,params[param])
    
    return config

def parse_args():
    params, unparsed = override_args()

    if ("resume",) in params:
        with open(f"stelr_eval_run_{params[('resume',)]}/config.json","r") as input:
            config = json.load(input)
    else: config = config_from_file(unparsed[0])
    #print(params)
    config = update_config(config, params)
    
    return config

def setup_run(config):
    run_id = random.randint(1000000,9999999) #generate a random run ID
    while len(glob.glob(f"stelr_eval_run_{run_id}/")) > 0:
        run_id = random.randint(1000000,9999999) #generate a random run ID
    tmp_dir = f"{config['out']}/stelr_eval_run_{run_id}"
    mkdir(tmp_dir)

    input_names, default_inputs = find_default_inputs()
    for reference in input_names:
        for input_file in input_names[reference]:
            maplist = ["input options",reference,input_file]
            link = f"{tmp_dir}/{getdict(input_names, maplist[1:])}"
            try:
                if maplist[-1] in getdict(default_inputs, maplist[1:-1]):
                    source = getdict(default_inputs, maplist[1:])
                if input_file in getdict(config, maplist[:-1]):
                    print(getdict(config, maplist[:-1]))
                    source = getdict(config, maplist)
                if check_exist(source):
                    symlink(source, link)
                    print(source)
                else: 
                    config["input options"][reference][input_file] = source
                    continue
            except: traceback.print_exc()
            config["input options"][reference][input_file] = link
    
    config["output"] = []
    for sim in config["simulation parameters"]:
        if not sim == "seed":
            sim_dict = getdict(config,["simulation parameters",sim])
            if "file" in sim_dict:
                input_file = sim_dict["file"]
                if sim == "model":
                    sim_dict["file"] = link = abs_path(f"{tmp_dir}/{sim_dict['model name']}.model")
                    if check_exist(input_file):
                        symlink(input_file,link)
                else:
                    sim_dir = abs_path(f"{tmp_dir}/{sim}")
                    mkdir(sim_dir, False)
                    sim_dict["file"] = link = abs_path(f"{sim_dir}/simulated_reads.fq")
                    if check_exist(input_file):
                        symlink(input_file,link)
                    config["output"] += [
                        f"{sim_dir}/liftover_eval/annotation.filter.bed",
                        f"{sim_dir}/af_eval/telr_eval_af.json",
                        f"{sim_dir}/seq_eval/seq_eval.json"
                    ]
    
    config["run_id"] = run_id
    config["tmp_dir"] = tmp_dir

    if "use slurm" in config["slurm options"]:
        slurm_optionals={
            "--mail-user=":["slurm options","(optional) mail user"],
            "--mail-type=":["slurm options","(optional) mail type"]
        }
        slurm_optionals = " \\\n                ".join([f"{o}{getdict(config,slurm_optionals[o])}" for o in slurm_optionals if slurm_optionals[o][-1] in getdict(config,slurm_optionals[o][:-1])])
        slurm_command = f"""
        --jobs 200
        --use-conda
        --latency-wait 999
        --keep-going
        --cluster-status 'bash {config["eval_src_dir"]}/slurm_status.sh'
        --cluster-cancel 'scancel'
        --cluster""".strip().split("'")
        for index in range(len(slurm_command)):
            if(index/2 == int(index/2)):
                slurm_command[index] = slurm_command[index].split()
            else:
                slurm_command[index] = [slurm_command[index]]
        slurm_command = [item for sublist in slurm_command for item in sublist] + [
            f"""sbatch \\
                --partition={config["slurm options"]["partition"]} \\
                --nodes=1 \\
                --tasks-per-node={{threads}} \\
                --mem={{resources.mem_mb}} \\
                --time={{resources.runtime}} \\
                --parsable \\
                {slurm_optionals}"""
                ]
        config["slurm options"]["use slurm"] = slurm_command

    config_path = f"{tmp_dir}/config.json" # path to config file
    config["config_file"] = config_path
    with open(config_path, "w") as conf:
        json.dump(config, conf) #write config file as json    

    return config

def find_default_inputs():
    test_dir = "/src/".join(abs_path(__file__).split("/src/")[:-1]) + "/test/evaluation"
    input_names = {
        "community reference":{
            "file":"community_reference.fasta",
            "te annotation file path":f"community_annotation.bed"
        },
        "mapping reference":{
            "file":"mapping_reference.fasta",
            "te annotation file path":f"mapping_annotation.bed",
            "regular recombination region file path":"regular_recomb.bed"
        },
        "te library":{"file":"library.fa"}
    }
    default_inputs = {
        "community reference":{
            "accession":"GCF_000001215.4",
            "te annotation file path":f"{test_dir}/community_annotation.bed"
        },
        "mapping reference":{
            "accession":"GCA_003401745.1",
            "te annotation file path":f"{test_dir}/mapping_annotation.bed",
            "regular recombination region file path":f"{test_dir}/regular_recomb.bed"
        },
        "te library":{"file":f"{test_dir}/library.fa"}
    }
    return [input_names, default_inputs]




def run_workflow(snakefile, config):
    print(f"\nSTELR_evaluation run ID {config['run_id']}\n")
    command = [
        "snakemake",
        "-s", snakefile,
        "--configfile", config["config_file"],
        "--use-conda"
    ]
    optional_args = {("resources","threads"):"--cores",("resources","estimated memory"):"--resources"}
    try:
        for arg in optional_args:
            if getdict(config,arg):
                command += [optional_args[arg], str(getdict(config,arg))]
    except: pass
    if "use slurm" in config["slurm options"]:
        if config["slurm options"]["use slurm"]:
            command += config["slurm options"]["use slurm"]
    try:
        if not "resume" in config:
            subprocess.run(command, cwd=config["tmp_dir"])
        else:
            subprocess.run(command + ["--unlock"], cwd=config["tmp_dir"])
            subprocess.run(command + ["--rerun-incomplete", "--rerun-triggers","mtime","--retries", "0"], cwd=config["tmp_dir"])
            #subprocess.run(command + ["--ignore-incomplete", "--rerun-triggers","mtime"], cwd=config["tmp_dir"])
        #for output_file in config["output"]:
        #    os.rename(f"{config['tmp_dir']}/{output_file}",f"{config['out']}/{output_file}")
        if not getdict(config,("TELR parameters","Additional options","Keep intermediate files")):
            rmtree(config['tmp_dir'])
    except Exception as e:
        print(e)
        print("STELR_evaluation failed!")
        sys.exit(1)
    print(f"STELR_evaluation run {config['run_id']} finished.")

if __name__ == '__main__':
    main()
