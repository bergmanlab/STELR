import subprocess
import json
import os
import sys

class paf_file:
    def __init__(self, file_path):
        def format_paf_line(line):
            if line.strip():
                line = line.split("\t")
                for n in [1,7,8,9,10,11]:
                    line[n] = int(line[n])
                return line
            else:
                return None

        if os.path.isfile(file_path):
            self.file_path = file_path            
            with open(file_path,"r") as paf_file:
                paf_data = [format_paf_line(line) for line in paf_file]
        
        else: paf_data = [format_paf_line(line) for line in file_path.split("\n")]
        
        self.data = [line for line in paf_data if line]
        
        max_len = max([line[10] for line in self.data])
        self.longest = [line for line in self.data if line[10] == max_len][0]
    
    def count(self):
        return len(self.data)
    
    def get(self, key):
        if not "characteristics" in self.__dict__:
            header = ["","query_len","","","","chrom","","start","end","matches","align_len","map_qual"]
            self.characteristics = {key:self.longest[header.index(key)] for key in header if key}
            self.characteristics["map_prop"] = self.characteristics["matches"]/self.characteristics["query_len"]
            self.characteristics["blast_id"] = self.characteristics["matches"]/self.characteristics["align_len"]
        return self.characteristics[key]
    
    def paftools_summary(self):
        if not "paftools_fields" in self.__dict__:
            self.paftools_fields = {
                "insertions":{
                    "1bp":0,
                    "2bp":0,
                    "[3,50)":0,
                    "[50,1000)":0,
                    "total":0
                },
                "deletions":{
                    "1bp":0,
                    "2bp":0,
                    "[3,50)":0,
                    "[50,1000)":0,
                    "total":0
                },
                "reference bases covered":0,
                "substitutions":0
            }
            
            paftools_summary = subprocess.run(["paftools.js","call","-l","100","-L","100",self.file_path], capture_output=True, text=True).stderr

            for line in paftools_summary.split("\n"):
                try:
                    value = int(line.split()[0])
                    field = [field for field in self.paftools_fields if field in line][0]
                    if type(self.paftools_fields[field]) is int: self.paftools_fields[field] = value
                    else:
                        subfield = [subfield for subfield in self.paftools_fields[field] if subfield in line][0]
                        self.paftools_fields[field][subfield] = value
                        self.paftools_fields[field]["total"] += value
                except: pass
        
        return self.paftools_fields

class bed_file:
    def __init__(self, bed_file):
        self.get_index = 0
        if type(bed_file) is list:
            if not type(bed_file[0]) is list: bed_file = [bed_file]
            self.data = bed_file
        else:
            def format_bed_line(line):
                if line.strip():
                    line = line.split("\t")
                    for n in [1,2]:
                        line[n] = int(line[n])
                    return line
                else: return [0,0,0,0,0,0]
            with open(bed_file,"r") as input_file:
                self.data = [format_bed_line(line) for line in input_file if line.strip()]
    
    def get(self, key, i=0):
        if not ("characteristics" in self.__dict__ and self.get_index == i):
            self.get_index = i
            header = ["chrom", "start", "end", "name", "score", "strand"]
            self.characteristics = {key:self.data[i][header.index(key)] for key in header if key}
        return self.characteristics[key]
    
    def write_out(self, file_path):
        if type(self.data[0]) is list:
            lines = ["\t".join([str(i) for i in line]) for line in self.data]
        else: lines = ["\t".join([str(i) for i in self.data])]
        with open(file_path,"w") as output_file:
            output_file.write("\n".join(lines))
    
    def faidx(self, reference_fasta, fasta_out = "", make_file = True, save_sequence = False, i=0):
        coord = f"{self.get('chrom',i)}:{self.get('start',i)+1}-{self.get('end',i)}"
        sequence = subprocess.run(["samtools","faidx",reference_fasta,coord], capture_output=True, text=True).stdout
        if make_file:
            with open(fasta_out,"w") as output_file:
                output_file.write(sequence)
            self.fasta = fasta_out
        if save_sequence:
            self.sequence = "".join(sequence.strip().split(">")[1].split("\n")[1:])

    def get_sequence(self, reference_fasta = None, i=0):
        if reference_fasta:
            self.faidx(reference_fasta, make_file = False, save_sequence = True, i = i)
        return self.sequence

class process:
    def __init__(self, command=[]):
        self.command = [str(i) for i in command]
    
    def add(self, appendment):
        if appendment:
            if type(appendment) is list:
                self.command = self.command + [str(i) for i in appendment]
            else: self.command.append(str(appendment))
    
    def output(self):
        return subprocess.run(self.command,capture_output=True,text=True).stdout
    def run(self):
        subprocess.run(self.command)    
    def write_to(self, output_file):
        with open(output_file,"w") as output_file:
            subprocess.run(self.command, stdout=output_file)

def minimap2(target, query, *args, **kwargs):

    class minimap_preset:
        def __init__(self):
            self.arg = "x"
            self.preset = "asm10"

        def add(self, a):
            self.arg = a + self.arg
        
        def set(self, a):
            self.preset = a

    def is_member(arg, l):
        a = arg.lower()
        if a in l:
            if type(l) is dict:
                return l[arg]
            return a
        if not a:
            return list[0]
        print(f"Error: provided arg {arg} is not in list of valid args {l}")
        sys.exit(1)    

    preset = minimap_preset()

    default_args = {
        "secondary":"--secondary=no",
        "verbose":["-v","0"],
    }

    commands = {
        "threads":str,
        "format":lambda a: is_member(a,{"paf":"","sam":"a"}),
        "cigar":lambda _: preset.add("c"),
        "preset":lambda x: preset.set(x),
        "cs":lambda _: "--cs"
    }

    p = process(["minimap2"])

    for arg in args:
        arg = arg.lower()
        if arg in default_args:
            default_args.pop(arg)
        if arg in commands:
            p.add(commands[arg](""))
    
    for kwarg in kwargs:
        key = kwarg.lower()
        if key in default_args:
            default_args.pop(key)
        if key in commands:
            p.add(commands[key](kwargs[kwarg]))
    
    for arg in default_args:
        p.add(default_args[arg])
    
    p.add([f"-{preset.arg}",preset.preset])    
    p.add([target,query])
    
    output_args = [arg for arg in ["-o","output"] if arg in kwargs]
    if len(output_args) == 0:
        return p.output()
    else:
        p.add(["-o",kwargs[output_args[0]]])
        p.run()
    

