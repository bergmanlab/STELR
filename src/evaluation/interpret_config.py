import sys
import json
import os
import traceback
telr_dir = f"{__file__.split('/evaluation')[0]}/telr"
sys.path.insert(0,telr_dir)
from STELR_utility import memory_format, check_exist, getdict, setdict

def config_from_file(file_path):
    f = open(file_path, "r")
    text = f.read()
    f.close()

    category = {}
    title = ""
    subtitle = ""
    for line in text.split("\n"):
        if line.strip():
            if line[0].strip():
                if ":" in line:
                    title = line.split(":")[0].strip().lower()
                    category[title] = {}
            else:
                if title == "input options":
                    if "community reference" in line.lower():
                        subtitle = "community reference"
                        category[title][subtitle] = {}
                    elif "mapping reference" in line.lower():
                        subtitle = "mapping reference"
                        category[title][subtitle] = {}
                    elif "te library" in line.lower():
                        subtitle = "te library"
                        category[title][subtitle] = {}
                    else:
                        line = checkline(line)
                        if line.check: 
                            category[title][subtitle][line.item] = line.content
                        elif line.item == "file path or genbank accession number":
                            if check_exist(line.content):
                                category[title][subtitle]["file"] = line.content
                            else: 
                                category[title][subtitle]["accession"] = line.content
                elif title == "simulation parameters":
                    if "seed" in line.lower():
                        category[title]["seed"] = int(line.split(":")[1].strip())
                    elif "Model:" in line:
                        subtitle = "model"
                        category[title][subtitle] = {}
                        category[title][subtitle]["model name"] = model_name = line.split(":")[1].strip()
                        category[title][subtitle]["file"] = f"{model_name}.model"
                    elif subtitle == "model":
                        line = checkline(line)
                        if line.check: 
                            category[title][subtitle][line.item] = line.content
                        subtitle = None
                    elif "check all" in line.lower():
                        line = checkline(line)
                        if line.check:
                            for simulation in [
                                "50x_diploid_homozygous",
                                "50x_diploid_heterozygous",
                                "50x_tetraploid_simplex",
                                "50x_tetraploid_duplex",
                                "50x_tetraploid_triplex",
                                "50x_tetraploid_quadruplex",
                                "100x_diploid_homozygous",
                                "100x_diploid_heterozygous",
                                "100x_tetraploid_simplex",
                                "100x_tetraploid_duplex",
                                "100x_tetraploid_triplex",
                                "100x_tetraploid_quadruplex",
                                "150x_diploid_homozygous",
                                "150x_diploid_heterozygous",
                                "150x_tetraploid_simplex",
                                "150x_tetraploid_duplex",
                                "150x_tetraploid_triplex",
                                "150x_tetraploid_quadruplex",
                                "200x_diploid_homozygous",
                                "200x_diploid_heterozygous",
                                "200x_tetraploid_simplex",
                                "200x_tetraploid_duplex",
                                "200x_tetraploid_triplex",
                                "200x_tetraploid_quadruplex"]:
                                category[title][simulation] = {}
                    else:
                        line = checkline(line)
                        try:
                            if line.item == "file":
                                try: category[title][subtitle][line.item] = line.content
                                except: pass
                            else:
                                subtitle = line.item.replace(" ","_")
                                if line.check:
                                    category[title][subtitle] = {}
                        except: traceback.print_exc()
                else:
                    try:
                        line = checkline(line)
                        if line.content:
                            category[title][line.item] = line.content
                        elif line.check == None:
                            subtitle = line.item
                        elif line.check:
                            if line.item in ["telr version 1.x", "use slurm"]:
                                category[title][line.item] = True
                            else: category[title][subtitle] = line.item
                    except: pass
    format_dict = {
        ("resources","estimated memory"):memory_format,
        ("resources","threads"):int
    }
    for param in format_dict:
        setdict(category,param,format_dict[param](getdict(category,param)))
    return category

class checkline:
    def __init__(self, line):
        self.check = None
        self.item = None
        self.content = None
        try: 
            self.check = line.split("[")[1].split("]")[0].strip()
            self.item = line.split("]")[1].split("#")[0].strip()
        except:
            self.item = line.split("#")[0].strip()
        if(self.item): self.item = self.item.lower()
        if ":" in self.item:
            self.item = self.item.split(":")[0].strip()
            self.content = ":".join(line.split(":")[1:]).split("#")[0].strip()
        if "use file" in self.item:
            self.item = "file"
            if not os.path.exists(self.content):
                self.content = None
        #self.print()
    
    def check_true(self):
        if self.check:
            return True
        return False
    
    def print(self):
        print(f"{self.check}\t{self.item}\t{self.content}")

if __name__ == "__main__":
    print(config_from_file(sys.argv[1]))