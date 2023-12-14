#this is just an idea, and an incomplete one, not in use
import sys
import json
import os
import traceback
telr_dir = f"{__file__.split('/evaluation')[0]}/stelr"
sys.path.insert(0,telr_dir)
from STELR_utility import memory_format, check_exist, getdict, setdict

def jprint(text):
    print(json.dumps(text,indent=4))

def decompose_dict(d):
    d_list = []
    for key in d:
        if type(d[key]) is dict and d[key]: d_list.append([key,decompose_dict(d[key])])
        else: d_list.append(key)
    return d_list

def flatten(l):
    if type(l) is list:
        flat_list = []
        for x in range(len(l)):
            flat_x = flatten(l[x])
            if flat_x[0]:
                flat_x = [[[x] + item[0], item[1]] for item in flat_x]
                flat_list += flat_x
            else:
                flat_x[0].insert(0,x)
                flat_list += [flat_x]
        return flat_list
    else: return [[],l]

def config_from_file(file_path):
    groups = {}
    depth = 0
    m = []
    last_m = {0:[]}
    with open(file_path,"r") as input:
        for line in input:
            line_depth = line.index(line.strip())
            strip = line.strip()
            if(line_depth > depth):
                last_m[line_depth] = m.copy()
                m += [strip]
                setdict(groups,m,{})
                depth = line_depth
            elif strip:
                m = last_m[line_depth] + [strip]
                setdict(groups,m,{})
            depth = line_depth

    num_levels = len(last_m)

    #jprint(groups)

    return decompose_dict(groups)
    

if __name__ == "__main__":
    args = config_from_file(sys.argv[1])
    #jprint(args)
    list_coords = {item[1]:item[0] for item in flatten(args)}
    jprint(flat_args)
    temp = args.copy()
    for number in [1,1,2,0][:-1]:
        temp = temp[number]
    print(temp)

    args_dict = {
        "input options":{
            "community reference"
        }
    }