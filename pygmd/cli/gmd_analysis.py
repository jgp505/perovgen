# coding : utf-8
# Cppyright (c) Green Materials Designs Team.

import os
import sys

import numpy as np
import pandas as pd

import yaml
from collections import defaultdict
from shutil import copy

from pymatgen.io.vasp.outputs import Vasprun
from perovgen.pygmd.shell import ShellPath
from perovgen.pygmd.input_structure import load_structure,GMDStructure
from perovgen.pygmd.analysis.energy import GMDPhaseDiagram

def graphyaml(string) :
    path = os.path.abspath(__file__).split(os.sep)[:-2]
    path.append("analysis")
    path.append("graph.yaml")
    path = os.sep.join(path)

    with open(path) as fi :
        graph = yaml.load(fi,Loader=yaml.FullLoader)
    return graph[string]

def pdospath() :
    path = os.path.abspath(__file__).split(os.sep)[:-3]
    path.append("pdos")
    path.append("pdos")
    path = os.sep.join(path)
    return path

def analysis(args) :
    if args.graph :
        dic = graphyaml(args.graph)
        with open('graph.yaml','w') as fi :
            yaml.dump(dic,fi,default_flow_style=False)
        print("%s generate the graph.yaml file"%(os.getcwd()))

    elif args.input :
        shell, numberlist = ShellPath().check()
        while True :
            number = int(input("Please enter the number >> "))
            if number in numberlist :
                break
        mn = str(input("Please enter the mode ex) RCDB >> "))

        with open("input.gmd",'w') as fi :
            fi.write("### Perovgen input.gmd\n")
            fi.write("KPOINTS:\nA\nCONSTK=30\n\n")
            fi.write("INCAR:\nMPJ\n\n")
            fi.write("SHELL:\n%s\n\n"%(shell[number-1][-1]))
            fi.write("CALMODE:\n%s\n"%(mn))
        fi.close()

    elif args.pdos :
        path = pdospath()
        s = load_structure(os.getcwd())[0][0]
        if not s :
            print("Structure file doesn't exist!\n") 
            sys.exit(1)
        # Total DOS
        print("Total DOS")
        os.system("%s dos width=0.03"%(path))

        # partial DOS
        dic = defaultdict(list)
        print("pDOS")
        for e,i in enumerate(s.species):
            dic[i.symbol].append(e+1)
        for k,v in dic.items():
            for s in ['tot','s','p'] :
                with open("LIST",'w') as fi :
                    fi.write("{}_{}\n".format(k,s))
                    for i in v :
                        if s == "tot" :
                            fi.write("%i tot tot\n"%(i))
                        elif s == "s" :
                            fi.write("%i s tot\n"%(i))
                        elif s == "p" :
                            fi.write("%i px tot\n"%(i))
                            fi.write("%i py tot\n"%(i))
                            fi.write("%i pz tot\n"%(i))
                os.system("{} pdos width=0.03".format(path))

    elif args.convert:
        files = load_structure(args.convert)[0]
        for e,f in enumerate(files) : 
            GMDStructure(f).convert_to_cif()