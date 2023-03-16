# coding : utf-8
# Cppyright (c) Green Materials Designs Team.

import os
import sys

import numpy as np
import pandas as pd

import yaml
from collections import defaultdict
from shutil import copy

from pymatgen.io.vasp.outputs import Vasprun, BSVasprun
from perovgen.pygmd.shell import ShellPath
from perovgen.pygmd.input_structure import load_structure,GMDStructure
import perovgen.pygmd.analysis.electronic as electronic


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

def write(args) :
    if args.graph :
        path = os.path.abspath(__file__).split(os.sep)[:-2]
        path.append("analysis")
        path.append("graph.yaml")
        path = os.sep.join(path)
        copy(path,os.getcwd())
        print("%s generate the graph.yaml file"%(os.getcwd()))
        
    elif args.shell :
       path = os.path.split(os.path.abspath(__file__))[0]
       path+='/shell.yaml'
       copy(path,os.getcwd())
       print("%s generate the shell.yaml file"%(os.getcwd()))

    elif args.input :
        shell, numberlist = ShellPath().check()
        while True :
            number = int(input("Please enter the number >> "))
            if number in numberlist :
                break
        modecheck = ['bulkmodulus','relaxation','chgcar','dos','band',
                    'effective_mass','dielectric','HSE','absorption']

        with open("input.gmd",'w') as fi :
            fi.write("### Perovgen input.gmd\n")
            fi.write("KPOINTS:\nA\nCONSTK=30\n\n")
            fi.write("INCAR:\nPBE\n\n")
            fi.write("SHELL:\n%s\n\n"%(shell[number-1][-1]))
            fi.write("CALMODE:\n")
            for m in modecheck :
                fi.write("%s = T\n"%(m))
        fi.close()
        print("%s generate the input.gmd file"%(os.getcwd()))

    elif args.pdos :
        width = args.pdos[0]
        s = load_structure(os.getcwd())[0][0]
        if not s :
            print("Structure file doesn't exist!\n") 
            sys.exit(1)
        electronic.GMDAnalysis().partialDOS(width=width,structure=s)    

    elif args.convert:
        files = load_structure(args.convert)[0]
        for e,f in enumerate(files) : 
            GMDStructure(f).convert_to_cif()
            
    elif args.hse :
        if "HSE" in os.listdir(os.getcwd()) :
            vrun = BSVasprun("HSE/vasprun.xml")
            bandgap = vrun.as_dict()['output']['bandgap']
        else :
            bandgap = None
            
        if "HSE_SOC" in os.listdir(os.getcwd()) :
            vrun = BSVasprun("HSE_SOC/vasprun.xml")
            bandgap_soc = vrun.as_dict()['output']['bandgap']
        else :
            bandgap_soc = None

        with open("hse_bandgap.gmd",'w') as fi :
            if bandgap == None :
                fi.write("HSE : None\n")
            else :
                fi.write("HSE : %.5f\n"%(bandgap))
            if bandgap_soc == None :
                fi.write("HSE_SOC : None\n")
            else :
                fi.write("HSE_SOC : %.5f\n"%(bandgap_soc))
            fi.close()
        