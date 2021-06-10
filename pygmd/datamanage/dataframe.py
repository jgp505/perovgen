import os
import sys

import collections 
import shutil

import numpy as np
import pandas as pd 

from pymatgen.core import Structure, Element
from pymatgen.io.vasp.outputs import Vasprun, BSVasprun

from perovgen.pygmd.base import createFolder

def average_bond(structure, atom1, atom2):
    dm = structure.distance_matrix
    index1 = [i for i in range(len(structure.species)) if Element(atom1) == structure.species[i]]
    index2 = [i for i in range(len(structure.species)) if Element(atom2) == structure.species[i]]
    avg_length=[]
    for i in index1:
        list1=[]
        for j in index2 :
            if dm[i,j] < 4 :
                list1.append(dm[i,j])
        avg_length.append(np.array(list1).mean())
    return avg_length

def OrganizeData(path,readonly=True):
    pwd = os.getcwd()

    modecheck = collections.defaultdict(list)
    formula = path.split(os.path.sep)[-1]

    # create first folder 
    firstname = "{}/{}".format(pwd, formula)
    if not os.path.isdir(firstname) :
        os.makedirs(firstname)
    else :
        pass

    for f in os.listdir(path) :
        if "mode" in f :
            name, number, mode= f.split("_")[:-1]
            savefile_list = ["POSCAR","POTCAR","INCAR","KPOINTS", # input files
                    "CONTCAR","OUTCAR","OSZICAR","vasprun.xml"] # output files
            if mode == "C" :
                savefile_list.append("CHGCAR")
            elif mode == "E" or mode == "H" :
                savefile_list.append("EM")

            # create second folder
            secondname = "{}/{}_{}".format(firstname, name, number)
            if not os.path.isdir(secondname):
                os.makedirs(secondname)
            else :
                # create third folder
                thirdfolder = "{}/{}_mode".format(secondname,mode)
                if not os.path.isdir(thirdfolder) :
                    os.makedirs(thirdfolder)
                    for s in savefile_list :
                        shutil.copy("{}/{}/{}".format(path,f,s),"{}/{}".format(thirdfolder,s))
                        if readonly :
                            os.system("chmod 744 {}/{}".format(thirdfolder,s))
                else :
                    print("{} is already exists".format(thirdfolder))