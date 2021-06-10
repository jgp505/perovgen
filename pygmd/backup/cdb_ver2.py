import os
import sys

import numpy as np
import pandas as pd
from collections import defaultdict

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.bandstructure import HighSymmKpath

from perovgen.pygmd.config import ShellPath
from perovgen.pygmd.base import GMDStructure, load_structure

def read_input(inputpath):
    f = open(inputpath)
    ff = f.readlines()
    dic=defaultdict(list)
    
    for i in ff :
        n = i.split("\n")[0] # classficiation 
        if "#" in n :
            pass
        else :
            if ":" in n :
                name = n.split(":")[0]
            else :
                if n == '' : # delete the blank
                    pass
                else :
                    try :
                        if type(int(n)) == int :
                            u = int(n)
                        elif type(float(n)) == float :
                            u = float(n)
                        elif type(bool(n)) == bool :
                            if n == "Ture" or n == "ture" :
                                u = True
                            else :
                                u = False
                        dic[name].append(u)
                    except :
                        if name == "METHOD":
                            dic[name].extend(n)
                        else :
                            dic[name].append(n)

    # KPOINTS                
    if not "KPOINTS" in [*dic] :
        dic["KPOINTS"]=40
    else :
        if len(dic["KPOINTS"]) != 1 :
            print("KPOINTS value have to be one!")
            sys.exit(0)
    # SHELL 
    if "SHELL" in [*dic] :
        if len(dic["SHELL"]) != 1 :
            print("SHELL value have to be one!")
            sys.exit(0)

    # Check the values
    nonvalue=[]
    for i in [*dic]:
        if not i in ["KPOINTS","INCAR","METHOD","SHELL", "STRUC"] :
            if i == "STRUC" :
                pass
            else :
                nonvalue.append(i)
    if len(nonvalue) != 0 :
        print(nonvalue, "doesn't exist! \nPlease enter the mode\n mode is",["KPOINTS","INCAR","METHOD","SHELL", "STRUC"])
        sys.exit(0)
    return dic

class RunningShell :
    def __init__(self, shell, name):
        self.shell_path = "{}".format(ShellPath().shellpath)
        print(self.shell_path)
        try :
            vaspsh = open("{}/{}".format(self.shell_path,shell),"r")
        except :
            print("There isn't {} file".format(shell))
            sys.exit(1)
        self.vaspsh_list = vaspsh.readlines()
        for i in self.vaspsh_list :
            if "-N" in i :
                nameline = i.split()
                name_index = self.vaspsh_list.index(i)
                del self.vaspsh_list[name_index]
                del nameline[-1]
        lines2=''
        for i in nameline :
            lines2 += i+""
        lines2 += name + '\n'
        self.vaspsh_list.insert(name_index,lines2)

    def write_vaspsh(self):
        with open("vasp.sh","w") as fi :
            for i in self.vaspsh_list :
                fi.write(i)
        fi.close()
        
    def SOC_read_vaspsh(self):
        for i in range(len(self.vaspsh_list)):
            if 'std' in self.vaspsh_list[i] :
                std_line = self.vaspsh_list[i].split("#")
                if len(std_line) == 1 :
                    a = list(std_line[0])
                    a.insert(0,"#")
                    aa=""
                    for j in a :
                        aa += j
                    self.vaspsh_list[i] = aa
                else :
                    pass
            elif 'ncl' in self.vaspsh_list[i] :
                ncl_line = self.vaspsh_list[i].split("#")
                if len(ncl_line) == 1 :
                    a = list(std_line[0])
                    a.insert(0,"#")
                    aa=""
                    for j in a :
                        aa += j
                    self.vaspsh_list[i] = aa
                else :
                    self.vaspsh_list[i] = ncl_line[-1]
            if 'gam' in self.vaspsh_list[i] :
                std_line = self.vaspsh_list[i].split("#")
                if len(std_line) == 1 :
                    a = list(std_line[0])
                    a.insert(0,"#")
                    aa=""
                    for j in a :
                        aa += j
                    self.vaspsh_list[i] = aa
                else :
                    pass
        with open("vasp.sh","w") as fi :
            for i in self.vaspsh_list :
                fi.write(i)
        fi.close()

class PerovInputs :
    def __init__(self,inputfiles):
        self.inputfiles = read_input(inputfiles)
        self.poscar = load_structure(self.inputfiles["STRUC"])
        self.incar = self.inputfiles["INCAR"]
        self.kpoints = self.inputfiles["KPOINTS"]
        self.shell = self.inputfiles["SHELL"]
        self.method = self.inputfiles["METHOD"]
    
    def _poscar(self,delete_selective_dynamics=False):
        for e,p in enumerate(self.poscar) :
            for i in range(p.num_sites):
                try : 
                    p.replace(i, p.species[i].element)
                except :
                    pass
                if delete_selective_dynamics :
                    try :
                        p.replace(i,p.species[i].element,properties=None)
                    except :
                        p.replace(i,p.species[i],properties=None)
                self.poscar[e]=p
        return self.poscar
    
    def _incar(self):
        user_incar_settings={"SYSTEM":"Structure Optimization","PREC":"Accurate","ISTART":0,"ISPIN":1,"LREAL":"A",
                                 "ENCUT":520,"IBRION":2,"ISIF":3,"EDIFF":1E-4,"EDIFFG":-1E-2,"NSW":500,"ALGO":"Normal",
                                 "LCHARG":False,"LWAVE":False,"SIGMA":0.05,'ISMEAR':0}
        if self.incar[0] == "default":
            while True :
                method = input(str("Please Enter the conditions >> "))
                if method not in ["PBE","PBESol","VDW","SCAN","MP"] :
                    print("Please Retry\n You only enter the PBE,PBESol,VDW,SCAN, and MP")
                else :
                    break
            if method == "PBE":
                pass
            elif method == "PBESol":
                user_incar_settings["GGA"]="PS"
            elif method == "IVDW" :
                user_incar_settings["IVDW"]=21
            elif method == "SCAN" :
                user_incar_settings["METAGGA"]="SCAN"
            elif method == "MP" :
                user_incar_settings["MP"]=None
                
        elif self.incar[0] == "PBE":
            pass
        elif self.incar[0]  == "PBESol":
            user_incar_settings["GGA"]="PS"
        elif self.incar[0] == "IVDW" :
            user_incar_settings["IVDW"]=21
        elif self.incar[0]  == "SCAN" :
            user_incar_settings["METAGGA"]="SCAN"
        else :
            user_incar_settings={"EDIFFG":1E-2,"EDIFF":1E-4,"NSW":500,"LCHARG":False,"ICHARG":0}
            for inc in self.incar :
                n = inc.split("=")
                try :
                    if type(int(n[1])) == int :
                        u = int(n[1])
                    elif type(float(n[1])) == float :
                        u = float(n[1])
                    user_incar_settings[n[0]]=u
                except :
                    user_incar_settings[n[0]]=n[1]
        return user_incar_settings
    
    def calculation_process(self,method):
        user_incar_settings=self._incar()
        if method == "R":
            del user_incar_settings["EDIFFG"]
        elif method == "C":
            user_incar_settings["NSW"] = 0
            user_incar_settings["LCHARG"] = True
            user_incar_settings["EDIFF"] = 1E-6
        elif method == "B" or method == "E":
            user_incar_settings["EDIFF"]=1E-6
            user_incar_settings["NSW"]=0
            user_incar_settings["LCHARG"] = False
            user_incar_settings["ICHARG"] = 11
        return user_incar_settings
        
    def CopyCHGCAR(path):
        b = os.path.abspath(path)
        if os.path.isfile(b):
            ops = os.path.split(b)[0]
            if "CHGCAR" in os.listdir(ops) :
                chg_path = "{}/CHGCAR".format(ops)
            else :
               print("It isn't have CHGCAR file")
               sys.exit(0)
        else :
            if "CHGCAR" in os.listdir(b) :
                chg_path = "{}/CHGCAR".format(chg_path)
            else :
               print("It isn't have CHGCAR file")
               sys.exit(0)
        return chg_path
            
    def MakingKpointBand(structure,path):
        kpath_info = HighSymmKpath(structure).kpath['kpoints']
        kpath = HighSymmKpath(structure).kpath['path'][0]
        with open(path,'w') as fi :
            fi.write("kpoints\n21\nL\nR\n")
            for i in range(len(kpath)-1):
                fi.write("%.3f %.3f %.3f !%s\n"%(kpath_info[kpath[i]][0],kpath_info[kpath[i]][1],kpath_info[kpath[i]][2],kpath[i]))
                fi.write("%.3f %.3f %.3f !%s\n"%(kpath_info[kpath[i+1]][0],kpath_info[kpath[i+1]][1],kpath_info[kpath[i+1]][2],kpath[i+1]))
                fi.write("\n")        