
# coding : utf-8
# Cppyright (c) Green Materials Designs Team.
'''
This module provides crystal data base
'''
import os
import sys
import time
import subprocess
import glob
from shutil import copyfile
import ast 

import yaml
import numpy as np
import pandas as pd
from collections import defaultdict, Counter

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.bandstructure import HighSymmKpath

from perovgen.pygmd.shell import ShellPath
from perovgen.pygmd.input_structure import load_structure

def read_input(inputpath):
    f = open(inputpath)
    ff = f.readlines()
    dic=defaultdict(list)
    repetit=[]

    for i in ff :
        n = i.split("\n")[0] # classficiation 
        if "#" in n :
            pass
        else :
            if ":" in n :
                name = n.split(":")[0]
                repetit.append(name)
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

    # Check the values
    if [*dic] :
        nonvalue=[i for i in [*dic] if not i in ["STRUC","INCAR","METHOD","SHELL","KPOINTS"]]
        if len(nonvalue) >= 3  :
            print(nonvalue, "doens't exist!")
            sys.exit(1)
    else :
        print(["METHOD", "INCAR", "SHELL"],"doens't exist!")
        sys.exit(1)
    
    # Remove repeating classes
    c=Counter(repetit)
    for k,v in c.items() :
        if v != 1 :
            print("{} overlapped. Please check input file".format(k))
            sys.exit(0)
    
    else :
        if len(dic["KPOINTS"]) != 1 :
            print("KPOINTS value have to be one!")
            sys.exit(0)
    # SHELL 
    if "SHELL" in [*dic] :
        if len(dic["SHELL"]) != 1 :
            print("SHELL value have to be one!")
            sys.exit(0)

    # MODE 
    if "METHOD" in [*dic] :
        modecheck = ["M","R","C","B","D","A","U","E",'G', # one mode
        "RC","CB","CD","DB","BD","BE","EU","CA",'CG', # two mode
        "RCB","RCD","CBD","CDB","DBE","CBE","BEU","RCA",'CBG', # THREE MODE
        "RCBD","RCDB","CBEU","DBEU","RCBE","CDBE",'RCBG', # four mode
        "RCBEU","RCDBE","CDBEU","RCBGE", # FIVE MODE
        "RCDBEU",'RCDBGE',# SIX MODE
        "RCADBEU"] # SEVEN MODE
        mn = ''
        for i in dic["METHOD"] :
            mn+=i
        if not mn in modecheck :
            print("There isn't {} mode\nPlease Check your mode".format(dic["METHOD"]))
            print("mode only :", modecheck)
            sys.exit(0)
        
    if len(nonvalue) != 0 :
        print(nonvalue, "doesn't exist! \nPlease enter the mode\n mode is",["KPOINTS","INCAR","METHOD","SHELL", "STRUC"])
        sys.exit(0)
    return dic

class GMDIncar :
    def __init__(self, gmdincar) :
        if gmdincar["INCAR"][0] in ["PBESol","VDW","SCAN","PBE"]:
            stream = open("%s/GMDincar.yaml"%(os.path.dirname(__file__)),'r')
            incar = yaml.load(stream)
            if gmdincar["INCAR"][0] == "PBE": 
                pass
            elif gmdincar["INCAR"][0] == "PBESol" :
                incar["GGA"]='Ps'
            elif gmdincar["INCAR"][0] == "VDW" :
                incar["IVDW"]=21
            elif gmdincar["INCAR"][0] == "SCAN" :
                incar["METAGGA"] = "SCAN"
                incar["LUSE_VDW"] = True
                incar["BPARAM"]=15.7
                incar["LASPH"]=True

            if len(gmdincar["INCAR"]) != 1 :
                incar1 = self._Incarmethod(gmdincar["INCAR"][1:])
                for k,v in incar1.items() :
                    incar[k]=v
            else :
                pass
        else :
            print("Please Enter the PBESol, PBE, VDW, SCAN, PBE, MPJ")
            sys.exit(1)
        self.incar = incar

    def _Incarmethod(self,incarstring) :
        incar = {}
        for inc in incarstring :
            n = inc.split("=")
            try :
                if type(int(n[1])) == int :
                    u = int(n[1])
                elif type(float(n[1])) == float :
                    u = float(n[1])
                incar[n[0]]=u
            except :
                if n[0] == "MAGMOM" :
                    incar['MAGMOM']=ast.literal_eval(n[1])
                else :
                    incar[n[0]]=n[1]
        return incar

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
            chg_path = "{}/CHGCAR".format(path)
        else :
            print("It isn't have CHGCAR file")
            sys.exit(0)
    return chg_path
            
def MakingKpointBand(structure,path):
    kpath_info = HighSymmKpath(structure).kpath['kpoints']
    kpath = HighSymmKpath(structure).kpath['path'][0]
    if len(kpath) <= 5 :
        kpath = HighSymmKpath(structure).kpath['path'][0]
        kpath.extend(HighSymmKpath(structure).kpath['path'][1])
    with open(path,'w') as fi :
        fi.write("kpoints\n21\nL\nR\n")
        for i in range(len(kpath)-1):
            fi.write("%.3f %.3f %.3f !%s\n"%(kpath_info[kpath[i]][0],kpath_info[kpath[i]][1],kpath_info[kpath[i]][2],kpath[i]))
            fi.write("%.3f %.3f %.3f !%s\n"%(kpath_info[kpath[i+1]][0],kpath_info[kpath[i+1]][1],kpath_info[kpath[i+1]][2],kpath[i+1]))
            fi.write("\n")        

def MakingInpcar(structure, path, nelect, kpoints):
    with open(path,'w') as fi :
        for i in kpoints :
            fi.write("%.3f "%(i))
        fi.write("\n0.01\n%i\nV\n"%(nelect))
        for i in structure.lattice.matrix :
            fi.write("%.5f %.5f %.5f\n"%(i[0],i[1],i[2]))
    fi.close()

class PerovInputs :
    def __init__(self, structure, is_selective=False):
        self.is_selective=is_selective
        for i in range(structure.num_sites):
            try : 
                structure.replace(i, structure.species[i].element)
            except : 
                structure.replace(i, structure.species[i])
        self.poscar = structure
        if self.is_selective :
            for i in range(self.poscar.num_sites):
                try :
                    self.poscar.replace(i,self.poscar.species[i].element, properties=None)
                except :
                    self.poscar.replace(i,self.poscar.species[i],properties=None)

    def _incarmode(incar, method):
        if method == "R" :
            pass
            #incar["EDIFF"] = 1E-6
        elif method == "C":
            incar["NSW"] = 0
            incar["LCHARG"] = True
            incar["EDIFF"] = 1E-6
        elif method == "B" or method == "E" or method == "D" or method == 'G':
            incar["NSW"] = 0
            incar["EDIFF"]=1E-6
            incar["LCHARG"] = False
            incar["ICHARG"] = 11
            incar["SIGMA"]=0.01
            incar["ISMEAR"] = 0
            if method == "D" :
                incar["SIGMA"]=0.01
                incar["ISMEAR"] = -5
        elif method == "U" :
            incar["LREAL"] = False
            incar["NSW"] = 0 
            incar["EDIFF"]=1E-8
            incar["ISMEAR"]=0
            incar["IBRION"]=8
            incar["LEPSILON"]=True
            incar["SIGMA"]=0.01
            if "NCORE" in [*incar] :
                del incar["NCORE"]
        elif method == "A" :
            incar["EDIFF"]=1E-8
            incar["LPEAD"]=True
            incar["NEDPS"]=2000
            incar["LOPTICS"]=True
            incar["CSHIFT"]=0.100
            incar["SIGMA"]=0.01
            incar["LSORBIT"]=True
        elif method == 'G' :
            incar['NKRED']=2
            incar['LHFCALC']=True
            incar['ALGO']='All'
            incar['ISYM']=0
            incar['SYMPREC']=1E-8
            incar['AEXX']=0.25
            incar['HFSCREEN']=0.2
            incar['TIME']=0.4
            incar['PRECFOCK']='FAST'
        #elif method == "M" :
        #    incar["EDIFF"]=1E-5
        #    incar["ENCUT"]=400
            
        return incar

    def inputfolder(self, incar, number):
        if "MAGMOM" in [*incar] :
            magmom = incar['MAGMOM']
            del incar['MAGMOM']
        else :
            magmom = None
        mpr = MPRelaxSet(self.poscar, user_incar_settings=incar,user_potcar_functional=None)
        vi = mpr.get_vasp_input()
        if number == 'G' :
            # number = 25 and even number
            lattice = []
            for l in ['a','b','c'] :
                kl = self.poscar.as_dict()['lattice'][l]
                if kl < 25 :
                    if round(25/kl)%2 == 0 :
                        lattice.append(round(25/kl))
                    else :
                        lattice.append(round(25/kl)+1)
                else :
                    lattice.append(1)
                    if not vi['INCAR']['NKREDX'] and vi['INCAR']['NKREDY'] :
                        vi['INCAR']['NKREDX']=2
                        vi['INCAR']['NKREDY']=2
            vi["KPOINTS"].kpts[0] = lattice
        elif number != None : 
            lattice = [self.poscar.as_dict()['lattice'][l] for l in ['a','b','c']]
            a,b,c = [i if i < number else number for i in lattice]
            vi["KPOINTS"].kpts[0] = [round(number/a),round(number/b),round(number/c)]
        if magmom :
            vi['INCAR']['MAGMOM'] = magmom
        return vi

class RunningShell :
    def __init__(self, shell, name, path):
        self.shell_path = ShellPath().shellpath
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
        self.path = path

    def write_vaspsh(self):
        with open("{}/vasp.sh".format(self.path),"w") as fi :
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
        with open("{}/vasp.sh".format(self.path),"w") as fi :
            for i in self.vaspsh_list :
                fi.write(i)
        fi.close()

    def running_mode(self,soc=False, run=True):
        pwd = os.getcwd()
        if soc : 
            #self.write_vaspsh()
            self.SOC_read_vaspsh()
        else :
            self.write_vaspsh()
        if run :
            os.chdir(self.path)
            subprocess.check_call(['qsub','vasp.sh'])
            #os.system("qsub < vasp.sh")
            os.chdir(pwd)