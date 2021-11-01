
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
from collections import defaultdict

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Kpoints, Incar, Kpoints_supported_modes

from perovgen.pygmd.shell import ShellPath
from perovgen.pygmd.input_structure import load_structure, GMDStructure

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
            
def MakingInpcar(structure, path, nelect, kpoints):
    with open(path,'w') as fi :
        for i in kpoints :
            fi.write("%.3f "%(i))
        fi.write("\n0.01\n%i\nV\n"%(nelect))
        for i in structure.lattice.matrix :
            fi.write("%.5f %.5f %.5f\n"%(i[0],i[1],i[2]))
    fi.close()

class inputgmd :
    def GMDKpoints(self,kpoints) :
        kpoint={"POINT":list(kpoints[0])[0]}
            
        for k in kpoints[1:] :
            n,v = k.split("=")
            n = n.replace(" ","")
            v = v.replace(" ","")
            if n == 'CONSTK':
                kpoint[n]=int(v)
            elif n == 'KPTS' :
                kpoint[n]=ast.literal_eval(v)
                    
        if not "CONSTK" in [*kpoint] and "KPTS" in [*kpoint] :
            kpoint['CONSTK']=False
        elif "CONSTK" in [*kpoint] and not "KPTS" in [*kpoint]:
            kpoint['KPTS']=False
        return kpoint
    
    def string_to_dict(self,incars) :
        incar = {}
        for inc in incars :
            n,v = inc.split("=")
            n = n.replace(" ","")
            v = v.replace(" ","")
            
            try :
                if type(int(v)) == int :
                    v = int(v)
                elif type(float(v)) == float :
                    v = float(v)
                incar[n]=v
                
            except :
                if n == "MAGMOM" : # list type
                    incar['MAGMOM']=ast.literal_eval(v)
                else :
                    incar[n]=v
        return incar
    
    def GMDIncar(self,incar) :
        if incar[0] == "MPJ" :   
            incars = {"ISMEAR":1}
        else :
            stream = open("%s/GMDincar.yaml"%(os.path.dirname(__file__)),'r')
            incars = yaml.load(stream)
            if incar[0] == "PBE" : 
                pass
            elif incar[0] == "PBESol" :
                incars["GGA"]='Ps'
            elif incar[0] == "VDW" :
                incars["IVDW"]=21
            elif incar[0] == "SCAN" :
                incars["METAGGA"] = "SCAN"
                incars["LUSE_VDW"] = True
                incars["BPARAM"]=15.7
                incars["LASPH"]=True
                
        if len(incar[1:]) != 0 :
            incar1 = self.string_to_dict(incar[1:])
            for k,v in incar1.items() :
                incars[k]=v
        return incars
    
    def __init__(self, path) :
        self.path = path
        f = open(path,'r')
        ff = f.readlines()
        self.inputgmd = defaultdict(list)
    
        self.class_type_list = ['KPOINTS','INCAR','SHELL','CALMODE']
        self.modecheck = ["M","R","C","B","D","A","U","E",'G', # one mode
            "RC","CB","CD","DB","BD","BE","EU","CA",'CG','RG', # two mode
            "RCB","RCD","CBD","CDB","DBE","CBE","BEU","RCA",'CBG', # THREE MODE
            "RCBD","RCDB","CBEU","DBEU","RCBE","CDBE",'RCBG', # four mode
            "RCBEU","RCDBE","CDBEU","RCBGE", # FIVE MODE
            "RCDBEU",'RCDBGE','RGCDBE',# SIX MODE
            "RCADBEU"] # SEVEN MODE
    
        for i in ff :
            line = i.split("\n")[0]
            if not '#' in line :
                if len(line.split(":")) == 2 :
                    if line.split(":")[0] in self.class_type_list :
                        classname = line.split(":")[0]
                    else :
                        print(line.split(":")[0],"doesn't exist class in the input.gmd")
                        sys.exit(1)
                else :
                    if not line == '' :
                        self.inputgmd[classname].append(line)

        # Check the KPOINTS class
        if not list(self.inputgmd['KPOINTS'][0])[0] in ["A", "G", "M"] :
            print("KPOINTS first line have to type A(Auto) or G(Gamma) or M(Monkhorst)!")
            sys.exit(1)
        elif len(self.inputgmd['KPOINTS']) > 4 :
            print("KPOINTS class can only have 3 or less values.")
            sys.exit(1)
        else :
            self.kpoints = self.GMDKpoints(self.inputgmd['KPOINTS'])
        
        # Check the SHELL class
        if len(self.inputgmd['SHELL']) != 1:
            print("SHELL class must have only one value!")
            sys.exit(1)
        else :
            self.shell = self.inputgmd['SHELL'][0]
        
        # Check the CALMODE class 
        if not self.inputgmd['CALMODE'][0] in self.modecheck :
            print("CALMODE calss is wrong mode!")
            sys.exit(1)
        else :
            self.calmode = list(self.inputgmd['CALMODE'][0])

        # Check the INCAR class
        if not self.inputgmd['INCAR'][0] in ['PBE','PBESol','VDW','SCAN','MPJ'] :
            print("INCAR first line have to type PBE or PBESol or VDW or SCAN or MPJ")
            sys.exit(1)
        else :  
            self.exchange_corr = self.inputgmd['INCAR'][0]
            self.incar = self.GMDIncar(self.inputgmd['INCAR'])
            
        self.inputgmd['KPOINTS'] = self.kpoints
        self.inputgmd['INCAR'] = self.incar
        self.inputgmd['SHELL'] = self.shell
        self.inputgmd['CALMODE'] = self.calmode

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
        full_formula = GMDStructure(self.poscar).formula(reduced=False)
        try :
            symmetry, groupnumber = self.poscar.get_space_group_info()
        except :
            groupnumber = 0
        self.naming = "{0}_{1:03d}".format(full_formula, groupnumber)

    def incarmode(self,incar, method):
        if method == "R" :
            pass
        elif method == "C":
            incar["NSW"] = 0
            incar["LCHARG"] = True
            incar["EDIFF"] = 1E-6
        elif method == "B" or method == "E" or method == "D":
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
            incar['ALGO']='D' # ver 3.6.8 
            incar['NSW']=1 # one-shot hybrid 
            incar['ISYM']=0
            incar['SYMPREC']=1E-8
            incar['AEXX']=0.25
            incar['HFSCREEN']=0.2
            incar['TIME']=0.4
            incar['PRECFOCK']='FAST'
        '''
        elif method == "M" :
            incar["EDIFF"]=1E-5
            incar["ENCUT"]=400
        '''
        return incar

    def inputfolder(self, inputs, method, soc=False):
        incar = self.incarmode(incar=inputs.incar,method=method) ; kpoints = inputs.kpoints
        mpr = MPRelaxSet(self.poscar, user_incar_settings=incar,user_potcar_functional=None) #"PBE_52")
        vi = mpr.get_vasp_input()
        if soc :
            del vi['INCAR']['MAGMOM']
            vi['INCAR']['LSORBIT'] = True
        vi['KPOINTS'].comment = "{0}_{1}".format(method, self.naming)

        if inputs.exchange_corr != 'MPJ' :
            warning = '''
                    [Warning] both CONSTK and KPTS exist in KPOINTS class in input.gmd.\n
                    It is reflected as CONSTK.
                    '''
            if kpoints['CONSTK'] and kpoints['KPTS'] :
                print(warning)

            # KPOINTS mode
            if kpoints['POINT'] == 'G' :
                vi['KPOINTS'].style = Kpoints_supported_modes.Gamma
            elif kpoints['POINT'] == 'M' :
                vi['KPOINTS'].style = Kpoints_supported_modes.Monkhorst
            elif kpoints['POINT'] == 'A' :
                pass

            # HSE06 and DOS calculation setting Gamma and kpoints 
            if method == "G" :
                vi['KPOINTS'].style = Kpoints_supported_modes.Gamma
                # CONSTK = 25 and even number
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
                        vi['INCAR']['NKREDX']=2
                        vi['INCAR']['NKREDY']=2
                vi["KPOINTS"].kpts[0] = lattice
            elif method == 'D' :
                vi['KPOINTS'].style = Kpoints_supported_modes.Gamma
                if (kpoints['CONSTK'] and kpoints['KPTS']) or kpoints['CONSTK']:
                    number = kpoints['CONSTK']*2
                    lattice = [self.poscar.as_dict()['lattice'][l] for l in ['a','b','c']]
                    a,b,c = [i if i < number else number for i in lattice]
                    vi["KPOINTS"].kpts[0] = [round(number/a),round(number/b),round(number/c)]
                elif kpoints['KPTS'] :
                    vi['KPOINTS'].kpts[0] = np.array(kpoints['KPTS'])*2
            elif method == 'B' :
                bandinfo = HighSymmKpath(self.poscar)
                vi['KPOINTS']=vi['KPOINTS'].automatic_linemode(divisions=21,ibz=bandinfo)
            else :
                if (kpoints['CONSTK'] and kpoints['KPTS']) or kpoints['CONSTK']:
                    number = kpoints['CONSTK']
                    lattice = [self.poscar.as_dict()['lattice'][l] for l in ['a','b','c']]
                    a,b,c = [i if i < number else number for i in lattice]
                    vi["KPOINTS"].kpts[0] = [round(number/a),round(number/b),round(number/c)]
                elif kpoints['KPTS'] :
                    vi['KPOINTS'].kpts[0] = kpoints['KPTS']
        else :
            # HSE06 and DOS calculation setting Gamma and kpoints 
            if method == "G" :
                vi['KPOINTS'].style = Kpoints_supported_modes.Gamma
                # CONSTK = 25 and even number
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
                        vi['INCAR']['NKREDX']=2
                        vi['INCAR']['NKREDY']=2
                vi["KPOINTS"].kpts[0] = lattice
            elif method == 'D' :
                vi['KPOINTS'].style = Kpoints_supported_modes.Gamma
            elif method == 'B' :
                bandinfo = HighSymmKpath(self.poscar)
                vi['KPOINTS']=vi['KPOINTS'].automatic_linemode(divisions=21,ibz=bandinfo)
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
            lines2 += i+" "
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
            self.SOC_read_vaspsh()
        else :
            self.write_vaspsh()

        if run :
            os.chdir(self.path)
            subprocess.check_call(['qsub','vasp.sh'])
            #os.system("qsub < vasp.sh")
            os.chdir(pwd)
