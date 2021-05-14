# coding : utf-8
# Cppyright (c) Green Materials Designs Team.

import os
import sys

import numpy as np
import random as rd
import subprocess
from collections import Counter
import yaml
from collections import defaultdict

from pymatgen.core import Structure, Element, Composition
from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.outputs import Vasprun, BSVasprun
from pymatgen.ext.matproj import MPRester
from shutil import copyfile

def load_structure(path, sformat=True):
    '''
    The structure file that exists in the PATH is read and return
    in the form of a list 
    
    Args :
        path(str) : OS.PATH
    '''
    struclist=[] ; spath=[]
    if type(path) == str :
        path = [path]
    for p in path :
        path1 = os.path.abspath(p)
        if os.path.isfile(path1):
            try :
                struclist.append(Structure.from_file(path1))
                spath.append(path1)
            except :
                pass
        else :
            for j in os.listdir(path1):
                try :
                    struclist.append(Structure.from_file("%s/%s"%(p,j)))
                    spath.append(j)
                except :
                    pass
    if sformat :
        return struclist
    else :
        return spath

def graphyaml(parameter):
    '''
    Auxiliary file for drawing BAND and DOS
    '''
    BD = dict(bs_projection = 'elements', 
			dos_projection =  'elements', 
			vb_energy_range= 4, 
			cb_energy_range= 4, 
			fixed_cb_energy= False, 
			egrid_interval= 2, 
			font= 'Arial', 
			axis_fontsize= 20, 
			tick_fontsize= 15, 
			legend_fontsize= 14, 
			bs_legend= 'best', 
			dos_legend= 'best', 
			rgb_legend= True, 
			fig_size= (11,8.5))
    D = dict(zero_to_efermi=True,
			stack=True,
			fig_size= (12,8),
			xlim= (-6,4), 
			ylim= None, 
			font_size= 30,
			color="r")
    B = dict(fig_size= (12,8), 
			zero_to_efermi= True, 
			fontsize= 20, 
			xlim= None,
			ylim= (-4,4), 
			color = 'b',
			vbm_cbm_marker= True, 
			linewidth= 1)
    if parameter == "BD" :
        return BD
    elif parameter == "B":
        return B
    elif parameter == "D":
        return D
    elif parameter == "All":
        print("Only BD, B and D options exist.\nPlease retry")

def createFolder(directory):
    '''
    Create the directory 
    '''
    try :
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError :
        print("Error : Creating directory. " + directory)

class MPJClass :
    def __init__(self):
        self.mpr=MPRester()
        
    def createStructure(mplist):
        s = []
        for mp in mplist :
            ss = self.mpr.get_structure_by_material_id(mp)
            name = ss.composition.get_reduced_formula_and_factor()[0]
            ss.to(filename="%s_%s.cif"%(name,mp))
        return s


class ShellPath :
    '''
    Designate the position of shell scripts and 
    check the shell scripts registered with gmd.
    '''
    def __init__(self) :
        self.shellpath = "%s/shell"%(os.path.split(os.path.dirname(__file__))[0])
        self.files = [f for f in os.listdir(self.shellpath)]

    def register_shell(self,shellpath):
        '''
        Save shell file
        '''
        if os.path.isfile(shellpath):
            filename = os.path.basename(shellpath)
            filepath = os.path.abspath("%s/%s"%(self.shellpath,filename))
            if os.path.isfile(filepath):
                print("%s is already enrolled\n"%(filename))
            else :
                copyfile(shellpath,filepath)
                print("%s Successfully Save "%(filepath))
        else :
            print("%s isn't file"%(shellpath))

    def check(self):
        '''
        Check shell file registered 
        '''
        shell=[f for f in self.files if f.split(".")[-1]== "sh"]
        if len(shell) != 1 :
            print("shell scripts are",shell)
        else:
            print("shell scirpt is",shell[0])

    def remove(self) :
        '''
        Remove the shell file registered
        '''
        shell=[f for f in self.files if f.split(".")[-1]== "sh"]
        print("list of registered shell script : ",shell)
        while True :
            r = str(input("Please enter the script name for removal >> "))
            if r in shell :
                break
            else :
                print("{} doens't exist in shell script".format(r))
        os.remove("{}/{}".format(self.shellpath, r))
        print("{} scirpt is deleted".format(r))

class GMDStructure :
    '''
    This module provides transforming 
    structure files(ex. cif, POSCAR).
    '''
    def __init__(self, structure):
        self.structure = structure
        self.coords = np.dot(structure.frac_coords, structure.lattice.matrix)
        self.species = structure.species

    def _options(self, delete_charge=True, delete_sd=False):
        for i in range(self.structure.num_sites):
            if delete_charge :
                try :
                    self.structure.replace(i, self.structure.species[i].element)
                except :
                    pass
            if delete_sd :
                try :
                    self.structure.replace(i,self.structure.species[i].element, properties=None)
                except :
                    self.structure.replace(i,self.structure.species[i],properties=None)
        return self.structure

    def _split_molecule(self):
        d = self.structure.distance_matrix
        hcn_coords = defaultdict(list)
        for i in range(len(self.coords)) :
            if self.species[i].symbol == "C" :
                hcn_coords["C"].append(i)
            elif self.species[i].symbol == "H":
                hcn_coords["H"].append(i)
            elif self.species[i].symbol == "N" :
                hcn_coords['N'].append(i)
        # the number H and N of the round C
        molecule = defaultdict(list)
        for c in hcn_coords['C'] :
            chbonding = np.where(d[c,hcn_coords['H'][0]:hcn_coords['H'][-1]+1] < 1.5)[0]
            cnbonding = np.where(d[c,hcn_coords['N'][0]:hcn_coords['N'][-1]+1] < 1.5)[0]
            if len(chbonding) == 0 and len(cnbonding) == 3 :
                molecule["GUA"].append(1)
            elif len(chbonding) == 1 and len(cnbonding) == 2 :
                nhbonding = np.where(d[hcn_coords['N'][0]+cnbonding[0],hcn_coords['H'][0]:hcn_coords['H'][-1]+1]<1.5)[0]
                if len(nhbonding) == 2 :
                    molecule["FA"].append(1)
                else :
                    molecule['Zolium'].append(1)
            elif len(chbonding) == 3 and len(cnbonding) == 1:
                nhbonding = np.where(d[hcn_coords['N'][0]+cnbonding[0],hcn_coords['H'][0]:hcn_coords['H'][-1]+1]<1.5)[0]
                if len(nhbonding) == 3 :
                    molecule["MA"].append(1)
                elif len(nhbonding) == 2 :
                    molecule['DMA'].append(1)
                elif len(nhbonding) == 1 :
                    molecule['triMA'].append(1)
                elif len(nhbonding) == 0 :
                    molecule['tetraMA'].append(1)
            else :
                molecule['H'].append(1)
                molecule['C'].append(1)
                molecule['N'].append(1)

        for k,v in molecule.items() :
            if k == "DMA" :
                molecule['DMA']=int(len(v)/2)
            elif k == 'triMA' :
                molecule['triMA']=int(len(v)/3)
            elif k == 'tetraMA':
                molecule['tetraMA']=int(len(v)/4)
            else :
                molecule[k]=len(v)
        return molecule
        
    def vaspname(self):
        sn = self.structure.composition.get_el_amt_dict()
        vaspname=''
        if 'C' in sn and 'H' in sn and 'N' in sn : 
            mole=self._split_molecule()
            for k,v in mole.items() :
                if v == 1 :
                    vaspname += k
                else :
                    vaspname += "%s%i"%(k,v)

                if k == 'FA' :
                    sn['C']-=v
                    sn['H']-=v*5
                    sn['N']-=v*2
                elif k == 'MA' :
                    sn['C']-=v
                    sn['N']-=v
                    sn['H']-=v*6
                elif k == 'GUA' :
                    sn['C']-=v
                    sn['N']-=v*3
                    sn['H']-=v*6
                elif k == 'DMA' :
                    sn['C']-=v*2
                    sn['N']-=v
                    sn['H']-=v*8
                elif k == 'triMA' :
                    sn['C']-=v*3
                    sn['N']-=v
                    sn['H']-=v*10
                elif k == 'tetraMA':
                    sn['C']-=v*4
                    sn['N']-=v
                    sn['H']-=v*12
                elif k == 'Zolium' :
                    sn['C']-=v*3
                    sn['N']-=v*2
                    sn['H']-=v*5
        for k,v in sn.items() :
            if int(v) == 1 :
                vaspname += "%s"%(k)
            elif int(v) <= 0 :
                pass
            else :
                vaspname += "%s%i"%(k,v) 
        return vaspname 

def _inputgmd(path):
    while True :
        shell = str(input("Please enter the shell name >> "))
        if shell in ShellPath().files :
            break
        else :
            print("There isn't %s in shell script"%(shell))
        
    modecheck = ["R","C","B","D","U",
        "RC","CB","CD","RU","BD","CU","BU","DU", # Two mode
        "RCB","CBD","RCU","RCD","BEH","BHE","CBD","CDB", # Three mode
        "RCBD","RCDB","CBEH","CBHE","BEHU","BHEU", # Four mode
        "RCBEH","RCBHE","CDBHE","CDBEH","RCDBU","CBEHU","CBHEU","DBEHU","DBHEU",# Five mode
        "RCDBEH","RCDBHE","CDBHEU", "CDBEHU",# Six mode
        "RCDBHEU","RCDBEHU"]
    while True :
        mn = str(input("Please enter the mode >> "))
        if not mn in modecheck :
            print("There isn't {} mode\nPlease Check your mode".format(mn))
            print("mode only :", modecheck)
        else :
            break
    with open("input.gmd","w") as fi :
        fi.write("#TEST FOR INPUTGMD\n")
        fi.write("STRUC:\n")
        for i in load_structure(path,sformat=False) :
            fi.write(i)
            fi.write("\n")
        fi.write("KPOINTS:\n40\n")
        fi.write("INCAR:\nMPJ\n")
        fi.write("SHELL:\n%s\n"%(shell))
        fi.write("METHOD:\n%s\n"%(mn))
    fi.close()

class pdos :
    def __init__(self):
        self.path = "%s/pdos/pdos"%(os.path.split(os.path.dirname(__file__))[0])
