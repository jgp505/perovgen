
import os
import sys

import numpy as np
import pandas as pd
from collections import Counter

from pymatgen.core import Structure, Element
from perovgen.pygmd.base import GMDStructure

'''
본 모듈은 mc 구조파일을 제공해주는 모듈입니다.
'''

def rotation(the, ph, ps):
    theta = np.radians(the)
    phi = np.radians(ph)
    psi = np.radians(ps)
    ct, st = np.cos(theta), np.sin(theta)
    cph, sph = np.cos(phi), np.sin(phi)
    cps, sps = np.cos(psi), np.sin(psi)
    R = np.array(((ct*cps, -ct*sps, st), (sph*st*cps+cph*sps, -sph*st*sps +cph*cps, -sph*ct), (-cph*st*cps+sph*sps, cph*st*sps+sph*cps, cph*ct)))
    return R

def _selective_bool(species, type_num):
    length = len(species)
    if type_num == 1:
        sp = [[True, True, True] for i in range(length)]
    elif type_num == 2:
        sp = [[False, False, False] for i in range(length)]
    elif type_num == 3:
        sp = []
        for i in range(length):
            if species[i].symbol in ["H", "C", "N"]:
                sp.append([False, False, False])
            else:
                sp.append([True, True, True])
    elif type_num == 4:
        sp = []
        for i in range(length):
            if species[i].symbol in ["H", "C", "N"]:
                sp.append([True, True, True])
            else:
                sp.append([False, False, False])
    return sp

class RandomPosition:
    def _loadmolecule(name):
        path = "%s/zero_coords"%(os.path.split(__file__)[0])
        if name == "MA":
            ms = Structure.from_file("%s/MA.json"%(path))
        elif name == "FA":
            ms = Structure.from_file("%s/FA.json"%(path))
        elif name == "GUA":
            ms = Structure.from_file("%s/GUA.json"%(path))
        elif name == "diMA":
            ms = Structure.from_file("%s/DiMA.json"%(path))
        elif name == "triMA":
            ms = Structure.from_file("%s/TriMA.json"%(path))
        elif name == "tetraMA":
            ms = Structure.from_file("%s/TetraMA.json"%(path))
        elif name == "Zolium":
            ms = Structure.from_file("%s/Imidazolium.json"%(path))
        return ms
    
    def __init__(self, random_coord=True, random_degree=True):
        self.random_coord = random_coord
        self.random_degree = random_degree
        if self.random_coord == False and self.random_degree == False :
            if os.path.isfile("%s/MOLE"%(os.getcwd())) :
                rmole = open("%s/MOLE"%(os.getcwd()),"r").readlines()
                self.cord = [int(i.split()[0]) for i in rmole]
                self.degree = [list(map(int,i.split()[1:])) for i in rmole]
            else :
                print("Please make the 'MOLE' file!\n")
                sys.exit(0)
        elif self.random_coord == False and self.random_degree == True :
            if os.path.isfile("%s/MOLE"%(os.getcwd())) :
                rmole = open("%s/MOLE"%(os.getcwd()),"r").readlines()
                self.cord = [int(i.split()[0])-1 for i in rmole]
                self.degree = None
            else :
                print("Please make the 'MOLE' file!\n")
                sys.exit(0)
        else :
            self.cord = None
            self.degree = None


    def _inputinform(self):
        while True :
            inputatom = str(input("Enter the elements >> "))
            try :
                Element(inputatom)
                break
            except :
                print("\nThere is no %s in structure\nPlease one more enter\n"%(inputatom))
        if self.cord != None :
            change = 1
        else :
            change = float(input("Enter the change ratio(ex.0.7) >> "))
            if change > 1 :
                change = 1

        while True :
            inputfix = str(input("Do you want to apply selective dynamic?(Y or N) >> "))
            if inputfix == "Y":
                print("\n#########################################################################################################")
                print("#\t1. Do you want to apply T to all elements? Please enter 1\t\t\t\t\t#")
                print("#\t2. Do you want to apply F to all elements? Please enter 2\t\t\t\t\t#")
                print("#\t3. Do you want to apply T to the inorganic element and F to the organic element? Please enter 3 #")
                print("#\t4. Do you want to apply F to the inorganic element and T to the organic element? Please enter 4 #")
                print("#########################################################################################################\n")
                while True :
                    fixcalc = int(input("Please refer to the options above >> "))
                    if fixcalc in [1,2,3,4]:
                        break
                    else :
                        print("There isn't in option")
                break
            elif inputfix == "N" :
                fixcalc = None
                break
            else :
                print("Please enter the Y or N\n")
        if self.random_coord == False and self.random_degree == False :
            multiple = 1
        elif self.random_coord == True and self.random_degree == False :
            self.degree = np.array(list(map(int,input("Please enter the degree >> ").split())))
            if change != 1 :
                multiple = int(input("Enter the number of the times to repeat >> "))
            else :
                multiple = 1
        else :
            multiple = int(input("Enter the number of the times to repeat >> "))
        return inputatom, change, fixcalc, multiple

    def tiltingmolecule(self, s, ms, inputatom, changenum, fixcalc=None):
        s_matrix = s.lattice.matrix
        s_species = s.species
        ms_species = ms.species
        s_coord = np.dot(s.frac_coords, s_matrix).tolist()
        ms_coord = np.dot(ms.frac_coords, ms.lattice.matrix)

        name = GMDStructure(s).vaspname()
        # pick the index of inputatom and shuffle 
        s_index = np.array([e for e,c in enumerate(s_species) if c.symbol == Element(inputatom).symbol])
        if self.random_coord :
            np.random.shuffle(s_index)
            self.cord = s_index[:round(changenum*len(s_index))]
        else :
            for r in self.cord :
                if not r in s_index :
                    print("%i isn't %s index\nPlease revise the index number\n" % (r, inputatom))
                    sys.exit(0)
        r_coord = [s_coord[f] for f in self.cord]

        # make array of degree 
        if self.random_degree :
            self.degree = np.random.randint(360, size=(len(r_coord),3))
            new_ = self.degree
        else :
            if not self.random_coord :
                self.degree = np.array(self.degree)
                new_ = self.degree
            else:
                new_ = [[self.degree[0], self.degree[1], self.degree[2]] for i in range(len(r_coord))]
        # add coordination
        lthe=[];lph=[];lps=[]
        for coord, degree in zip(r_coord,new_):
            index = s_coord.index(coord)
            del s_coord[index]
            del s_species[index]

            the, ph, ps = degree[0], degree[1], degree[2]
            lthe.append(degree[0])
            lph.append(degree[1])
            lps.append(degree[2])
            R = rotation(the, ph, ps)
            dot = np.dot(ms_coord, R)
            s_coord.extend(np.add(dot, coord).tolist())
            s_species.extend(ms_species)
        # csv 
        df = pd.DataFrame({'index':self.cord,'Theta':lthe,'Phi':lph,'Psi':lps})
        df = df.set_index('index')

        # add the properties of selective_dynamics 
        if fixcalc != None :
            sp = _selective_bool(species=s_species, type_num=fixcalc)
            new_s = Structure(s_matrix, s_species, s_coord, coords_are_cartesian=True,site_properties={"selective_dynamics":sp})
        else :
            new_s = Structure(s_matrix, s_species, s_coord, coords_are_cartesian=True)
        new_s.sort()
        return new_s, df