# coding: utf-8
# Cppyright (c) Green Materials Designs Team.
# Hanbat National University, Korea 

import os
import sys

import numpy as np
import pandas as pd

from collections import defaultdict
from mendeleev.fetch import fetch_table
from pymatgen.core import Structure

def load_structure(path):
    '''
    The structure file that exists in the PATH is read and return
    in the form of a list 
    
    Args :
        path(str) : OS.PATH
    '''
    struclist=[] ; spath=[]
    if type(path) == str :
        path = [path]

    if path :
        for p in path :
            path1 = os.path.abspath(p)
            if os.path.isfile(path1):
                try :
                    s = Structure.from_file(path1)
                    struclist.append(s)
                    spath.append(path1)
                except :
                    pass
            else :
                for j in os.listdir(path1):
                    try :
                        s = Structure.from_file("%s/%s"%(p,j))
                        struclist.append(s)
                        spath.append(j)
                    except :
                        pass
    return struclist, spath

class GMDStructure :
    """
    This module provides structural files(ex. cif, POSCAR) 
    needed for automatic calculation
        """
    def __init__(self, structure):
        self.structure = structure
        self.structure.sort()
        self.coords = np.dot(self.structure.frac_coords, self.structure.lattice.matrix)
        self.species = self.structure.species
        
    def organizing(self, sort=None) :
        new_species = []
        for e in self.structure.species :
            if not e.symbol in new_species :
                new_species.append(e.symbol)
        element = fetch_table("elements")[['atomic_number','symbol']].set_index("symbol").to_dict()['atomic_number'] 
        if sort != None :
            element = dict()
            for i,e in enumerate(sort) :
                if not e in new_species :
                    return False
                else :
                    element[e]=i
        new_species.sort(key=lambda x : element[x])
        indexes1 = defaultdict(list)
        for i in range(self.structure.num_sites) :
            indexes1[self.structure.species[i].symbol].append(i)
    
        indexes = dict()
        for k1 in new_species :
            indexes[k1]=indexes1[k1]
    
        struc_data = defaultdict(list)
        for k in new_species :
            for v in indexes[k] :
                struc_data['species'].append(self.structure[v].specie)
                struc_data['coords'].append(self.structure[v].coords.tolist())
                if self.structure[v].properties :
                    for k1 in [*self.structure[v].properties] :
                        struc_data[k1].append(self.structure[v].properties[k1])
        if len([*struc_data]) == 2 :
            struc = Structure(self.structure.lattice.matrix, struc_data['species'], struc_data['coords'], coords_are_cartesian=True)
        else :
            struc = Structure(self.structure.lattice.matrix, struc_data['species'], struc_data['coords'], 
                          site_properties={'selective_dynamics':struc_data['selective_dynamics']}, coords_are_cartesian=True)
        return struc

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
    
    def formula_dict(self,reduced=True) :
        if reduced :
            sn = self.structure.composition.to_reduced_dict
        else :
            sn = self.structure.composition.get_el_amt_dict()
        sn_dict = dict()
        if 'C' in sn and 'H' in sn and 'N' in sn : 
            mole=self._split_molecule()
            for k,v in mole.items() :
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
                sn_dict[k]=v
            if sn['C'] <= 0 :
                del sn['C']
            else :
                sn_dict['C'] = sn['C']
            if sn['H'] <= 0 :
                del sn['H']
            else :
                sn_dict['H'] = sn['H']
            if sn['N'] <= 0 :
                del sn['N']
            else :
                sn_dict['N'] = sn['N']
        for k,v in sn.items():
            sn_dict[k]=v
        return sn_dict
    
    def formula(self, reduced=True):
        if reduced :
            sn = self.structure.composition.to_reduced_dict
        else :
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

    def convert_to_cif(self, name=False) :
        '''
        convert POSCAR to cif. if name is false, file name is [formula]_[symmetry number].cif.
        
        name (str) = designate the name
        '''
        if not name :
            try :
                symmetry, groupnumber = self.structure.get_space_group_info()
            except :
                groupnumber = 0
            formula = self.formula(reduced=True)
            name = "{0}_{1:03d}".format(formula, groupnumber)
        return self.structure.to(filename="{}.cif".format(name))