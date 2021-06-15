import os
import sys

import collections 
import shutil

import numpy as np
import pandas as pd 

from pymatgen.core import Structure, Element
from pymatgen.io.vasp.outputs import Vasprun, BSVasprun
import xml.etree.ElementTree as ET

from perovgen.pygmd.base import GMDStructure, createFolder
from perovgen.pygmd.analysis.electronic import BSPlotting
from perovgen.pygmd.analysis.energy import GMDPhaseDiagram, GMDExcitonbinding

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
            elif mode == "D" :
                savefile_list.append("PROCAR")

            # create second folder
            secondname = "{}/{}_{}".format(firstname, name, number)
            if not os.path.isdir(secondname):
                os.makedirs(secondname)
                
            # create third folder
            thirdfolder = "{}/{}_mode".format(secondname,mode)
            if not os.path.isdir(thirdfolder) :
                os.makedirs(thirdfolder)
                for s in savefile_list :
                    try :
                        shutil.copy("{}/{}/{}".format(path,f,s),"{}/{}".format(thirdfolder,s))
                    except FileNotFoundError :
                        pass
                    if readonly :
                        os.system("chmod 744 {}/{}".format(thirdfolder,s))
            else :
                print("{} is already exists".format(thirdfolder))

class GMDDataFrame :
    def __init__(self,path) :
        self.path = path
        self.vrun = Vasprun("{}/vasprun.xml".format(self.path),parse_potcar_file=True)
        self.dataframe = collections.defaultdict(list)

    def UnitCellInfo(self) :
        s = self.vrun.final_structure 
        self.dataframe["Formula"].append(GMDStructure(s).vaspname())
        self.dataframe["Formula"].append(GMDStructure(s).vaspname())
        self.dataframe["a"].append(s.lattice.a)
        self.dataframe["b"].append(s.lattice.b)
        self.dataframe["c"].append(s.lattice.c)
        self.dataframe["alpha"].append(s.lattice.alpha)
        self.dataframe["beta"].append(s.lattice.beta)
        self.dataframe["gamma"].append(s.lattice.gamma)
        self.dataframe["volume"].append(s.volume)
        space,number = s.get_space_group_info()
        self.dataframe["Symmetry"].append(space)
        self.dataframe["SPN"].append(number)
        self.dataframe["Energy"].append("%.5f"%(self.vrun.final_energy))
        entries,pd1 = GMDPhaseDiagram(["{}/vasprun.xml".format(self.path)]).get_phasediagram()
        h = pd1.get_e_above_hull(entries[0])
        d = pd1.get_form_energy_per_atom(entries[0])
        self.dataframe["E_hull"].append(h)
        self.dataframe["E_form"].append(d)

        element = list(s.composition.get_el_amt_dict().items())
        oxistate = s.composition.oxi_state_guesses()[0]
        self.dataframe["Asite"].append(element[0][0])
        self.dataframe["A_num"].append(element[0][1])
        self.dataframe["Bsite"].append(element[1][0])
        self.dataframe["B_num"].append(element[1][1])
        self.dataframe["Xsite"].append(element[2][0])
        self.dataframe["X_num"].append(element[2][1])
        self.dataframe["Asite_oxistate"].append(oxistate[element[0][0]])
        self.dataframe["Bsite_oxistate"].append(oxistate[element[1][0]])
        self.dataframe["Xsite_oxistate"].append(oxistate[element[2][0]])
        self.dataframe["W"].append("%.5f"%(s.composition.weight))
        self.dataframe["total_electrons"].append(s.composition.total_electrons)
        self.dataframe["electronegativity"].append(s.composition.average_electroneg)
        return self.dataframe

    def BandInfo(self) :
        bsp = BSPlotting(vasprun="{}/vasprun.xml".format(self.path),
            kpoints = "{}/KPOINTS".format(self.path))
        bi = bsp._bandinform()
        self.dataframe["Ef"].append(bsp.bs.efermi)
        self.dataframe["Eg_direct"].append(bi['E_g']['Direct'])
        if bsp.bsdict['band_gap']['direct'] :
            self.dataframe["Eg_indirect"].append(bi['E_g']['Direct'])
        else :
            self.dataframe["Eg_indirect"].append(bi['E_g']['Indirect'])
        self.dataframe["Eg"].append(bsp.bsdict['band_gap']['direct'])
            
        cbm = bsp.bsdict['cbm']
        vbm = bsp.bsdict['vbm']
        cbm_kindex=cbm['kpoint_index'] ; vbm_kindex = vbm['kpoint_index']
        cbm_bindex =cbm['band_index'] ; vbm_bindex = vbm['band_index']

        self.dataframe["CBM_Kindex"].append(cbm_kindex)
        self.dataframe["CBM_Bindex"].append(cbm_bindex)
        self.dataframe["VBM_Kindex"].append(vbm_kindex)
        self.dataframe["VBM_Bindex"].append(vbm_bindex)
        return self.dataframe

    def EMInfo(self) :
        em_e = open("{}/EM".format(self.path),'r').readlines()[-12:]
        E = GMDExcitonbinding.harm_mean_em(em_e)
        if self.path.split(os.path.sep)[-2] == "E" :
            self.dataframe["me*"].append(E)
        else :
            self.dataframe["mh*"].append(E)
        return self.dataframe

    def DielectricInfo(self,e_path, h_path) :
        # effective mass
        em_e = open("{}/EM".format(e_path),'r').readlines()[-12:]
        em_h = open("{}/EM".format(h_path),'r').readlines()[-12:]

        E, H = GMDExcitonbinding.harm_mean_em(em_e), GMDExcitonbinding.harm_mean_em(em_h)

        # Outcar
        outcar = open("{}/OUTCAR".format(self.path),"r").readlines()
        dielec = GMDExcitonbinding.dielectricconst(outcar)
        self.dataframe["Dielectric Const"].append(dielec)
                        
        # reduced effecitve mass
        mu = (E*H)/(E+H)
        exciton = (mu*13.605692)/(dielec**2)
        self.dataframe["Exciton Binding"].append(exciton*1E+3)
        return self.dataframe