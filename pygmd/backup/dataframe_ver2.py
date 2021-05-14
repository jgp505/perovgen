import os
import sys

import collections 
import numpy as np
import pandas as pd 

from pymatgen.core import Structure, Element
from pymatgen.io.vasp.outputs import Vasprun, BSVasprun

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

class DataFrame :
    def _make_name(list):
        vaspname = ''
        sn = collections.Counter(list)
        for i in sn :
            try :
                if sn[i] == 1 :
                    vaspname += i.symbol
                else :
                    vaspname += "{}{}".format(i.symbol,sn[i])
            except :
                if sn[i] == 1:
                    vaspname += i
                else :
                    vaspname += "{}{}".format(i,sn[i])
        return vaspname

    def make_csv(path,bandinfo=False):
        data = collections.defaultdict(list)
        for p in path :
            try :
                vrun = Vasprun("{}/vasprun.xml".format(p))
            except :
                continue

            s = vrun.final_structure
            data["PATH"].append(p)
            data["Formula"].append(DataFrame._make_name(s.species))
            data["a"].append(s.lattice.a)
            data["b"].append(s.lattice.b)
            data["c"].append(s.lattice.c)
            data["alpha"].append(s.lattice.alpha)
            data["beta"].append(s.lattice.beta)
            data["gamma"].append(s.lattice.gamma)
            data["volume"].append(s.volume)
            space,number = s.get_space_group_info()
            data["Space Group"].append("{}[{}]".format(space, number))
            data["Energy"].append(vrun.final_energy)

            '''
            dm = s.distance_matrix

            for i in s.species :
                if i.group == 1 :
                    asite= i
                elif i.group == 15 :
                    msite = i
                elif i.group == 17 :
                    xsite = i

            index = average_bond(structure=s,atom1=msite,atom2=xsite)
            dm = s.distance_matrix
            if len(index) == 2 :
                data["M1"].append(index[0])
                data["M2"].append(index[1])
                data["M3"].append(0)
                data["M4"].append(0)
            else :
                data["M1"].append(index[0])
                data["M2"].append(index[1])
                data["M3"].append(index[2])
                data["M4"].append(index[3])
            '''

            if bandinfo:
                run = BSVasprun("{}/vasprun.xml".format(p),parse_projected_eigen=True)
                bs = run.get_band_structure("KPOINTS")
                data["Ef"].append(vrun.efermi)
                bandgap = "%.3f(Indirect)"%(run.as_dict()['output']['bandgap'])
                if run.as_dict()['output']['is_gap_direct']:
                    bandgap = "%.3f(Direct)"%(run.as_dict()['output']['bandgap'])
                data["Eg"].append(bandgap)

        if bandinfo :
            df = pd.DataFrame(data,columns=["PATH","Formula","a","b","c","alpha","beta","gamma",
                "volume","Space Group","Energy",
                "Ef","Eg"])
        else :
            df = pd.DataFrame(data,columns=["PATH","Formula","a","b","c","alpha","beta","gamma",
            "volume","Space Group","Energy"])
        return df 