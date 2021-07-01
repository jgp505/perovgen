import os
import sys

import pandas as pd
import numpy as np

from collections import defaultdict
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core import Element, Structure
from perovgen.pygmd.base import GMDStructure
from perovgen.pygmd.analysis.electronic import BSPlotting
from perovgen.pygmd.analysis.energy import GMDPhaseDiagram, GMDExcitonbinding
data = defaultdict(list)
for f1 in os.listdir(os.getcwd()):
    if os.path.isdir(f1) :
        files=[]
        for f2 in os.listdir(f1) :
            if os.path.isdir("{}/{}".format(f1,f2)) :
                for f3 in os.listdir("{}/{}".format(f1,f2)) :
                    data[os.path.abspath("{}/{}".format(f1,f2))].append(f3)
            else :
                pass

result=[]
index1 = 1
for k,v in data.items():
    print(k.split(os.path.sep)[-1],end="\n\n")
    dataframe=defaultdict(list)
    dataframe['path'].append(k)
    dataframe['Index'].append(index1)
    index1+=1
    for v1 in v :
        if v1 == "R_mode":
            s = Structure.from_file("{}/{}/POSCAR".format(k,v1))
            s.make_supercell([[2,0,0],[0,2,0],[0,0,2]])
            index = [*s.composition.get_el_amt_dict()]
            species = [ele.symbol for ele in s.species]
            
            distance = s.distance_matrix
            bsite , xsite =[], []
            for e,u in enumerate(species) :
                if u == index[1] :
                    bsite.append(e)
                elif u == index[2]:
                    xsite.append(e)
            if len(distance[bsite[0],xsite][distance[bsite[0],xsite]<=3.6]) == 6 :
                dataframe['Perovskite'].append(True)
            else :
                dataframe['Perovskite'].append(False)

        elif v1 == "C_mode" :
            vrun = Vasprun("{}/{}/vasprun.xml".format(k,v1),parse_potcar_file=False)
            s = vrun.final_structure
            dataframe["formula"].append(GMDStructure(s).vaspname())
            dataframe["reduced_formula"].append(s.composition.reduced_formula)


            dataframe['vpa'].append(s.volume/len(s))
            total_rad = 0
            for site in s.species :
                total_rad += site.atomic_radius_calculated ** 3
            dataframe['packing fraction'].append(4*np.pi*total_rad/(3*s.volume))

            dataframe["volume"].append(s.volume)
            dataframe["density"].append(s.density)
            space,number = s.get_space_group_info()
            dataframe["Space Group"].append(number)
            dataframe["Energy"].append("%.5f"%(vrun.final_energy))
            
            entries,pd1 = GMDPhaseDiagram(["{}/{}/vasprun.xml".format(k,v1)]).get_phasediagram()
            h = pd1.get_e_above_hull(entries[0])
            d = pd1.get_form_energy_per_atom(entries[0])
            dataframe["E_hull"].append(h)
            dataframe["E_form"].append(d)
            dataframe["total_electrons"].append(s.composition.total_electrons)
            dataframe["electronegativity"].append(s.composition.average_electroneg)

        elif v1 == "B_mode" :
            bsp = BSPlotting(vasprun="{}/{}/vasprun.xml".format(k,v1), 
            kpoints = "{}/{}/KPOINTS".format(k,v1))
            bi = bsp._bandinform()
            #dataframe["Ef"].append(bsp.bs.efermi)
            #dataframe['Eg'].append(bsp.bsdict['band_gap']['energy'])
            dataframe["Eg_direct"].append(bi['E_g']['Direct'])
            if bsp.bsdict['band_gap']['direct'] :
                dataframe["Eg_indirect"].append(bi['E_g']['Direct'])
            else :
                dataframe["Eg_indirect"].append(bi['E_g']['Indirect'])
            dataframe["Eg"].append(bsp.bsdict['band_gap']['direct'])
            
            cbm = bsp.bsdict['cbm']
            vbm = bsp.bsdict['vbm']

            cbm_kindex=cbm['kpoint_index'] ; vbm_kindex = vbm['kpoint_index']
            cbm_bindex =cbm['band_index'] ; vbm_bindex = vbm['band_index']

            dataframe["CBM_Kindex"].append(cbm_kindex)
            dataframe["CBM_Bindex"].append(cbm_bindex)
            dataframe["VBM_Kindex"].append(vbm_kindex)
            dataframe["VBM_Bindex"].append(vbm_bindex)
            
            
    df = pd.DataFrame(dataframe,columns=["Index","path","formula","reduced_formula","total_electrons","electronegativity",
                        "volume","density","vpa","packing fraction","Space Group","Energy","E_hull","E_form",
                        "Eg","Eg_indirect","Eg_direct","CBM_Kindex","CBM_Bindex","VBM_Kindex","VBM_Bindex","Perovskite"],
                        dtype=object)
    result.append(df)
df1 = pd.concat(result)
df1 = df1.set_index("Index")
df1.to_csv("Result/Result_ver2.csv")
