import os
import sys

import pandas as pd
from collections import defaultdict
from pymatgen.io.vasp.outputs import Vasprun
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
checklist = open("checklist.txt",'w')
for k,v in data.items():
    print(k.split(os.path.sep)[-1],end="\n\n")
    dataframe=defaultdict(list)
    for v1 in v :
        if v1 == "C_mode" :
            vrun = Vasprun("{}/{}/vasprun.xml".format(k,v1),parse_potcar_file=True)
            s = vrun.final_structure
            dataframe["Formula"].append(GMDStructure(s).vaspname())
            dataframe["a"].append(s.lattice.a)
            dataframe["b"].append(s.lattice.b)
            dataframe["c"].append(s.lattice.c)
            dataframe["alpha"].append(s.lattice.alpha)
            dataframe["beta"].append(s.lattice.beta)
            dataframe["gamma"].append(s.lattice.gamma)
            dataframe["volume"].append(s.volume)
            space,number = s.get_space_group_info()
            dataframe["Space Group"].append("{}[{}]".format(space, number))
            dataframe["Energy"].append("%.5f"%(vrun.final_energy))
            entries,pd1 = GMDPhaseDiagram(["{}/{}/vasprun.xml".format(k,v1)]).get_phasediagram()
            h = pd1.get_e_above_hull(entries[0])
            d = pd1.get_form_energy_per_atom(entries[0])
            dataframe["E_hull"].append(h)
            dataframe["E_form"].append(d)

            element = list(s.composition.get_el_amt_dict().items())
            oxistate = s.composition.oxi_state_guesses()[0]
            dataframe["Asite"].append(element[0][0])
            dataframe["A_num"].append(element[0][1])
            dataframe["Bsite"].append(element[1][0])
            dataframe["B_num"].append(element[1][1])
            dataframe["Xsite"].append(element[2][0])
            dataframe["X_num"].append(element[2][1])
            dataframe["Asite_oxistate"].append(oxistate[element[0][0]])
            dataframe["Bsite_oxistate"].append(oxistate[element[1][0]])
            dataframe["Xsite_oxistate"].append(oxistate[element[2][0]])
            dataframe["W"].append("%.5f"%(s.composition.weight))
            dataframe["total_electrons"].append(s.composition.total_electrons)
            dataframe["electronegativity"].append(s.composition.average_electroneg)

        elif v1 == "B_mode" :
            bsp = BSPlotting(vasprun="{}/{}/vasprun.xml".format(k,v1), 
            kpoints = "{}/{}/KPOINTS".format(k,v1))
            bi = bsp._bandinform()
            dataframe["Ef"].append(bsp.bs.efermi)
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

        elif v1 == "E_mode" or v1 == "H_mode" :
            em_e = open("{}/{}/EM".format(k,v1),'r').readlines()[-12:]
            E = GMDExcitonbinding.harm_mean_em(em_e)
            if v1 == "E_mode" :
                dataframe["me*"].append(E)
            else :
                dataframe["mh*"].append(E)

        elif v1 == "U_mode" :
            # effective mass
            em_e = open("{}/E_mode/EM".format(k),'r').readlines()[-12:]
            em_h = open("{}/H_mode/EM".format(k),'r').readlines()[-12:]
            E, H = GMDExcitonbinding.harm_mean_em(em_e), GMDExcitonbinding.harm_mean_em(em_h)

            # Outcar
            outcar = open("{}/{}/OUTCAR".format(k,v1)).readlines()
            dielec = GMDExcitonbinding.dielectricconst(outcar)
            dataframe["Dielectric Const"].append(dielec)
                        
            # reduced effecitve mass
            mu = (E*H)/(E+H)
            exciton = (mu*13.605692)/(dielec**2)
            dataframe["Exciton Binding"].append(exciton*1E+3)
        
    df = pd.DataFrame(dataframe,columns=["Formula","Asite","Bsite","Xsite","A_num","B_num","X_num","Asite_oxistate","Bsite_oxistate",
                                        "Xsite_oxistate","W","total_electrons","electronegativity",
                                        "a","b","c","alpha","beta","gamma","volume","Space Group","Energy","E_hull","E_form",
                                        "Ef","Eg","Eg_indirect","Eg_direct","CBM_Kindex","CBM_Bindex","VBM_Kindex","VBM_Bindex","me*","mh*",
                                        "Dielectric Const","Exciton Binding"],index=[0],dtype=object)
    checklist.write("%s\n"%(k.split(os.path.sep)[-1]))
    result.append(df)
df1 = pd.concat(result)
df1.to_csv("Result.csv")
