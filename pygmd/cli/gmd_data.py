import os
import sys

import pandas as pd
from collections import defaultdict

from pymatgen.io.vasp.outputs import Vasprun
from perovgen.pygmd.base import GMDStructure
from perovgen.pygmd.analysis.electronic import BSPlotting
from perovgen.pygmd.analysis.energy import GMDPhaseDiagram, GMDExcitonbinding


def dataframe(args):
    result = []
    for p1 in args.path :
        path = [f for f in os.listdir(p1) if "mode" in f]
        dict1 = defaultdict(list)
        for p in path :
            name, number ,mode = p.split("_")[:-1]
            dict1["{}_{}".format(name,number)].append(mode)

        for k,v in dict1.items():
            data =defaultdict(list)
            if len(v) == 6 or len(v) == 7 or len(v) == 5 :
                for v1 in v :
                    if v1 == "C" :
                        vrun = Vasprun("{}/{}_{}_mode/vasprun.xml".format(os.path.abspath(p1),k,v1))
                        s = vrun.final_structure
                        data["Formula"].append(GMDStructure(s).vaspname())
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
                        entries,pd = GMDPhaseDiagram(["{}/{}_{}_mode/vasprun.xml".format(os.path.abspath(p1),k,v1)]).get_phasediagram()
                        h = pd.get_e_above_hull(entries[0])
                        d = pd.get_form_energy_per_atom(entries[0])
                        data["E_hull"].append(h)
                        data["E_form"].append(d)

                    elif v1 == "B" :
                        bsp = BSPlotting(vasprun="{}/{}_{}_mode/vasprun.xml".format(os.path.abspath(p1),k,v1), 
                            kpoints = "{}/{}_{}_mode/KPOINTS".format(os.path.abspath(p1),k,v1))
                        data["Ef"].append(bsp.bs.efermi)
                        data["Eg"].append(bsp.bsdict['band_gap']['energy'])
                        cbm = bsp.bsdict['cbm']
                        vbm = bsp.bsdict['vbm']
                        cbm_kindex=cbm['kpoint_index'] ; vbm_kindex = vbm['kpoint_index']
                        cbm_bindex =cbm['band_index'] ; vbm_bindex = vbm['band_index']

                        data["CBM_Kindex"].append(cbm_kindex)
                        data["CBM_Bindex"].append(cbm_bindex)
                        data["VBM_Kindex"].append(vbm_kindex)
                        data["VBM_Bindex"].append(vbm_bindex)
                        print("{}_{}_mode".format(k,v1))
                        bsp.printinform()

                    elif v1 == "E" or v1 == "H" :
                        em_e = open("{}/{}_{}_mode/EM".format(os.path.abspath(p1),k,v1),'r').readlines()[-12:]
                        E = GMDExcitonbinding.harm_mean_em(em_e)
                        if v1 == "E" :
                            data["me*"].append(E)
                            print("\n{}_{}_mode".format(k,v1))
                            print("\nConduction Effective Mass : ", E)
                        else :
                            data["mh*"].append(E)
                            print("\n{}_{}_mode".format(k,v1))
                            print("\nValence Effective Mass : ", E)

                    elif v1 == "U" :
                        print("\n{}_{}_mode".format(k,v1))
                        # effective mass
                        em_e = open("{}/{}_E_mode/EM".format(os.path.abspath(p1),k),'r').readlines()[-12:]
                        em_h = open("{}/{}_H_mode/EM".format(os.path.abspath(p1),k),'r').readlines()[-12:]
                        E, H = GMDExcitonbinding.harm_mean_em(em_e), GMDExcitonbinding.harm_mean_em(em_h)

                        # Outcar
                        outcar = open("{}/{}_U_mode/OUTCAR".format(os.path.abspath(p1),k)).readlines()
                        dielec = GMDExcitonbinding.dielectricconst(outcar)
                        data["Dielectric Const"].append(dielec)
                        print("Dielectric Constant : ", dielec)

                        # reduced effecitve mass
                        mu = (E*H)/(E+H)
                        exciton = (mu*13.605692)/(dielec**2)
                        print("Exciton binding Energy : ",exciton*1E+3,"(meV)")
                        data["Exciton Binding"].append(exciton*1E+3)
                    elif v1 == "R" :
                        pass
                    print()
                df = pd.DataFrame(data,columns=["Formula","a","b","c","alpha","beta","gamma","volume","Space Group","Energy",
                                               "Ef","Eg","CBM_Kindex","CBM_Bindex","VBM_Kindex","VBM_Bindex",
                                                "me*","mh*","Dielectric Const","Exciton Binding"])
                result.append(df)
    df1 = pd.concat(result)
    df1.to_csv("Result.csv")
            



