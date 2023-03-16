import os
import sys
from tkinter import W

from pymatgen.core import Structure
import numpy as np
import pandas as pd

from shutil import copy, move
from collections import defaultdict

from perovgen.pygmd.input_structure import load_structure
from perovgen.pygmd.analysis.electronic import BSPlotting, GMDAnalysis
from perovgen.pygmd.analysis.energy import GMDExcitonbinding
from perovgen.pygmd.analysis.energy import GMDAbsorption
from perovgen.pygmd.cli.gmd_write import pdospath

def analyze_data(args):
    if args.dos or args.pdos:
        path = pdospath()
        if args.dos:
            os.system("%s dos width=%f"%(path,args.dos[0]))
        elif args.pdos :
            s = load_structure(os.getcwd())[0][0]
            if not s :
                print("Structure file doesn't exist!\n")
                sys.exit(1)
            dic = defaultdict(list)
            for e,i in enumerate(s.composition.elements) :
                dic[i.symbol].append(e+1)
                for k,v in dic.items():
                    for s in ['tot','s','p'] :
                        with open("LIST",'w') as fi :
                            fi.write("{}_{}\n".format(k,s))
                            for i in v :
                                if s == "tot" :
                                    fi.write("%i tot tot\n"%(i))
                                elif s == "s" :
                                    fi.write("%i s tot\n"%(i))
                                elif s == "p" :
                                    fi.write("%i px tot\n"%(i))
                                    fi.write("%i py tot\n"%(i))
                                    fi.write("%i pz tot\n"%(i))
                        os.system("%s pdos width=%f"%(path,args.pdos[0]))
    elif args.band :
        bsp = BSPlotting() 
        bsp.write_to_json(os.getcwd())
        
    elif args.absorp:
        if len(args.absorp[0]) == 1 :
            df = GMDAbsorption.make_csv(args.absorp[0])
            df.to_csv("%s/absorp.csv"%(os.getcwd()))
            print("The file generated in {}/absorp.csv".format(os.getcwd()))
        else :
            print("Please enter the path")
        
    elif args.em :
        if len(args.em) == 1 :
            pwd = os.getcwd() ; abspath = os.path.abspath(args.em[0])
            emc = GMDAnalysis()
            if "electron" in os.listdir(abspath) and "hole" in os.listdir(abspath):
                os.chdir("{}/electron".format(abspath))
                if not "EM" in os.listdir(".") :
                    emc.effectivemass(".", secondstep=True)
                electron = "{}/EM".format(os.getcwd())

                os.chdir("{}/hole".format(abspath))
                if not "EM" in os.listdir(".") :
                    print(os.getcwd())
                    emc.effectivemass(".", secondstep=True)
                hole = "{}/EM".format(os.getcwd())

            elif "electron" in os.listdir(abspath) :
                os.chdir("{}/electron".format(abspath))
                if not "EM" in os.listdir(".") :
                    emc.effectivemass(".", secondstep=True)
                electron = "{}/EM".format(os.getcwd())
                hole = None

            elif "hole" in os.listdir(abspath) :
                os.chdir("{}/hole".format(abspath))
                if not "EM" in os.listdir(".") :
                    emc.effectivemass(".", secondstep=True)
                hole = "{}/EM".format(os.getcwd())
                electron = None
            else :
                print("Please calculation effective mass")
                
            os.chdir(pwd)
            with open("effectivemass.gmd",'w') as fi :
                if electron == None and hole == None :
                    fi.write("None")
                elif electron == None :
                    hole = GMDExcitonbinding.harm_mean_em(hole)
                    fi.write("VBM\n")
                    fi.write("x : %.3f\n"%(hole["x"]))
                    fi.write("y : %.3f\n"%(hole["y"]))
                    fi.write("z : %.3f\n"%(hole["z"]))
                    fi.write("harm_mean_hole : %.3f\n"%(hole["HarmMean"]))
                elif hole == None :
                    electron = GMDExcitonbinding.harm_mean_em(electron)
                    fi.write("CBM\n")
                    fi.write("x : %.3f\n"%(electron["x"]))
                    fi.write("y : %.3f\n"%(electron["y"]))
                    fi.write("z : %.3f\n"%(electron["z"]))
                    fi.write("harm_mean_electron : %.3f\n"%(electron["HarmMean"]))
                else :
                    electron = GMDExcitonbinding.harm_mean_em(electron)
                    hole = GMDExcitonbinding.harm_mean_em(hole)
                    harm_mean = 2/((1/electron['HarmMean']) + (1/hole['HarmMean']))
                    
                    fi.write("CBM\n")
                    fi.write("x : %.3f\n"%(electron["x"]))
                    fi.write("y : %.3f\n"%(electron["y"]))
                    fi.write("z : %.3f\n"%(electron["z"]))
                    fi.write("harm_mean_electron : %.3f\n"%(electron["HarmMean"]))
                    fi.write("\n\n")
                    fi.write("VBM\n")
                    fi.write("x : %.3f\n"%(hole["x"]))
                    fi.write("y : %.3f\n"%(hole["y"]))
                    fi.write("z : %.3f\n"%(hole["z"]))
                    fi.write("harm_mean_hole : %.3f\n"%(hole["HarmMean"]))
                    fi.write("mu value : %.5f"%(harm_mean))
            fi.close()
            print("The file generated in {}/effectivemass.gmd".format(pwd))
        else :
            print("Please enter the path for effective mass")
            
    elif args.de :
        if len(args.de) == 1 :
            dielectric = GMDExcitonbinding.dielectricconst("{}/OUTCAR".format(os.path.abspath(args.de[0])))

            array1 = dielectric['Array']['high_frequency_dielectric']
            array2 = dielectric['Array']['ionic_contribution']
            high_mean = dielectric['geom_mean']
            
            with open("dielectric.gmd",'w') as fi :
                fi.write("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT\n")
                fi.write("====================================\n")
                for i in array1 :
                    for j in i :
                        fi.write("%.5f   "%(j))
                    fi.write("\n")
                fi.write("====================================\n\n")

                fi.write("MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION\n")
                fi.write("====================================\n\n")
                for i in array2 :
                    for j in i :
                        fi.write("%.5f   "%(j))
                    fi.write("\n")
                fi.write("====================================\n\n")

                fi.write("static_dielectric\n")
                fi.write("====================================\n")
                for i in range(3) :
                    for j in range(3) :
                        fi.write("%.5f   "%(array1[i][j]+array2[i][j]))
                    fi.write("\n")
                fi.write("====================================\n\n")

                fi.write("high_frequency_dielectric\n")
                fi.write("====================================\n")
                for i in array1 :
                    for j in i :
                        fi.write("%.5f   "%(j))
                    fi.write("\n")
                fi.write("====================================\n\n")

                fi.write("Geom. Mean : %.5f"%(high_mean))
            fi.close()
            print("The file generated in {}/dielectric.gmd".format(os.path.abspath(args.de[0])))
            
    elif args.exciton :
        if len(args.exciton) == 2 :
            # Harm Mean effective mass
            pwd = os.getcwd() ; harm_abspath = os.path.abspath(args.exciton[0])
            geom_abspath = os.path.abspath(args.exciton[1])
            emc = GMDAnalysis()
            if "electron" in os.listdir(harm_abspath) and "hole" in os.listdir(harm_abspath):
                os.chdir("{}/electron".format(harm_abspath))
                if not "EM" in os.listdir(".") :
                    emc.effectivemass(".", secondstep=True)
                electron = GMDExcitonbinding.harm_mean_em("{}/EM".format(os.getcwd()))

                os.chdir("{}/hole".format(harm_abspath))
                if not "EM" in os.listdir(".") :
                    emc.effectivemass(".", secondstep=True)
                hole = GMDExcitonbinding.harm_mean_em("{}/EM".format(os.getcwd()))
            else :
                print("Please calculation effective mass")
                sys.exit(1)
            os.chdir(pwd)
            harm_mean = 2/((1/electron['HarmMean']) + (1/hole['HarmMean']))
            
            # Geom Mean Dielectric Constant 
            dielectric = GMDExcitonbinding.dielectricconst("{}/OUTCAR".format(geom_abspath))
            geom_mean = dielectric['geom_mean']

            exciton = GMDExcitonbinding.excitonbinding(harm_mean,geom_mean)

            with open("excitonbinding.gmd",'w') as fi :
                fi.write("%.5f meV"%(exciton*100))
            fi.close()
            print("The file generated in {}/excitonbinding.gmd".format(pwd))
                
        else :
            print("gmd3 data --exciton [effective_mass_path] [dielectric_path]")