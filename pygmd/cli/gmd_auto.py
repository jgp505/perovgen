import os
import sys
import time
import subprocess

from shutil import copyfile
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun

from perovgen.pygmd.input_structure import load_structure, GMDStructure
from perovgen.pygmd.autocal.inputset import *
from perovgen.pygmd.autocal.algorithm import openingphrase
from perovgen.pygmd.analysis.electronic import BSPlotting, DOSPlotting, GMDAnalysis
from perovgen.pygmd.analysis.energy import GMDExcitonbinding

def analyze_calculation(args):
    pwd = os.getcwd()
    if args.inputpath :
        inputs = inputgmd(args.inputpath)
    else :
        print("\nUserError : gmd3 auto [path of input.gmd] [path of structure] [optional]\n")
        sys.exit(1) 
    
    if not args.strucpath :
        print("\nUserError : gmd3 auto [path of input.gmd] [path of structure] [optional]\n")
        sys.exit(1) 
    else :
        if not load_structure(args.strucpath) :
            print("\n[Warning] the structure file does not exist.\n")
            sys.exit(1) 
        else :
            strucpath = load_structure(args.strucpath)

    openingphrase(inputs,strucpath[-1])
    if args.directory :
        for struc, filename in zip(strucpath[0], strucpath[-1]) :
            fn = os.path.basename(filename).split(".")[0]
            for mt in inputs.calmode:
                runfolder = []
                p = PerovInputs(structure=struc, is_selective=args.deleteselect)
                try :
                    symmetry, groupnumber = struc.get_space_group_info()
                except :
                    groupnumber = 0

                # Revised INPUT FILES
                vi = p.inputfolder(inputs=inputs, method=mt,soc=args.soc)
                inputs = inputgmd(args.inputpath)

                # Writing the Input directory
                full_formula = GMDStructure(struc).formula(reduced=False)
                pretty_formula = GMDStructure(struc).formula()
                root_dir = "{0}/{1}_{2:03d}".format(pwd, pretty_formula, groupnumber)

                folder_name = "{0}/{1}_{2}_{3}".format(root_dir,mt,full_formula,fn) 
                vi.write_input(output_dir=folder_name) ; runfolder.append(folder_name)
                if mt == "E" :
                    folder_name_H = "{0}/H_{1}_{2}".format(root_dir,full_formula,fn) 
                    vi.write_input(output_dir=folder_name_H); runfolder.append(folder_name_H)

                for runf in runfolder :
                    naming = os.path.basename(runf)
                    rs = RunningShell(shell=inputs.shell, name=naming, path=runf)
                    rs.running_mode(soc=args.soc, run=False)

    else :
        for struc, filename in zip(strucpath[0], strucpath[-1]) :
            fn = os.path.basename(filename).split(".")[0]

            # directory
            f = open("{}/calnohup.py".format(os.path.dirname(__file__)),"r")
            ff = f.readlines()
            ff.insert(0,"inputpath='{}'\n".format(os.path.abspath(args.inputpath)))
            ff.insert(1,"strucpath='{}'\n".format(filename))
            if args.deleteselect :
                ff.insert(2,"ds=True\n")
            else :
                ff.insert(2,"ds=False\n")
            if args.soc :
                ff.insert(3,"orbit=True\n")
            else :
                ff.insert(3,"orbit=False\n")
            ff.insert(4,"mole=False\n")
            ff.insert(5,"\n")

            # directory name
            try :
                symmetry, groupnumber = struc.get_space_group_info()
            except :
                groupnumber = 0
            pretty_formula = GMDStructure(struc).formula()
            root_dir = "{0}/{1}_{2:03d}".format(pwd, pretty_formula, groupnumber)

            os.makedirs(root_dir,exist_ok=True)
            with open("{}/run_{}.py".format(root_dir,fn),"w") as fi :
                for i in ff :
                    fi.write(i)
            fi.close()
            #os.system("python {}/run_{}.py".format(root_dir,fn))
            os.system("nohup python -u {}/run_{}.py > {}/gmd_{}.out &".format(root_dir,fn,root_dir,fn))