import os
import sys
import time
import subprocess

from shutil import copyfile

from sympy import N
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun

from perovgen.pygmd.shell import ShellPath
from perovgen.pygmd.input_structure import load_structure, GMDStructure
from perovgen.pygmd.autocal.inputset import *
from perovgen.pygmd.autocal.algorithm import openingphrase
from shutil import move, copy

def analyze_calculation(args):
    pwd = os.getcwd()
    if args.inputpath :
        inputs = inputgmd(args.inputpath)
    else :
        shell, numberlist = ShellPath().check()
        while True :
            number = int(input("Please enter the number >> "))
            if number in numberlist :
                break
        modecheck = ['bulkmodulus','relaxation','chgcar','dos','band',
                    'effective_mass','dielectric','HSE','absorption']

        with open("input.gmd",'w') as fi :
            fi.write("### Perovgen input.gmd\n")
            fi.write("KPOINTS:\nA\nCONSTK=30\n\n")
            fi.write("INCAR:\nPBE\n\n")
            fi.write("SHELL:\n%s\n\n"%(shell[number-1][-1]))
            fi.write("CALMODE:\n")
            for m in modecheck :
                fi.write("%s = T\n"%(m))
        fi.close()
        print("Generate input.gmd")
        print("\nUserError : gmd3 auto [path of input.gmd] [path of structure] [optional]\n")
        print("[path of input.gmd] ex. -i input.gmd")
        print("[path of structure] ex. -p *.cif POSCAR*")
        print("[optional] --soc (spin orbit coupling)\n\t -d (direcotry, no calculation)\n\t --de (remove the selective dynamics")
        sys.exit(1) 
    
    if not args.strucpath :
        print("\nUserError : gmd3 auto [path of input.gmd] [path of structure] [optional]\n")
        print("optional : --soc (spin orbit coupling)\n\t -d (direcotry, no calculation)\n\t --de (remove the selective dynamics")
        sys.exit(1) 
    else :
        if not load_structure(args.strucpath) :
            print("\n[Warning] the structure file does not exist.\n")
            sys.exit(1) 
        else :
            strucpath = load_structure(args.strucpath)
            for struc,filename in zip(strucpath[0],strucpath[1]) :
                for e in struc.composition.elements :
                    if e.Z >= 56 and "HSE" in [*inputs.calmode]:
                        print(filename, "will be calculated considering SOC\n\n")
                        break
                
    openingphrase(inputs,strucpath[-1])
    for struc, filename in zip(strucpath[0], strucpath[-1]) :
        fn = os.path.basename(filename).split(".")[0]
        formula = PerovInputs(struc).naming
        root_dir = "{0}/{1}_{2}".format(pwd,formula,fn)
        os.makedirs(root_dir,exist_ok=True)
        if args.directory :
            for k,v in inputs.calmode.items():
                poscar = PerovInputs(struc)
                inputs = inputgmd(args.inputpath)
                naming = poscar.naming+"_"+fn
                runfolder = []
                if k == 'bulkmodulus' :
                    volumes = v*poscar.poscar.volume
                    for v1, vol in zip(v, volumes) :
                        struc.scale_lattice(vol)
                        PI = PerovInputs(struc)
                        vi = PI.inputfolder(inputs=inputs, method=k, soc=args.soc)
                        folder_name = "%s/%s/%.2f"%(root_dir,k,v1)
                        vi.write_input(folder_name)
                        runfolder.append(folder_name)
                
                elif k == 'dos' or k == 'band' or k == 'effective_mass' :
                    vi = poscar.inputfolder(inputs=inputs, method=k, soc=args.soc)
                    if k == 'effective_mass' :
                        folder_name_elec = "%s/chgcar/%s/electron"%(root_dir,k)
                        folder_name_hole = "%s/chgcar/%s/hole"%(root_dir,k)
                
                        # CBM effective mass
                        vi.write_input(folder_name_elec)
                        runfolder.append(folder_name_elec)
                
                        # VBM effective mass
                        vi.write_input(folder_name_hole)
                        runfolder.append(folder_name_hole)
                    else :
                        folder_name = "%s/chgcar/%s"%(root_dir,k)
                        vi.write_input(folder_name)
                        runfolder.append(folder_name)
                        
                elif k == 'HSE' :
                    hse_soc = False
                    for e in struc.composition.elements :
                        if e.Z >= 56 :
                            hse_soc = True
                            break
                    if hse_soc == True :
                        vi1 = poscar.inputfolder(inputs=inputs, method=k, soc=False)
                        vi2 = poscar.inputfolder(inputs=inputs, method=k, soc=True)
                        folder_name_nsoc = "%s/HSE"%(root_dir)
                        folder_name_soc = "%s/HSE_SOC"%(root_dir)
                        vi1.write_input(folder_name_nsoc)
                        vi2.write_input(folder_name_soc)
                        runfolder.append(folder_name_nsoc)
                        runfolder.append(folder_name_soc)
                    else :
                        folder_name = "%s/%s"%(root_dir,k)
                        vi = poscar.inputfolder(inputs=inputs, method=k, soc=args.soc)
                        vi.write_input(folder_name)
                        runfolder.append(folder_name)
                    
                else :
                    folder_name = "%s/%s"%(root_dir,k)
                    vi = poscar.inputfolder(inputs=inputs, method=k, soc=args.soc)
                    vi.write_input(folder_name)
                    runfolder.append(folder_name)
            
                e = 0
                for runf in runfolder :
                    if k == 'bulkmodulus' :
                        rs = RunningShell(shell=inputs.shell, name="%s_%s_%.2f"%(naming,k,v[e]),path=runf)
                        e += 1
                        rs.running_mode(soc=args.soc, run=False) # If run is true, running the calculation
                    elif k == 'HSE' :
                        runf1 = os.path.basename(runf)
                        if "HSE_SOC" == runf1 :
                            rs = RunningShell(shell=inputs.shell, name="%s_HSE_SOC"%(naming),path=runf)
                            rs.running_mode(soc=True, run=False) # If run is true, running the calculation
                        elif "HSE" == runf1 :
                            rs = RunningShell(shell=inputs.shell, name="%s_HSE_SOC"%(naming),path=runf)
                            rs.running_mode(soc=False, run=False) # If run is true, running the calculation
                            
                    else :
                        rs = RunningShell(shell=inputs.shell, name="%s_%s"%(naming,k), path = runf)
                        rs.running_mode(soc=args.soc, run=False) # If run is true, running the calculation
                
        else :
            # directory
            f = open("{}/calnohup.py".format(os.path.dirname(__file__)),"r")
            ff = f.readlines()
            # copy the input.gmd file 
            copy(os.path.abspath(args.inputpath), "{}/input.gmd".format(root_dir))
            #ff.insert(0,"inputpath='{}'\n".format(os.path.abspath(args.inputpath)))
            ff.insert(0,"inputpath='{}'\n".format("{}/input.gmd".format(root_dir)))
            ff.insert(1,"strucpath='{}'\n".format(filename))
            if args.deleteselect :
                ff.insert(2,"ds=True\n")
            else :
                ff.insert(2,"ds=False\n")
            if args.soc :
                ff.insert(3,"orbit=True\n")
            else :
                ff.insert(3,"orbit=False\n")
            #ff.insert(4,"mole=False\n")
            ff.insert(4,"\n")
            with open("{}/run_{}.py".format(root_dir,fn),"w") as fi :
                for i in ff :
                    fi.write(i)
            fi.close()

            #os.system("python {}/run_{}.py".format(root_dir,fn))
            os.system("nohup python -u {}/run_{}.py > {}/gmd_{}.out &".format(root_dir,fn,root_dir,fn))
