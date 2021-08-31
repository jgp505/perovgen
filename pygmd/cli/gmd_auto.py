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

pwd = os.getcwd()

def analyze_calculation(args):
    inputs = read_input(args.input)
    
    # initial structural path
    if args.path :
        strucpath = args.path
    else :
        strucpath = inputs["STRUC"]

    if not strucpath :
        print("Please check the structural files!!")
        sys.exit(1)
    else :
        try :
            openingphrase(inputs,strucpath)
        except FileNotFoundError :
            print("the structure file does not exist.\n")
            sys.exit(1)

    # start calculation
    if args.directory :
        for mt in inputs["METHOD"] :
            os.chdir(pwd)
            runfolder = []

            # Make Folder
            struclist, filenames = load_structure(strucpath)
            if not struclist : 
                print("\n[Warning] the structure file does not exist.\n")
                sys.exit(0)

            # copy where the CHGCAR is located
            if mt == "B" or mt == "D" or mt=="E" :
                path = [os.path.abspath(b) for b in strucpath]
                kpath_list = [CopyCHGCAR(i) for i in path]

                # Data for making INPCAR using emc-master modules
                if mt == "E" :
                    nelect=[];kpoints=[]
                    for b in path :
                        if os.path.isfile(b) :
                            b = os.path.split(b)[0]
                        a = subprocess.check_output(['grep','NELECT','%s/OUTCAR'%(b)])
                        bsp = BSPlotting(vasprun=os.path.abspath("{}/vasprun.xml".format(b)), kpoints=os.path.abspath("{}/KPOINTS".format(b)))

                        try : 
                            vbm = bsp.bs.kpoints[bsp.bsdict['vbm']['kpoint_index'][0]].as_dict()['fcoords']
                        except : 
                            vbm = [0.000, 0.000, 0.000]
                        try :
                            cbm = bsp.bs.kpoints[bsp.bsdict['cbm']['kpoint_index'][0]].as_dict()['fcoords']
                        except : 
                            cbm = [0.000, 0.000, 0.000]

                        if args.soc : 
                            nelect.append((int(float(a.split()[2])),int(float(a.split()[2]))+1))
                        else :
                            nelect.append((int(float(a.split()[2])/2),int(float(a.split()[2])/2)+1))
                        kpoints.append((vbm, cbm))
                else :
                    pass

            for e,s in enumerate(struclist) :

                if args.deleteselect :
                    p = PerovInputs(structure=s,is_selective=True)
                else :
                    p = PerovInputs(structure=s)

                # Control incar 
                incar = GMDIncar(inputs).incar
                modeincar = PerovInputs._incarmode(incar=incar, method=mt)
                if args.soc :
                    modeincar["LSORBIT"]=True

                # Writing the folder 
                if not "None" in inputs["INCAR"] :
                    vi = p.inputfolder(incar=modeincar, number=inputs["KPOINTS"][0])
                else :
                    vi = p.inputfolder(incar=modeincar, number=None)

                symmetry, groupnumber = s.get_space_group_info()
                symmetry = symmetry.split("/")[0]
                fn = filenames[e].split("/")[-1].split(".")[0]

                root_dir = "{0}/{1}_{2:03d}/".format(pwd, pretty_formula, groupnumber)
        
                fn = 1
                if os.path.exists(root_dir) :
                    for exist in os.listdir(root_dir) :
                        if mt in exist and full_formula in exist :
                            fn += 1 
            
                folder_name = "{0}/{1}_{2}_{3:02d}".format(root_dir,mt,full_formula,fn) 
                vi.write_input(output_dir=folder_name)
                runfolder.append(folder_name)

                if mt == "E" :
                    folder_name_H = "{0}/H_{1}_{2:02d}".format(root_dir,full_formula,fn) 
                    vi.write_input(output_dir=folder_name_H)
                    runfolder.append(folder_name_H)

                # Copy the other files to generated folder
                if mt == "D" or mt=="B" or mt == "E" :
                    copyfile(kpath_list[e],"{}/{}/CHGCAR".format(pwd,folder_name))
                    if mt == "B" :
                        MakingKpointBand(s,"{}/{}/KPOINTS".format(pwd,folder_name))
                    elif mt == "E" :
                        copyfile(kpath_list[e],"{}/{}/CHGCAR".format(pwd,folder_name))
                        MakingInpcar(s,"{}/{}/INPCAR".format(pwd,folder_name),nelect[e][1],kpoints[e][1])
                        copyfile(kpath_list[e],"{}/{}/CHGCAR".format(pwd,folder_name_H))
                        MakingInpcar(s,"{}/{}/INPCAR".format(pwd,folder_name_H),nelect[e][1],kpoints[e][1])

            for runf in runfolder :
                rs = RunningShell(shell = inputs["SHELL"][0],name=runf, path=os.path.join(pwd,runf))
                if mt == "E" :
                    emc = GMDAnalysis()
                    emc.effectivemass(path="{}/{}".format(pwd,runf),secondstep=False)
                    os.chdir(pwd)

                if args.soc :
                    incar = open("{}/{}/INCAR".format(pwd,runf),'r').readlines()
                    index = [e for e,inc in enumerate(incar) if "MAGMOM" in inc]
                    del incar[index[0]]
                    with open("{}/{}/INCAR".format(pwd,runf),'w') as fi :
                        for f in incar :
                            fi.write(f)
                    fi.close()
                rs.running_mode(soc=args.soc, run=False)
                print("Generate folder {}/{}".format(pwd,runf))
    else :
        for stru in strucpath :
            f = open("{}/calnohup.py".format(os.path.dirname(__file__)),"r")
            ff = f.readlines()
            ff.insert(0,"inputs={}\n".format(dict(inputs)))
            ff.insert(1,"path='{}'\n".format(str(stru)))
            if args.deleteselect :
                ff.insert(2,"ds=True\n")
            else :
                ff.insert(2,"ds=False\n")
            if args.soc :
                ff.insert(3,"orbit=True\n")
            else :
                ff.insert(3,"orbit=False\n")
            ff.insert(4,"mole=False")
            ff.insert(5,"\n")

            with open("{}/running_calnohup.py".format(os.path.dirname(__file__)),"w") as fi :
                for i in ff :
                    fi.write(i)
            fi.close()
            #os.system("python {}/running_calnohup.py".format(os.path.dirname(__file__)))
            name = stru.split("/")[-1].split(".")[0]
            os.system("nohup python -u {}/running_calnohup.py > output_{}.gmd &".format(os.path.dirname(__file__),name))
            time.sleep(1)