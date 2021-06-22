import os
import sys
import time
import subprocess

from shutil import copyfile
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun

from perovgen.pygmd.base import load_structure, GMDStructure
from perovgen.pygmd.autocal.cdb import *
from perovgen.pygmd.autocal.auto_process import openingphrase
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
            print()
            print("="*10)
            print("We can't find file or directory", strucpath)
            print("="*10)
            sys.exit(1)

    # start calculation
    if args.directory :
        for mt in inputs["METHOD"] :
            os.chdir(pwd)
            runfolder =[]
            folder_list=[folder for folder in os.listdir(".") if "{}_mode".format(mt) in folder]
            nf = len(folder_list)
            if mt == "B" or mt == "D" or mt=="E" :
                    path = [os.path.abspath(b) for b in strucpath]
                    kpath_list = [CopyCHGCAR(i) for i in path]
                    if mt == "E" :
                        nelect=[];kpoints=[]
                        for b in path :
                            if os.path.isfile(b) :
                                b = os.path.split(b)[0]
                            a = subprocess.check_output(['grep','NELECT','%s/OUTCAR'%(b)])
                            nelect.append(int(float(a.split()[2])/2))
                            bsp = BSPlotting(vasprun=os.path.abspath("{}/vasprun.xml".format(b)), kpoints=os.path.abspath("{}/KPOINTS".format(b)))
                            try : 
                                vbm = bsp.bs.kpoints[bsp.bsdict['vbm']['kpoint_index'][0]].as_dict()['fcoords']
                            except : 
                                vbm = [0.000, 0.000, 0.000]
                            try :
                                cbm = bsp.bs.kpoints[bsp.bsdict['cbm']['kpoint_index'][0]].as_dict()['fcoords']
                            except : 
                                cbm = [0.000, 0.000, 0.000]
                            nelect.append((int(float(a.split()[2])/2),int(float(a.split()[2])/2)+1))
                            kpoints.append((vbm, cbm))

            # Make Folder
            struclist = [load_structure(strucpath)][0]
            if not struclist : 
                print("Please Enter the strucfile")
                sys.exit(0)

            for e,s in enumerate(struclist) :
                nf+=1
                v = GMDStructure(s).vaspname()
                p = PerovInputs(structure=s)
                if args.deleteselect :
                    p = p._poscar(delete_selective=True)

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
                folder_name ="%s_%i_%s_mode"%(v,nf,mt) 
                vi.write_input(output_dir=folder_name)
                runfolder.append(folder_name)
                if mt == "E" :
                    vi.write_input(output_dir="%s_%i_H_mode"%(v,nf))
                    runfolder.append("%s_%i_H_mode"%(v,nf))

                if mt == "D" or mt == "B" :
                    copyfile(kpath_list[e],"{}/{}/CHGCAR".format(pwd,folder_name))
                    if mt == "B" :
                        MakingKpointBand(s,"{}/{}/KPOINTS".format(pwd,folder_name))
                elif mt == "E" :
                    copyfile(kpath_list[e],"{}/{}/CHGCAR".format(pwd,folder_name))
                    copyfile(kpath_list[e],"%s/%s_%i_H_mode/CHGCAR"%(pwd,v,nf))
                    MakingInpcar(s,"{}/{}/INPCAR".format(pwd,folder_name),nelect[e][1],kpoints[e][1])
                    MakingInpcar(s,"%s/%s_%i_H_mode/INPCAR"%(pwd,v,nf),nelect[e][0],kpoints[e][0])

            for runf in runfolder :
                rs = RunningShell(shell = inputs["SHELL"][0],name=runf, path=os.path.join(pwd,runf))
                if mt == "E" :
                    emc = GMDAnalysis()
                    emc.effectivemass(path="{}/{}".format(pwd,runf),secondstep=False)
                    os.chdir(pwd)
                rs.running_mode(soc=args.soc, run=False)
                print("Generate folder {}/{}".format(pwd,runf))
    else :
        f = open("{}/calnohup.py".format(os.path.dirname(__file__)),"r")
        ff = f.readlines()
        ff.insert(0,"inputs={}\n".format(dict(inputs)))
        ff.insert(1,"path={}\n".format(str(strucpath)))
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
        os.system("nohup python {}/running_calnohup.py > output.gmd &".format(os.path.dirname(__file__)))