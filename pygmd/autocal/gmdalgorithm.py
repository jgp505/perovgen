# coding : utf-8
# Cppyright (c) Green Materials Designs Team.
'''
This module provides 
'''
import os
import sys
import time
import subprocess
import glob
from shutil import copyfile
from tabulate import tabulate

import yaml
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.io.vasp.inputs import Kpoints, Kpoints_supported_modes
from pymatgen.core import Structure, Element

from perovgen.pygmd.input_structure import load_structure, GMDStructure
from perovgen.pygmd.autocal.inputset import *
from perovgen.pygmd.autocal.substitute import _selective_bool
from perovgen.pygmd.analysis.electronic import BSPlotting, DOSPlotting, GMDAnalysis
from perovgen.pygmd.analysis.energy import GMDPhaseDiagram, GMDExcitonbinding


pwd = os.getcwd()
def controlincar(namelist):
    """
    controlincar
    """
    # Control incar 
    if namelist["INCAR"][0] in ["PBESol","VDW","SCAN","PBE","MPJ"]:
        stream = open("%s/GMDincar.yaml"%(os.path.dirname(__file__)),'r')
        incar = yaml.load(stream)
        if namelist["INCAR"][0] == "PBESol" :
            incar["GGA"]='Ps'
        elif namelist["INCAR"][0] == "VDW" :
            incar["IVDW"]=21
        elif namelist["INCAR"][0] == "SCAN" :
            incar["METAGGA"] = "SCAN"
            incar["LUSE_VDW"] = True
            incar["BPARAM"]=15.7
            incar["LASPH"]=True
        elif namelist["INCAR"][0] == "MPJ" :
            incar = {"SIGMA":0.2,'ISMEAR':0,"GGA_COMPAT":False,"LASPH":True}
        else : 
            print("Please Enter the PBESol, PBE, VDW, SCAN, PBE, MPJ")

        if len(namelist["INCAR"]) != 1 :
            incar1 = Incarmethod(namelist["INCAR"][1:])
            for k,v in incar1.items() :
                incar[k]=v
        else :
            pass
    else :
        incar = Incarmethod(namelist["INCAR"])
    return incar

def openingphrase(inputs,strucpath):
    text = """
    version 3.6.0
         
                 W E C O M E 
                GMD AUTO MODE

    copyright @ Hanbat National University, Korea
    created by Jong Goo Park
            """

    table = [[text]]
    output = tabulate(table, tablefmt='grid')
    print(output)
    print("Method of Exchange Correlation : ",inputs["INCAR"][0])
    print("The number of the calculated structures : ",len(load_structure(strucpath)[0]))
    print("Calculation mode : ", inputs["METHOD"])
    print("Running Shell name : ", inputs["SHELL"][0])
    print("Running start time : ",time.strftime("%c",time.localtime(time.time())))

def fileopen(data, target):
    with open(data,'r',encoding='UTF8') as file:
        text = file.read()
        if target in text   :
            flag = True
            splitdata = text.split()
        else :
            flag = False
            splitdata = None
    return flag, splitdata

def Recalculate():
    subprocess.call(['rm','vasprun.xml'])
    subprocess.call(['rm','OUTCAR'])
    folder_list=[folder for folder in os.listdir(".") if folder.endswith("initial")]
    filelist = [u for u in os.listdir(".") if "energy" == u]
    if not folder_list :
        subprocess.call(['cp','POSCAR','POSCAR_initial'])
    subprocess.call(['cp','CONTCAR','POSCAR'])
    subprocess.check_call(['qsub','vasp.sh'])

def Process(inputs, strucpath, ds=False, orbit=False):
    openingphrase(inputs,strucpath)

    # start calculation
    pwd = os.getcwd()
    e = 0
    s,filenames = load_structure(strucpath)
    s = s[0] ; filenames = filenames[0]
    if not s : 
        print("\n[Warning] the structure file does not exist.\n")
        sys.exit(0)

    for mt in inputs["METHOD"] :
        os.chdir(pwd)
        runfolder = []

        # copy where the CHGCAR is located
        if mt == "B" or mt == "D" or mt=="E" :
            path = [os.path.abspath(strucpath[0])]
            kpath_list = [CopyCHGCAR(path[0])]
            print(path, kpath_list)

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

                    if orbit : 
                        nelect.append((int(float(a.split()[2])),int(float(a.split()[2]))+1))
                    else :
                        nelect.append((int(float(a.split()[2])/2),int(float(a.split()[2])/2)+1))
                    kpoints.append((vbm, cbm))
	
        full_formula = GMDStructure(s).formula(reduced=False)
        pretty_formula = GMDStructure(s).formula()
        
        if ds :
            p = PerovInputs(structure=s,is_selective=True)
        else :
            p = PerovInputs(structure=s)

        # Control incar
        incar = GMDIncar(inputs).incar
        modeincar = PerovInputs._incarmode(incar=incar, method=mt)

        print(modeincar)
        if orbit :
            modeincar["LSORBIT"]=True
    
        # Writing the folder 
        if not "None" in inputs["INCAR"] :
            vi = p.inputfolder(incar=modeincar, number=inputs["KPOINTS"][0])
        else :
            vi = p.inputfolder(incar=modeincar, number=None)

        try :
            symmetry, groupnumber = s.get_space_group_info()
        except :
            groupnumber = 0
        fn = filenames.split("/")[-1].split(".")[0]

        # [pretty formula]_[symmetry number]/[mode]_[full_formula]_[filename]
        folder_name = "{0}_{1}_{2:02d}".format(mt,full_formula,groupnumber) 
        vi.write_input(output_dir=folder_name)
        runfolder.append(folder_name)

        if mt == "E" :
            folder_name_H = "H_{0}_{1:02d}".format(full_formula,groupnumber) 
            vi.write_input(output_dir=folder_name_H)
            runfolder.append(folder_name_H)

        # Copy the other files to generated folder
        if mt == "D" or mt=="B" or mt == "E" :
            copyfile(kpath_list[e],"{}/CHGCAR".format(folder_name))
        if mt == "B" :
            # Revise the KPOINTS for calculating band
            MakingKpointBand(s,"{}/KPOINTS".format(folder_name))
        elif mt == "E" :
            # generate INPCAR file to execute emc-master modules
            copyfile(kpath_list[e],"{}/CHGCAR".format(folder_name))
            MakingInpcar(s,"{}/INPCAR".format(folder_name),nelect[e][1],kpoints[e][1])
            copyfile(kpath_list[e],"{}/CHGCAR".format(folder_name_H))
            MakingInpcar(s,"{}/INPCAR".format(folder_name_H),nelect[e][1],kpoints[e][1])
        elif mt == 'D' :
            # change the mode
            kpoints = Kpoints.from_file("{}/KPOINTS".format(folder_name))
            kpoints.style = Kpoints_supported_modes.Gamma
            kpoints.write_file("{}/KPOINTS".format(folder_name))


        for runf in runfolder :
            naming = runf.split("/")[-1]
            rs = RunningShell(shell = inputs["SHELL"][0],name=naming, path=runf)
            if mt == "E" :
                emc = GMDAnalysis()
                emc.effectivemass(path="{}".format(runf),secondstep=False)
                os.chdir(pwd)
            
            if orbit :
                incar = open("{}/INCAR".format(runf),'r').readlines()
                index = [e for e,inc in enumerate(incar) if "MAGMOM" in inc]
                del incar[index[0]]
                with open("{}/INCAR".format(runf),'w') as fi :
                    for f in incar :
                        fi.write(f)
                fi.close()

            rs.running_mode(soc=orbit, run=True)

        # Running Check
        while True :
            time.sleep(10)
            path1 = [] 
            for j in runfolder :
                os.chdir(os.path.join(pwd,j))
                try :
                    vrun = Vasprun("%s/vasprun.xml"%(os.path.join(pwd,j)),parse_potcar_file=True)
                    ionicsteps = vrun.nionic_steps
                    nsw = vrun.incar['NSW']
                    if nsw == 1 :
                        path1.append("%s/CONTCAR"%(os.path.join(pwd,j)))
                    elif nsw == ionicsteps :
                        print("[Notice] Realculation because ionic step is same NSW value")
                        Recalculate()
                        time.sleep(10)
                    else :
                        path1.append("%s/CONTCAR"%(os.path.join(pwd,j)))
                except ET.ParseError:
                    # Error Check in R-mode
                    if mt == "R" :
                        boolen, targetlist = fileopen("OUTCAR","accuracy")
                        if not boolen : 
                            try :
                                number = subprocess.check_output(['tail','-n','1','OSZICAR']).decode('utf-8')
                                number1 = number.split()[0]
                                time.sleep(180)
                                number2 = subprocess.check_output(['tail','-n','1','OSZICAR']).decode('utf-8')
                                number3 = number2.split()[0]
                                if int(number1) == int(number3) :
                                    print("[Notice] Realculation because it has not yet obtained a stabilizing structure.")
                                    Recalculate()
                                    time.sleep(10)
                                else :
                                    pass
                            except :
                                pass
                        else :
                            path1.append("%s/CONTCAR"%(os.path.join(pwd,j)))
                    else :
                        pass
                except FileNotFoundError :
                    pass
                except AttributeError :
                    pass

            if len(runfolder) == len(path1) :
                os.chdir(pwd)
                break

        # Properties for DOS and effective mass
        if mt == "E" :
            for i in runfolder :
                i = os.path.abspath(i)
                analysis_path = GMDAnalysis()
                analysis_path.effectivemass(path=i, secondstep=True)
                em_e = open("{}/EM".format(i),'r').readlines()[-12:]
                E = GMDExcitonbinding.harm_mean_em(em_e)
                print("{} of effective mass is ".format(os.path.split(i)[-1]), E)
            os.chdir(pwd)

        index = [path1[j] for j in range(len(path1)) if "H_mode" in path1[j]]
        if index :
            for j in index : 
                u = path1.index(j)
                del path1[u]

        # CONTCOAR TO POSCAR 
        strucpath = path1
        print("\n%s mode finished time : "%(mt),end="")
        print(time.strftime("%c\n",time.localtime(time.time())))