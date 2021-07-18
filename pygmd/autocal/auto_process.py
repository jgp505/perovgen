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
from pymatgen.core import Structure, Element

from perovgen.pygmd.base import load_structure, GMDStructure, pdos, fileopen, count_word
from perovgen.pygmd.autocal.cdb import *
from perovgen.pygmd.autocal.substitute import _selective_bool
from perovgen.pygmd.analysis.electronic import BSPlotting, DOSPlotting, GMDAnalysis
from perovgen.pygmd.analysis.energy import GMDPhaseDiagram, GMDExcitonbinding


pwd = os.getcwd()
def controlincar(namelist):
    """
    docstring
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
    version 3.5.1
         
                 W E C O M E 
                GMD AUTO MODE

    copyright @ Hanbat National University, Korea
    created by Jong Goo Park
            """

    table = [[text]]
    output = tabulate(table, tablefmt='grid')
    print(output)
    print("Method of Exchange Correlation : ",inputs["INCAR"][0])
    print("The number of the calculated structures : ",len(load_structure(strucpath)))
    print("Calculation mode : ", inputs["METHOD"])
    print("Running Shell name : ", inputs["SHELL"][0])
    print("Running start time : ",time.strftime("%c",time.localtime(time.time())))

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
            else :
                pass

        # Make Folder
        struclist = load_structure(strucpath)
        if not struclist : 
            print("Please Enter the strucfile")
            sys.exit(0)

        for e,s in enumerate(struclist) :
            nf+=1
            v = GMDStructure(s).vaspname()
            if ds :
                p = PerovInputs(structure=s,is_selective=True)
            else :
                p = PerovInputs(structure=s)

            # Control incar 
            incar = GMDIncar(inputs).incar
            modeincar = PerovInputs._incarmode(incar=incar, method=mt)
            if orbit :
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

            if mt == "D" or mt=="B" :
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
            
            if orbit :
                incar = open("{}/{}/INCAR".format(pwd,runf),'r').readlines()
                index = [e for e,inc in enumerate(incar) if "MAGMOM" in inc]
                del incar[index[0]]
                with open("{}/{}/INCAR".format(pwd,runf),'w') as fi :
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
                        #print("Realculation")
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
                                    print("Recalculation")
                                    Recalculate()
                                    time.sleep(10)
                                else :
                                    pass
                            except ValueError :
                                pass
                            except IndexError:
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
        if mt == "D" or mt == "E" :
            for i in runfolder :
                i = os.path.abspath(i)
                analysis_path = GMDAnalysis()
                if mt == "D" :
                    os.chdir(i)
                    os.system("%s dos width=0.03"%(analysis_path.pdos))
                    print("## Total DOS is Done ##")
                    s = load_structure(".")[0]
                    path = analysis_path.partialDOS(structure=s)
                    print("\n## Partial DOS is Done ##")
                    #dos = open("{}/dos".format(i),'r').readlines()[1:]
                    #fig = plt.figure(figsize=(12,8))
                    #DOSPlotting(vasprun="vasprun.xml",dos=dos).get_plot().savefig("{}.pdf".format(i.split("/")[-1]))
                else :
                    analysis_path.effectivemass(path=i, secondstep=True)
                    em_e = open("{}/EM".format(i),'r').readlines()[-12:]
                    E = GMDExcitonbinding.harm_mean_em(em_e)
                    print("{} of effective mass is ".format(os.path.split(i)[-1]), E)
                os.chdir(pwd)

        # CONTCOAR TO POSCAR 
        index = [path1[j] for j in range(len(path1)) if "H_mode" in path1[j]]
        if index :
            for j in index : 
                u = path1.index(j)
                del path1[u]

        strucpath = path1
        print("\n%s mode finished time : "%(mt),end="")
        print(time.strftime("%c\n",time.localtime(time.time())))

def AutoMolOpt(strucpath, inputs) :
    pwd = os.getcwd()
    for isif in ["ISIF7", "ISIF2", "ISIF3"] :
        runfolder = [];struclist2=[]
        for path in strucpath :
            #os.chdir(path)
            if isif == "ISIF7" :
                struclist=load_structure(path)
            elif isif == "ISIF2" :
                path1 = ["{}/ISIF7/{}/CONTCAR".format(path,f) for f in os.listdir("{}/ISIF7".format(path)) if "M_mode" in f]
                struclist = load_structure(path1)
            elif isif == "ISIF3" :
                path1 = ["{}/ISIF2/{}/CONTCAR".format(path,f) for f in os.listdir("{}/ISIF2".format(path)) if "M_mode" in f]
                struclist = load_structure(path1)
            makepath = "{}/{}".format(path,isif)
            os.makedirs(makepath)

            nf = 0
            for e,s in enumerate(struclist) :
                nf+=1
                g = GMDStructure(s)
                v = g.vaspname()

                # revised poscar
                if isif == "ISIF2" :
                    sp = _selective_bool(species=g.species,type_num=4)
                    new_s = Structure(s.lattice.matrix, g.species, g.coords, coords_are_cartesian=True,site_properties={"selective_dynamics":sp})
                    p = PerovInputs(structure=new_s)
                else :
                    p = PerovInputs(structure=s,is_selective=True)
                    
                # Control incar 
                incar = GMDIncar(inputs).incar
                modeincar = PerovInputs._incarmode(incar=incar, method="M")
                if isif == "ISIF2" : 
                    modeincar["ISIF"]=2
                elif isif == "ISIF7" : 
                    modeincar["ISIF"]=7
                    modeincar["NSW"]=21
                    del modeincar["EDIFFG"]
                elif isif == "ISIF3" :
                    modeincar["ISIF"]=3
                    modeincar["EDIFF"]=1E-6
                    modeincar["EDIFFG"]=-0.01
                    modeincar["ENCUT"]=520

                # Writing the folder 
                if isif == "ISIF7" or isif == "ISIF2" :
                    vi = p.inputfolder(incar=modeincar, number=20)
                else :
                    if inputs["INCAR"][0] != "MPJ" :
                        vi = p.inputfolder(incar=incar, number=40)
                    else :
                        vi = p.inputfolder(incar=incar, number=None)
                folder_name ="%s_%i_M_mode"%(v,nf)
                vi.write_input(output_dir="{}/{}".format(makepath,folder_name))
                runfolder.append("{}/{}".format(makepath,folder_name))

        for runf in runfolder :
            rs = RunningShell(shell = inputs["SHELL"][0],name=os.path.split(runf)[-1], path=runf)
            rs.running_mode(soc=False, run=True)
        # Running Check
        while True :
            time.sleep(10)
            path1 = [] 
            for j in runfolder :
                os.chdir(j)
                try :
                    vrun = Vasprun("vasprun.xml",parse_potcar_file=True)
                    ionicsteps = vrun.nionic_steps
                    nsw = vrun.incar['NSW']
                    if nsw == ionicsteps :
                        Recalculate()
                    else :
                        path1.append("%s/CONTCAR"%(os.path.join(pwd,j)))
                except ET.ParseError:
                    pass
            if len(runfolder) == len(path1) :
                os.chdir(pwd)
                break
        
        print("\n%s finished time : "%(isif),end="")
        print(time.strftime("%c\n",time.localtime(time.time())))

    for path in strucpath : 
        path1 = ["{}/ISIF3/{}/vasprun.xml".format(path,f) for f in os.listdir("{}/ISIF3".format(path)) if "M_mode" in f]
        energy = [Vasprun(p1,parse_potcar_file=True).final_energy for p1 in path1] 
        index = energy.index(min(energy))
        direct = os.path.dirname(path1[index])
        os.system("mv {} {}/MINIMUM".format(direct, os.path.dirname(direct)))