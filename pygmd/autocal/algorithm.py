
# coding: utf-8
# Cppyright (c) Green Materials Designs Team.
# Hanbat National University, Korea 

import os
import sys
import time
import subprocess
import glob
from shutil import copyfile, copy
from tabulate import tabulate
import logging

import yaml
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

from pymatgen.io.vasp.outputs import Vasprun, Eigenval, BSVasprun
from pymatgen.io.vasp.inputs import Incar
from pymatgen.core import Structure
from pymatgen.analysis.eos import EOS

from perovgen.pygmd.input_structure import load_structure
from perovgen.pygmd.autocal.inputset import *
from perovgen.pygmd.analysis.electronic import BSPlotting, GMDAnalysis
from perovgen.pygmd.analysis.energy import GMDAbsorption, GMDExcitonbinding

def openingphrase(inputs, strucpath):
    
    text = """
    version 3.8.5
         
                 W E C O M E 
                   PEROVGEN

    copyright @ Hanbat National University, Korea
    created by Jong Goo Park
            """

    table = [[text]]
    output = tabulate(table, tablefmt='grid')
    print(output)
    print("Method of Exchange Correlation : ",inputs.exchange_corr)
    print("The number of the calculated structures : ",len(strucpath))
    print("Running Shell name : ", inputs.shell)
    print("Calculation mode : ", [*inputs.calmode])
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

def Process(calpath,inputspath, strucpath, ds=False, soc=False):
    struc, filename = load_structure(strucpath)[0][0], load_structure(strucpath)[1][0]
    inputs = inputgmd(inputspath)
    fn = os.path.basename(filename).split(".")[0]

    #Generate Folder
    for k,v in inputs.calmode.items():
        poscar = PerovInputs(struc)
        inputs = inputgmd(inputspath)
        naming = poscar.naming+"_"+fn
        runfolder = []
        if k == 'bulkmodulus' :
            volumes = v*poscar.poscar.volume
            for v1, vol in zip(v, volumes) :
                struc.scale_lattice(vol)
                PI = PerovInputs(struc)
                vi = PI.inputfolder(inputs=inputs, method=k, soc=soc)
                folder_name = "%s/%s/%.2f"%(calpath,k,v1)
                if soc :
                    folder_name = "%s/%s/%.2f_SOC"%(calpath,k,v1)
                vi.write_input(folder_name)
                runfolder.append(folder_name)
                
        elif k == 'dos' or k == 'band' or k == 'effective_mass' :
            vi = poscar.inputfolder(inputs=inputs, method=k, soc=soc)
            if k == 'effective_mass' :
                folder_name_elec = "%s/chgcar/%s/electron"%(calpath,k)
                if soc :
                    folder_name_elec = "%s/chgcar_SOC/%s_SOC/electron"%(calpath,k)
                folder_name_hole = "%s/chgcar/%s/hole"%(calpath,k)
                if soc :
                    folder_name_hole = "%s/chgcar_SOC/%s_SOC/hole"%(calpath,k)
                
                # CBM effective mass
                vi.write_input(folder_name_elec)
                if soc :
                    copy("%s/chgcar_SOC/CHGCAR"%(calpath), folder_name_elec)
                    copy("%s/chgcar_SOC/INPCAR_ele"%(calpath), "%s/INPCAR"%(folder_name_elec))
                else :
                    copy("%s/chgcar/CHGCAR"%(calpath), folder_name_elec)
                    copy("%s/chgcar/INPCAR_ele"%(calpath), "%s/INPCAR"%(folder_name_elec))
                runfolder.append(folder_name_elec)
                
                # VBM effective mass
                vi.write_input(folder_name_hole)
                if soc :
                    copy("%s/chgcar_SOC/CHGCAR"%(calpath), folder_name_hole)
                    copy("%s/chgcar_SOC/INPCAR_hole"%(calpath), "%s/INPCAR"%(folder_name_hole))
                else :
                    copy("%s/chgcar/CHGCAR"%(calpath), folder_name_hole)
                    copy("%s/chgcar/INPCAR_hole"%(calpath), "%s/INPCAR"%(folder_name_hole))
                runfolder.append(folder_name_hole)
            else :
                folder_name = "%s/chgcar/%s"%(calpath,k)
                if soc :
                    folder_name = "%s/chgcar_SOC/%s_SOC"%(calpath,k)
                vi.write_input(folder_name)
                if soc :
                    copy("%s/chgcar_SOC/CHGCAR"%(calpath),folder_name)
                else :
                    copy("%s/chgcar/CHGCAR"%(calpath),folder_name)
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
                folder_name_nsoc = "%s/HSE"%(calpath)
                folder_name_soc = "%s/HSE_SOC"%(calpath)
                vi1.write_input(folder_name_nsoc)
                vi2.write_input(folder_name_soc)
                runfolder.append(folder_name_nsoc)
                runfolder.append(folder_name_soc)
            else :
                folder_name = "%s/%s"%(calpath,k)
                if soc :
                    folder_name = "%s/%s_SOC"%(calpath,k)
                vi = poscar.inputfolder(inputs=inputs, method=k, soc=soc)
                vi.write_input(folder_name)
                runfolder.append(folder_name)
        else :
            folder_name = "%s/%s"%(calpath,k)
            if soc :
                folder_name = "%s/%s_SOC"%(calpath,k)
            vi = poscar.inputfolder(inputs=inputs, method=k, soc=soc)
            vi.write_input(folder_name)
            runfolder.append(folder_name)
            
        # Running
        e = 0
        for runf in runfolder :
            if k == 'bulkmodulus' :
                rs = RunningShell(shell=inputs.shell, name="%s_%s_%.2f"%(naming,k,v[e]),path=runf)
                e += 1
                rs.running_mode(soc=soc, run=True) # If run is true, running the calculation
            elif k == 'effective_mass' :
                emc = GMDAnalysis()
                emc.effectivemass(runf, secondstep=False)
                rs = RunningShell(shell=inputs.shell, name="%s_%s"%(naming,k), path = runf)
                rs.running_mode(soc=soc, run=True) # If run is true, running the calculation
            elif k == 'HSE' :
                runf1 = os.path.basename(runf)
                if "HSE" == runf1 :
                    rs = RunningShell(shell=inputs.shell, name="%s_HSE"%(naming),path=runf)
                    rs.running_mode(soc=False, run=True) # If run is true, running the calculation
                elif "HSE_SOC" == runf1 :
                    rs = RunningShell(shell=inputs.shell, name="%s_HSE_SOC"%(naming),path=runf)
                    rs.running_mode(soc=True, run=True) # If run is true, running the calculation
            else :
                rs = RunningShell(shell=inputs.shell, name="%s_%s"%(naming,k), path = runf)
                rs.running_mode(soc=soc, run=True) # If run is true, running the calculation
            
        # Running Check
        while True :
            time.sleep(10)
            path1 = []
            for realpath in runfolder :
                os.chdir(realpath)
                try : 
                    vrun = Vasprun("%s/vasprun.xml"%(realpath),parse_potcar_file=True)
                    ionicsteps = vrun.nionic_steps
                    nsw = vrun.incar['NSW']
                    if nsw == 1 :
                        path1.append("%s/CONTCAR"%(realpath))
                        os.remove("%s/CHG"%(realpath))
                    elif nsw == ionicsteps :
                        print("[NOTICE] Recalculation")
                        Recalculate()
                        time.sleep(10)
                    else :
                        path1.append("%s/CONTCAR"%(realpath))
                        os.remove("%s/CHG"%(realpath))
                except ET.ParseError :
                    # Error check in R-mode
                    incar = Incar.from_file("%s/INCAR"%(realpath))
                    if k == 'relaxation' and incar['ISIF'] == 3 : # molecular calculation is traped by ISIF
                        boolen, targetlist = fileopen("OUTCAR","accuracy")
                        if not boolen :
                            try :
                                number = subprocess.check_output(['tail','-n','1','OSZICAR']).decode('utf-8')
                                number1 = number.split()[0]
                                time.sleep(180)
                                number2 = subprocess.check_output(['tail','-n','1','OSZICAR']).decode('utf-8')
                                number2 = number2.split()[0]
                                if int(number1) == int(number2) :
                                    print("[NOTICE] Recalculation")
                                    Recalculate()
                                    time.sleep(10)
                                else :
                                    pass
                            except :
                                pass
                        else :
                            path1.append("%s/CONTCAR"%(realpath))
                            os.remove("%s/CHG"%(realpath))
                except :
                    pass
            if len(runfolder) == len(path1) :
                os.chdir(os.getcwd())
                break
            
        # Ananlysis
        if k == 'bulkmodulus' :
            os.chdir("%s/%s"%(calpath, k))
            direct = [f for f in os.listdir(".") if os.path.isdir(f)]
            direct.sort()
            energy, volume = [], []
            for d in direct :
                vrun = Vasprun("%s/vasprun.xml"%(d))
                energy.append(vrun.final_energy)
                volume.append(vrun.final_structure.volume)
                struc = vrun.final_structure        

            with open("bulkmodulus.gmd", 'w') as fi :
                for e,v in zip(energy,volume) :
                    fi.write("%.4f %.4f\n"%(v,e))
            try :
                eos = EOS(eos_name='birch_murnaghan')
                eos_fit = eos.fit(volume, energy)
                struc.scale_lattice(eos_fit.v0)
            except :
                indexs = energy.index(min(energy))
                min_volume = volume[indexs]
                struc.scale_lattice(min_volume)
            struc.to(filename="POSCAR_optimized")
            
        elif k == 'dos' or k == 'band' :
            analysis_path = GMDAnalysis()
            if soc :
                os.chdir("%s/chgcar_SOC/%s_SOC"%(calpath, k))
                os.remove("%s/chgcar_SOC/%s_SOC/CHGCAR"%(calpath,k))
            else :
                os.chdir("%s/chgcar/%s"%(calpath, k))
                os.remove("%s/chgcar/%s/CHGCAR"%(calpath,k))
            if k == 'dos' :
                os.system("%s dos width=0.03"%(analysis_path.pdos))
                analysis_path.partialDOS(structure=struc)
            elif k == 'band' :
                bsp = BSPlotting(vasprun='vasprun.xml',kpoints='KPOINTS')
                bsp.write_to_json(path=calpath)
                nelect = Eigenval("EIGENVAL").nelect
                if soc :
                    nelect = (int(nelect),int(nelect)+1)
                else :
                    nelect = (int(nelect/2),int(nelect/2)+1)
                try :
                    vbm = bsp.bs.kpoints[bsp.bsdict['vbm']['kpoint_index'][0]].as_dict()['fcoords']
                    cbm = bsp.bs.kpoints[bsp.bsdict['cbm']['kpoint_index'][0]].as_dict()['fcoords']
                except :
                    vbm = [0.000, 0.000, 0.000]
                    cbm = [0.000, 0.000, 0.000]
                if soc :
                    MakingInpcar(struc, "%s/chgcar_SOC/INPCAR_ele"%(calpath),nelect[1],cbm)
                    MakingInpcar(struc, "%s/chgcar_SOC/INPCAR_hole"%(calpath),nelect[0],vbm)
                else :
                    MakingInpcar(struc, "%s/chgcar/INPCAR_ele"%(calpath),nelect[1],cbm)
                    MakingInpcar(struc, "%s/chgcar/INPCAR_hole"%(calpath),nelect[0],vbm)
            
        elif k == 'effective_mass' :
            emc = GMDAnalysis()
            for runf in runfolder :
                emc.effectivemass(runf, secondstep=True)
            if soc :
                os.remove("%s/chgcar_SOC/%s_SOC/electron/CHGCAR"%(calpath,k))
                os.remove("%s/chgcar_SOC/%s_SOC/hole/CHGCAR"%(calpath,k))
            else :
                os.remove("%s/chgcar/%s/electron/CHGCAR"%(calpath,k))
                os.remove("%s/chgcar/%s/hole/CHGCAR"%(calpath,k))
            os.chdir(calpath)
            hole = GMDExcitonbinding.harm_mean_em("%s/chgcar/%s/hole/EM"%(calpath,k))
            electron = GMDExcitonbinding.harm_mean_em("%s/chgcar/%s/electron/EM"%(calpath,k))
            harm_mean = 2/((1/electron['HarmMean']) + (1/hole['HarmMean']))
            with open("effectivemass.gmd",'w') as fi :
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
            
        elif k == 'dielectric' :
            os.chdir(calpath)
            dielectric = GMDExcitonbinding.dielectricconst("{}/{}/OUTCAR".format(calpath,k))
            array1 = dielectric['Array']['high_frequency_dielectric']
            array2 = dielectric['Array']['ionic_contribution']
            high_mean = dielectric['geom_mean']
            with open ("dielectric.gmd",'w') as fi :
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
            
        elif k == 'HSE' :
            os.chdir(calpath)
            vrun = BSVasprun("HSE/vasprun.xml")
            bandgap = vrun.as_dict()['output']['bandgap']
            if hse_soc :
                vrun = BSVasprun("HSE_SOC/vasprun.xml")
                bandgap_soc = vrun.as_dict()['output']['bandgap']

            with open("hse_bandgap.gmd",'w') as fi :
                fi.write("HSE : %.5f\n"%(bandgap))
                if hse_soc :
                    fi.write("HSE_SOC : %.5f\n"%(bandgap_soc))
            fi.close()

        elif k == 'absorption' :
            pass
            #df = GMDAbsorption.make_csv("{}/{}".format(calpath,k))
            #df.to_csv("%s/absorp.csv"%(calpath))
        
        # relaxed structure
        if k == 'bulkmodulus' : 
            struc = Structure.from_file("%s/bulkmodulus/POSCAR_optimized"%(calpath))
        elif k == 'relaxation' or k == 'substitute_relax' :
            struc = Structure.from_file("%s/%s/CONTCAR"%(calpath,k))
        else :
            if 'chgcar' in [*inputs.calmode] :
                if soc :
                    struc = Structure.from_file("%s/chgcar_SOC/CONTCAR"%(calpath))
                else :
                    struc = Structure.from_file("%s/chgcar/CONTCAR"%(calpath))
            elif 'relaxation' in [inputs.calmode] :
                struc = Structure.from_file("%s/relaxation/CONTCAR"%(calpath))
            elif 'bulkmodulus' in [*inputs.calmode] :
                struc = Structure.from_file("%s/relaxation/CONTCAR"%(calpath))
            else :
                struc = load_structure(strucpath)[0][0]
        os.chdir(calpath)           