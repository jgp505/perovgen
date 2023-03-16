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
import logging

import yaml
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun, Outcar, Eigenval
from pymatgen.io.vasp.inputs import Kpoints, Kpoints_supported_modes
from pymatgen.core import Structure, Element

from perovgen.pygmd.input_structure import load_structure, GMDStructure
from perovgen.pygmd.autocal.inputset import *
from perovgen.pygmd.autocal.substitute import _selective_bool
from perovgen.pygmd.analysis.electronic import BSPlotting, DOSPlotting, GMDAnalysis
from perovgen.pygmd.analysis.energy import GMDPhaseDiagram, GMDExcitonbinding
from perovgen.pygmd.autocal.algorithm import openingphrase, Recalculate, fileopen

def Process(calpath,inputspath, strucpath, ds=False, soc=False):
    
    struc, filename = load_structure(strucpath)[0][0], load_structure(strucpath)[1][0]
    inputs = inputgmd(inputspath)
    fn = os.path.basename(filename).split(".")[0]

    for k,v in inputs.calmode.items() :
        poscar = PerovInputs(struc)
        inputs = inputgmd(inputspath)
        naming = poscar.naming + "_" + fn
        
        runfolder = []
        if k == 'bulkmodulus' :
            volumes = v*poscar.poscar.volume
            for v1, vol in zip(v, volumes) :
                struc.scale_lattice(vol)
                PI = PerivInputs(struc)
                vi = PI.inputfolder(inputs=inputs, method=k, soc=soc)
                folder_name = "%s/%s/%.2f"%(calpath,k,v1)
                vi.write_input(folder_name)
                runfolder.append(folder_name)

        elif k == 'dos' or k == 'band' or k == 'effective_mass' :
            vi = poscar.inputfolder(inputs=inputs, method=k, soc=soc)
            if k == 'effective_mass' :
                folder_name_elec = "%s/chgcar/%s/electron"%(calpath,k)
                folder_name_hole = "%s/chgcar/%s/hole"%(calpath,k)
                
                # CBM effective mass
                vi.write_input(folder_name_elec)
                copy("%s/chgcar/CHGCAR"%(calpath), folder_name_elec)
                copy("%s/chgcar/INPCAR_ele"%(calpath), "%s/INPCAR"%(folder_name_elec))
                runfolder.append(folder_name_elec)
                
                # VBM effective mass
                vi.write_input(folder_name_hole)
                copy("%s/chgcar/CHGCAR"%(calpath), folder_name_hole)
                copy("%s/chgcar/INPCAR_hole"%(calpath), "%s/INPCAR"%(folder_name_hole))
                runfolder.append(folder_name_hole)
            else :
                folder_name = "%s/chgcar/%s"%(calpath,k)
                vi.write_input(folder_name)
                copy("%s/chgcar/CHGCAR"%(calpath),folder_name)
                runfolder.append(folder_name)
        else :
            folder_name = "%s/%s"%(calpath,k)
            vi = poscar.inputfolder(inputs=inputs, method=k, soc=soc)
            vi.write_input(folder_name)
            runfolder.append(folder_name)
            
        # Running
        e = 0
        for runf in runfolder :
            if k == 'bulkmodulus' :
                rs = RunningShell(shell=inputs.shell, name="%s_%s_%.2f"%(naming,k,v[e]),path=runf)
                e += 1
            elif k == 'effective_mass' :
                emc = GMDAnalysis()
                emc.effectivemass(runf, secondstep=False)
                rs = RunningShell(shell=inputs.shell, name="%s_%s"%(naming,k), path = runf)
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
                    if k == 'relaxation' :
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
            os.chdir("%s/chgcar/%s"%(calpath, k))
            if k == 'dos' :
                os.system("%s dos width=0.03"%(analysis_path.pdos))
                analysis_path.partialDOS(structure=struc)
            elif k == 'band' :
                bsp = BSPlotting(vasprun='vasprun.xml',kpoints='KPOINTS')
                bsp.write_to_json()
                nelect = Eigenval("EIGENVAL").nelect
                nelect = (int(nelect/2),int(nelect/2)+1)
                try :
                    vbm = bsp.bs.kpoints[bsp.bsdict['vbm']['kpoint_index'][0]].as_dict()['fcoords']
                    cbm = bsp.bs.kpoints[bsp.bsdict['cbm']['kpoint_index'][0]].as_dict()['fcoords']
                except :
                    vbm = [0.000, 0.000, 0.000]
                    cbm = [0.000, 0.000, 0.000]
                MakingInpcar(struc, "%s/chgcar/INPCAR_ele"%(calpath),nelect[1],cbm)
                MakingInpcar(struc, "%s/chgcar/INPCAR_hole"%(calpath),nelect[0],vbm)
            os.remove("%s/chgcar/%s/CHGCAR"%(calpath,k))
            
        elif k == 'effective_mass' :
            emc = GMDAnalysis()
            for runf in runfolder :
                emc.effectivemass(runf, secondstep=True)
            struc = Structure.from_file("%s/chgcar/CONTCAR"%(calpath))             
            os.remove("%s/chgcar/%s/electron/CHGCAR"%(calpath,k))
            os.remove("%s/chgcar/%s/hole/CHGCAR"%(calpath,k))
            
        elif k == 'HSE' :
            vrun = BSVasprun("vasprun.xml")
            bandgap = vrun.as_dict()['output']['bandgap']
            with open("hse_bandgap.gmd",'w') as fi :
                fi.write("HSE mode")
                fi.write("%.5f\n"%(bandgap))
            fi.close()
        
        #elif k == 'absorption' :
        
        # relaxed structure
        if k == 'bulkmodulus' : 
            struc = Structure.from_file("%s/bulkmodulus/POSCAR_optimized"%(calpath))
        elif k == 'relaxation' :
            struc = Structure.from_file("%s/relaxation/CONTCAR"%(calpath))
        elif 'chgcar' in [*inputs.calmode] :
            if k == 'dos' or k == 'band' or k == 'effective_mass' :
                struc = Structure.from_file("%s/chgcar/CONTCAR"%(calpath))
        else :
            struc = load_structure(strucpath)[0][0]
        os.chdir(calpath)           
                