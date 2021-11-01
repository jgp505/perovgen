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
    inputs = inputgmd(inputspath); struc, filename = load_structure(strucpath)[0][0], load_structure(strucpath)[1][0]
    fn = os.path.basename(filename).split(".")[0]
    chgpath = None ; bandpath = None

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    file_handler = logging.FileHandler('{}/{}.log'.format(calpath,fn))
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    if not "C" in inputs.calmode :
        for d in os.listdir(calpath) :
            if os.path.isdir("{}/{}".format(calpath,d)):
                mode = d.split("_")[0]
                fn1 = '_'.join(d.split("_")[2:])
                if mode == 'C' and fn == fn1 :
                    chgpath = "{}/{}/CHGCAR".format(calpath,d)
        if chgpath == None :
            logger.info("C moode doens't exist")
            sys.exit(1)

    if not "B" in inputs.calmode :
        for d in os.listdir(calpath) :
            if os.path.isdir("{}/{}".format(calpath,d)):
                mode = d.split("_")[0]
                fn1 = '_'.join(d.split("_")[2:])
                if mode == 'B' and fn == fn1 :
                    bandpath = "{}/{}".format(calpath,d)
        if bandpath == None :
            logger.info("B moode doens't exist")
            sys.exit(1)

    for mt in inputs.calmode :
        os.chdir(calpath)
        runfolder = []

        if mt == 'E' :
            nelect1 = Eigenval("{}/EIGENVAL".format(bandpath)).nelect
            if soc :
                nelect = (int(nelect1),int(nelect1)+1)
            else :
                nelect = (int(nelect1/2),int(nelect1/2)+1)

            bsp = BSPlotting(vasprun=os.path.abspath("{}/vasprun.xml".format(bandpath)), kpoints=os.path.abspath("{}/KPOINTS".format(bandpath)))
            try : 
                vbm = bsp.bs.kpoints[bsp.bsdict['vbm']['kpoint_index'][0]].as_dict()['fcoords']
            except : 
                vbm = [0.000, 0.000, 0.000]
            try :
                cbm = bsp.bs.kpoints[bsp.bsdict['cbm']['kpoint_index'][0]].as_dict()['fcoords']
            except : 
                cbm = [0.000, 0.000, 0.000]
            kpoints = (vbm, cbm)

        p = PerovInputs(structure=struc,is_selective=ds)

        # Revised INPUT FILES
        vi = p.inputfolder(inputs=inputs, method=mt,soc=soc)
        inputs = inputgmd(inputspath)

        try :
            symmetry, groupnumber = struc.get_space_group_info()
        except :
            groupnumber = 0

        # Designate the folder name
        full_formula = GMDStructure(struc).formula(reduced=False)
        pretty_formula = GMDStructure(struc).formula()
            
        if "{0}_{1}_{2:02d}".format(mt,full_formula,groupnumber) in os.listdir(calpath) :
            logging.info("{0}_{1}_{2:02d} is already exists!".format(mt, full_formula, groupnumber))
            sys.exit(1)
        folder_name = "{0}/{1}_{2}_{3:02d}".format(calpath,mt,full_formula,groupnumber) 

        vi.write_input(output_dir=folder_name)
        runfolder.append(folder_name)

        if mt == "E" :
            folder_name_H = "{0}/H_{1}_{2:02d}".format(calpath,full_formula,groupnumber) 
            vi.write_input(output_dir=folder_name_H)
            runfolder.append(folder_name_H)

        # Copy the other files to generated folder
        if mt == "D" or mt=="B" :
            copyfile(chgpath,"{}/CHGCAR".format(folder_name))
        elif mt == "E" :
            copyfile(chgpath,"{}/CHGCAR".format(folder_name))
            copyfile(chgpath,"{}/CHGCAR".format(folder_name_H))
            MakingInpcar(struc,"{}/INPCAR".format(folder_name),nelect[1], kpoints[1])
            MakingInpcar(struc,"{}/INPCAR".format(folder_name_H),nelect[0], kpoints[0])
            logging.info("INPCAR is generated in {} mode ".format(mt))

        for runf in runfolder :
            if mt == 'E' :
                emc = GMDAnalysis() 
                emc.effectivemass(path="{}".format(runf),secondstep=False)
                logging.info("KPOINTS is fixed in {} mode ".format(mt))
            naming = os.path.basename(runf)
            rs = RunningShell(shell=inputs.shell, name=naming, path=runf)
            rs.running_mode(soc=soc, run=True)
            logging.info("{} mode calculation is being started".format(mt))

        # Running Check
        while True :
            time.sleep(10)
            path1 = [] 
            for j in runfolder :
                os.chdir(os.path.join(calpath,j))
                try :
                    vrun = Vasprun("%s/vasprun.xml"%(os.path.join(calpath,j)),parse_potcar_file=True)
                    ionicsteps = vrun.nionic_steps
                    nsw = vrun.incar['NSW']
                    if nsw == 1 :
                        path1.append("%s/CONTCAR"%(os.path.join(calpath,j)))
                    elif nsw == ionicsteps :
                        print("[Notice] Realculation because ionic step is same NSW value")
                        logging.info("{} mode calculation is recalulated (ionic step is over)".format(mt))
                        Recalculate()
                        time.sleep(10)
                    else :
                        path1.append("%s/CONTCAR"%(os.path.join(calpath,j)))
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
                                    print("[Notice] Reculation because it has not yet obtained a stabilizing structure.")
                                    logging.info("{} mode calculation is recalculated (unstbilized tructure)".format(mt))
                                    Recalculate()
                                    time.sleep(10)
                                else :
                                    pass
                            except :
                                pass
                        else :
                            path1.append("%s/CONTCAR"%(os.path.join(calpath,j)))
                    else :
                        pass
                except FileNotFoundError :
                    pass
                except AttributeError :
                    pass

            if len(runfolder) == len(path1) :
                os.chdir(calpath)
                logger.info("{} mode is calculated".format(mt))
                break

        # Properties for DOS and effective mass
        if mt == "C" :
            chgpath = CopyCHGCAR(path1[0])
        elif mt == 'B' :
            bandpath = os.path.split(path1[0])[0]
            bsp = BSPlotting(vasprun=os.path.abspath("{}/vasprun.xml".format(bandpath)), kpoints=os.path.abspath("{}/KPOINTS".format(bandpath)))
            #bsp.get_plot().savefig("{}/BS.pdf".format(calpath))
            bsp.printinform(path=calpath)

        elif mt == "D" or mt == "E" :
            for i in runfolder :
                i = os.path.abspath(i)
                analysis_path = GMDAnalysis()
                if mt == "D" :
                    os.chdir(i)
                    os.system("%s dos width=0.03"%(analysis_path.pdos))
                    analysis_path.partialDOS(structure=struc)
                    logger.info("D mode folder is generated pdos files")
                elif mt == 'E' :
                    analysis_path.effectivemass(path=i, secondstep=True)
                    em_e = open("{}/EM".format(i),'r').readlines()[-12:]
                    E = GMDExcitonbinding.harm_mean_em(em_e)
                    logger.info("E mode folder is generated EM files")
                    del runfolder[-1]
                os.chdir(calpath)
        # CONTCOAR TO POSCAR 
        print("%s mode finished time : "%(mt),end=" ")
        print(time.strftime("%c\n",time.localtime(time.time())))
