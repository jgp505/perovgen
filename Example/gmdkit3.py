#!/usr/local/anaconda3-2019.10/bin/python
# coding : utf-8
# Pymtagen Citation

import os
import sys
import time

import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import palettable
import yaml
from shutil import copyfile
from tabulate import tabulate

from pymatgen.io.vasp.outputs import Vasprun, BSVasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter

from perovgen.pygmd.base import graphyaml, _inputgmd, load_structure,ShellPath,GMDStructure,createFolder
from perovgen.pygmd.auto_calculation.crystal_data_base.cdb import *
from perovgen.pygmd.auto_calculation import auto_process
from perovgen.pygmd.auto_calculation.molecular_cation import molecule
from perovgen.pygmd.properties_analysis.electronic_properties import _load_yaml, DOSPlotting, BSPlotting

def opening_phrase(dictionary):
    text = ''
    for k,v in dictionary.items() :
        text+="%s : %s\n"%(k,v)
    table = [[text]]
    output = tabulate(table, tablefmt='grid')
    print(output)

def yesorno(phrase):
    while True :
        i = str(input(phrase))
        if i == "Y" :
            return True
        elif i == "N" :
            return False
        else :
            print("Please enter the Y or N")
    
print("+"+("-"*63)+"+")
print("|\t\tGMDKIT Version: 3.5.1 (30 Jun, 2021)\t\t|")
print("|\t\t %s \t\t\t|"%(time.strftime("%c", time.localtime(time.time()))))
print("+"+("-"*63)+"+")

# Select the argument 
while True :
    mode_dic = {"AUTO" : "auto calculation",
    "Replace" : "the inorganic material is substituted with an organic material(e.g. MA(Methylammonium),FA(Formamidinium))",
    "PLOT" : "Plotting for electronic structure(e.g. Band Structure, DOS)",
    "FILE" : "Structure file for VASP calculation (e.g. POSCAR, CONTCAR) and file for electronic structure plotting",
    "CONFIG" : "Manage shell script for VASP calculation",
    "0" : "quit program"}
    opening_phrase(mode_dic)
    mode = str(input("Please enter the mode >> "))
    if not mode in [*mode_dic] :
        print("{} doesn't exist".format(mode))

    # AUTO MODE    
    if mode == 'AUTO' :
        while True :
            # read inputfile
            if not "input.gmd" in os.listdir("."):
                path = str(input("Please enter the input.gmd path >> "))
                path = os.path.abspath(path)
                inputs = read_input("{}/input.gmd".format(path))
            else :
                break
        inputs = read_input("input.gmd")

        if not inputs['STRUC'] :
            strucpath = load_structure(os.getcwd(),sformat=False)
            if len(strucpath) == 0 :
                path = str(input("Please enter the structural path >> "))
                print(os.path.abspath(path))
                strucpath = load_structure(os.path.abspath(path),sformat=False)
        else :
            strucpath = inputs['STRUC']
        ds = yesorno("Delete selective dynamics? (Y or N) >> ")
        soc = yesorno("SOC mode (Y or N)>> ")

        # write the running_calnohup for autocalculation
        autopath = "{}/pygmd/cli".format(os.path.dirname(ShellPath().shellpath))
        auto_process.openingphrase(inputs["METHOD"],strucpath,inputs["SHELL"][0])
        f = open("{}/calnohup.py".format(autopath),"r")
        ff = f.readlines()
        ff.insert(0,"inputs={}\n".format(dict(inputs)))
        ff.insert(1,"path={}\n".format(str(strucpath)))
        if ds :
            ff.insert(2,"ds=True\n")
        else :
            ff.insert(2,"ds=False\n")
        if soc :
            ff.insert(3,"orbit=True\n")
        else :
            ff.insert(3,"orbit=False\n")
            ff.insert(4,"\n")
            with open("{}/running_calnohup.py".format(autopath),"w") as fi :
                for i in ff :
                    fi.write(i)
            fi.close()
        #os.system("python {}/running_calnohup.py".format(os.path.dirname(__file__)))
        os.system("nohup python {}/running_calnohup.py &".format(autopath))
        break

    # MOLE MODE
    elif mode == "MOLE" :
        # Select the molecule 
        while True :
            mole_dic = {'MA':'Methylammonium','FA':'Formamidinium','GUA':'Guanidium','diMA':'dimethylammonium',
                    'triMA':'trimethylammonium','tetraMA':'tetramethylammonium','Zolium':'Imidazolium'}
            opening_phrase(mole_dic)
            molemode = str(input("Please enter the mode >> "))
            if not molemode in [*mole_dic]:
                print("{} doesn't exist.".format(molemode))
            else :
                break
        rp = yesorno("Random Position (Y or N) >> ")
        rd = yesorno("Random Degree (Y or N) >> ")
        print(rp)

        if not rp :
            if not 'MOLE' in os.listdir(".") :
                print("MOLE file required")
                sys.exit(0)
        else :
            pass

        ms = molecule.RandomPosition._loadmolecule(molemode)
        pwd = os.getcwd()
        cur = pwd.split(os.path.sep)[-1]

        l=load_structure(pwd,sformat=True)
        if not l :
            path = str(input("Please enter the structural path >> "))
            l = load_structure(os.path.abspath(path))
                    
        mole = molecule.RandomPosition(random_coord=rp,random_degree=rd)
        inputatom, changenum, fixcalc, multiple = mole._inputinform()

        if fixcalc == None :
            cur1 = "%s_ration%.2f%%_multiple%i_%s"%(cur,changenum*100,multiple,molemode)
        else :
            cur1 = "%s_ration%.2f%%_multiple%i_selective%i_%s"%(cur,(changenum*100),multiple,fixcalc,molemode)

        for s in l :
            vaspname = GMDStructure(s).vaspname()
            for i in range(multiple) :
                m = mole.tiltingmolecule(s, ms, inputatom=inputatom, changenum=changenum,
                fixcalc = fixcalc)
                m[0].to(filename="POSCAR_%s_multiple%i"%(vaspname,i+1))
                m[1].to_csv("%s_%i.csv"%(vaspname,i+1))
        createFolder(os.path.join(cur1))
        os.system("mv %s/POSCAR_*_multiple* %s/%s"%(pwd,pwd,cur1))
        os.system("mv %s/*.csv %s/%s"%(pwd,pwd,cur1))
        print("Genreate %s/%s"%(pwd,cur1))
        break

    # PLOT MODE
    elif mode == 'PLOT' :
        while True :
            plot_dic={"B":"Plotting band structure","D":"Plotting DOS","BD":"Plotting band struture & DOS",
            "P":"Plotting partial DOS","C":"check the band information"}
            opening_phrase(plot_dic)
            plotmode = str(input("Please enter the mode >> "))
            if not plotmode in [*plot_dic]:
                print("{} doesn't exist.".format(plotmode))
            else :
                outfile1 = yesorno("Do you want to save the file (Y or N) >> ")
                if outfile1 : 
                    filename = str(input("Please enter the name saved >> "))
                break

        if plotmode == "B" or plotmode == "C" :
            bsp = BSPlotting(vasprun="vasprun.xml")
            if plotmode == "C" :
                bsp.band_inform()
            else :
                para = _load_yaml("B")
                plt = bsp.get_plot(figsize=para["fig_size"],zero_to_efermi=para["zero_to_efermi"],
			    fontsize=para["fontsize"],ylim=para["ylim"],color=para["color"],
			    vbm_cbm_marker=para["vbm_cbm_marker"],linewidth=para["linewidth"])

        elif plotmode == "D" or plotmode == "P" :
            para = _load_yaml("D")
            # check the dos file
            if plotmode == "D" :
                if not "dos" in os.listdir(".") :
                    path = str(input("Please enter the dos file path >> "))
                    path = os.path.abspath(path)
                    path = "{}/dos".format(path)
                else :
                    path = "{}/dos"
                    f = open(path,'r')
                    filelist = f.readlines()[1:]
                    dp = DOSPlotting(vasprun=path,dos=filelist, zero_to_efermi=para["zero_to_efermi"],stack=para["stack"]) 
                    plt = dp.get_plot(figsize=para["fig_size"],xlim=para["xlim"],ylim=para["ylim"],
                                    fontsize=para["font_size"],color=para["color"])
            else :
                path =[p for p in os.listdir(".") if ".dos" in p]
                filelist = [];ylim=[]
                color = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
                for e,i in enumerate(path) :
                    f = open(i,"r")
                    filelist = f.readlines()[1:]
                    dp = DOSPlotting(vasprun=path,dos=filelist, zero_to_efermi=para["zero_to_efermi"],stack=para["stack"]) 
                    plt = dp.get_plot(figsize=para["fig_size"],xlim=para["xlim"],ylim=para["ylim"],
                        fontsize=para["font_size"],color=color[e%9],label=i.split(".")[0])
                    ax = plt.gca()
                    ylim.append(ax.get_ylim()[-1])
                    ax.set_ylim(0,max(ylim))
                    plt.legend(fontsize=para["font_size"]/2, loc="upper right")
        elif plotmode == "BD":
            run = Vasprun(path,parse_dos=True)
            dos = run.complete_dos
            vrun = BSVasprun(path,parse_projected_eigen=True)
            bs =vrun.get_band_structure('KPOINTS',efermi=dos.efermi)

            bdpara=_load_yaml("BD")
            bsdosplot = BSDOSPlotter(bs_projection=bdpara['bs_projection'],
            dos_projection=bdpara['dos_projection'],
			vb_energy_range=bdpara['vb_energy_range'],
			cb_energy_range=bdpara['cb_energy_range'],
			fixed_cb_energy=bdpara['fixed_cb_energy'],
			egrid_interval=bdpara['egrid_interval'], 
			font=bdpara['font'], 
			axis_fontsize=bdpara['axis_fontsize'],
			tick_fontsize=bdpara['tick_fontsize'], 
			legend_fontsize=bdpara['legend_fontsize'],
			bs_legend=bdpara['bs_legend'],
			dos_legend=bdpara['dos_legend'],
			rgb_legend=bdpara['rgb_legend'], 
			fig_size=bdpara['fig_size'])
            plt = bsdosplot.get_plot(bs, dos=dos)

        if plt :
            if outfile1 : # 수정
                plt.savefig("%s.pdf"%(filename))
                break
            else :
                plt.show()
                break
        break

    elif mode == "FILE" :
        while True :
            file_dic={"G":"garaph.ymal","S":"substitution",
            "I":"input.gmd"}
            opening_phrase(file_dic)
            filemode = str(input("Please enter the mode >> "))
            if not filemode in [*file_dic]:
                print("{} doesn't exist.".format(filemode))
            else :
                if filemode == "G" :
                    while True :
                        fileg = str(input("Please enter the style of graph.yaml (B, D, BD) >> "))
                        if not fileg in ["B", "D", "BD"] :
                            print("{} doesn't exist!".format(fileg))
                        else :
                            break
                    data = graphyaml(fileg)
                    data["name"]=fileg
                    with open('graph.yaml','w') as out :
                        yaml.dump(data,out,default_flow_style=False)
                    print("%s generate the graph.yaml file"%(os.getcwd()))
                    break
                elif filemode == "S" :
                    # get path 
                    path = os.getcwd()
                    strucpath = load_structure(path)
                    if not strucpath :
                        path = str(input("Please enter the path >> "))
                        path = os.path.abspath(path)
                        strucpath = load_structure(path)

                    print("FILE MODE - S")
                    inputatom = input("Enter the atom to be substituted >> ")
                    substitute = input("Enter the atom to substitute >> ")
                    change = float(input("Enter the ratio >> "))
                    if change > 1 :
                        change = 1.0
                    if change == 1 :
                        multiple = 1 
                    else :
                        multiple = int(input("Enter the number of the times to repeat >> "))
                    name = os.path.split(os.path.abspath(path))[-1]
                    directory = "%s_ratio_%.2f%%_to_%s_from_%s"%(name,(change*100),inputatom,substitute)
                    createFolder(directory)

                    for i in range(multiple) :
                        for s in strucpath :
                            gs = GMDStructure(s)
                            vs = gs.vaspname()
                            struc = gs.substitution(atom1=inputatom, atom2=substitute, ratio=change)
                            struc.sort()
                            filename="POSCAR_%s_multiple_%i_from_%s"%(GMDStructure(struc).vaspname(), i+1, vs)
                            struc.to(filename="%s/%s"%(directory, filename))

                    print("generate the %s/%s"%(os.getcwd(),directory))
                    break

                elif filemode == "I" :
                    # get path 
                    path = os.getcwd()
                    strucpath = load_structure(path,sformat=False)
                    if not strucpath :
                        path = str(input("Please enter the path >> "))
                        path = os.path.abspath(path)
                        strucpath = load_structure(path,sformat=False)
                    _inputgmd(strucpath)

                    print("Generate {}/input.gmd".format(os.getcwd()))
                    break

        break
    elif mode == "CONFIG" :
        while True :
            config_dic={"S":"shell","C":"check","R":"remove"}
            opening_phrase(config_dic)
            configmode = str(input("Please enter the mode >> "))
            if not configmode in [*config_dic]:
                print("{} doesn't exist.".format(configmode))
            else :
                sp = ShellPath()
                if configmode == "S" :
                    while True :
                        path = str(input("Please enter the path >> "))
                        path = os.path.abspath(path)
                        if os.path.basename(path) in os.listdir(os.path.dirname(path)) :
                            break
                        else :
                            print("{} doens't exist file".format(path))
                    sp.register_shell(path)
                    break
                elif configmode == "C" :
                    sp.check()
                    break
                elif configmode == "R" :
                    sp.remove()
                    break
        break
    elif mode == '0' :
        break
        