import os
import sys

import pandas as pd
import numpy as np

import palettable
from collections import defaultdict

from pymatgen.io.vasp.outputs import Vasprun, BSVasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter

from perovgen.pygmd.analysis.electronic import _load_yaml, BSPlotting, DOSPlotting
from perovgen.pygmd.analysis.energy import GMDAbsorption

def plot(args):
    plt=None
    if args.band or args.bandcheck :
        bsp = BSPlotting(vasprun='vasprun.xml',kpoints='KPOINTS') #line the path of vasprun.xml using the argument
        if args.bandcheck :
            dic = bsp._bandinform()
            for k,v in dic.items() :
                if k == 'energies' : 
                    pass
                elif k == 'distances' :
                    pass
                elif k == 'E_g' or k == 'VBM' or k == 'CBM':
                    print(k)
                    for k1, v1 in dic['E_g'].items() :
                        print('    ',k1,":",v1)
                else:
                    print(k,":",v)
     
        elif args.band :
            para = _load_yaml("B")
            plt = bsp.get_plot(figsize=para["fig_size"],zero_to_efermi=para["zero_to_efermi"],
            fontsize=para["fontsize"],spindownoff=True,color=para["color"],ylim=para["ylim"],vbm_cbm_marker=para["vbm_cbm_marker"])

    elif args.dos or args.partial :
        para = _load_yaml("D")
        if args.dos :
            f = open(args.name,"r")
            filelist = f.readlines()[1:]
            dp = DOSPlotting(vasprun='vasprun.xml',dos=filelist, zero_to_efermi=para["zero_to_efermi"],stack=para["stack"]) 
            plt = dp.get_plot(figsize=para["fig_size"],xlim=para["xlim"],ylim=para["ylim"],
                fontsize=para["font_size"],color=para["color"])
            plt.legend(frameon=False)
        if args.partial :
            if args.name :
                filelist = [];ylim=[]
                color = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
                for e,i in enumerate(args.name) :
                    f = open(i,"r")
                    filelist = f.readlines()[1:]
                    dp = DOSPlotting(vasprun='vasprun.xml',dos=filelist, 
                                     zero_to_efermi=para["zero_to_efermi"],stack=para["stack"]) 
                    plt = dp.get_plot(figsize=para["fig_size"],xlim=para["xlim"],ylim=para["ylim"],
                        fontsize=para["font_size"],color=color[e%9],label=i.split(".")[0])
                    ax = plt.gca()
                    ylim.append(ax.get_ylim()[-1])
                    ax.set_ylim(0,max(ylim))
                    plt.legend(fontsize=para["font_size"]/2, loc="upper right",frameon=False)

                # writing script
                files = open("%s/pdos.py"%(os.path.dirname(__file__)),'r').readlines()
                files.insert(0,"name=%s\n"%(str(args.name)))
                with open("pdos.py",'w') as fi :
                    for f in files :
                        fi.write(f)
                fi.close() 
            else :
                print("Please use the -n argument")
                sys.exit(0)

    elif args.bdos :
        run = Vasprun('vasprun.xml',parse_dos=True)
        dos = run.complete_dos
        vrun = BSVasprun('vasprun.xml',parse_projected_eigen=True)
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
  
    elif args.bulk :
        if len(args.name) == 1 :
            bulkdata = pd.read_csv(args.name[0],sep=' ',header=None)
            import matplotlib.pyplot as plt
            from pymatgen.analysis.eos import EOS
            plt.rcParams['font.size'] = 14
            try :
                eos = EOS(eos_name='birch_murnaghan')
                eos_fit = eos.fit(bulkdata[0],bulkdata[1])
                 
                plt.figure(figsize=(8,6))
                plt.plot(bulkdata[0],bulkdata[1],'o',color='r')
                plt.plot(bulkdata[0],eos_fit.func(bulkdata[0]),'--',color='r')
                plt.xlabel("Volume (Å$^{3}$)")
                plt.ylabel("Energy (eV)")
            except :
                plt.figure(figsize=(8,6))
                plt.plot(bulkdata[0],bulkdata[1],'x-',color='r')
                plt.xlabel("Volume (Å$^{3}$)")
                plt.ylabel("Energy (eV)")
        else :
        	print("Please enter the one argument")
   
    elif args.absorp :
        if len(args.name) == 1 :
            apara=_load_yaml("A")
            print(apara)
            plt = GMDAbsorption.plotting(path=os.path.abspath(args.name[0]),fig_size = apara['fig_size'],
                    fontsize = apara['fontsize'],
                    xlim = apara['xlim'],
                    ylim = apara['ylim'],
                    sharex = apara['sharex'])
            
    if plt :
        if args.out_file :
            plt.savefig("%s.%s"%(args.out_file,args.format))
        else :
            plt.show()
