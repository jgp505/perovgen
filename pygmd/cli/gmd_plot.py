import os
import sys

import palettable

from pymatgen.io.vasp.outputs import Vasprun, BSVasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter

from perovgen.pygmd.analysis.electronic import _load_yaml, BSPlotting, DOSPlotting

def plot(args):
	plt=None
	if args.band or args.bandcheck :
		bsp = BSPlotting(vasprun='vasprun.xml',kpoints='KPOINTS') #line the path of vasprun.xml using the argument
		if args.bandcheck :
			bsp.printinform()
		elif args.band :
			para = _load_yaml("B")
			plt = bsp.get_plot(figsize=para["fig_size"],zero_to_efermi=para["zero_to_efermi"],
			fontsize=para["fontsize"],spindownoff=True,color=para["color"],ylim=para["ylim"],
			vbm_cbm_marker=para["vbm_cbm_marker"])

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
					dp = DOSPlotting(vasprun='vasprun.xml',dos=filelist, zero_to_efermi=para["zero_to_efermi"],stack=para["stack"]) 
					plt = dp.get_plot(figsize=para["fig_size"],xlim=para["xlim"],ylim=para["ylim"],
						fontsize=para["font_size"],color=color[e%9],label=i.split(".")[0])
					ax = plt.gca()
					ylim.append(ax.get_ylim()[-1])
				ax.set_ylim(0,max(ylim))
				plt.legend(fontsize=para["font_size"]/2, loc="upper right",frameon=False)
			else :
				print("Please use the -n argument")
				sys.exit(0)

			filelist=f.readlines()[1:]


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

	if plt :
		if args.out_file :
			plt.savefig("%s.%s"%(args.out_file,args.format))
		else :
			plt.show()
