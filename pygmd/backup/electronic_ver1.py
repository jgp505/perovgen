import os
import sys

import numpy as np
import pandas as pd
import scipy.interpolate as sci
import yaml
import palettable
import matplotlib.pyplot as plt

from perovgen.pygmd.base import graphyaml

from pymatgen.io.vasp.outputs import Vasprun, BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter, BSDOSPlotter
from pymatgen.electronic_structure.core import Spin, OrbitalType

def _load_yaml(paraname):
	pwd = os.getcwd()
	graph = "{}/graph.yaml".format(pwd)
	if not os.path.isfile(graph):
		loading = graphyaml(paraname)
	else :
		stream = open(graph,"r")
		loading = yaml.load(stream)
	return loading

class DOSPlotting :
	def __init__(self, vasprun="vasprun.xml",dos=None, zero_to_efermi=True, stack=True) :
		self.zero_to_efermi = zero_to_efermi 
		self.stack = stack
		self.efermi = Vasprun("vasprun.xml").efermi
		self.energies = [float(i.split()[0])-self.efermi if self.zero_to_efermi else float(i.split()[0]) for i in dos]
		self.densities = [float(i.split()[1]) for i in dos]
	
	def get_plot(self, figsize=(12,8), xlim=None, ylim=None, fontsize=30,color="r",label="Total DOS") :
		plt.rcParams['font.family']='Arial'
		plt.rcParams['figure.figsize']=figsize
		plt.rcParams['font.size']=fontsize
		plt.rcParams['lines.linewidth'] = 3
		plt.rcParams['axes.linewidth'] = 3

		allpts = []
		allpts.extend(list(zip(self.energies, self.densities)))

		if xlim :
			plt.xlim(xlim)
		if ylim :
			plt.ylim(ylim)
		else :
			xlim = plt.xlim()
			relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
		plt.ylim((min(relevanty),max(relevanty)))

		if self.stack :
			plt.fill_between(self.energies,0,self.densities,color=color,alpha=.1)
		xlabel = "Energies(eV)"
		if self.zero_to_efermi :
			xlabel = "Energies-E$_{F}$(eV)"
			ylim = plt.ylim()
			plt.plot([0,0], ylim, 'k--', linewidth=2)
		
		plt.xlabel(xlabel)
		plt.ylabel("Density of States")

		plt.plot(self.energies, self.densities, color=color, label=label)
		plt.tight_layout()
		return plt

class BSPlotting :
	def __init__(self, vasprun="vasprun.xml"):
		self.vrun = Vasprun(vasprun,parse_dos=True)
		self.bsrun = BSVasprun(vasprun,parse_projected_eigen=True)
		self.bs = self.bsrun.get_band_structure("KPOINTS",efermi=self.vrun.efermi)
		self.bsp = BSPlotter(self.bs)
		self.data =self.bsp.bs_plot_data()

	def band_inform(self):
		ad = self.bsrun.as_dict()['output']

		bandgap = "%.3f(Indirect)"%(ad['bandgap'])
		if ad['is_gap_direct'] :
			bandgap = "%.3f(Direct)"%(ad['bandgap'])

		print("\nnumber of bands : {}".format(self.bs.nb_bands))
		print("number of kpoints :", len(self.bs.kpoints))
		print("fermi energy : %.3f"%(self.bs.efermi))
		print("band gap : {}\n".format(bandgap))

	def get_plot(self,figsize=(12,8),zero_to_efermi=True,fontsize=20,ylim=(-4,6),color='b',vbm_cbm_marker=True,linewidth=1):
		import scipy.interpolate as scint
		from pymatgen.util.plotting import pretty_plot

		ad = self.bsrun.as_dict()['output']
		bandgap = "%.3f(Indirect)"%(ad['bandgap'])
		if ad['is_gap_direct'] :
			bandgap = "%.3f(Direct)"%(ad['bandgap'])
		print("band gap : {}".format(bandgap))

		plt = pretty_plot(figsize[0],figsize[1])

		plt.rcParams['lines.linewidth'] = linewidth
		plt.rcParams['font.size'] = fontsize

		for d in range(len(self.data['distances'])):
			for i in range(self.bs.nb_bands):
				plt.plot(self.data['distances'][d],
					[self.data['energy'][d][str(Spin.up)][i][j] 
					for j in range(len(self.data['distances'][d]))],
					color=color, ls='-')
				if self.bsp._bs.is_spin_polarized :
					plt.plot(self.data['distances'[d],
						[self.data['energy'][d][str(Spin.up)][i][j] 
						for j in range(len(self.data['distances'][d]))]],
						color='r', ls='-')
		self.bsp._maketicks(plt)
		plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
		ylabel = r'$\mathrm{Energy\ (eV)}$'
		if zero_to_efermi :
			ylabel = r'$\mathrm{E\ -\ E_{VBM}\ (eV)}$'
		plt.ylabel(ylabel)
		plt.xlim(0, self.data['distances'][-1][-1])


		if ylim is None :
			emin = -10 ; emax=10
			if self.bsp._bs.is_metal():
				if zero_to_efermi:
					plt.ylim(emin, emax)
				else :
					plt.ylim(self.vrun.efermi+emin, self.vrun.efermi+emax)
			else :
				if vbm_cbm_marker :
					for cbm in self.data['cbm']:
						plt.scatter(cbm[0], cbm[1], color='r', marker='o', s=100)
					for vbm in self.data['vbm']:
						plt.scatter(vbm[0], vbm[1], color='g', marker='o', s=100)
				plt.ylim(self.data['vbm'][0][1] + emin, self.data['cbm'][0][1]+emax)
		else :
			plt.ylim(ylim)
			if not self.bsp._bs.is_metal() and vbm_cbm_marker:
				for cbm in self.data['cbm']:
					plt.scatter(cbm[0], cbm[1], color='r', marker='o',s=100)
				for vbm in self.data['vbm']:
					plt.scatter(vbm[0], vbm[1], color='g', marker='o',s=100)
			plt.tight_layout()

		if not zero_to_efermi :
			ef = self.vrun.efermi
			plt.axhline(ef, linewidth=linewidth, color='k', ls="--")
		else :
			ax = plt.gca()
			xlim = ax.get_xlim()
			ax.hlines(0, xlim[0], xlim[1], linestyles="dashed", color='k')

		if self.bsp._bs.is_spin_polarized :
			ax.plot((),(),color=color,ls='-',label="spin up")
			ax.plot((),(),color='r',ls='-',label="spin dwon")
			ax.legend(loc="upper left")
		return plt
