import os
import sys

import palettable
from perovgen.pygmd.base import load_structure
from perovgen.pygmd.analysis.electronic import DOSPlotting, GMDAnalysis

color = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

# Make the LIST and load the .dos file

s = load_structure(".")[0]
path = GMDAnalysis().partialDOS(structure=s)
name =[f for f in os.listdir(".") if ".dos" in f]

ylim=[]
for e,i in enumerate(name) :
	f = open(i,'r')
	filelist = f.readlines()[1:]
	dp = DOSPlotting(vasprun="vasprun.xml",dos=filelist)
	plt = dp.get_plot(color=color[e%9],label=i.split(".")[0])
	ax = plt.gca()
	ylim.append(ax.get_ylim()[-1])
	ax.set_ylim(0,max(ylim))
	plt.legend(fontsize=15)
plt.show()