import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from perovgen.pygmd.analysis.electronic import _load_yaml, DOSPlotting
import palettable

para = _load_yaml("D")
color = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

ylim = []
for e,i in enumerate(name) :
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
plt.show()