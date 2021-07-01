import os
import sys

from perovgen.pygmd.analysis.electronic import BSPlotting, _load_yaml

# load the path of vasprun.xml and KPOINTS
bsp = BSPlotting(vasprun='Cs3Bi2I9_B/vasprun.xml', kpoints='Cs3Bi2I9_B/KPOINTS')

# write the band_inform.log 
bsp.printinform()

# load the graph.yaml and plot BS
plt = bsp.get_plot()
plt.show()