import os
import sys

from perovgen.pygmd.analysis.energy import GMDPhaseDiagram

# make the GMDPhaseDiagram
base = GMDPhaseDiagram(vasprunpath=['Cs3Bi2I9_C/vasprun.xml'],extractcomp=None)

# make the phase diagram and entries for hull energy and formation energy
entries, pd = base.get_phasediagram()
df = GMDPhaseDiagram.get_data(entries, pd)
# save the csv file
df.to_csv("PhaseDiagram.csv")