inputs={'KPOINTS': [40], 'INCAR': ['PBE'], 'SHELL': ['full2.sh'], 'METHOD': ['R', 'C', 'D', 'B', 'E']}
path=['../EntryWithCollCode26058 (1).cif', '../EntryWithCollCode39823 (1).cif', '../EntryWithCollCode431322 (1).cif']
ds=False
orbit=False
mole=False
import os
import sys
import time

from perovgen.pygmd.autocal.auto_process import Process, AutoMolOpt

if mole :
    AutoMolOpt(strucpath,inputs)
else :
    Process(inputs=inputs, strucpath=path,ds=ds,orbit=orbit)
