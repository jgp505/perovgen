inputs={'KPOINTS': [40], 'INCAR': ['MPJ'], 'SHELL': ['vasp_36.sh'], 'METHOD': ['E', 'U']}
path=['Cs20In20Br60_1_B_mode/CONTCAR']
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
