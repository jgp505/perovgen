inputs={'KPOINTS': [25], 'INCAR': ['PBE'], 'METHOD': ['G'], 'SHELL': ['vasp.sh']}
path='CsSnI3_221/HSE/B_CsSnI3_01/CONTCAR'
ds=False
orbit=False
mole=False
import os
import sys
import time

from perovgen.pygmd.autocal.algorithm import Process, AutoMolOpt
from perovgen.pygmd.input_structure import load_structure

if mole :
    AutoMolOpt(strucpath,inputs)
else :
    struclist = load_structure([path]) 
    if not struclist :
        print("\n[Warning] the structure file does not exist.\n")
        sys.exit(0)
    Process(inputs=inputs, strucpath=[path],ds=ds,orbit=orbit)

