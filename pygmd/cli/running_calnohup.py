strucpath=['/home/jgp505/Lead-Free/Calculation/FASnBr/FASnBr3']
mole=True
inputs={'KPOINTS': [40], 'INCAR': ['MPJ'], 'SHELL': ['full2.sh'], 'METHOD': ['M']}

import os
import sys
import time

from perovgen.pygmd.autocal.auto_process import Process, AutoMolOpt

if mole :
    AutoMolOpt(strucpath,inputs)
else :
    Process(inputs=inputs, strucpath=path,ds=ds,orbit=orbit)
