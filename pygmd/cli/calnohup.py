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

