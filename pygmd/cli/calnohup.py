import os
import sys
import time

from perovgen.pygmd.autocal.algorithm import Process, AutoMolOpt
from perovgen.pygmd.input_structure import load_structure
from perovgen.pygmd.autocal.inputset import inputgmd

if mole :
    AutoMolOpt(strucpath,inputs)
else :
    Process(calpath=os.path.dirname(__file__),inputspath=inputpath, strucpath=strucpath,ds=ds,soc=orbit)

