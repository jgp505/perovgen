import os
import sys

import numpy as np
import pandas as pd
import json
from collections import defaultdict

from pymatgen.core import Structure,composition
from pymatgen.analysis.eos import EOS
from pymatgen.io.vasp.outputs import Vasprun,BSVasprun,Procar

from perovgen.pygmd.analysis.electronic import BSPlotting
from perovgen.pygmd.autocal.inputset import inputgmd

def GMDDataFrame(calpath, inputpath) :
    dataframe = defaultdict(list)
    inputs = inputgmd(inputpath).calmode
    id1, name = os.path.absname(calpath).split("_")[-2:]
    dataframe['material_id'].append(id1)
    dataframe['composition'].append(name)
   
    # bulkmodulus
    calculation = dict()
    if "bulkmodulus" in [*inputs]:
       	calculation['bulkmodulus'] = True