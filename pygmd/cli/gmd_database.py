import os
import sys

from pymatgen.core import Structure
import numpy as np
import pandas as pd

from shutil import copy, move
from collections import defaultdict

from perovgen.pygmd.input_structure import load_structure,GMDStructure

input_abspath = '/home/jgp505/Database/01_input'
def database(args) :
    if args.inputcif :
        files = [f for f in os.listdir(".") if "cif" in f]
        for s in files :
            struc = Structure.from_file(s)
            for i in range(struc.num_sites) :
                try :
                    struc.replace(i, struc.species[i].element)
                except :
                    struc.replace(i, struc.species[i])
            is1 = GMDStructure(struc).savefile()