import os
import sys

from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Vasprun

from perovgen.pygmd.input_structure import GMDStructure, load_structure
from perovgen.pygmd.shell import ShellPath
from perovgen.pygmd.autocal.substitute import RandomMolecule, InputInform, RandomAtom
from perovgen.pygmd.autocal.algorithm import controlincar, openingphrase
from perovgen.pygmd.autocal.inputset import *

import time
import subprocess

def molecule(args) :
    if args.ma :
        mafa = "MA"
    elif args.fa :
        mafa = "FA"
    elif args.gua :
        mafa = "GUA"
    elif args.dima :
        mafa = "diMA"
    elif args.trima :
        mafa = "triMA"
    elif args.tetrama :
        mafa = "tetraMA"
    elif args.zolium :
        mafa = "Zolium"
    ms = RandomMolecule._loadmolecule(mafa)
    inputatom, changenum, fixcalc, multiple = InputInform(random_coord=args.position, random_degree=args.degree).molecule_input()
    mole = RandomMolecule(random_coord=args.position, random_degree=args.degree)
    # Input  
    strucpath = []
    structure = load_structure(args.path)[0]
    e = 0
    for s in structure :
        for i in range(multiple) : 
            molestruc, csv = mole.tiltingmolecule(s, ms, inputatom=inputatom, changenum=changenum,fixcalc = fixcalc)
            vaspname = GMDStructure(molestruc).formula(reduced=False)
            if not os.path.exists(vaspname):
                os.makedirs(vaspname)
            molestruc.to(filename="{0}/POSCAR_{1}_{2:04d}".format(vaspname,vaspname,e+1))
            if args.csv :
                if not os.path.exists("{}/CSV".format(vaspname)):
                    os.makedirs("{}/CSV".format(vaspname))
                csv.to_csv("{0}/CSV/{1}_{2:04d}.csv".format(vaspname,vaspname,e+1))
            e+=1
    return strucpath

def element(args):
    r_atom = RandomAtom(random_coord=args.position, random_degree=args.degree)
    atom1, atom2, change, multiple = InputInform(random_coord=args.position, random_degree=args.degree).atom_input()
    #Input 
    strucpath = []
    structure = load_structure(args.path)[0]
    e = 0
    for s in structure :
        for i in range(multiple) :
            s1, s = r_atom.substitution(structure=s,atom1=atom1, atom2=atom2, ratio=change)
            vaspname = GMDStructure(s1).formula(reduced=False)
            if not os.path.exists(vaspname) :
                os.makedirs(vaspname)
            s1.to(filename="{0}/POSCAR_{1}_{2:04d}".format(vaspname,vaspname,e+1))
            e+=1
            #if args.input : 
            #    strucpath.append(os.path.abspath("POSCAR_{0}_times{1:04d}".format(vaspname,i+1)))
    return strucpath

def randomreplace(args) :
    if args.sub :
        strucpath =element(args)
    else :
        strucpath=molecule(args)