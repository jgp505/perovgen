import os
import sys

from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Vasprun

from perovgen.pygmd.base import GMDStructure, load_structure, ShellPath
from perovgen.pygmd.autocal.substitute import RandomMolecule, InputInform, RandomAtom
from perovgen.pygmd.autocal.auto_process import controlincar, openingphrase
from perovgen.pygmd.autocal.cdb import *

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
    structure = load_structure(args.path)
    for s in structure :
        for i in range(multiple) : 
            molestruc, csv = mole.tiltingmolecule(s, ms, inputatom=inputatom, changenum=changenum,fixcalc = fixcalc)
            vaspname = GMDStructure(molestruc).vaspname()
            if not os.path.exists(vaspname):
                os.makedirs(vaspname)
                if args.input : 
                    strucpath.append(os.path.abspath(vaspname))
            molestruc.to(filename="{0}/POSCAR_times{1:04d}".format(vaspname,i+1))
            if args.csv :
                if not os.path.exists("{}/CSV".format(vaspname)):
                    os.makedirs("{}/CSV".format(vaspname))
                csv.to_csv("{0}/CSV/{1}_{2:04d}.csv".format(vaspname,vaspname,i+1))
    return strucpath

def element(args):
    r_atom = RandomAtom(random_coord=args.position, random_degree=args.degree)
    atom1, atom2, change, multiple = InputInform(random_coord=args.position, random_degree=args.degree).atom_input()
    #Input 
    strucpath = []
    structure = load_structure(args.path)
    for s in structure : 
        for i in range(multiple) :
            s1 = r_atom.substitution(structure=s,atom1=atom1, atom2=atom2, ratio=change)
            vaspname = GMDStructure(s1).vaspname()
            s1.to(filename="POSCAR_{0}_times{1:04d}".format(vaspname,i+1))
            if args.input : 
                strucpath.append(os.path.abspath("POSCAR_{0}_times{1:04d}".format(vaspname,i+1)))
    return strucpath

def randomreplace(args) :
    if args.sub :
        if args.input :
            inputs = read_input(args.input)
        strucpath = element(args)
    else :
        if args.input :
            inputs = read_input(args.input)
            if not inputs["METHOD"] == ["M"] :
                print("Please enter the M in input.gmd")
                sys.exit(1)
        strucpath=molecule(args)

    if args.input :
        openingphrase(inputs, strucpath)
        if args.sub :
            f = open("{}/calnohup.py".format(os.path.dirname(__file__)),"r")
            ff = f.readlines()
            ff.insert(0,"inputs={}\n".format(dict(inputs)))
            ff.insert(1,"path={}\n".format(str(strucpath)))
            ff.insert(2,"ds=False\n")
            ff.insert(3,"orbit=False\n")
            ff.insert(4,"mole=False")
            ff.insert(5,"\n")
            with open("{}/running_calnohup.py".format(os.path.dirname(__file__)),"w") as fi :
                for i in ff :
                    fi.write(i)
            fi.close()
        else :
            f = open("{}/calnohup.py".format(os.path.dirname(__file__)),"r")
            ff = f.readlines()
            ff.insert(0,"strucpath={}\n".format(str(strucpath)))
            ff.insert(1,"mole=True\n")
            ff.insert(2,"inputs={}\n".format(dict(inputs)))
            ff.insert(3,"\n")
            with open("{}/running_calnohup.py".format(os.path.dirname(__file__)),"w") as fi :
                for i in ff :
                    fi.write(i)
            fi.close()
        #os.system("python {}/running_calnohup.py".format(os.path.dirname(__file__)))
        os.system("nohup python {}/running_calnohup.py > output.gmd &".format(os.path.dirname(__file__)))
    else :
        print("Done")