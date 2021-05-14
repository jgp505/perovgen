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

def mole(args) :
    # relax option
    if args.relax :
        inputs = read_input(args.relax)
        if "R" in inputs["METHOD"] :
            inputs["INCAR"].append("ISIF=2")
        else :
            print("Please enter the R in input.gmd")
            sys.exit(1)
    
    if not args.sub : 
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
        inputatom, changenum, fixcalc, multiple = InputInform(random_coord=args.position, random_degre=args.degree).molecule_input()
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
                    if args.relax : 
                        strucpath.append(os.path.abspath(vaspname))
                else :
                    pass
                molestruc.to(filename="{0}/POSCAR_times{1:04d}".format(vaspname,i+1))
                if args.csv : 
                    csv.to_csv("{0}/{1}_{2:04d}.csv".format(vaspname,vaspname,i+1))
        if args.relax :
            openingphrase(inputs,strucpath)
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
            os.system("nohup python {}/running_calnohup.py &".format(os.path.dirname(__file__)))
    else :
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
                if args.relax : 
                    strucpath.append(os.path.abspath("POSCAR_{0}_times{1:04d}".format(vaspname,i+1)))
        if args.relax :
            openingphrase(inputs,strucpath)
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
            #os.system("python {}/running_calnohup.py".format(os.path.dirname(__file__)))
            os.system("nohup python {}/running_calnohup.py &".format(os.path.dirname(__file__)))