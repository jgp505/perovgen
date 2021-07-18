# coding : utf-8
# Cppyright (c) Green Materials Designs Team.

import os
import sys

import yaml
from shutil import move
from perovgen.pygmd.base import graphyaml, load_structure, createFolder, GMDStructure,_inputgmd,pdos,MPJClass
from collections import defaultdict

def graph_yaml(args):
    '''
    Generate graph.yaml
    '''
    data = graphyaml(args.graph)
    data["name"]=args.graph
    with open('graph.yaml','w') as out :
        yaml.dump(data,out,default_flow_style=False)
    print("%s generate the graph.yaml file"%(os.getcwd()))

def strain(args):
    '''
    Transform the structure
    '''
    path = os.path.abspath(args[0])
    strain = float(args[1])

    struc = load_structure(path)
    filename = []

    for e,s in enumerate(struc) :
        gs = GMDStructure(s)
        n = gs.vaspname()
        s.apply_strain(strain)
        filename.append("POSCAR_%i_%s"%(e,n))
        s.to(filename="POSCAR_%i_%s"%(e,n))

    createFolder("GMD_file_strain")
    for m in filename :
        move(m,"GMD_file_strain")

def files(args):
    if args.graph :
        graph_yaml(args)
    elif args.strain :
        try : 
            float(args.strain[1])
        except :
            print("Wrongg Argument!")
            sys.exit(0)
        strain(args.strain)
    elif args.input :
        _inputgmd(path=args.input)
    elif args.pdos :
        path = pdos().path
        s = load_structure(os.getcwd())[0]

        # Total DOS
        print("Total DOS")
        os.system("%s dos width=0.03"%(path))

        # p-DOS
        print("pdos")
        dic = defaultdict(list)

        for e,i in enumerate(s.species):
            dic[i.symbol].append(e+1)
        for k,v in dic.items():
            for s in ['tot','s','p'] :
                with open("LIST",'w') as fi :
                    fi.write("{}_{}\n".format(k,s))
                    for i in v :
                        if s == "tot" :
                            fi.write("%i tot tot\n"%(i))
                        elif s == "s" :
                            fi.write("%i s tot\n"%(i))
                        elif s == "p" :
                            fi.write("%i px tot\n"%(i))
                            fi.write("%i py tot\n"%(i))
                            fi.write("%i pz tot\n"%(i))
                os.system("%s pdos width=0.03"%(path))
    elif args.mp :
        if len(args.mp) != 1 :
            print("Please write the mplist")
            sys.exit(0)
        else :
            mp = open(os.path.abspath(args.mp[0]),'r') 
            f = mp.readlines();mplist=[]
            for i in f :
                ff = i.split("\n")
                mplist.append(ff[0])
        mpr = MPJClass()
        s = []
        for mp in mplist :
            ss = mpr.mpr.get_structure_by_material_id(mp)
            name = ss.composition.get_reduced_formula_and_factor()[0]
            ss.to(filename="%s_%s.cif"%(name,mp))

    elif args.transition :
        files = load_structure(args.transition)
        for e,f in enumerate(files) : 
            
            name = GMDStructure(f).vaspname()
            f.to(filename="{}_GMD_{}.cif".format(name,e))