# coding : utf-8
# Cppyright (c) Green Materials Designs Team.
import os
import sys

import pandas as panda
import numpy as np
from collections import defaultdict

from pymatgen.core import Element, Composition
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.ext.matproj import MPRester
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram

mpr = MPRester()
class GMDPhaseDiagram :
    def __init__ (self, vasprunpath) :
        mpr = MPRester()
        self.vasprunpath = vasprunpath
        self.compatibility = MaterialsProject2020Compatibility()
        
    def get_entries(self) :
        entries = []
        '''
        for v in self.vasprunpath : 
            vrun = Vasprun(v,parse_potcar_file=False)
            entry = vrun.get_computed_entry(inc_structure=True)
            entry = self.compatibility.process_entry(entry)
            entries.append(entry)
        '''
        return entries
    
    def entries_element(self) :
        entries = self.get_entries()
        element = []
        for e in entries :
            for c in e.composition.elements :
                if not c in element : 
                    element.append(c.symbol)
        element = list(set(element))
        return element
    
    def get_phasediagram(self, extractcomp=False) :
        element = self.entries_element()
        mp_entries = mpr.get_entries_in_chemsys(element)
        if extractcomp :   
            if Element("Cl").symbol in element : 
                element.append("Cl2")
            if Element("H").symbol in element :
                element.append("H2")
            if Element("N").symbol in element :
                element.append("N2")
            self.exceptcomp.extend(element)
            except_entries = [e for e in mp_entries if e.composition.reduced_formula in extractcomp]
            entries = self.get_entries() + except_entries
        else :
            entries = self.get_entries() + mp_entries  
        pd = PhaseDiagram(entries)
        return entries, pd
    
    def get_data(entries, pd) :
        data =defaultdict(list) ; number = 1
        for e in entries:
            decomp, ehull = pd.get_decomp_and_e_above_hull(e)
            comp = e.composition 
            formula_num = comp.get_integer_formula_and_factor()[1]
            if e.entry_id == None : 
                data["Materials ID"].append("GMD-{}".format(number))
                number += 1
            else :
                data["Materials ID"].append(e.entry_id)
            data["Composition"].append(e.composition.reduced_formula)
            data["Formula Num"].append(formula_num)
            data["Total E(eV/FU)"].append(e.energy/formula_num)
            data["Ehull"].append(ehull)    
            data["FormE"].append(pd.get_form_energy(e))
            data["Decomposition"].append(" + ".join(["%.2f %s" % (v, k.composition.formula) for k, v in decomp.items()]))
            data['Corrections'].append(e.correction/formula_num)
        df = panda.DataFrame(data, columns=["Materials ID","Composition","Ehull", "FormE","Formula Num","Total E(eV/FU)","Decomposition","Corrections"]).set_index("Materials ID")
        return df

class GMDExcitonbinding :
    def harm_mean_em(em):
        xe = float(em[0].split(":")[-1].split("\n")[0])
        ye = float(em[4].split(":")[-1].split("\n")[0])
        ze = float(em[8].split(":")[-1].split("\n")[0])
        hm = (3*xe*ye*ze)/((xe*ye)+(xe*ze)+(ye*ze))
        return hm

    def dielectricconst(outcar) :
        index = [e for e,o in enumerate(outcar) if "DIELECTRIC" in o]
    
        static = outcar[index[-2]+2:index[-2]+5]
        ionic = outcar[index[-1]+2:index[-1]+5]
        static_n, ionic_n = [],[]
        for s,i in zip(static, ionic) : 
            static_n.append(s.split())
            ionic_n.append(i.split())
        array1 = np.array(static_n).astype(np.float)
        array2 = np.array(ionic_n).astype(np.float)
        dielec = (np.sum(array1)+np.sum(array2))/3
        return dielec

#class AbsorptionCoefficient : 
#    def __init__(self, shell) :
#        self.opticsh = os.system('plotopticsAsymmetriy.sh')

#    def 