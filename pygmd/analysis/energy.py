# coding : utf-8
# Cppyright (c) Green Materials Designs Team.
import os
import sys

import pandas as panda
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from collections import defaultdict

from pymatgen.core import Element, Composition
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.ext.matproj import MPRester
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram

class GMDPhaseDiagram :
    def __init__ (self, vasprunpath) :
        #mpr = MPRester()
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
                data["Materials ID"].append("None-{}".format(number))
                number += 1
            else :
                data["Materials ID"].append(e.entry_id)
            data["Composition"].append(e.composition.reduced_formula)
            data["Formula Num"].append(formula_num)
            data['energy'].append(e.energy)
            data['energy_per_atom'].append(e.energy_per_atom)
            data["formation_energy"].append(pd.get_form_energy(e))
            data['formation_energy_per_atom'].append(pd.get_form_energy_per_atom(e))
            data["Ehull"].append(ehull)    
            data["Decomposition"].append(" + ".join(["%.2f %s" % (v, k.composition.formula) for k, v in decomp.items()]))
            data['Corrections'].append(e.correction/formula_num)
        df = panda.DataFrame(data).set_index("Materials ID")
        return df

class GMDExcitonbinding :
    def harm_mean_em(em_path) :
        EM = open(em_path,'r',encoding='utf-8').readlines()
        for f in EM :
            if "(0)" in f :
                xe = float(f.split(":")[-1].split("\n")[0])
            elif "(1)" in f :
                ye = float(f.split(":")[-1].split("\n")[0])
            elif "(2)" in f :
                ze = float(f.split(":")[-1].split("\n")[0])
        hm = (3*xe*ye*ze)/((xe*ye)+(xe*ze)+(ye*ze))
        return {"x":xe,"y":ye,"z":ze,"HarmMean":abs(hm)}

    def dielectricconst(outcar_path) :
        outcar = open(outcar_path,'r').readlines()
        index = [e for e,o in enumerate(outcar) if "DIELECTRIC" in o]

        static = outcar[index[0]+2:index[0]+5]
        ionic = outcar[index[-1]+2:index[-1]+5]
        
        static_n, ionic_n = [],[]
        for s,i in zip(static, ionic) : 
            static_n.append(s.split())
            ionic_n.append(i.split())
            
        array1 = np.array(static_n).astype(np.float)
        array2 = np.array(ionic_n).astype(np.float)
        high_mean = 3/((1/array1[0].sum())+(1/array1[1].sum())+(1/array1[2].sum())) 

        data = {"Array":{"high_frequency_dielectric":array1,"ionic_contribution":array2,"static_dielectric":array1+array2},
                "geom_mean":high_mean}
        return data
    
    def excitonbinding(harm_mean, geom_mean) :
        return (harm_mean*13.605692/(geom_mean**2))

class GMDAbsorption :
    # version 1
    # created date : 2019/01/29
    def make_csv(path) :
        pwd = os.getcwd()
        os.chdir(path)
        if not "optics.dat" in os.listdir(".") :
            sh_path = "%s/plotopticsAsymmetriy.sh"%(os.path.dirname(__file__))
            os.system("sh %s"%(sh_path))
        
        hbar = 6.58E-16 ; c = 3.0E+8
        fi = open("optics.dat",'r').readlines()
        dataframe = defaultdict(list)
        for f in fi[:1000] :
            energy,ximg,xreal,yimg,yreal,zimg,zreal = f.split()
            absor_x = ((np.sqrt(2)*float(energy)/hbar)*np.sqrt((np.sqrt(float(xreal)**2 + float(ximg)**2)-float(xreal))))/100/c
            absor_y = ((np.sqrt(2)*float(energy)/hbar)*np.sqrt((np.sqrt(float(yreal)**2 + float(yimg)**2)-float(yreal))))/100/c
            absor_z = ((np.sqrt(2)*float(energy)/hbar)*np.sqrt((np.sqrt(float(zreal)**2 + float(zimg)**2)-float(zreal))))/100/c
        
            dataframe['energy'].append(float(energy))
            dataframe['X'].append(absor_x)
            dataframe['Y'].append(absor_y)
            dataframe['Z'].append(absor_z)
        os.chdir(pwd)
        df = panda.DataFrame(dataframe).set_index("energy")
        return df 
    
    def plotting(path, fig_size=(8,6),fontsize=14, xlim=(0,5),ylim=None,sharex=True,unit='eV') :
        df = panda.read_csv(path,index_col=0)
        plt.rcParams['font.size'] = 14
        absorp_data = df.values
        if unit == 'nm' :
            hbar = 6.58E-16 ; c = 3.0E+8
            x = df.index
            x /= hbar*c
        elif unit == 'eV' :
            pass
        else :
            print("unit is only nm or eV")
            sys.exit(1)
        
        fig, axes = plt.subplots(2,1,figsize=fig_size,sharex=sharex)
        axes[0].plot(x,absorp_data[:,0],'r--',label='X')
        axes[0].plot(x,absorp_data[:,1],'g--',label='Y',alpha=.7)
        axes[0].plot(x,absorp_data[:,2],'b--',label='Z',alpha=.5)
        axes[0].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2e"))
        axes[0].legend(frameon=False)

        axes[1].semilogy(x,absorp_data[:,0],'r--',label='X')
        axes[1].semilogy(x,absorp_data[:,1],'g--',label='Y',alpha=.7)
        axes[1].semilogy(x,absorp_data[:,2],'b--',label='Z',alpha=.5)
        axes[1].legend(frameon=False)
            
        if xlim :
            allpts = []
            allpts.extend(list(zip(x, absorp_data[:,0])))
            allpts.extend(list(zip(x, absorp_data[:,1])))
            allpts.extend(list(zip(x, absorp_data[:,2])))
            relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
            axes[0].set_xlim(xlim)
            axes[0].set_ylim((min(relevanty),max(relevanty)))
        
        if ylim : 
            axes[0].set_ylim(ylim) 

        axes[0].set_title(r"Energy - Absorption Coefficient ($\alpha$)")
        axes[1].set_xlabel("Energy (%s)"%(unit))
        axes[0].set_ylabel(r"$\alpha$ (cm$^{-1}$)")
        axes[1].set_ylabel(r"$ln$($\alpha$) (cm$^{-1}$)")
        fig.tight_layout()
        return plt