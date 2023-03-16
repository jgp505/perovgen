import os
import pstats
import sys
from tkinter import W

import numpy as np
import pandas as pd
import json
from collections import defaultdict

from pymatgen.core import Structure
from pymatgen.analysis.eos import EOS
from pymatgen.io.vasp.outputs import Vasprun, BSVasprun

from perovgen.pygmd.input_structure import GMDStructure
from perovgen.pygmd.analysis.electronic import BSPlotting, DOSPlotting, GMDAnalysis
from perovgen.pygmd.analysis.energy import GMDExcitonbinding
from perovgen.pygmd.autocal.inputset import inputgmd

def checkcal(path) :
    try :
        vrun = Vasprun("{}/vasprun.xml".format(path))
        return True
    except :
        return False

class GMDExtractData :
    def __init__(self, root_path, inputgmd_path) :
        self.root_path = os.path.abspath(root_path)
        self.inputgmd = inputgmd(inputgmd_path)
        self.dataframe = dict()
        self.columns = {'bulkmodulus':['b0','bulkdata_path'], 
                        'chgcar':['optimized_structure','full_formula','pretty_formula','spacegroup','energy','energy_per_atom','elements'],
                        'band':['NB','NK','E_g','CBM','VBM','xlabel','band_path','E_f'], 
                        'dos':['dos_path'], 
                        'HSE':['HSE_gap'],
                        'effective_mass':['me_x','me_y','me_z','mh_x','mh_y','mh_z','electron_mean','hole_mean','reudced_effective_mass'],
                        'dielectric':['high_frequency_dielectric','ionic_contribution','static_dielectric','geom_mean'], 
                        'exciton':['exciton_binding_energy']}
        
    def input_version(self) :
        relaxtion_path = "{}/relaxation".format(self.root_path)
        vrun = Vasprun("%s/vasprun.xml"%(relaxtion_path))
        version = vrun.vasp_version
        potcar = vrun.potcar_symbols
        kpoint = vrun.kpoints.as_dict()
        incar = vrun.incar
        self.dataframe['inputs'] = potcar
        self.dataframe['kpoints'] = kpoint
        self.dataframe['incar'] = incar
        self.dataframe['VASP-version'] = version
        return self.dataframe
        
    def bulkmodulus(self) :
        bulkpath = "{}/bulkmodulus".format(self.root_path)
        if 'bulkmodulus.gmd' in os.listdir(bulkpath) and os.path.getsize("{}/bulkmodulus.gmd".format(bulkpath)) != 0:
            bulkdata = pd.read_csv("{}/bulkmodulus.gmd".format(bulkpath),sep=' ', header=None).values
        else :
            direct = [f for f in os.listdir(bulkpath) if os.path.isdir(f)]
            direct.sort()
            energy, volume = [], []
            for d in direct :
                vrun = Vasprun("%s/%s/vasprun.xml"%(bulkpath,d))
                energy.append(vrun.final_energy)
                volume.append(vrun.final_structure.volume)

            with open("%s/bulkmodulus.gmd"%(bulkpath), 'w') as fi :
                for e,v in zip(energy,volume) :
                    fi.write("%.4f %.4f\n"%(v,e))
        try :
            eos = EOS(eos_name='birch_murnaghan').fit(bulkdata[:,0],bulkdata[:,1])
            self.dataframe['b0'] = "%.5f"%(eos.b0)
        except :
            self.dataframe['b0'] = None
        self.dataframe['bulkdata_path'] = "{}/bulkmodulus.gmd".format(bulkpath)
        return self.dataframe
    
    def chgcar(self, soc=False) :
        chgcar_path = "{}/chgcar".format(self.root_path)
        if soc :
            chgcar_path = "{}/chgcar_SOC".format(self.root_path)
        if checkcal(chgcar_path) :
            vrun = Vasprun("{}/vasprun.xml".format(chgcar_path))
            struc = vrun.final_structure 
            gmd_struc =  GMDStructure(struc)

            #structure        
            struc_json = struc.as_dict()
            self.dataframe['optimized_structure'] = struc_json  
            self.dataframe['full_formula'] = gmd_struc.formula(reduced=False)
            self.dataframe['pretty_formula'] = gmd_struc.formula(reduced=True)
            self.dataframe['spacegroup'] = struc.get_space_group_info()[1]
            element_name = ''
            for e,ele in enumerate(struc.composition.elements) :
                if e+1 == len(struc.composition.elements) :
                    element_name += ele.symbol
                else :
                    element_name += '{}-'.format(ele.symbol)
            self.dataframe['elements'] = element_name
            
            #energy
            self.dataframe['energy'] = vrun.final_energy
            self.dataframe['energy_per_atom'] = vrun.final_energy / struc.composition.num_atoms
        else :
            for col in self.columns['chgcar'] :
                self.dataframe[col] = None

    def band(self, soc = False) :
        band_path = "{}/chgcar/band".format(self.root_path)
        if soc :
            band_path = "{}/chgcar_SOC/band_SOC".format(self.root_path)
            
        if checkcal(band_path) :
            bsp = BSPlotting(vasprun='{}/vasprun.xml'.format(band_path),kpoints='{}/KPOINTS'.format(band_path))
            band_inform = bsp._bandinform()
            energies = band_inform['energies']
            distances = band_inform['distances']
        
            #del bi['energies']
            #del bi['distances'] 
       
            data = ''
            for ib in range(bsp.bs.nb_bands) :
                    for sp in bsp.bs.bands.keys():
                        for xpath, epath in zip(distances, energies[str(sp)]):
                            for x1,y1 in zip(xpath, epath[ib]) :
                                data += '%i %.5f %.5f %i\n'%(ib, x1, y1, sp)
            for k,v in band_inform.items():
                if k == 'energies' or k == 'distances' :
                    pass
                elif k == 'E_g' or k == 'CBM' or k == 'VBM' :
                    for k1,v1 in v.items():
                        self.dataframe["{}_{}".format(k,k1)] = v1
                else :
                    self.dataframe[k] = v
            self.dataframe['band_data'] = data
        return self.dataframe
   
    def dos(self, soc=False, zero_to_efermi=False) :
        path =  '{}/chgcar/dos'.format(self.root_path)
        if soc :
            path = "{}/chgcar_SOC/dos_SOC".format(self.root_path)
        self.dataframe['dos_path']=[]
        if checkcal(path):
            for i in os.listdir(path) :
                if "dos" in i and "out" != i.split(".")[-1] :
                    dos = open("{}/{}".format(path,i),'r').readlines()[1:]
                    #dp = DOSPlotting(vasprun = '{}/vasprun.xml'.format(path),dos=dos, zero_to_efermi=zero_to_efermi)
                    energies = [float(i.split()[0])-self.efermi if zero_to_efermi else float(i.split()[0]) for i in dos]
                    densities = [float(i.split()[1]) for i in dos]
                    self.dataframe['dos_path'].append({"{}_energies".format(i):energies,
                                                       "{}_densities".format(i):densities})
        return self.dataframe
    
    def effectivemass(self,soc=False) :
        abspath = '{}/chgcar/effective_mass'.format(self.root_path)
        if soc :
            abspath = "{}/chgcar_SOC/effective_mass_SOC".format(self.root_path)
        pwd = os.getcwd() 
        emc = GMDAnalysis()
        if "electron" in os.listdir(abspath) and "hole" in os.listdir(abspath):
            os.chdir("{}/electron".format(abspath))
            if not "EM" in os.listdir(".") :
                emc.effectivemass(".", secondstep=True)
            electron = "{}/EM".format(os.getcwd())

            os.chdir("{}/hole".format(abspath))
            if not "EM" in os.listdir(".") :
                emc.effectivemass(".", secondstep=True)
            hole = "{}/EM".format(os.getcwd())
        os.chdir(pwd)
                
        electron = GMDExcitonbinding.harm_mean_em(electron)
        hole = GMDExcitonbinding.harm_mean_em(hole)
        harm_mean = 2/((1/electron['HarmMean']) + (1/hole['HarmMean']))
        self.dataframe['me_x'] = electron['x']
        self.dataframe['me_y'] = electron['y']
        self.dataframe['me_z'] = electron['z']
        self.dataframe['electron_mean'] = electron['HarmMean']
        self.dataframe['mh_x'] = hole['x']
        self.dataframe['mh_y'] = hole['y']
        self.dataframe['mh_z'] = hole['z']
        self.dataframe['hole_mean'] = hole['HarmMean']
        self.dataframe['reudced_effective_mass'] = harm_mean
        return self.dataframe

    def dielectricconst(self) :
        path = "{}/dielectric".format(self.root_path)
        if checkcal(path) :
            array = GMDExcitonbinding.dielectricconst("{}/dielectric/OUTCAR".format(self.root_path))
            for k,v in array.items():
                if k == 'Array' :
                    for k1, v1 in v.items() :
                        self.dataframe[k1] = v1.tolist()
                else :
                    self.dataframe[k] = v
        else :
            self.dataframe["high_frequency_dielectric"] = None 
            self.dataframe["ionic_contribution"]= None 
            self.dataframe["static_dielectric"] = None 
            self.dataframe["geom_mean"] = None
        return self.dataframe
    
    def excitonbinding(self) :
        abspath = '{}/chgcar/effective_mass'.format(self.root_path)
        pwd = os.getcwd() 
        emc = GMDAnalysis()
        if "electron" in os.listdir(abspath) and "hole" in os.listdir(abspath):
            os.chdir("{}/electron".format(abspath))
            if not "EM" in os.listdir(".") :
                emc.effectivemass(".", secondstep=True)
            electron = "{}/EM".format(os.getcwd())

            os.chdir("{}/hole".format(abspath))
            if not "EM" in os.listdir(".") :
                emc.effectivemass(".", secondstep=True)
            hole = "{}/EM".format(os.getcwd())
        os.chdir(pwd)
                
        electron = GMDExcitonbinding.harm_mean_em(electron)
        hole = GMDExcitonbinding.harm_mean_em(hole)
        harm_mean = 2/((1/electron['HarmMean']) + (1/hole['HarmMean']))
        if checkcal("{}/dielectric".format(self.root_path)) :
            geom_mean = GMDExcitonbinding.dielectricconst("{}/dielectric/OUTCAR".format(self.root_path))['geom_mean']
            self.dataframe['exciton_binding_energy'] = (harm_mean*13.605692*100)/(geom_mean**2)
        else :
            self.dataframe['exciton_binding_energy'] = None
        return self.dataframe
    
    def HSE(self) :
        global hse_gap
        global hse_soc_gap
        
        hsepath = '{}/HSE'.format(self.root_path)
        hsesocpath = '{}/HSE_SOC'.format(self.root_path)
        if 'hse_bandgap.gmd' in os.listdir(self.root_path) :
            hsedata = open("{}/hse_bandgap.gmd".format(self.root_path),'r').readlines()
            if len(hsedata) == 1 :
                if 'None' == hsedata[0].split(":")[-1].split()[0] :
                    hse_gap = None
                else :
                    hse_gap = float(hsedata[0].split(":")[-1].split()[0])
                hse_soc_gap = None
            else :
                if 'None' == hsedata[0].split(":")[-1].split()[0] :
                    hse_gap = None
                else :
                    hse_gap = float(hsedata[0].split(":")[-1].split()[0])
                if 'None' == hsedata[1].split(":")[-1].split()[0] :
                    hse_soc_gap = None
                else :
                    hse_soc_gap = float(hsedata[1].split(":")[-1].split()[0])
                
            self.dataframe['HSE_gap'] = {"HSE_bandgap":hse_gap,'HSE_SOC_bandgap':hse_soc_gap}
        else :
            if checkcal(hsepath) :
                vrun = BSVasprun("{}/vasprun.xml".format(hsepath))
                hse_gap = vrun.as_dict()['output']['bandgap']
            else :
                hse_gap = None
            if checkcal(hsesocpath) :
                vrun = BSVasprun("{}/vasprun.xml".format(hsesocpath))
                hse_soc_gap = vrun.as_dict()['output']['bandgap']
            else :
                hse_soc_gap = None
        self.dataframe['HSE_gap'] = {"HSE_bandgap":hse_gap,'HSE_SOC_bandgap':hse_soc_gap}
        return self.dataframe 