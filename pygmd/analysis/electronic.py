import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import yaml
import json

from perovgen.pygmd.cli.gmd_write import graphyaml
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun, BSVasprun

def _load_yaml(paraname):
    pwd = os.getcwd()
    graph = "{}/graph.yaml".format(pwd)
    if os.path.isfile(graph):
        stream = open(graph,"r")
        loading = yaml.load(stream)[paraname]
        #loading = {paraname:loading}
    else :
        loading = graphyaml(paraname)
    return loading

class DOSPlotting :
    """
    Class for plotting DOS. Use the pdos program
    """
    def __init__(self, vasprun="vasprun.xml",dos=None, zero_to_efermi=True, stack=True) :
        """
        Args:
            vasprun (str, optional): Path the vasprun.xml . Defaults to "vasprun.xml".
            dos ([type], optional): data with pdos. Defaults to None.
            zero_to_efermi (bool, optional): if to shift all DOS to have zero energy at the fermi energy. Defaults to True.
            stack (bool, optional): If to plot the DOS as a astacked area graph. Defaults to True.
        """
        self.zero_to_efermi = zero_to_efermi 
        self.stack = stack
        self.efermi = Vasprun(vasprun,parse_potcar_file=False).efermi
        self.energies = [float(i.split()[0])-self.efermi if self.zero_to_efermi else float(i.split()[0]) for i in dos]
        self.densities = [float(i.split()[1]) for i in dos]
        
    def get_plot(self, figsize=(12,8), xlim=(-4,6), ylim=None, fontsize=30,color="r",label="Total DOS") :
        """
        Args:
            figsize (tuple, optional): specifies figure size. Defaults to (12,8).
            xlim (tuple, optional): specifies the x-axis limits. Defaults to (-4,6).
            ylim ([type], optional): specifies the y-axis limits. Defaults to None.
            fontsize (int, optional): specifies font size in figure. Defaults to 30.
            color (str, optional): specifies graph color. Defaults to "r".
            label (str, optional): specifies legend label. Defaults to "Total DOS".

        Returns:
            plt (matplotlib.pyplot) 
        """
        plt.rcParams['font.family']='Arial' # font familiy is Arial is fixed
        plt.rcParams['figure.figsize']=figsize
        plt.rcParams['font.size']=fontsize
        plt.rcParams['lines.linewidth'] = 3 # figure line width is 3 is fixed
        plt.rcParams['axes.linewidth'] = 3 # figure box width is 3 is fixed

        allpts = []
        allpts.extend(list(zip(self.energies, self.densities)))

        if xlim :
            plt.xlim(xlim)
        if ylim :
            plt.ylim(ylim)
        else :
            xlim = plt.xlim()
            relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
            plt.ylim((min(relevanty),max(relevanty)))

        if self.stack :
            plt.fill_between(self.energies,0,self.densities,color=color,alpha=.1)
        xlabel = "Energies(eV)"
        if self.zero_to_efermi :
            xlabel = "Energies-E$_{F}$(eV)"
            ylim = plt.ylim()
            plt.plot([0,0], ylim, 'k--', linewidth=2)
		
        plt.xlabel(xlabel)
        plt.ylabel("Density of States")

        plt.plot(self.energies, self.densities, color=color, label=label)
        plt.tight_layout()
        return plt

class BSPlotting :
    def __init__(self, vasprun='vasprun.xml',kpoints='KPOINTS'):
        """
        Args:
            vasprun (str, optional): path the vasprun.xml . Defaults to "vasprun.xml".
            kpoints (str, optional): path the KPOINTS. Defaults to 'KPOINTS'.
        """
        self.bsrun = BSVasprun(vasprun, parse_potcar_file=False,parse_projected_eigen=True)
        self.bs = self.bsrun.get_band_structure(kpoints)
        self.bsdict = self.bs.as_dict()

    def _xlabels(self):
        """
        initialize all the k-point labels and k-point x-distances 
        Returns:
            uniq(dict): Set x-label according to eigenvalues 
        """
        steps=[];uniq_d=[];uniq_l=[]

        for br in self.bs.branches :
            s, e = br['start_index'], br['end_index']
            labels = br['name'].split("-")

            steps.append(e+1)

            if labels[0] == labels[1] :
                continue

            for i,l in enumerate(labels) :
                if l.startswith("\\") or "_" in l :
                    labels[i] = "$"+l+"$"

            if uniq_d != [] and labels[0] != uniq_l[-1] :
                uniq_l[-1] += "$\\mid$" + labels[0]
                uniq_l.append(labels[1])
                uniq_d.append(self.bs.distance[e])
            else :
                uniq_l.extend(labels)
                uniq_d.extend([self.bs.distance[s], self.bs.distance[e]])
        del steps[-1]

        uniq=defaultdict(list)
        uniq['steps'].extend(steps)
        uniq['distance'].extend(uniq_d)
        uniq['labels'].extend(uniq_l)

        return uniq

    def _bandinform(self):
        """
        Band structure information
        
        Returns:
            band_inform (dict): Dict of band data
            {"NB","NK","E_g":{"Direct":...,"Indirect":...},"CBM":{"Direct":...,"Indirect":...},"VBM":{"Direct":...,"Indirect":...}}
        """
        band_inform = dict()

        # the number of the bands and kpoints
        band_inform['NB'] = self.bs.nb_bands
        band_inform['NK'] = len(self.bs.kpoints)

        # Fermi energy & band gap & CBM and VBM

        band_inform['E_f'] = self.bs.efermi
        eg = self.bsdict['band_gap']['energy']
        cbm = self.bsdict['cbm']
        vbm = self.bsdict['vbm']

        cbm_kindex=cbm['kpoint_index'] ; vbm_kindex = vbm['kpoint_index']
        cbm_bindex =cbm['band_index'] ; vbm_bindex = vbm['band_index']

        if eg != 0 :
            cbm1 = [(self.bs.distance[index], cbm['energy']) for index in cbm_kindex]
            vbm1 = [(self.bs.distance[index], vbm['energy']) for index in vbm_kindex]
            
            if self.bsdict['band_gap']['direct'] :
                direct_eg = eg
                indirect_eg = eg
                
                cbm2 = cbm1 ; vbm2 = vbm1
            else :
                direct_dict = self.bs.get_direct_band_gap_dict()[Spin.up]
                indirect_eg = eg
                direct_eg = direct_dict['value']
                direct_kindex = direct_dict['kpoint_index']
                
                vbm2 = [(self.bs.distance[direct_kindex], self.bs.bands[Spin.up][direct_dict['band_indices'][0],direct_kindex])]
                cbm2 = [(self.bs.distance[direct_kindex], self.bs.bands[Spin.up][direct_dict['band_indices'][1],direct_kindex])]
                
            band_inform['E_g'] = {"Direct":direct_eg,"Indirect":indirect_eg}
            band_inform['CBM'] = {"Direct":cbm2, "Indirect": cbm1}
            band_inform['VBM'] = {"Direct":vbm2, "Indirect": vbm1}
        else :
            band_inform['E_g'] = {"Direct":0, "Indirect" : 0}
            band_inform['CBM'] = {"Direct":None,"Indirect":None}
            band_inform['VBM'] = {"Direct":None,"Indirect":None}

        # the position of xlabel
        sum1 = 0 ; name , distance = '', '' 
        xlabel = []
        for d,l in zip(self._xlabels()['distance'],self._xlabels()['labels']):
            xlabel.append(['%.5f'%(d),l])
        band_inform['xlabel'] = xlabel 
        
        # Energies and distances
        steps = [br["end_index"] + 1 for br in self.bs.branches][:-1]
        energies={}
        for sp in self.bs.bands.keys():
            energies[str(sp)]=np.hsplit(self.bs.bands[sp], steps)
        distances = np.split(self.bs.distance, steps)

        band_inform['energies'] = energies
        band_inform['distances'] = distances

        return band_inform
    
    def write_to_json(self,path=os.getcwd()):
        """
        write the band information from _bandinform function.
        Args:
            path ([os], optional): specifies the path to save the band_inform.gmd. Defaults to os.getcwd().
        """
        
        bi = self._bandinform()
        energies = bi['energies']
        distances = bi['distances']
        
        del bi['energies']
        del bi['distances'] 
       
        with open('%s/band_inform.json'%(path),'w') as fi :
            json.dump(bi,fi,indent=4)
        fi.close()

        data = ''
        for ib in range(self.bs.nb_bands) :
                for sp in self.bs.bands.keys():
                    for xpath, epath in zip(distances, energies[str(sp)]):
                        for x1,y1 in zip(xpath, epath[ib]) :
                            data += '%i %.5f %.5f %i\n'%(ib, x1, y1, sp)

        with open("%s/band.gmd"%(path),'w') as fi :
            fi.write(data)
        fi.close()
        
    def get_plot(self, figsize=(12,8), zero_to_efermi=True,color='b',ylim=(-4,6), 
                 fontsize=32, spindownoff=True, vbm_cbm_marker=True, hse_bandgap = False):
        """
        Args:
            figsize (tuple, optional): specifies the figure size. Defaults to (12,8).
            zero_to_efermi (bool, optional): if to shift all DOS to have zero energy at the fermi energy. Defaults to True.
            color (str, optional): specifies the bnad data color. Defaults to 'b'.
            ylim (tuple, optional): specifies the y-axis limits. Defaults to (-4,6).
            fontsize (int, optional): specifies the font size in figure. Defaults to 32.
            spindownoff (bool, optional): if ture, spin down band is not shown. Defaults to True.
            vbm_cbm_marker (bool, optional): if true, vbm and cbm has been marker by red and green, respectively. Defaults to True.

        Returns:
            plt (matplotlib.pyplot) 
        """
        # Figure 
        plt.rcParams['figure.figsize'] = figsize
        plt.rcParams['font.size']=fontsize
        plt.rcParams['font.family'] = 'Arial'

        plt.figure(figsize=figsize)
        # get information from def
        bi = self._bandinform()
        label = self._xlabels()

        # consider the vbm energy
        zero_energy = 0
        if zero_to_efermi :
            if self.bsdict['vbm']['energy'] == None :
                zero_energy = 0
            else :
                zero_energy = self.bsdict['vbm']['energy']
                plt.axhline(0,color='k',lw=1,ls='--')

        # Plotting energies
        bandgap = self.bsdict['band_gap']['energy']
        if hse_bandgap :
            for ib in range(self.bs.nb_bands) :
                    for sp in self.bs.bands.keys():
                        for xpath, epath in zip(bi['distances'], bi['energies'][str(sp)]):
                                if str(sp) == '-1' and spindownoff == False :
                                    plt.plot(xpath, epath[ib] - zero_energy,color='r')
                                    plt.plot(xpath, epath[ib] - zero_energy + (hse_bandgap-bandgap),color='r--')
                                else :
                                    plt.plot(xpath, epath[ib] - zero_energy,color=color)
                                    if len(epath[ib][epath[ib] - zero_energy <= 0]-zero_energy) == len(epath[ib] - zero_energy):
                                        plt.plot(xpath, epath[ib] - zero_energy ,color='r',ls='-.',lw=2)
                                    else :
                                        plt.plot(xpath, epath[ib] - zero_energy +(hse_bandgap-bandgap),color='r',ls='-.',lw=2)
                                    
        else :
            for ib in range(self.bs.nb_bands) :
                    for sp in self.bs.bands.keys():
                        for xpath, epath in zip(bi['distances'], bi['energies'][str(sp)]):
                                if str(sp) == '-1' and spindownoff == False :
                                    plt.plot(xpath, epath[ib] - zero_energy,color='r')
                                else :
                                    plt.plot(xpath, epath[ib] - zero_energy,color=color)

        # decorating the plot
        plt.xticks(label['distance'],label['labels'])
        plt.xlim(min(label['distance']),max(label['distance']))
        plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
        if zero_to_efermi : 
            ylabel = r'$\mathrm{E\ -\ E_{VBM}\ (eV)}$' 
            if self.bsdict['vbm']['energy'] == None : 
                ylabel = r'$\mathrm{Energy\ (eV)}$'
        else :
            ylabel = r'$\mathrm{Energy\ (eV)}$'
        plt.ylabel(ylabel, fontsize=30)
        for i in range(len(label['distance'])):
            plt.axvline(label['distance'][i],color='k',lw=1)
        plt.ylim(ylim)
        
        # cbm and vbm
        if bandgap != 0 and vbm_cbm_marker :
            if self.bsdict['band_gap']['direct'] :
                for c in bi['CBM']['Direct'] :
                    plt.scatter(c[0],c[1]-zero_energy,color='g',s=(fontsize*5))
                    if hse_bandgap :
                        plt.scatter(c[0],c[1]-zero_energy +(hse_bandgap-bandgap) ,color='g',s=(fontsize*5),marker='D')
                for v in bi['VBM']['Direct'] :
                    plt.scatter(v[0], v[1]-zero_energy,color='#FF0000',s=(fontsize*5))
            else :
                for c in bi['CBM']['Indirect'] :
                    plt.scatter(c[0],c[1]-zero_energy,color='g',s=(fontsize*5))
                    if hse_bandgap :
                        plt.scatter(c[0],c[1]-zero_energy +(hse_bandgap-bandgap) ,color='g',s=(fontsize*5),marker='D')
                for v in bi['VBM']['Indirect'] :
                    plt.scatter(v[0], v[1]-zero_energy,color='#FF0000',s=(fontsize*5))

                for c in bi['CBM']['Direct'] :
                    plt.scatter(c[0],c[1]-zero_energy,color='purple',s=(fontsize*5)/2)
                    if hse_bandgap :
                        plt.scatter(c[0],c[1]-zero_energy +(hse_bandgap-bandgap) ,color='purple',s=(fontsize*5)/2,marker='D')
                for v in bi['VBM']['Direct'] :
                    plt.scatter(v[0], v[1]-zero_energy,color='y',s=(fontsize*5)/2)
        plt.tight_layout()
        return plt

class GMDAnalysis :
    def __init__(self) :
        sep = os.path.abspath(__file__).split(os.sep)
        pwd = ''
        for j in sep :
            pwd += '{}/'.format(j)
            if j == 'perovgen' :
                break
        self.pdos = "{}pdos/pdos".format(pwd)
        self.emc = "{}/emc-master/emc.py".format(pwd)
        
    def effectivemass(self,path,secondstep) :
        os.chdir(path)
        if secondstep :
            before_filelist = os.listdir(".")
            os.system("python %s INPCAR EIGENVAL"%(self.emc))
            after_filelist = os.listdir(".")
            for i in after_filelist : 
                if not i in before_filelist :
                    os.system("mv {}/{} {}/EM".format(path,i,path))
                    break
        else :
            os.system("python %s INPCAR"%(self.emc))

    def partialDOS(self,width,structure):
        dic = defaultdict(list)
        for e,j in enumerate(structure.species) :
            dic[j.symbol].append(e+1)
        for k,v in dic.items():
            for s in ['tot','s','p','d'] :
                with open("LIST","w") as fi :
                    fi.write("{}_{}\n".format(k,s))
                    for u in v :
                        if s == "tot" :
                            fi.write("%i tot tot\n"%(u))
                        elif s == "s" :
                            fi.write("%i s tot \n"%(u))
                        elif s == "p" :
                            fi.write("%i px tot\n"%(u))
                            fi.write("%i py tot\n"%(u))
                            fi.write("%i pz tot\n"%(u))
                        elif s == 'd' :
                            fi.write("%i dxy tot\n"%(u))
                            fi.write("%i dyz tot\n"%(u))
                            fi.write("%i dz2 tot\n"%(u))
                            fi.write("%i dxz tot\n"%(u))
                            fi.write("%i x2-y2 tot\n"%(u))
                fi.close()
                os.system("%s pdos width=%f"%(self.pdos,width))