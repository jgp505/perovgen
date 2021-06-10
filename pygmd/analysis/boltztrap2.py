import os
import sys

from pymatgen.electronic_structure.boltztrap2 import *
from monty.serialization import loadfn
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter

vrun = Vasprun('vasprun.xml')
data = VasprunBSLoader(vrun)

bs = vrun.get_band_structure()
nele = vrun.parameters['NELECT']
st = vrun.final_structure
data = VasprunBSLoader(bs,structure=st,nelect=nele)
data = VasprunBSLoader.from_file('vasprun.xml')

# Interpolate Band
bztInterp = BztInterpolator(data,lpfac=10,energy_range=1.5,curvature=True)
sbs = bztInterp.get_band_structure()

# Compute and plot DOS
#tot_dos = bztInterp.get_dos()
#tot_proj_dos = bztInterp.get_dos(partial_dos=True,progress=False)

# save and load coefficients
# set fname argument to specify a different file name
bztInterp = BztInterpolator(data,lpfac=10,energy_range=1.5,curvature=True,
                            save_bztInterp=True,fname='bztInterp.json')

# set fname argument to specify a different file name
bztTransp = BztTransportProperties(bztInterp,temp_r = np.arange(300,1300,300), doping=10.**np.arange(16,23),
                                   save_bztTranspProps=True,fname='bztTranspProps.json')
bztTransp = BztTransportProperties(bztInterp,load_bztTranspProps=True,fname='bztTranspProps.json')

bztPlotter = BztPlotter(bztTransp,bztInterp)
bztPlotter.plot_props('C','mu','temp').show()