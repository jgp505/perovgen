import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pymatgen.analysis.diffraction.xrd import XRDCalculator

class GMDXRDPattern :
    # Ver 1
    # 2021-12-23
    # Made by Kyoung Ha Min
    
    def __init__(self, structure) :
        self.struc = structure
    
    def dataframe(self,wavelength='CuKa',sigma=0.03) :
        xrd = XRDCalculator(wavelength=wavelength).get_pattern(self.struc)
        xrd_x = xrd.as_dict()['x']
        xrd_y = xrd.as_dict()['y']
        dataframe = dict()
        maximum = 0
        x = np.linspace(1,120,11900) # same length with VESTA 
        y = 0 
        for i,j in zip(xrd_x,xrd_y):
            mu = i
            y1 = (1 / np.sqrt(2 * np.pi * sigma**2)) * np.exp(-(x-mu)**2 / (2 * sigma**2))
            y += y1*j
        y = y*(100/y.max())
        dataframe['2theta'] = x
        dataframe['Intensity'] = y 
        
        df = pd.DataFrame(dataframe).set_index("2theta")
        return df
    
    def get_plot(self,figsize=(8,6),xlim=(0,50),color='b',linewidth=1) :
        dataframe = self.dataframe()
        x = dataframe.index
        y = dataframe['Intensity']
        plt.figure(figsize=figsize)
        plt.plot(x, y, lw=linewidth, color=color, label = self.struc.composition.reduced_formula)
        plt.xlim(xlim)
        plt.legend(frameon=False)
        return plt
