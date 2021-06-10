import os
import sys

import numpy as np

from pymatgen.core import Structure

class Octahedron :
    def __init__ (self, initial_structure, final_structure) :
        self.s1 = initial_structure
        self.s2 = final_structure

    def qelonagation(self, index) :
        s1_species = self.s1.species ; s2_species = self.s2.species
        s1_matrix = self.s1.distance_matrix ; s2_matrix = self.s2.distance_matrix

        index1 = s1_matrix[index].argsort()[1:7]
        index2 = s2_matrix[index].argsort()[1:7]

        poscar = s1_matrix[3][index1]
        contcar = s2_matrix[3][index1]

        elonagation= 0 
        for p,c in zip(poscar, contcar):
            elonagation+=(p/c)**2
        elonagation /= len(poscar)

        return elonagation

    def angle_variance(self):
        print("Prepare")