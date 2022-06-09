# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 21:30:48 2021

@author: Sushant Mahat

contains work of Jan Jaeken 2016
"""
"""
the code is supposed to take a matrix of shape 6X6 and convert it to an appropriate tensor of 3x3x3x3 shape. The code was adapted from Jan Jaeken's work regarding elastic constants. They were working in python 2 so some changes to the language (indexing was abit off for python 3, 00 was treated as 0, so the counting had to start at 11)
had to be made to make the code work.

It is a seperated code becuase I recon I mgiht have need of converting between voigt_not and back again in the future for other programs unrelated to photoelasticity.
"""
# TLDR: coverts back and forth between voigt matrix and tensorial notation

import numpy as np

#I define a indexed system for voigt notation
Voigt = {11: 0, 22: 1, 33: 2, 23: 3, 32: 3, 13: 4, 31: 4, 12: 5, 21: 5}
#this basically means if I call 00 item from list VOIGT, it will return a value of 0
#similarly an item 01 will return value of 5
#tis is from christoffel.py from authors Jaeken

def de_voigt(C_ij):
    "Turn a 6x6 matrix into a 3x3x3x3 tensor according to Voigt notation."
    C_ijkl = [[[[C_ij[Voigt[10*(i+1)+j+1]][Voigt[10*(1+k)+l+1]]
                 for i in range(3)] for j in range(3)]
                 for k in range(3)] for l in range(3)]
    return C_ijkl

def de_voigt2(vec):
    """Turn a 6-dim vector into a 3x3 tensor according to Voigt notation."""
    T_ij = [[vec[Voigt[10*(i+1)+j+1]] for i in range(3)] for j in range(3)]
    return T_ij

def voigt(C_ijkl):
    "Turn a 3x3x3x3 tensor into a 6x6 matrix"
    # Note from Jaeken original:  Divide by 2 because symmetrization will double the main diagonal
    
    C_ij = np.zeros((6,6))
    
    C_ij[0,0] = 0.5*C_ijkl[0][0][0][0]
    C_ij[1,1] = 0.5*C_ijkl[1][1][1][1]
    C_ij[2,2] = 0.5*C_ijkl[2][2][2][2]
    C_ij[3,3] = 0.5*C_ijkl[1][2][1][2]
    C_ij[4,4] = 0.5*C_ijkl[0][2][0][2]
    C_ij[5,5] = 0.5*C_ijkl[0][1][0][1]

    C_ij[0,1] = C_ijkl[0][0][1][1]
    C_ij[0,2] = C_ijkl[0][0][2][2]
    C_ij[0,3] = C_ijkl[0][0][1][2]
    C_ij[0,4] = C_ijkl[0][0][0][2]
    C_ij[0,5] = C_ijkl[0][0][0][1]

    C_ij[1,2] = C_ijkl[1][1][2][2]
    C_ij[1,3] = C_ijkl[1][1][1][2]
    C_ij[1,4] = C_ijkl[1][1][0][2]
    C_ij[1,5] = C_ijkl[1][1][0][1]

    C_ij[2,3] = C_ijkl[2][2][1][2]
    C_ij[2,4] = C_ijkl[2][2][0][2]
    C_ij[2,5] = C_ijkl[2][2][0][1]

    C_ij[3,4] = C_ijkl[1][2][0][2]
    C_ij[3,5] = C_ijkl[1][2][0][1]

    C_ij[4,5] = C_ijkl[0][2][0][1]
    
    
    return C_ij + C_ij.T

def voigt2(C_ij):
    "Turn a 3x3 tensor into a 6x1 matrix"
    C_i = np.zeros(6)
     
    C_i[0] = C_ij[0,0]
    C_i[1] = C_ij[1,1]
    C_i[2] = C_ij[2,2]
 
    C_i[3] = C_ij[1][2]
    C_i[4] = C_ij[0][2]
    C_i[5] = C_ij[0][1]
     
    return C_i    
     
     
# following line sof codes are for testing purposes.

"""
uses makeMatrix or atleast needs a file created by makeMAtrix.   
create an input file (.ini) with the appropriate file name to run this test.
"""
'''
import configparser  

config = configparser.ConfigParser()
config.read('tempMatrix.ini')

M_tensor = list(map(float, config.get('SCAN', 'Matrix').split()))

M_matrix = np.reshape(M_tensor, (6,6))

M_tensor = de_voigt(M_matrix)

M_test   = voigt(M_tensor)

print(M_matrix)
print(M_test)
'''
