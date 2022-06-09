# -*- coding: utf-8 -*-
"""
Modified on Tue Jan 25 21:44:09 2022

This code can be used to rotate second and fourth order tensors 
(represented in 6x6 Voigt notation).

work partly based on paper by: slezak and mocek 2020
@author: Sushant Mahat
"""
#required modules: 
import sympy as sp
import numpy as np
np.set_printoptions(precision=3, suppress = True)

#Definition of invariants:
#the following matrices ensure symmetry when 
#changing from voigt to tensor notation
Iinv = sp.diag(1,1,1,2,2,2)
I    = sp.diag(1,1,1,0.5,0.5,0.5)


#Inputs:

#define terms of the 6x6 photoleastic constant matrix;
#
#examples of symbolic values:
p11 = sp.symbols('p11')
p12 = sp.symbols('p12')
#p13 = sp.symbols('p13')
p44 = sp.symbols('p44') 


##one can also use numerical values
#p11 = 0.15
#p12 = 0.095
#p44 = 0.072

#define rows of an initial 6x6 photoelastic constant (before rotation) for 
# cubic cyrtsals:

Prow1 = [p11, p12, p12, 0, 0, 0]
Prow2 = [p12, p11, p12, 0, 0, 0]
Prow3 = [p12, p12, p11, 0, 0, 0]
Prow4 = [0, 0, 0, p44, 0, 0]
Prow5 = [0, 0, 0, 0, p44, 0]
Prow6 = [0, 0, 0, 0, 0, p44]


#define the 6x6 photoelastic constant for hexagon:

#Prow1 = [p11, p12, p13, 0, 0, 0]
#Prow2 = [p12, p11, p13, 0, 0, 0]
#Prow3 = [p13, p13, p11, 0, 0, 0]
#Prow4 = [0, 0, 0, p44, 0, 0]
#Prow5 = [0, 0, 0, 0, p44, 0]
#Prow6 = [0, 0, 0, 0, 0, 0.5 * (p11-p12)]

#create fourth order tensor in Voigt notation from its row components:
matrix_P = sp.Matrix([Prow1, Prow2, Prow3, Prow4, Prow5, Prow6])


#define strain tensor in Voigt notation:
strain_M = sp.Matrix([0,0,1,0,0,0])


#define terms of the 3x3 rotation tensor:
##Symbolic rotation terms:
#a11 = sp.symbols('a11')
#a12 = sp.symbols('a12')
#a13 = sp.symbols('a13')
#
#a21 = sp.symbols('a21')
#a22 = sp.symbols('a2')
#a23 = sp.symbols('a23')
#
#a31 = sp.symbols('a31')
#a32 = sp.symbols('a32')
#a33 = sp.symbols('a33')

#we can also define actual numerical values to the rotation matrix:
#look into substituion in sympy notes to refine this section

##below is defined for (110) [1-10]
#a11 = 0.7071 
#a12 = -0.7071
#a13 = 0
#
#a21 = 0
#a22 = 0
#a23 = -1
#
#a31 = 0.7071
#a32 = 0.7071
#a33 = 0
####

#below is the code for (111) [0-11]
#
#a11 = -0.408 
#a12 = 0.816 
#a13 = 0.408
#
#a21 = 0.707 
#a22 = 0
#a23 = -0.707
#
#a31 = 0.577
#a32 = 0.577
#a33 = 0.577

#below is the code for (111) [10-1]
##
#a11 = -0.70711
#a12 = 0
#a13 = 0.70711
#
#a21 = 0.40825 
#a22 = -0.8165
#a23 = 0.40825
#
#a31 = 0.57735
#a32 = 0.57735
#a33 = 0.57735

#below is the code for (311) [0-11]:

#a11 = 0 
#a12 = -0.70711
#a13 =  0.70711
#
#a21 = 0.4264
#a22 =  -0.6396
#a23 = -0.6396
#
#a31 = 0.90453
#a32 = 0.30151
#a33 = 0.30151


#below is the code for (-2 20 5):

a11 = 0.99249   
a12 = 0.07511 
a13 =  0.09656

a21 =  0.07511
a22 =   0.24891 
a23 = -0.96561

a31 = -0.09656 
a32 =  0.96561  
a33 = 0.2414

#below is the code for (322):

#a11 = 0.64351 
#a12 = -0.23766 
#a13 =  -0.72761
#
#a21 = -0.23766   
#a22 =    0.84156  
#a23 =  -0.48507
#
#a31 = 0.72761 
#a32 =   0.48507  
#a33 =  0.48507
#


#define rows of a 6x6 rotation matrix in voigt notation:

Arow1 = [a11*a11, a12*a12, a13*a13, 2*a12*a13, 2*a13*a11, 2*a11*a12]
Arow2 = [a21*a21, a22*a22, a23*a23, 2*a22*a23, 2*a23*a21, 2*a21*a22]
Arow3 = [a31*a31, a32*a32, a33*a33, 2*a32*a33, 2*a31*a33, 2*a31*a32]

Arow4 = [a21*a31, a22*a32, a23*a33, a23*a32 + a33*a22, a21*a33 + a31*a23, 
         a22*a31 + a32*a21]

Arow5 = [a31*a11, a32*a12, a33*a13, a33*a12 + a13*a32, a13*a31 + a11*a33,
         a32*a11 + a12*a31]

Arow6 = [a11*a21, a12*a22, a13*a23, a13*a22 + a23*a12, a11*a23 + a21*a13, 
         a12*a21 + a11*a22]
#create rotation matrix from the rows
matrix_A = sp.Matrix([Arow1, Arow2, Arow3, Arow4, Arow5, Arow6])




#calculate transverse rotation matrix:
T_matrix_A = matrix_A.T

#calculate rotated fourth order tensor
rot_P = matrix_A * matrix_P * T_matrix_A #* Iinv

#calculate change in second order tensor when coupled to strain:
del_D = rot_P * strain_M

sp.pprint(rot_P.evalf(n=3))
sp.pprint(del_D.evalf(n=3))
#need to optimize this:
#sp.init_printing(use_unicode=True, wrap_line=True, no_global=True)
#sp.pprint(rot_P,full_prec = False)



















