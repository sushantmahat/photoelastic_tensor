# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 04:47:18 2020

main input file routine. Imports all necessary files and runs the material
routine to get photoelastic constants at different rotations

@author: WorkRel
"""

import configparser  
import numpy as np

# the packages below were homemade
import voigt_not
import material


"""
shape2tensor = voigt_not.de_voigt
shape2matrix = voigt_not.voigt
"""

# reading the input file

config = configparser.ConfigParser()
config.read('input.ini')

#read in main photoelastic constant and other inputs throught input.ini file:
photoelastic_Matrix = list(map(float, config.get('SCAN', 'photoelasticMatrix').split()))
photoelastic_Matrix = np.reshape(photoelastic_Matrix, (6,6))
N_temp              = float(config.get('SCAN', 'refractiveIndexX'))
#N                   = N_temp*np.identity(3)

#other parameters:
P44Step  = 0.0003   # step size of P44 guesses
P44Range = 0.0100   # max difference between initial guess of P44 and all tested P44 values
totP44   = 2 * int(P44Range/P44Step)

P11Step  = 0.0003
P11Range = 0.0100 
totP11   = 2 * int(P11Range/P11Step)

P12Step  = 0.0003
P12Range = 0.0100 
totP12   = 2 * int(P12Range/P12Step)




#read in any orienation information available

zdir = None
xdir = None

try:
    zdir = list(map(float, config.get('SCAN', 'zdir').split()))
    zdir = np.reshape(zdir, (int(len(zdir)/3),3))
except:
    zdir = None
    
try:
    xdir = map(float, config.get('SCAN', 'xdir').split())
    xdir = np.reshape(xdir, (len(xdir)/3,3))
except:
    xdir = None



mat = material.Material(photoelastic_Matrix, N_temp)
#initiate value

'''
the section below cacluates the sum of absolute difference for different photoelastic 
constants used
'''

P11 = mat.matrix[1][1]
P12 = mat.matrix[0][1]
P44 = mat.matrix[3][3]




#for i in range(0,totP11,1):
#    '''
#    nP12 = P12+P12Range-i*P12Step
#    mat.ChangeP12(nP12)
#    #for j in range(len(zdir)):
#    #    mat.rotate_tensor(z_dir = zdir[j])
#    mat.rotate_tensor(z_dir = [5,6,-4])
#        #print(mat.rmatrix)
#    print(str(nP12),str(mat.FindA()))
#    
#    '''
#    nP11 = P11+P11Range-i*P11Step
#    mat.ChangeP11(nP11)
#    #for j in range(len(zdir)):
#    #    mat.rotate_tensor(z_dir = zdir[j])
#    mat.rotate_tensor(z_dir = [7,25,-15])
#    
#      
#    '''
#    nP44 = P44+P44Range-i*P44Step
#    mat.ChangeP44(nP44)
#    #for j in range(len(zdir)):
#    #    mat.rotate_tensor(z_dir = zdir[j])
#    mat.rotate_tensor(z_dir = [5,6,-4])
#        #print(mat.rmatrix)
#    print(str(nP44),str(mat.FindA()))
#    '''






'''commented lines should help in testing; uncomment to see different values at
ech step of the matrix manipulation'''
#print(mat.matrix)
#print(mat.matrix)
mat.rotate_tensor(z_dir = [1,1,0])

mat.strained(st3 = 1, st4 = 0, st5 = 0)
listA = []
listAng = []
np.set_printoptions(precision=4,suppress=True)
#print(mat.DNew)
#print(mat.testMatrix)
print(mat.matrix)
print(mat.rMatrix)
print(mat.diTensor)


for i,d in enumerate(range(0,360,5)):
    listA.append(mat.inLight(polAng = d))
    listAng.append(d)
    #print(mat.diTensor)
    
    
    
#print(listA)
#for i,d in enumerate(listA):
   #print (listAng[i],d/np.amin(listA))
print (np.max(listA)/np.min(listA))
for i,d in enumerate(range(0,360,2)):
    listA.append(mat.inLight(polAng = i))

  






#print(mat.matrix)
#print(mat.N2matrix)
#print(mat.FindA())
#print(mat.inLight())
#print(mat.diTensor)

#mat.ChangeP11(1)
#print(mat.matrix)
#mat.ChangeP12(2)
#print(mat.matrix)
#print(mat.oldMatrix)

#print(id(mat.oldMatrix))
#mat.rotate_tensor(z_dir = [0,0,1])
#mat.strained(st3 = 1, st4 = 0, st5 = 0)

#normalA = mat.inLight(polAng = 0)
#for i in range(0,360,2):
#    print(i, mat.inLight(polAng = i)/normalA)

#print(mat.inLight(polAng = 0))

