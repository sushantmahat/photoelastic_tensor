# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 04:47:18 2020

main input file routine. Imports all necessary files and runs the material
routine to get photoelastic constants at different rotations

@author: WorkRel
"""

import configparser  
import numpy as np
np.set_printoptions(precision=5, suppress = True)
# the packages below were homemade
import voigt_not
import material
import matplotlib.pyplot as plt

"""
shape2tensor = voigt_not.de_voigt
shape2matrix = voigt_not.voigt
"""

# reading the input file


config = configparser.ConfigParser()
config.read('tempMatrix.ini')

#read in main photoelastic constant:

photoelastic_Matrix = list(map(float, config.get('SCAN', 'Matrix').split()))
photoelastic_Matrix = np.reshape(photoelastic_Matrix, (6,6))
N_temp = 3.5
#N_temp              = float(config.get('SCAN', 'refractiveIndexX'))
#N                   = N_temp*np.identity(3)



#read in any orienation information available
try:
    zdir = map(float, config.get('SCAN', 'zdir').split())
except:
    zdir = None
try:
    xdir = map(float, config.get('SCAN', 'xdir').split())
except:
    xdir = None


mat = material.Material(photoelastic_Matrix, N_temp)

mat.rotate_tensor(z_dir = [-2, 20, 5]) #, x_dir = [0, -1, 1])

mat.strained(st3 = 1, st5 = 0, st4 = 0)
#print(mat.matrix)
print(mat.rMatrix)
print(mat.rot_mat)
#print(mat.Dtensor)
#print(mat.matrix)
#print(mat.rot_mat)

print(mat.diTensor)

normal = mat.inLight(polAng = 0)
print(normal)


datx = []
daty = []
for i in range(0,180,5):
    #print(i,'\t', mat.inLight(polAng = i)/normal)
    datx.append(i)
    daty.append(mat.inLight(polAng = i)/normal)
    print(i,'\t', mat.inLight(polAng = i)/normal)
    #print(mat.diTensor)
print(mat.inLight(polAng = 90))


ax = plt.scatter(datx,daty)
plt.text(0,max(daty),str(max(daty)))
plt.axhline(y = max(daty), label = 0)








