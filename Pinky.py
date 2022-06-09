# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 09:40:59 2022
Pinky et al. Rotation in cubic crystals; codes to calculate photoelastic constant
of uniaxial texture
@author: Mahat
"""

import numpy as np
np.set_printoptions(precision=3, suppress = True)
#definition of euler angles 1 (eu1) and euler angle 2 (eu2)
eu1 = np.radians(0)
eu2 = np.radians(45)


#defininf old system:
P11 = 0.15
P12 = 0.095
P44 = 0.072

#definition of cosine and sine terms to simply any future equations:
cos1 = np.cos(eu1)
sin1 = np.sin(eu1)

cos2 = np.cos(eu2)
sin2 = np.cos(eu2)

#computing squares

sqcos1 = cos1**2
sqsin1 = sin1**2


sqcos2 = cos2**2
sqsin2 = sin2**2



qcos1 = cos1**4
qsin1 = sin1**4


qcos2 = cos2**4
qsin2 = sin2**4

#define rotation matrix; doesn't agree with my code or expected results from hkl rotations:

Rot_M = np.array([[cos2*cos1, - sin1, -sin2*cos1], 
         [cos2*sin1, cos1, -sin2*sin1],
         [sin2, 0, cos2]])

print(Rot_M)

#defining photoelastic constants:

NP11 =  P11 * (qcos2 * (qcos1 + qsin1) + qsin2) + 2 * P12 * (qcos2 * sqcos1 * sqsin1 + sqcos2 * sqsin2 ) + \
 P44 * (qcos2 * sqcos1 * sqsin1 + sqcos2 * sqsin2)

NP22 = P11 * (qcos1 + qsin1) + 2 * P12 * (sqcos1 * sqsin1) + P44 * (sqcos1 * sqcos1)

NP33 = P11 * (qsin2 * (qcos1 + qsin1) + qcos2) + 2 * P12 * (qsin2 * sqcos1 * sqsin1 + sqcos2 * sqsin2 ) + \
 P44 * (qsin2 * sqcos1 * sqsin1 + sqcos2 * sqsin2)

NP12 = 2 * P11 * (sqcos1 * sqsin1) + P12 * (sqcos2 * (qcos1 + qsin1) + sqsin1) - P44 * (sqcos1 * sqsin1 * sqcos2)
NP21 = NP12

NP13 = P11 * (sqsin2 * sqcos2 * (qcos1 + qsin1 + 1)) + P12 * (2 * sqcos2 * sqsin2 * sqcos1 * sqsin1 + qcos2 + qsin2) - \
 P44 * (sqcos2 * sqsin2 - sqcos1 * sqsin1 * sqcos2 * sqsin2)
NP31 = NP13

NP23 = 2 * P11 * (sqcos2 * sqcos1 * sqsin1) + P12 * (sqsin2 * (qcos1 + qsin1) + sqcos2) - P44 * (sqcos1 * sqsin1 * sqsin2) 
NP32 = NP23

NP66 = 2 * P11 * (sqcos2 * sqcos1 * sqsin1) - 2 * P12 * (sqcos2 * sqcos1 * sqsin1) + P44 * (sqsin2 * sqsin1 + qcos1 * sqcos2)

NP44 = 2 * P11 * (sqcos2 * sqcos1 * sqsin1) - 2 * P12 * (sqsin2 * sqcos1 * sqsin1) + P44 * (sqsin2 * qsin1 + qcos1 * sqcos2)

NP55 = P11 * (sqsin2 * sqcos2 * (qcos1 + qsin1 + 1)) + 2 * P12 * (sqcos2 * sqsin2 * sqcos1 * sqsin1 - sqcos2 * sqsin2) \
 + P44 * (qcos2 * sqcos1 + sqcos1 * sqsin1 * sqcos2 * sqsin2 + qsin2 * sqsin1)




print(np.array([[NP11,NP21,NP31,0,0,0], [NP12,NP22,NP32,0,0,0],[NP13,NP23,NP33,0,0,0],[0,0,0,NP44,0,0],[0,0,0,0,NP55,0],[0,0,0,0,0,NP66]]))