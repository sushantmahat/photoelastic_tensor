# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 16:33:22 2021

test code; ismple arithematic. calculated C3 and C4 as described in Physucal 
review B 75, 2007, Pezeral et al.
@author: Mahat
"""
import numpy as np

#
#P11 = 5
#P33 = 4
#P12 = 2
#P13 = 1
#P44 = 6
#
#thetaD = 26.56
#
#
#T = np.deg2rad(thetaD)
#
#c  = np.cos(T)
#s  = np.sin(T)
#ss = np.sin(T)**2
#cc = np.cos(T)**2
#
#C31 = P12* ss + P13* cc
#C32 = (P11*cc+P13*ss)*ss + (P13*cc+P33*ss)*cc - 4*P44*cc*ss
#C33 = (P11*ss+P13*cc)*ss + (P13*ss+P33*cc)*cc + 4*P44*cc*ss
#C34 = c*s*(-1*P11*ss+P13*ss-1*P13*cc+P33*cc-2*P44*np.cos(2*T))
#
#C41 = c*s*(P13-P12)
#C42 = c*s*(cc*(P13-P11)+ss*(P33-P13)+2*P44*np.cos(2*T))
#C43 = c*s*(ss*(P13-P11)+cc*(P33-P13)-2*P44*np.cos(2*T))
#C44 = ((c*s)**2)*((P11-P13)+(P33-P13)+P44*np.cos(2*T)**2) 
#
#print(C31,C32,C33,C34)
#print(C41,C42,C43,C44)s = np.sin(np.deg2rad(-1* polAng/2))




''' test get_rot_mat functon from the material.py main body below.'''

idmat = np.identity(3)

def norm(v):
    """Return the Pythagorean norm of a 3-dim vector."""
    return np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])


def get_rot_mat(vector1, vector2):
    """Return a rotation matrix that rotates vector2 towards vector1."""
    vector1 = np.array(vector1)/norm(vector1)
    vector2 = np.array(vector2)/norm(vector2)
    rotvec = np.cross(vector2, vector1)

    sin_angle = norm(rotvec)
    cos_angle = np.sqrt(1.0 - sin_angle*sin_angle)
    if sin_angle > 1e-10:
        dir_vec = rotvec/sin_angle
    else:
        return idmat

    ddt = np.outer(dir_vec, dir_vec)
    skew = np.array([[        0.0, -dir_vec[2],  dir_vec[1]],
                     [ dir_vec[2],         0.0, -dir_vec[0]],
                     [-dir_vec[1],  dir_vec[0],        0.0]])

    mtx = ddt + cos_angle * (idmat - ddt) - sin_angle * skew
    return mtx    


print(get_rot_mat([1,-1,0], [0,0,1]))














        
        