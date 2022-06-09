# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 13:29:26 2021


defining a class called material

contains all relevant info on the materuial like stiffness constants, density, and photoelastic constant


@author: Mahat
"""

import numpy as np


import voigt_not

shape2tensor = voigt_not.de_voigt
shape2matrix = voigt_not.voigt
shapeV2tensor= voigt_not.de_voigt2

idmat = np.identity(3)
#identity marix for various purposes

"""
Naming convention:
    old*         = any value before changing P11, P12 or P44. 
    r*           = any value after rotation
    *matrix*     = Voigt notation
    *tensor*     = tensorial notation
"""


class Material:
    def __init__(self, photo_elastic, refractive_ind):
        self.matrix = photo_elastic
        
        self.tensor = shape2tensor(photo_elastic)
        
        self.N      = refractive_ind

        self.N2matrix= self.N * self.N * self.matrix
        #the term N2matrix is equivalent to the N matrix defined in equation 54
        #by Pezeril et al. 2007, but for cubic symmetry. For complex symmetries,
        #the refractive index will need to be adjusted.
        
        #value placeholders and initiations below:
        self.N2tensor= shape2tensor(self.N2matrix)

        
        #for posterity
        self.oldMatrix = photo_elastic
        self.clear_direction()
        
        
    def clear_direction(self):
        """Clear all direction-dependent data"""
        self.strain    = None
        self.amp       = None #amplitude of del R vs polarization angle function
        self.diTensor  = None #change in dielectric tensor
        self.strain    = None #strain tesnor
        self.rMatrix   = None #rotated matrix
        self.rTensor   = None #rotated tensor
        self.rN2tensor = None #rotated refractive index
        self.rN2matrix = None #rotated refractive index in matrix form (1x6)

    def rotate_tensor(self, rot_mat = None, x_dir = None, z_dir=None):
        #Applies rotation given a matrix or any comibantion of z and x directions
        self.clear_direction()              #remove all existing rotated values
        
        #initilize all values that need to be rotated
        self.rN2tensor = self.N2tensor * 1  
        self.rtensor   = self.tensor   * 1
        
        #define new orientation of crystal and prepare rotation tensor:
        if rot_mat is None:
            rot = idmat
            if z_dir is not None and x_dir is None:
                z_dir = z_dir / norm(z_dir)
                rot = get_rot_mat(z_dir, [0.0, 0.0, 1.0])
            if x_dir is not None and z_dir is None:
                x_dir = x_dir / norm(x_dir)
                rot = get_rot_mat(x_dir, [1.0, 0.0, 0.0])
            if x_dir is not None and z_dir is not None:
                x_dir = x_dir / norm(x_dir)
                z_dir = z_dir / norm(z_dir)

                x_dir -= np.dot(x_dir, z_dir) * z_dir # Gram-Schmidt
                x_dir = x_dir / norm(x_dir)
                y_dir = np.cross(z_dir, x_dir)

                rot = np.array([x_dir, y_dir, z_dir])
        else:
            rot = rot_mat
        
        #start rotating to new orientation:
        for i in range(4):
            self.rtensor  = np.tensordot(rot, self.rtensor, (1,i))
            self.rmatrix  = shape2matrix(self.rtensor)
            self.rN2tensor= np.tensordot(rot, self.rN2tensor, (1,i))
            self.rN2matrix= shape2matrix(self.rN2tensor)
    
    '''defining photoelasticity specific funcitons now
    more specifically deining the 6x1 strain matrix that will be multiplied into 
    the photoelastic tensor. by default, all values will be set to 0, with only on
    longitudinal velocity in the z direction set to 1'''
    
    def strained(self, st1 = 0, st2 = 0, st3 = 1, st4 = 0, st5 = 0, st6 = 0):
        
        #if the sample is not rotated, set values to unrotated crystal
        if self.rN2tensor is None:
            self.rN2tensor = self.N2tensor * 1 
        
        
        strain = np.zeros((3,3))
        #define all elements of the strain tensor
        strain[0][0] = st1
        strain[1][1] = st2
        strain[2][2] = st3
        strain[0][1] = st6
        strain[1][0] = st6
        strain[0][2] = st5
        strain[2][0] = st5
        strain[1][2] = st4
        strain[2][1] = st4
        
        #set the current strain value as tensor
        self.strain  = strain
        
        #calculate new tensor for dielectric constant:
        self.diTensor = np.tensordot(self.rN2tensor, self.strain, 2)
        
  
    
    def inLight(self, polAng = 0):
        
        #if the crystal is not strained, it will need to be strained
        if self.strain is None or self.diTensor is None:
            self.strained()
            
        #check how the polarization angle is specified       
        self.LightP = [np.cos(np.deg2rad(polAng)),np.sin(np.deg2rad(polAng)),0]
        
        
        #polarization information needs to be normalized
        self.LightP = np.array(self.LightP)/norm(self.LightP)
        #light interacts with the changed dielectric tensor
        #interact   = np.matmul(self.diTensor,self.LightP)
        interact2d = []
        #interact2d.append(interact[0] * self.LightP[0])
        #interact2d.append(interact[1] * self.LightP[1])
        
        interact2d.append((self.LightP[0] * self.diTensor[0][0]+self.LightP[1] * self.diTensor[0][1]) * self.LightP[0])
        interact2d.append((self.LightP[0] * self.diTensor[1][0]+self.LightP[1] * self.diTensor[1][1]) * self.LightP[1])
        #self.interact2d = interact2d #uncomment for testing
        
         
        
        #amplitude of reflected electric field from dn:
        self.amp1      = interact2d[0] + interact2d[1]  #+ interact3[1]**2 ) 
        #self.amp      = self.interact[0]*self.LightP[0] + self.interact[1]*self.LightP[1]
        return self.amp1
         
        
    def FindA(self, AngleRange = 180, stepS = 2 ):
        #find the amplitude of the reflecticted intesnity with probe polarization
        normalA  = self.inLight(polAng = 0)
        Aup      = 0
        Adown    = 100
        
        for i in range(0,AngleRange, stepS):
            Acur = self.inLight(polAng = i)/normalA
            if Aup < Acur:
                Aup = Acur
            if Adown > Acur:
                Adown = Acur
        return (Aup - Adown)/2
        
    
    def ChangeP11(self, P11):
        self.matrix = self.oldMatrix *1
        #the * 1 helps because otherwise python will save both oldmatrix and matrix as tags for same value
        #adding a mathematical function forces python to distinguish between the two instances
        
        self.matrix[0][0] = P11
        self.matrix[1][1] = P11
        self.matrix[2][2] = P11
        self.tensor = shape2tensor(self.matrix)
        self.N2matrix= self.N * self.N * self.matrix
        self.N2tensor= shape2tensor(self.N2matrix)
        
    def ChangeP12(self, P12):
        self.matrix = self.oldMatrix *1 
        self.matrix[0][1] = P12
        self.matrix[0][2] = P12
        self.matrix[1][0] = P12
        self.matrix[1][2] = P12
        self.matrix[2][0] = P12
        self.matrix[2][1] = P12
        self.tensor = shape2tensor(self.matrix)
        self.N2matrix= self.N * self.N * self.matrix
        self.N2tensor= shape2tensor(self.N2matrix)
        
    def ChangeP44(self, P44):
        self.matrix = self.oldMatrix *1
        self.matrix[3][3] = P44
        self.matrix[4][4] = P44
        self.matrix[5][5] = P44
        self.tensor = shape2tensor(self.matrix)
        self.N2matrix= self.N * self.N * self.matrix
        self.N2tensor= shape2tensor(self.N2matrix)
    
  


    
    ''' not implremented completely yet:    
    def set_direction_cartesian(self, direction):
        
        self.clear_direction()
        
        
        self.direction = direction / norm(direction)
        q = self.direction

        x = q[0]
        y = q[1]
        z = q[2]
        if z >= 1.0 or z <= -1.0:
            if z > 0.0:
                self.theta = 0.0
            else:
                self.theta = np.pi
            self.phi = 0.0
        else:
            self.theta = np.arccos(z)
            sin_theta = np.sqrt(1 - z**2)

            cos_phi = x/sin_theta

            self.phi = np.arccos(cos_phi)
            if y < 0.0:
                self.phi = 2.0*np.pi - self.phi

        self.material = np.dot(q, np.dot(q, self.tensor))
    '''

                
                
""" some definitions to aid with tensorial calculations """

def determinant(m):
    """Return the determinant of a 3x3 matrix."""
    return (m[0][0] * m[1][1] * m[2][2] -
            m[0][0] * m[1][2] * m[2][1] +
            m[0][1] * m[1][2] * m[2][0] -
            m[0][1] * m[1][0] * m[2][2] +
            m[0][2] * m[1][0] * m[2][1] -
            m[0][2] * m[1][1] * m[2][0])

def norm(v):
    """Return the Pythagorean norm of a 3-dim vector."""
    return np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def cofactor(m):
    """
    Return the cofactor matrix of a 3x3 matrix.
    """
    cof = np.empty((3, 3))

    cof[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1]
    cof[0][1] = m[1][2]*m[2][0] - m[1][0]*m[2][2]
    cof[0][2] = m[1][0]*m[2][1] - m[1][1]*m[2][0]

    cof[1][0] = m[0][2]*m[2][1] - m[0][1]*m[2][2]
    cof[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0]
    cof[1][2] = m[0][1]*m[2][0] - m[0][0]*m[2][1]

    cof[2][0] = m[0][1]*m[1][2] - m[0][2]*m[1][1]
    cof[2][1] = m[0][2]*m[1][0] - m[0][0]*m[1][2]
    cof[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0]

    return cof
                
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
                
                
                