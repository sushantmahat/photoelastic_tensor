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
        self.amp     = None
        self.diTensor= None
        self.strain  = None
        
        #for posterity
        self.oldMatrix = photo_elastic
      
        
        
    def clear_direction(self):
        """Clear all direction-dependent data"""
        self.direction = None

        self.theta = None
        self.phi = None
        self.strain = None

    def rotate_tensor(self, rot_mat = None, x_dir = None, z_dir=None):
        #Applies rotation given a matrix or any comibantion of z and x directions
        
        self.clear_direction()
        #for posterity
        
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
        
        for i in range(4):
            self.tensor = np.tensordot(rot, self.tensor, (1,i))
            self.matrix = shape2matrix(self.tensor)
            self.N2tensor= np.tensordot(rot, self.N2tensor, (1,i))
            self.N2matrix= shape2matrix(self.N2tensor)
    
    '''defining photoelasticity specific funcitons now
    more specifically deining the 6x1 strain matrix that will be multiplied into 
    the photoelastic tensor. by default, all values will be set to 0, with only on
    longitudinal velocity in the z direction set to 1'''
    
    def strained(self, st1 = 0, st2 = 0, st3 = 1, st4 = 0, st5 = 0, st6 = 0):
        
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
        self.diTensor = np.tensordot(self.N2tensor, self.strain, 2)
        
  
    
    def inLight(self, pol1 = None, pol2 = None, polAng = None):
        if self.strain is None or self.diTensor is None:
            self.strained()
            
        #in case someone jumps the gun and runs inLight before straining the sample
        
        if polAng is None and pol1 is None and pol2 is None:
            self.LightP = [1,0,0]
        elif polAng is None:
            self.LightP = [pol1,pol2,0]
        else:
            self.LightP = [np.sin(np.deg2rad(polAng)),np.cos(np.deg2rad(polAng)), 0]
            
        #self.LightP = np.array(self.LightP)/norm(self.LightP)
        self.interact = np.matmul(self.diTensor,self.LightP)
        self.amp      = np.sqrt(self.interact[0]**2 + self.interact[1]**2 ) 
        return self.amp 
         
        
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
        
    def ChangeP12(self, P12):
        self.matrix = self.oldMatrix *1 
        self.matrix[0][1] = P12
        self.matrix[0][2] = P12
        self.matrix[1][0] = P12
        self.matrix[1][2] = P12
        self.matrix[2][0] = P12
        self.matrix[2][1] = P12
        self.tensor = shape2tensor(self.matrix)
        
    def ChangeP44(self, P44):
        self.matrix = self.oldMatrix *1
        self.matrix[3][3] = P44
        self.matrix[4][4] = P44
        self.matrix[5][5] = P44
        self.tensor = shape2tensor(self.matrix)
        
        
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
                
                
                