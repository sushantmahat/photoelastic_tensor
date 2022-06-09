# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 07:59:21 2022

@author: WorkRel
"""
import numpy as np 



a1 = 0
a2 = -1.5
a3 = -a1 - a2 + 1 #eqn of plane for z = 110, x = 001

print(a3)



#a1 = 1
#a2 = -1
#a3 = 0

b1 = 0
b2 = 0
b3 = 1


absA = (a1**2 + a2**2 + a3**2 ) ** 0.5
absB = (b1**2 + b2**2 + b3**2 ) ** 0.5

a_b = a1 * b1 + a2 *b2 + a3 * b3

theta = np.rad2deg(np.arccos(a_b/(absA * absB)))

print(theta)


