# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 17:14:37 2021

@author: WorkRel
"""

"""
Standalone code for making the input file for photoelastic code. 
Made in a versatile manner to allow porting to other codes.
 
The only fucntion is to write a 6X6 matrix into a temmporary file that can then
be copied or renamed to make an input.ini file. Takes in info. about the 
symmerty and the unique elastic constant values. Info is given in a seperate 
file. (i.e. for cubic crystal, only C11, C12, and C44 will be provided)


Made just so to avoid having to edit changging an enitre 6X6 matrix i input. 
More chances of getting a value incorrect that way, and not necssary for our 
problem as we do not work with general symmetries.

currently supports simple cubic and isotropic symmetry
If working with general symmetry, recommend chaning the output file of this code
directly as that feeds into the input file of the main code.

a 0 symmetry refers to isotropic and a 1 refers to cubic symmetry by default.

"""



import configparser  
import numpy as np

configMM = configparser.ConfigParser()

# specify output file name here:L

outptNm = 'tempMatrix.ini'

#will read input file to determine what type of matrix to build.

configMM.read('matrixVal.ini')

M11 = configMM.get('SCAN', 'M11')
M44 = configMM.get('SCAN', 'M44')

symmetry = configMM.getfloat('SCAN', 'symmetry')



if symmetry == 0:
    M12 = M11
    M13 = M12
    M66 = M11
elif symmetry == 1:
    M12 = configMM.get('SCAN', 'M12')
    M13 = M12
    M66 = M44
elif symmetry == 2:
    M12 = configMM.get('SCAN', 'M12')
    M13 = configMM.get('SCAN', 'M13')
    M66 = 0.5*(float(M11) - float(M12))

#pretty self explanatory as everything is bruteforced
    
MatrixFile = open(outptNm, 'w')


MatrixFile.write('[SCAN] \n \n' )
MatrixFile.write('Matrix = \n')

MatrixFile.write('\t' + M11 + '\t' + M12 + '\t' + M13 + '\t'+ '0'+ '\t'+ 
                 '0' + '\t' + '0 \n')
MatrixFile.write('\t' + M12 + '\t' + M11 + '\t' + M13 + '\t'+ '0'+ '\t'+ 
                 '0' + '\t' + '0 \n')
MatrixFile.write('\t' + M13 + '\t' + M13 + '\t' + M11 + '\t'+ '0'+ '\t'+ 
                 '0' + '\t' + '0 \n')
MatrixFile.write('\t' + '0' + '\t' + '0' + '\t' + '0' + '\t'+ M44+ '\t'+ '0' +
                 '\t' + '0 \n')
MatrixFile.write('\t' + '0' + '\t' + '0' + '\t' + '0' + '\t'+ '0'+ '\t'+ M44 +
                 '\t' + '0 \n')
MatrixFile.write('\t' + '0' + '\t' + '0' + '\t' + '0' + '\t'+ '0'+ '\t'+ '0' +
                 '\t' + str(M66) + '\n')

MatrixFile.close()
#symmtery = configMM.getfloat('SCAN', 'symmetry')


#you can also use makeMatrix to make a matrix:
# I will need to improve this as the matrix become more complex

def makeMatrix(C11, C12, C44, outputFileName):
    FxnFile = open(outputFileName, 'w')
    FxnFile.write('[SCAN] \n \n' )
    FxnFile.write('Matrix = \n')
    FxnFile.write('\t' + M11 + '\t' + M12 + '\t' + M12 + '\t'+ '0'+ '\t'+ '0' 
                  + '\t' + '0 \n')
    FxnFile.write('\t' + M12 + '\t' + M11 + '\t' + M12 + '\t'+ '0'+ '\t'+ '0' 
                  + '\t' + '0 \n')
    FxnFile.write('\t' + M12 + '\t' + M12 + '\t' + M11 + '\t'+ '0'+ '\t'+ '0' 
                  + '\t' + '0 \n')
    FxnFile.write('\t' + '0' + '\t' + '0' + '\t' + '0' + '\t'+ M44+ '\t'+ '0' 
                  + '\t' + '0 \n')
    FxnFile.write('\t' + '0' + '\t' + '0' + '\t' + '0' + '\t'+ '0'+ '\t'+ M44 
                  + '\t' + '0 \n')
    FxnFile.write('\t' + '0' + '\t' + '0' + '\t' + '0' + '\t'+ '0'+ '\t'+ '0' 
                  + '\t' + M44 + '\n')
    FxnFile.close()
    
    











