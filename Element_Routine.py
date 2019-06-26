''' Author: Bauskar Nikhil '''
''' eMail: nikhil.bauskar@student.tu-freiberg.de '''
#from Material_Routine.py import materialRoutine
import numpy as np


def elementRoutine(U_e, T_m, X_e):
    # X_e = np.array([[1],[3]])
    # U_e = np.array([[0.],[0.]])
    gauss_points = np.array([0])
    gauss_weights = np.array([2])

    
    for gp in gauss_points:
        N = np.array([[0.5*(1-gp)], [0.5*(1+gp)]])
        dN_e = np.array([[-0.5], [0.5]])
       
        jacobian = np.matmul(dN_e.T,X_e)
        det_jacobian = jacobian
        r = np.matmul(N.T,X_e)

        B = np.array([
                        [ dN_e[0,0]/det_jacobian[0,0], dN_e[1,0]/det_jacobian[0,0] ], 
                        [ N[0,0]/r[0,0], N[1,0]/r[0,0] ], 
                        [ N[0,0]/r[0,0], N[1,0]/r[0,0] ]
                    ])
        
        epsilon = np.matmul(B,U_e)
        epsilon_pl = 0

        #Ct, sigma = materialRoutine(epsilon,epsilon_pl)
        Ct = 10
        sigma = np.array([[3], [2], [1]])
        Kt_e = np.matmul(B.T,B)*r**2*det_jacobian*Ct

        F_int_e = gauss_weights[gp]*np.matmul(B.T,sigma)*(r**2)*det_jacobian

        F_ext_e = (r**2)*sigma[0]*det_jacobian*(N)

    return Kt_e, F_int_e, F_ext_e
        
    