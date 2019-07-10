''' Author: Bauskar Nikhil '''
''' eMail: nikhil.bauskar@student.tu-freiberg.de '''
from Material_Routine import materialRoutine
import numpy as np


def elementRoutine(U_e, T_m, X_e,E_pl):
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
        epsilon_pl = E_pl
        # Calling Material Routine
        Ct, sigma, E_pl = materialRoutine(epsilon,epsilon_pl,T_m)
        # Gauss integration for element stiffness and F_int
        Kt_e = gauss_weights[gp]*np.matmul(np.matmul(B.T,Ct),B)*r**2*det_jacobian
        F_int_e = gauss_weights[gp]*np.matmul(B.T,sigma)*(r**2)*det_jacobian
        F_ext_e = np.zeros_like(N)
        
    return Kt_e, F_int_e, F_ext_e, E_pl
        
    