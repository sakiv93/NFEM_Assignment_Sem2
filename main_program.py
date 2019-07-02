import numpy as np 
import matplotlib.pyplot as plt

#To import co-ordinates of nodes from mesh_generation file
from mesh_generation import *
#To import material parameters from material_parameters file
from material_parameters import *
#Assignment 
from assignment_matrix import ASSIGNMENT_MATRIX
#For calling element routine in main
from Element_Routine import elementRoutine

#solving parameters
number_of_steps = 100
initial_tau = 0
final_tau = 1
delta_t = (final_tau-initial_tau)/number_of_steps
tau= 0.0
E_pls = np.zeros([3,1])
nElem = number_of_elements

#initialization of time like parameter and displacement vector
time=np.array([0])
U_g_0=np.zeros((nElem+1,1))

#loop for iterating load from 0% to 100% with user defined start,stop,minimum,maximum
for i in range(number_of_steps+1):
    tau=delta_t*i          #see if tau values are generating properly
    #U_g=U_g_0[-1]
    U_g=U_g_0  
    #append U_g_0 at the end of iteration
    ##### dU_g=np.zeros((nElem+1,1))

    #DO Newton_Raphson_method
    k=1
    while 1:
        Kt_g=np.zeros([nElem+1,nElem+1])
        G_global=np.zeros((nElem+1,1))
        F_g_int=np.zeros((nElem+1,1))
        for j in range(nElem):
            A=ASSIGNMENT_MATRIX(j+1,nElem)
            U_e=np.matmul(A,U_g)
            # calling Element routine 
            K_e,F_e_int,F_e_ext,E_pl=elementRoutine(U_e,tau,rnodes[j:j+2],E_pls[-1])  #Eplsilon_plastic to be saved globally
            print(tau)
            Kt_g=Kt_g+np.matmul(np.transpose(A),np.matmul(K_e,A))
            G_global=G_global+np.matmul(np.transpose(A),(F_e_int-F_e_ext))
            F_g_int=F_g_int+np.matmul(np.transpose(A),(F_e_int))
        #Implementation of essential boundary conditions in K and G
        #G_global[0,0]=-10
        U_g[0,0]=1/3*E_v*tau*10    # have to chane U_g assignment
        #print(U_g)
        #Reduced system of equations
        K_rg=Kt_g
        K_rg=np.delete(K_rg,0,axis=0)
        K_rg=np.delete(K_rg,0,axis=1)
        reduced_G_global=G_global
        reduced_G_global=np.delete(reduced_G_global,0,axis=0)
        #print(Kt_g)
        dU_g=np.matmul(np.linalg.inv(K_rg),-reduced_G_global)
        U_g[1:]=U_g[1:]+dU_g
        k=k+1
        
        
        #Implementation of essential boundary conditions in U_g
        if (np.linalg.norm(G_global)<=0.005*np.linalg.norm(F_g_int) or np.linalg.norm(dU_g)<=0.005*np.linalg.norm(U_g)) or k>=5:
            break
    print(k) 
    U_g_0 = U_g
    E_pls = E_pl
print(U_g_0)
    #U_g_0=np.append(U_g_0,U_g,axis=0)              # Appending to be take care of
    #E_pls=np.append(E_pls,E_pl,axis=1)             # Appending to be take care of
    #print(U_g)


