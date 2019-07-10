import numpy as np 
import matplotlib.pyplot as plt

#f=open('Output_Nfem.txt','a')
#f100=open('Output_Nfem_100.txt','a')
#fk=open('Output_Nfem_k.txt','a')

#To import co-ordinates of nodes from mesh_generation file
from mesh_generation import *
#To import material parameters from material_parameters file
from material_parameters import *
#Assignment 
from assignment_matrix import ASSIGNMENT_MATRIX
#For calling element routine in main
from Element_Routine import elementRoutine


number_of_steps = 10
initial_tau = 0
final_tau = 1
delta_t = (final_tau-initial_tau)/number_of_steps
tau= 0.0
nElem = number_of_elements


#initialization of time like parameter and displacement vector
# tau = tau_start
# allTimes = np.arange(tau_start+delta_t, tau_end+delta_t, delta_t)
# nSteps = len(allTimes)

U_g_0=np.zeros((nElem+1,1))
E_pls = np.zeros([nElem,3,1])

# global_U = np.zeros([nSteps,nElem+1,1])
# global_E_pl = np.zeros([nSteps,nElem,3,1])
# global_sigma = np.zeros([nSteps,nElem,3,1])


#loop for iterating load from 0% to 100% with user defined start,stop,minimum,maximum
for i in range(number_of_steps):#len(allTimes)):
    tau = round(delta_t*(i+1),5) #allTimes[i] # round(delta_t*(i+1),5)
    u_g = U_g_0 # global_U[i]
    u_g[0,0] = 1/3*E_v*tau*rnodes[0]
    E_plss = np.zeros_like(E_pls)
    # current_E_pl = global_E_pl[i]
    #DO Newton_Raphson_method
    k=1
    while 1:
        Kt_g=np.zeros([nElem+1,nElem+1])
        G_global=np.zeros((nElem+1,1))
        F_g_int=np.zeros((nElem+1,1))
        for j in range(nElem):
            A=ASSIGNMENT_MATRIX(j+1,nElem)
            u_e=np.matmul(A,u_g)
            # calling Element routine
            #a = global_E_pl[i,j]
            K_e,F_e_int,F_e_ext,E_pl = elementRoutine(u_e,tau,rnodes[j:j+2],E_pls[j])  #Eplsilon_plastic to be saved globally
            #K_e,F_e_int,F_e_ext,E_pl = elementRoutine(u_e,tau,rnodes[j:j+2],global_E_pl[i-1,j])  #Eplsilon_plastic to be saved globally
            Kt_g = Kt_g+np.matmul(np.transpose(A),np.matmul(K_e,A))
            G_global = G_global+np.matmul(np.transpose(A),(F_e_int-F_e_ext))
            F_g_int = F_g_int+np.matmul(np.transpose(A),(F_e_int))
            E_plss[j] = E_pl
            #current_E_pl[j] = E_pl
            
        #Reduced system of equations
        K_rg=Kt_g
        K_rg=np.delete(K_rg,0,axis=0)
        K_rg=np.delete(K_rg,0,axis=1)
        reduced_G_global=G_global
        reduced_G_global=np.delete(reduced_G_global,0,axis=0)
        dU_g=np.matmul(np.linalg.inv(K_rg),-reduced_G_global)
        u_g[1:]=u_g[1:]+dU_g

        if (np.linalg.norm(reduced_G_global,np.inf)<0.005*np.linalg.norm(F_g_int,np.inf) or np.linalg.norm(dU_g,np.inf)<0.005*np.linalg.norm(u_g[1:],np.inf)) :
            break
        else:
            k=k+1
    U_g_0 = u_g
    #global_U[i] = u_g
    E_pls=E_plss
    #global_E_pl[i] = current_E_pl


#Analytical solution:
#disp_ana=(rnodes[0]**3*E_v)/(3*rnodes**2)

# print('Displacements Analytical:\n',disp_ana,file=f100)
# print('Displacements Fem:\n',U_g_0,file=f100)
# print('Difference:\n',disp_ana-U_g_0,file=f100)

# f.close()
# f100.close()
# fk.close()  

# fig,ax=plt.subplots()
# ax.plot(rnodes,disp_ana/1e-6) 
# ax.plot(rnodes,U_g_0/1e-6) 
# plt.show()


    
    
    
    
    
    
    
    
    #U_g_0=np.append(U_g_0,u_g,axis=0)              # Appending to be take care of
    #E_pls=np.append(E_pls,E_pl,axis=1)             # Appending to be take care of
    #print(u_g)


