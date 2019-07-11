import numpy as np 
import matplotlib.pyplot as plt
from pathlib import Path

#To import co-ordinates of nodes from mesh_generation file
from mesh_generation import *
#To import material parameters from material_parameters file
from material_parameters import *
#Assignment 
from assignment_matrix import ASSIGNMENT_MATRIX
#For calling element routine in main
from Element_Routine import elementRoutine


#initialization of time like parameter and displacement vector and state variables
allTimes = np.arange(tau_start+delta_t, tau_end+delta_t, delta_t)
allTimes[-1] = tau_end
nSteps = len(allTimes)

U_g_0=np.zeros((nElem+1,1))
E_pls = np.zeros([nElem,3,1])
global_sigma = np.zeros([nSteps,nElem,3,1])

#loop for iterating load from 0% to 100% with user defined start,stop,minimum,maximum
for i in range(nSteps):
    tau = allTimes[i] 
    u_g = U_g_0 
    u_g[0,0] = 1/3*E_v*tau*rnodes[0]
    current_E_pl = np.zeros_like(E_pls)
    current_sigma = np.zeros([nElem,3,1])

    #DO Newton_Raphson_method
    k=1
    while 1:
        gauss_loc = np.zeros(nElem)
        Kt_g=np.zeros([nElem+1,nElem+1])
        G_global=np.zeros((nElem+1,1))
        F_g_int=np.zeros((nElem+1,1))
        for j in range(nElem):
            A=ASSIGNMENT_MATRIX(j+1,nElem)
            u_e=np.matmul(A,u_g)
            # calling Element routine
            K_e,F_e_int,F_e_ext,E_pl, sigma, r_gp = elementRoutine(u_e,tau,rnodes[j:j+2],E_pls[j])
            Kt_g = Kt_g+np.matmul(np.transpose(A),np.matmul(K_e,A))
            G_global = G_global+np.matmul(np.transpose(A),(F_e_int-F_e_ext))
            F_g_int = F_g_int+np.matmul(np.transpose(A),(F_e_int))
            current_E_pl[j] = E_pl
            current_sigma[j] = sigma        
            gauss_loc[j] = r_gp
        #Reduced system of equations
        K_rg=Kt_g
        K_rg=np.delete(K_rg,0,axis=0)
        K_rg=np.delete(K_rg,0,axis=1)
        reduced_G_global=G_global
        reduced_G_global=np.delete(reduced_G_global,0,axis=0)
        dU_g=np.matmul(np.linalg.inv(K_rg),-reduced_G_global)
        u_g[1:]=u_g[1:]+dU_g

        if (np.linalg.norm(reduced_G_global,np.inf)<0.005*np.linalg.norm(F_g_int,np.inf) or np.linalg.norm(dU_g,np.inf)<0.005*np.linalg.norm(u_g[1:],np.inf)) or k > 5 :     
            break
        else:
            k=k+1
    if k>5:
        print("Convergence criterion is not met, at load (peercentage):",tau*100)
        break
    else:
        U_g_0 = u_g
        E_pls=current_E_pl
        global_sigma[i] = current_sigma


''' POST PROCESSING '''
# # 1. Mesh Plot
fig,ax=plt.subplots()
ax.plot(rnodes, np.zeros_like(rnodes), marker = '*') 
ax.set(xlabel = 'r [$mm$]', title= 'Mesh Distribution')
# output_file_path = Path("Results", "Mesh_Plot.png")
# output_file_path.parent.mkdir(exist_ok=True)
# fig.savefig(output_file_path)
plt.show() 

# #Analytical solution:
disp_analytical = (rnodes[0]**3*E_v)/(3*rnodes**2)
sigma_rr_analytical = (-2*youngs_modulus*E_v*rnodes[0]**3)/(3*(1+poissons_ratio)*np.array(gauss_loc)**3)
  
fig,ax=plt.subplots(ncols=2)
ax[0].plot(rnodes,disp_analytical, marker = 'd') 
ax[0].plot(rnodes,U_g_0, marker = 'o')
ax[0].set(xlabel = 'r [$mm$]', ylabel = 'Nodal Displacements [$\mu$m]', title= 'Mesh Distribution') 

ax[1].plot(gauss_loc, sigma_rr_analytical, marker = 'd') 
ax[1].plot(gauss_loc, global_sigma[-1,:,0], marker = 'o')
ax[1].set(xlabel = 'r [$mm$]', ylabel = 'Stress at gauss points [$MPa$]', title= 'Mesh Distribution')
# output_file_path = Path("Results", "Analytical_Solutions.png")
# output_file_path.parent.mkdir(exist_ok=True)
# fig.savefig(output_file_path)     
plt.show() 

#Results:
fig,ax=plt.subplots(ncols=3,figsize=(15,10)) 
ax[0].plot(rnodes,U_g_0, marker = 'o')
ax[0].set(xlabel = 'r [$mm$]', ylabel = '$U_r$ [$mm$]', title= 'Nodal Displacements') 

ax[1].plot(gauss_loc, global_sigma[-1,:,0], marker = 'o')
ax[1].set(xlabel = 'r [$mm$]', ylabel = '$\sigma_{rr}$ [$MPa$]', title= 'Radial Stress at gauss points')

ax[2].plot(gauss_loc, global_sigma[-1,:,1], marker = 'o')
ax[2].set(xlabel = 'r [$mm$]', ylabel = '$\sigma_{\phi\phi}$ [$MPa$]', title= 'Longitudinal Stress at gauss points')

fig.subplots_adjust(wspace=0.3)
fig.suptitle("Distribution Results (Full Load)")
# output_file_path = Path("Results", "Results_fullLoad.png")
# output_file_path.parent.mkdir(exist_ok=True)
# fig.savefig(output_file_path)  
plt.show()

fig,ax=plt.subplots() 
ax.plot(allTimes, global_sigma[:,0,0], marker = '.')
ax.set(xlabel = '$tau$', ylabel = '$\sigma_{rr}$ [$MPa$]', title= 'Radial Stress at near $r_i$')
# output_file_path = Path("Results", "Results_fullLoad.png")
# output_file_path.parent.mkdir(exist_ok=True)
# fig.savefig(output_file_path)
plt.show()

# Convergence Study






