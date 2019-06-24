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
number_of_steps=10
initial_tau=0
final_tau=1
delta_t=(final_tau-initial_tau)/number_of_steps
tau=0.0

#initialization of time like parameter and displacement vector
time=np.array([0])
displacement_matrix_global_initial=np.zeros((number_of_elements+1,1))

#loop for iterating load from 0% to 100% with user defined start,stop,minimum,maximum
for i in range(number_of_steps):
    tau=tau+delta_t           #see if tau values are generating properly
    displacement_matrix_global=displacement_matrix_global_initial[-1]  
    #append displacement_matrix_global_initial at the end of iteration
    delta_displacement_matrix_global=np.zeros((number_of_elements+1,1))

    #DO Newton_Raphson_method
        k=1
        while 1:
            stiffness_matrix_global=np.zeros([number_of_elements+1,number_of_elements+1])
            G_global=np.zeros((number_of_elements+1,1))
            Force_internal_global=np.zeros((number_of_elements+1,1))
            for j in range(number_of_elements):
                assignment_matrix=ASSIGNMENT_MATRIX(j,number_of_elements)
                displacement_matrix_element=np.matmul(assignment_matrix,displacement_matrix_global)
                # calling Element routine 
                stiffness_matrix_element,force_internal_element,force_external_element=elementRoutine(displacement_matrix_element,tau,rnodes[j:j+2])
                stiffness_matrix_global=stiffness_matrix_global+np.matmul(np.transpose(assignment_matrix),np.matmul(stiffness_matrix_element,assignment_matrix))
                G_global=G_global+np.matmul(np.transpose(assignment_matrix),(force_internal_element-force_external_element))
                Force_internal_global=Force_internal_global+np.matmul(np.transpose(assignment_matrix),(force_internal_element))
            #Implementation of essential boundary conditions in K and G
            G_global[0,0]=-10
            delta_displacement_matrix_global=-G_global/stiffness_matrix_global
            displacement_matrix_global=displacement_matrix_global+delta_displacement_matrix_global
            k=k+1 
            #Implementation of essential boundary conditions in displacement_matrix_global
            if (np.linalg.norm(G_global)<=0.005*np.linalg.norm(Force_internal_global) or np.linalg.norm(delta_displacement_matrix_global)<=0.005*np.linalg.norm(displacement_matrix_global)) and k>=5:
            break
        displacement_matrix_global=np.append(displacement_matrix_global,displacement_matrix_global)


