import numpy as np 
import matplotlib.pyplot as plt

#To import co-ordinates of nodes from mesh_generation file
from mesh_generation import *
#To import material parameters from material_parameters file
from material_parameters import *
#Assignment 
from assignment_matrix import ASSIGNMENT_MATRIX

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
    displacement_matrix_global=displacement_matrix_global_initial
    #delta_displacement_matrix_global???

    #DO Newton_Raphson_method
        k=1
        stiffness_matrix_global=np.zeros([number_of_elements+1,number_of_elements+1])
        #G=???
        for j in range(number_of_elements):
            assignment_matrix=ASSIGNMENT_MATRIX(j,number_of_elements)
            displacement_matrix_element=np.matmul(assignment_matrix,displacement_matrix_global)
            ###  call Element routine with 
            ###  input displacement_matrix_element,tau,co-ordinates of nodes
            ###  output stiffness_matrix_element,force_internal_element,force_external_element
            stiffness_matrix_global=stiffness_matrix_global+np.matmul(np.transpose(assignment_matrix),np.matmul(stiffness_matrix_element,assignment_matrix))
            G=G+np.matmul(np.transpose(assignment_matrix),(force_internal_element-force_external_element))
        Implementation of essential boundary conditions in K and G
        solve K.delta_u=-G 
        displacement_matrix_global=displacement_matrix_global+delta_u
        Implementation of essential boundary conditions in displacement_matrix_global
        k=k+1
        if k>k_max restart increment with smaller delta_t
    #While norm(G)>   , norm(delta_u)>   

