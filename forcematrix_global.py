import numpy as np

##Fucntion ASSIGNMENT_MATRIX will be in main program
from assignment_matrix import ASSIGNMENT_MATRIX


##Function definition for External force, Have to be taken from element routine##
#from ELEMENT_FORCE_MATRIX element_force_matrix




def force_external(element_number,number_of_elements,??????):
    force_external=np.zeros([number_of_elements+1,1])
    for element_number in range(1,number_of_elements+1):
        #ASSIGNMENT_MATRIX from main program
        assignment_matrix=ASSIGNMENT_MATRIX((element_number),number_of_elements)
        #ELEMENT_FORCE_MATRIX from element routine
        element_force_matrix=ELEMENT_FORCE_MATRIX(??????)
        force_external=force_external+(np.matmul(np.transpose(assignment_matrix),ef))
    return force_external