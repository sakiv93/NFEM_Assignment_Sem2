import numpy as np


##Function definition for Element stiffness matrix, Have to be taken from element routine##
#from element_stiffness_matrix import ELEMENT_STIFFNESS_MATRIX


##Fucntion ASSIGNMENT_MATRIX will be in main program
from assignment_matrix import ASSIGNMENT_MATRIX



def GLOBAL_STIFFNESS_MATRIX(number_of_elements,?????):
    global_stiffness_matrix=np.zeros([number_of_elements+1,number_of_elements+1])
    element_stiffness_matrix=ELEMENT_STIFFNESS_MATRIX(??????????)
    for element_number in range(1,number_of_elements+1):
        assignment_matrix=ASSIGNMENT_MATRIX((element_number),number_of_elements)
        global_stiffness_matrix=global_stiffness_matrix+np.matmul(np.transpose(assignment_matrix),np.matmul(element_stiffness_matrix,assignment_matrix))
    return global_stiffness_matrix