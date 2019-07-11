import numpy as np
from material_parameters import *

#Input parameters
radius_internal = ri
radius_external = ro
number_of_elements = nElem 
#Mesh refinement factor
meshrefinement_factor= meshRefinementFactor

q=meshrefinement_factor**(1/(number_of_elements-1))
first_element=(radius_external-radius_internal)*(1-q)/(1-meshrefinement_factor*q)
rnode = radius_internal

#Function to extract co-ordinates of nodes in global system.
def COORDINATES_NODES_GLOBAL(number_of_elements,rnode,first_element,q):
        rnodes=np.array([rnode])
        for i in range(number_of_elements):
                rnode=rnode+first_element
                rnodes=np.append(rnodes,rnode)
                first_element=first_element*q
        return rnodes

r=COORDINATES_NODES_GLOBAL(number_of_elements,rnode,first_element,q)
rnodes_Transpose=np.array([r])
rnodes=rnodes_Transpose.T