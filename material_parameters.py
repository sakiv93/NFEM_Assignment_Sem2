''' Input Parameters '''

# Mesh Generation Parameters
ri = 10                         # radius of inclusion (um)
ro = 40                         # outer radius (um)
nElem = 10                    # number of elements
meshRefinementFactor = 5        # ratio of element sizes at outer and inner radius

#Material properties
youngs_modulus = 70000e-6       # N/um^2
poissons_ratio = 0.25
yield_stress = 70e-6            # N/um^2
E_v = 0.01

# Time, Load Parameters
delta_t = 0.1
tau_start = 0.0
tau_end = 1.0
