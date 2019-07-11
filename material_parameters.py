''' Input Parameters '''

# Mesh Generation Parameters
ri = 10e-3                         # radius of inclusion (mm)
ro = 40e-3                         # outer radius (mm)
nElem = 10                         # number of elements
meshRefinementFactor = 5           # ratio of element sizes at outer and inner radius

#Material properties
youngs_modulus = 70000          # MPa
poissons_ratio = 0.25
yield_stress = 70               # MPa
E_v = 0.01

# Time, Load Parameters
delta_t = 0.05
tau_start = 0.0
tau_end = 1.0
