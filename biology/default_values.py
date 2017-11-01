"""
Default values used in all analysis and simulations.
User can change them here and they will be globally applied.
Only "Total Lipid Level" will be included into logs
"""

# Parameter related
default_v = 1
default_k = 10

min_v, max_v = 0.01, 10
min_k, max_k = 1, 200

rand_around_v = 1
rand_around_k = 10

alpha = 0.05  # PIP2
beta = 0.05  # PI4P
gamma = 0.008  # DAG
delta = 0.1677  # PA
epsilon = 0.001  # CDPDAG

outer_iterations = 50000
inner_iterations = 20000
para_skip_threshold = 200
save_cutoff = 0.15

total_lipid_concentration = 200

# ODE settings
ode_end_time = 50000
ode_slices = ode_end_time * 10
