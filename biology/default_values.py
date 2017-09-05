# Parameter related
v = 1
k = 20

min_v, max_v = 0.01, 10
min_k, max_k = 0.1, 100

rand_around_v = 1
rand_around_k = 10

alpha = 0.05  # PIP2
beta = 0.05  # PI4P
gamma = 0.008  # DAG
delta = 0.1677  # PA
epsilon = 0.001  # CDPDAG

outer_iterations = 5000000
inner_iterations = 2000
para_skip_threshold = 200
save_cutoff = 0.15

total_lipid_concentration = 10

# ODE settings
ode_end_time = 50000
ode_slices = 5000
