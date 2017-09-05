from utils.functions import update_progress
from .helper import *


# Creates n different types of parameter values
def create_population(size: int, kinetics: str) -> list:
    population = []
    for i in range(size):
        population.append(get_random_enzymes(kinetics))
    return population


# We will multiply error with -1 so that lowest error will have highest fitness
def calculate_fitness(system: str, individual: dict):
    initial_conditions = get_random_concentrations(total_lipid_concentration, system)
    output = get_concentration_profile(system, initial_conditions, individual, ode_end_time, ode_slices)
    return Error(total_lipid_concentration, list(output[-1]), individual).total_error


# Randomly swap parameters between two sets
def reproduce(parent1: dict, parent2: dict, mutation_rate) -> dict:
    for i in list(parent1.keys()):
        if np.random.choice([True, False]):
            parent1[i], parent2[i] = parent2[i], parent1[i]
    if np.random.choice([True, False], p=[mutation_rate, 1 - mutation_rate]):
        parent1[np.random.choice(list(parent1.keys()))].randomize()
    return parent1


# Create new parameter set
def renew_population(system: str, population: list, mutation_rate: float):
    fitness_matrix = []
    for p in population:
        fitness_matrix.append(calculate_fitness(system, p))
    min_error = min(fitness_matrix)
    fitness_matrix = [sum(fitness_matrix) - x for x in fitness_matrix]
    fitness_matrix = [x / sum(fitness_matrix) for x in fitness_matrix]

    # Pick 2 highest probable parents without replacement with probabilities as fitness matrix
    new_generation = []
    for i in range(len(population)):
        p1, p2 = np.random.choice(fitness_matrix, 2, False, p=fitness_matrix)
        new_generation.append(
            reproduce(population[fitness_matrix.index(p1)], population[fitness_matrix.index(p2)], mutation_rate))

    return new_generation, min_error


def do_sweep(system: str, kinetics: str, population: int, mutation_rate: float):
    pop = create_population(population, kinetics)
    lowest_error = 100000
    progress_counter = 0
    update_progress(progress_counter / outer_iterations,
                    "lowest_error: %s current_error : %s" % (lowest_error, lowest_error))
    for i in range(outer_iterations):
        pop, current_fitness = renew_population(system, pop, mutation_rate)
        if current_fitness < lowest_error:
            lowest_error = current_fitness

        progress_counter += 1
        update_progress(progress_counter / outer_iterations,
                        "lowest_error: %s current_error : %s" % (lowest_error, current_fitness))
