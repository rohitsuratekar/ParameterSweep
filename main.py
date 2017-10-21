from analysis.mutants.experiment import light_flash
from analysis.mutants.mutant_check import calculate_mutant
from analysis.mutants.plot_mutant import plot
from analysis.mutants.visualize import visualize_sensitivity
from constants.namespace import *

system = S_OPEN_2


def calculate():
    calculate_mutant('input.txt', system, MICHAELIS_MENTEN, [0.1, 0.3, 0.5,
                                                             0.7, 1])


def vis():
    visualize_sensitivity('output/output.log', system)


def plot_mutant():
    plot(system)


def exp():
    light_flash(system)


vis()
