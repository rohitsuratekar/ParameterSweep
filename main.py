from analysis.mutants.mutant_check import calculate_mutant
from analysis.mutants.visualize import visualize
from constants.namespace import *

system = S_ONLY_FORWARD


def calculate():
    calculate_mutant('input.txt', system, MICHAELIS_MENTEN)


def vis():
    visualize('output/output.log', system)


vis()
