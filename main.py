from analysis.mutants.experiment import light_flash
from analysis.mutants.mutant_check import calculate_mutant
from analysis.mutants.plot_mutant import plot
from analysis.mutants.visualize import visualize, single_para_sensitivity
from analysis.sweeps.improve_parameter import improve
from analysis.sweeps.mass_action_reverse_calculation import calculate
from constants.namespace import *
from test.convert_to_latex import convert_to_latex

system = S_OPEN_2
kinetics = MASS_ACTION
expression_level = 0.1


def cal():
    calculate_mutant('input.txt', system, kinetics, [expression_level])


def vis():
    visualize('output/output.log', system, expression_level)


def plot_mutant():
    plot(system, expression_level)


def exp():
    light_flash(system)


def single():
    single_para_sensitivity(system)


def imp():
    improve(system, kinetics, 'analysis/mutants/para.txt', expression_level)


def latex():
    convert_to_latex('top_para.txt')


def convert_mm():
    calculate(system, kinetics, 0.1)


convert_mm()
