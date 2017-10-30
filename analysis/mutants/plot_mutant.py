"""
Plots bar graph of RDGA and LAZA mutant given parameter set in "para.txtx" file
"""

import matplotlib.pylab as plt
from matplotlib import gridspec

from analysis.helper import *


def get_pa_ratio(wt, mt):
    return ((mt[4] + mt[5]) / (mt[0] + mt[7])) / (
        (wt[4] + wt[5]) / (wt[0] + wt[7]))


def get_dag_ratio(wt, mt):
    return (mt[3] / (mt[0] + mt[7])) / (wt[3] / (wt[0] + wt[7]))


def get_parameter_set() -> dict:
    with open("analysis/mutants/para.txt", "r") as f:
        for line in f:
            raw = extract_enz_from_log(line)
            para = {}
            for key in raw.keys():
                para[key] = Enzyme.make_with_values(key, raw[key])
            return para


def print_rdga(val):
    print("RDGA")
    print("PA Ratio: %.3f (error: %.3f)" % (val[0], (val[0] - 1) * 100 / 1))
    print("DAG Ratio: %.3f (error: %.3f)" % (val[1], (val[1] - 1) * 100 / 1))


def print_laza(val):
    print("LAZA")
    print("PA Ratio: %.3f (error: %.3f)" % (val[0], (val[0] - 2.5) * 100 / 1))
    print("DAG Ratio: %.3f (error: %.3f)" % (val[1], (val[1] - 1) * 100 / 1))


def plot(system: str):
    enz = get_parameter_set()
    initial_con = get_random_concentrations(total_lipid_concentration, system)
    wt_output = get_concentration_profile(system, initial_con, enz,
                                          ode_end_time, ode_slices)
    enz[E_DAGK].v *= 0.1
    rdga_output = get_concentration_profile(system, initial_con, enz,
                                            ode_end_time, ode_slices)
    enz[E_DAGK].v *= 10
    enz[E_LAZA].v *= 0.1
    laza_output = get_concentration_profile(system, initial_con, enz,
                                            ode_end_time, ode_slices)
    enz[E_LAZA].v *= 10

    fig = plt.figure()
    gs = gridspec.GridSpec(8, 5)
    wt = [1, 1]
    rdga = [get_pa_ratio(wt_output[-1], rdga_output[-1]),
            get_dag_ratio(wt_output[-1], rdga_output[-1])]
    laza = [get_pa_ratio(wt_output[-1], laza_output[-1]),
            get_dag_ratio(wt_output[-1], laza_output[-1])]

    print_rdga(rdga)
    print_laza(laza)

    ind = np.arange(2)
    bar_width = 0.35
    ax2 = fig.add_subplot(gs[0:3, 1:4])

    ax2.bar(ind, wt, bar_width,
            color='g',
            label='WT')

    ax2.bar(ind + bar_width, rdga, bar_width,
            color='r',
            label='$rdgA^3$')
    ax2.set_ylim(0.5, 1.3)
    ax2.axhline(1, linestyle="--")
    ax2.set_xticks(ind + bar_width / 2)
    ax2.tick_params(axis=u'both', which=u'both', length=0)
    ax2.set_xticklabels(('$PA_{total}/PI_{total}$', '$DAG/PI_{total}$'))
    ax2.set_title("$rdgA^3$")
    ax2.set_ylabel("Fold change")

    ax3 = fig.add_subplot(gs[4:7, 1:4])

    ax3.bar(ind, wt, bar_width,
            color='g',
            label='WT')

    ax3.bar(ind + bar_width, laza, bar_width,
            color='r',
            label='Mutant')
    ax3.set_ylim(0.5, 3)
    ax3.axhline(1, linestyle="--")
    ax3.set_xticks(ind + bar_width / 2)
    ax3.tick_params(axis=u'both', which=u'both', length=0)
    ax3.set_xticklabels(('$PA_{total}/PI_{total}$', '$DAG/PI_{total}$'))
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    ax3.set_ylabel("Fold change")
    ax3.set_title("$laza^{22}$")
    plt.show()
