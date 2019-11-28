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


def print_normalized_values(enzymes, wt):
    for k in enzymes:
        if k != E_SOURCE:
            print("%s\t%.3f\t%.3f" % (
                k, enzymes[k].v / enzymes[E_PLC].v, enzymes[k].k / wt))
        else:
            print("%s\t%.3f\t%.3f" % (
                k, enzymes[k].v, enzymes[k].k / enzymes[E_PLC].v))


def print_regular(enzymes):
    for k in enzymes:
        print("%s\t%.3f\t%.3f" % (
            k, enzymes[k].v, enzymes[k].k))


def print_rdga(val):
    print("RDGA")
    print("PA Ratio: %.3f (error: %.3f)" % (val[0], (val[0] - 1) * 100 / 1))
    print("DAG Ratio: %.3f (error: %.3f)" % (val[1], (val[1] - 1) * 100 / 1))


def print_laza(val):
    print("LAZA")
    print("PA Ratio: %.3f (error: %.3f)" % (val[0], (val[0] - 2.5) * 100 / 1))
    print("DAG Ratio: %.3f (error: %.3f)" % (val[1], (val[1] - 1) * 100 / 1))


def plot(system: str, mutant_expression: float):
    enz = get_parameter_set()
    # enz = get_open2(MICHAELIS_MENTEN)

    # exp_pi = 3.5868
    # exp_pa = 0.6016
    # exp_dag = 0.0288
    #
    # total_lipid_before = 207.2753
    #
    # pf = exp_pi / 0.7777
    # multi_factor = 1 / total_lipid_before
    # multi_factor = multi_factor * pf
    #
    # temp = enz[E_PLC].v
    # for k in enz:
    #     if k != E_SOURCE:
    #         enz[k].v = enz[k].v / temp
    #         enz[k].k = enz[k].k * 6.22
    #     else:
    #         enz[k].k = enz[k].k / temp

    initial_con = get_random_concentrations(1, system)
    # initial_con = [0.6633562, 0.03816765, 0.03842087, 0.006625219,
    # 0.01672982, 0.1172366, 0.0007857711, 0.1169003] initial_con = [3.737,
    # 0.109, 0.102, 0.076, 0.803, 0.106, 0.045, 1.299]

    wt_output = get_concentration_profile(system, initial_con, enz,
                                          ode_end_time, ode_slices)

    s = ""
    for i in wt_output[-1]:
        s += str(round(i, 3)) + "\t"

    print("PI4P", wt_output[-1][1] / (wt_output[-1][0] + wt_output[-1][7]))
    print("PIP2", wt_output[-1][2] / (wt_output[-1][0] + wt_output[-1][7]))
    print("DAG", wt_output[-1][3] / (wt_output[-1][0] + wt_output[-1][7]))
    print("PA", (wt_output[-1][4] + wt_output[-1][5]) / (
            wt_output[-1][0] + wt_output[-1][7]))

    print(["%.4f" % x for x in wt_output[-1]])
    # print(sum(wt_output[-1]), wt_output[-1][0] + wt_output[-1][7])
    # print(pf, sum(wt_output[-1]),
    #       wt_output[-1][0] + wt_output[-1][7])

    enz[E_DAGK].mutate(mutant_expression)
    rdga_output = get_concentration_profile(system, initial_con, enz,
                                            ode_end_time, ode_slices)
    enz[E_DAGK].mutate(1 / mutant_expression)
    enz[E_LAZA].mutate(mutant_expression)
    laza_output = get_concentration_profile(system, initial_con, enz,
                                            ode_end_time, ode_slices)
    enz[E_LAZA].mutate(1 / mutant_expression)

    fig = plt.figure()
    gs = gridspec.GridSpec(8, 5)
    wt = [1, 1]
    rdga = [get_pa_ratio(wt_output[-1], rdga_output[-1]),
            get_dag_ratio(wt_output[-1], rdga_output[-1])]
    laza = [get_pa_ratio(wt_output[-1], laza_output[-1]),
            get_dag_ratio(wt_output[-1], laza_output[-1])]

    print_rdga(rdga)
    print_laza(laza)

    # print_normalized_values(enz, sum(wt_output[-1]))
    # print_regular(enz)

    ind = np.arange(2)
    bar_width = 0.35
    ax2 = fig.add_subplot(gs[0:3, 1:4])

    ax2.bar(ind, wt, bar_width,
            # color='g',
            color="#9ECC3B",
            label='WT')

    ax2.bar(ind + bar_width, rdga, bar_width,
            # color='r',
            color="#4A2261",
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
            # color='g',
            color="#9ECC3B",
            label='WT')

    ax3.bar(ind + bar_width, laza, bar_width,
            # color='r',
            color="#4A2261",
            label='Mutant')
    ax3.set_ylim(0.5, 3)
    ax3.axhline(1, linestyle="--")
    ax3.set_xticks(ind + bar_width / 2)
    ax3.tick_params(axis=u'both', which=u'both', length=0)
    ax3.set_xticklabels(('$PA_{total}/PI_{total}$', '$DAG/PI_{total}$'))
    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    ax3.set_ylabel("Fold change")
    ax3.set_title("$laza^{22}$")
    plt.savefig("mutant_para.png", format='png', dpi=300,
                bbox_inches='tight', transparent=True)
    # plt.show()


def k_sensitivity(mutant_expression: float):
    system = S_OPEN_2
    enz = get_parameter_set()

    initial_con = get_random_concentrations(1, system)

    s_array = np.linspace(0.1, 10, 30).tolist()
    s_array.append(1)
    pa_level = []
    print(enz.get(E_SINK).v)
    for f in s_array:
        enz.get(E_SOURCE).k *= f

        if f == 1:
            print(enz.get(E_SOURCE).k)

        wt_output = get_concentration_profile(system, initial_con, enz,
                                              ode_end_time, ode_slices)

        enz[E_DAGK].mutate(mutant_expression)
        rdga_output = get_concentration_profile(system, initial_con, enz,
                                                ode_end_time, ode_slices)
        enz[E_DAGK].mutate(1 / mutant_expression)
        enz[E_LAZA].mutate(mutant_expression)
        laza_output = get_concentration_profile(system, initial_con, enz,
                                                ode_end_time, ode_slices)
        enz[E_LAZA].mutate(1 / mutant_expression)

        rdga = [get_pa_ratio(wt_output[-1], rdga_output[-1]),
                get_dag_ratio(wt_output[-1], rdga_output[-1])]
        laza = [get_pa_ratio(wt_output[-1], laza_output[-1]),
                get_dag_ratio(wt_output[-1], laza_output[-1])]

        pa_level.append(get_pa_ratio(wt_output[-1], laza_output[-1]))

        enz.get(E_SOURCE).k /= f

    plt.scatter(s_array, pa_level)
    plt.xlabel("Factor (x $k_{source}$)")
    plt.ylabel("PA/Pi ratio in LAZA")
    plt.title(system)
    plt.axvline(1)
    plt.axhline(2.5)
    plt.show()


def lipid_sensitivity(mutant_expression: float):
    system = TEST_SYSTEM
    enz = get_parameter_set()

    s_array = np.linspace(0.1, 10, 30).tolist()
    s_array.append(1)
    pa_level = []

    for f in s_array:
        initial_con = get_random_concentrations(f, system)
        wt_output = get_concentration_profile(system, initial_con, enz,
                                              ode_end_time, ode_slices)

        enz[E_DAGK].mutate(mutant_expression)
        rdga_output = get_concentration_profile(system, initial_con, enz,
                                                ode_end_time, ode_slices)
        enz[E_DAGK].mutate(1 / mutant_expression)
        enz[E_LAZA].mutate(mutant_expression)
        laza_output = get_concentration_profile(system, initial_con, enz,
                                                ode_end_time, ode_slices)
        enz[E_LAZA].mutate(1 / mutant_expression)

        rdga = [get_pa_ratio(wt_output[-1], rdga_output[-1]),
                get_dag_ratio(wt_output[-1], rdga_output[-1])]
        laza = [get_pa_ratio(wt_output[-1], laza_output[-1]),
                get_dag_ratio(wt_output[-1], laza_output[-1])]

        pa_level.append(get_pa_ratio(wt_output[-1], laza_output[-1]))

    plt.scatter(s_array, pa_level)
    plt.xlabel("Total Lipid content ")
    plt.ylabel("PA/Pi ratio in LAZA")
    plt.title(system)
    plt.axvline(1)
    plt.axhline(2.5)
    plt.show()


def test_check(mutant_expression: float):
    enz = get_parameter_set()

    for system in [S_OPEN_2, TEST_SYSTEM]:
        initial_con = get_random_concentrations(6.2, system)
        wt_output = get_concentration_profile(system, initial_con, enz,
                                              ode_end_time, ode_slices)

        enz[E_DAGK].mutate(mutant_expression)
        rdga_output = get_concentration_profile(system, initial_con, enz,
                                                ode_end_time, ode_slices)
        enz[E_DAGK].mutate(1 / mutant_expression)
        enz[E_LAZA].mutate(mutant_expression)
        laza_output = get_concentration_profile(system, initial_con, enz,
                                                ode_end_time, ode_slices)
        enz[E_LAZA].mutate(1 / mutant_expression)

        rdga = [get_pa_ratio(wt_output[-1], rdga_output[-1]),
                get_dag_ratio(wt_output[-1], rdga_output[-1])]
        laza = [get_pa_ratio(wt_output[-1], laza_output[-1]),
                get_dag_ratio(wt_output[-1], laza_output[-1])]

        # print(sum(laza_output[-1]))
        # print("WT " + system, wt_output[-1])
        print("WT " + system,
              (wt_output[-1][4] + wt_output[-1][5]) / (wt_output[-1][0] +
                                                       wt_output[-1][7]))
        # print("LAZA " + system, laza_output[-1][4] + laza_output[-1][5],
        #       laza_output[-1][0] + laza_output[-1][7])


def do():
    # k_sensitivity(0.1)
    test_check(0.1)
    # plot(TEST_SYSTEM, 0.1)
