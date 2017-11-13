"""
Do light flash experiment with given parameter set in "para.txt"
"""

from itertools import product

import matplotlib.pylab as plt

from analysis.best_parameters import get_open2
from analysis.helper import *

lipid_color_pellet = ["#9C27B0", "#F44336", "#FFEB3B", "#2196F3",
                      "#009688", "#4CAF50", "#795548", "#00BCD4"]

lipid_names = [L_PMPI, L_PI4P, L_PIP2, L_DAG, L_PMPA, L_ERPA, L_CDPDAG, L_ERPI]
recovery_time = 100


def get_parameter_set() -> dict:
    # with open("analysis/mutants/para.txt", "r") as f:
    #     for line in f:
    #         raw = extract_enz_from_log(line)
    #         para = {}
    #         for key in raw.keys():
    #             para[key] = Enzyme.make_with_values(key, raw[key])
    #         return para

    return get_open2(MICHAELIS_MENTEN)


def get_stimulus_profile(system, enz, stimulus):
    initial_time = 10
    stimulus_time = 1

    final_time = []
    initial_con = get_random_concentrations(total_lipid_concentration,
                                            system)

    wt_output = get_concentration_profile(system, initial_con, enz,
                                          ode_end_time, ode_slices)

    ss_con = wt_output[-1]

    ss_output = get_concentration_profile(system, ss_con, enz, initial_time,
                                          initial_time * 100)

    time = np.linspace(0, initial_time, initial_time * 100)

    final_time.extend(time)
    final_concentration = ss_output

    enz.get(E_PLC).v *= stimulus

    ss_stimulated = get_concentration_profile(system, ss_output[-1], enz,
                                              stimulus_time,
                                              stimulus_time * 100)

    ss_time = np.linspace(final_time[-1], final_time[-1] + stimulus_time,
                          stimulus_time * 100)

    final_time.extend(ss_time)
    final_concentration = np.concatenate((final_concentration, ss_stimulated),
                                         axis=0)
    enz.get(E_PLC).v *= 1 / stimulus

    ss_recovery = get_concentration_profile(system, ss_stimulated[-1], enz,
                                            recovery_time,
                                            recovery_time * 100)

    rec_time = np.linspace(final_time[-1], final_time[-1] + recovery_time,
                           recovery_time * 100)

    final_time.extend(rec_time)
    final_concentration = np.concatenate((final_concentration, ss_recovery),
                                         axis=0)

    return final_concentration, final_time, ss_output[-1], ss_recovery


def plot_lipids(final_concentration, final_time, lipids: list, ss_level,
                with_label=""):
    for i in lipids:
        if with_label == "":
            plt.plot(final_time,
                     final_concentration[:, lipid_names.index(i)] / ss_level,
                     label=i,
                     color=lipid_color_pellet[lipid_names.index(i)],
                     linewidth=2)
        else:
            if len(lipids) == 1:
                plt.plot(final_time,
                         final_concentration[:,
                         lipid_names.index(i)] / ss_level,
                         linewidth=2, label=with_label)
            else:
                plt.plot(final_time,
                         final_concentration[:,
                         lipid_names.index(i)] / ss_level,
                         linewidth=2, label=with_label + " " + i)


def calculate_area(con_wt, con_mt, lipid, wt_ss, mt_ss):
    wt = con_wt[:, lipid_names.index(lipid)] / wt_ss
    mt = con_mt[:, lipid_names.index(lipid)] / mt_ss
    return np.trapz(wt) - np.trapz(mt)


def light_flash(system: str, factor, lipid):
    enz = get_parameter_set()

    # Feedback Settings
    enz[E_PIP5K].feedback_substrate = lipid[0]
    enz[E_PIP5K].feedback_factor = factor[0]

    enz[E_PITP].feedback_substrate = lipid[1]
    enz[E_PITP].feedback_factor = factor[1]

    final_concentration, final_time, ss_level_wt, rec_wt = get_stimulus_profile(
        system,
        enz, 50)

    plot_lipids(final_concentration, final_time, [L_PIP2],
                ss_level_wt[lipid_names.index(L_PIP2)], with_label="WT")

    enz.get(E_PIP5K).v *= 0.1
    final_concentration, final_time, ss_level_mt, rec_mt = get_stimulus_profile(
        system,
        enz, 50)
    plot_lipids(final_concentration, final_time, [L_PIP2],
                ss_level_mt[lipid_names.index(L_PIP2)], with_label="MT")

    calculate_area(rec_wt, rec_mt, L_PIP2,
                   ss_level_wt[lipid_names.index(L_PIP2)],
                   ss_level_mt[lipid_names.index(L_PIP2)])
    plt.legend(loc=0)
    plt.show()


def find_lowest_para(system: str):
    lowest_area = 0
    for l in product(lipid_names, lipid_names):
        for factor in product(range(15), range(15)):
            # Feedback Settings
            if l[0] == l[1]:
                break

            enz = get_parameter_set()
            enz[E_PIP5K].feedback_substrate = l[0]
            enz[E_PIP5K].feedback_factor = factor[0]

            enz[E_PITP].feedback_substrate = l[1]
            enz[E_PITP].feedback_factor = factor[1]

            final_concentration, final_time, ss_level_wt, rec_wt = get_stimulus_profile(
                system,
                enz, 50)

            enz.get(E_PIP5K).v *= 0.1
            final_concentration, final_time, ss_level_mt, rec_mt = get_stimulus_profile(
                system,
                enz, 50)

            area = calculate_area(rec_wt, rec_mt, L_PIP2,
                                  ss_level_wt[lipid_names.index(L_PIP2)],
                                  ss_level_mt[lipid_names.index(L_PIP2)])

            if area > lowest_area and area > 0:
                lowest_area = area
                print(area, l, factor)


def do_analysis(system):
    light_flash(system, [13, 0], [L_PI4P, L_ERPI])
    #find_lowest_para(system)
