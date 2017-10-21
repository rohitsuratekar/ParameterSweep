"""
Do light flash experiment with given parameter set in "para.txt"
"""

import matplotlib.pylab as plt

from analysis.helper import *

lipid_color_pellet = ["#9C27B0", "#F44336", "#2196F3", "#00BCD4", "#009688",
                      "#4CAF50", "#FFEB3B", "#795548"]

lipid_names = ["PMPI", "PI4P", "PIP2", "DAG", "PMPA", "ERPA", "CDPDAG", "ERPI"]


def get_parameter_set() -> dict:
    with open("analysis/mutants/para.txt", "r") as f:
        for line in f:
            raw = extract_enz_from_log(line)
            para = {}
            for key in raw.keys():
                para[key] = Enzyme.make_with_values(key, raw[key])
            return para


def light_flash(system: str):
    initial_con = get_random_concentrations(total_lipid_concentration,
                                            system)
    enz = get_parameter_set()

    wt_output = get_concentration_profile(system, initial_con, enz,
                                          ode_end_time, ode_slices)
    time = np.linspace(0, ode_end_time, ode_slices)

    for i in range(len(lipid_names)):
        plt.plot(time, wt_output[:, i], label=lipid_names[i],
                 color=lipid_color_pellet[i], linewidth=2)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
    plt.show()
