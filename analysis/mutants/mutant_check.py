"""
Calculates mutant ratios from the parameter set obtained from parameter sweep
"""

from analysis.helper import *
from utils.functions import update_progress
from utils.log import LOG, CURRENT_JOB


def convert_to_enzyme(para) -> dict:
    converted_para = {}
    for key in para.keys():
        converted_para[key] = Enzyme.make_with_values(key, para[key])
    return converted_para


def get_parameter_set(filename) -> list:
    parameters = []
    with open(filename, "r") as f:
        for line in f:
            parameters.append(extract_enz_from_log(line))
    return parameters


def get_pa_ratio(wt, mt):
    return ((mt[4] + mt[5]) / (mt[0] + mt[7])) / ((wt[4] + wt[5]) / (wt[0] + wt[7]))


def get_dag_ratio(wt, mt):
    return (mt[3] / (mt[0] + mt[7])) / (wt[3] / (wt[0] + wt[7]))


def calculate_mutant(filename: str, system: str, kinetics: str):
    log_data = {
        "UID": CURRENT_JOB,
        "System": system,
        "Kinetics": kinetics,
        "TotalLipid": total_lipid_concentration,
        "Analysis": "Calculating Mutant",
        "ODEEndTime": ode_end_time,
        "version": "3.0"}
    LOG.info(json.dumps(log_data, sort_keys=True))

    all_para = get_parameter_set(filename)
    progress_counter = 0
    update_progress(progress_counter / len(all_para), "Calculating mutants")
    for para in all_para:
        enz = convert_to_enzyme(para)
        initial_con = get_random_concentrations(total_lipid_concentration, system)
        wt_output = get_concentration_profile(system, initial_con, enz, ode_end_time, ode_slices)
        enz[E_DAGK].v *= 0.1
        rdga_output = get_concentration_profile(system, initial_con, enz, ode_end_time, ode_slices)
        enz[E_DAGK].v *= 10
        enz[E_LAZA].v *= 0.1
        laza_output = get_concentration_profile(system, initial_con, enz, ode_end_time, ode_slices)
        enz[E_LAZA].v *= 10

        mutant_ratio = {
            "RDGA": {
                "PA": round(get_pa_ratio(wt_output[-1], rdga_output[-1]), 4),
                "DAG": round(get_dag_ratio(wt_output[-1], rdga_output[-1]), 4)
            },
            "LAZA": {
                "PA": round(get_pa_ratio(wt_output[-1], laza_output[-1]), 4),
                "DAG": round(get_dag_ratio(wt_output[-1], laza_output[-1]), 4)
            },
        }

        data = {}
        for value in enz.values():
            data[value.name] = {
                "v": round(value.v, 4),
                "k": round(value.k, 4),
                "kinetics": value.kinetics
            }
        save_values = {
            "Mutants": mutant_ratio,
            "Enzymes": data
        }
        OUTPUT.info(json.dumps(save_values, sort_keys=True))
        progress_counter += 1
        update_progress(progress_counter / len(all_para), "Calculating mutants")
