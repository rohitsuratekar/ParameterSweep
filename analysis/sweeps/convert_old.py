"""This scripts converts old parameter values (which were calculated based
on analytical expression" into exact ODE equations i.e. by using K_source as
source. Ideally this script will take old parameter, randomize ONLY source
parameters and check for error. """

from analysis.helper import *
from utils.functions import update_progress
from utils.log import LOG, CURRENT_JOB


def calculate_wild_type_error(parameters):
    pmpi, pi4p, pip2, dag, pmpa, erpa, cdpdag, erpi = parameters
    pi = pmpi + erpi
    e1 = abs((pip2 / pi) - alpha) / alpha
    e2 = abs((pi4p / pi) - beta) / beta
    e3 = abs((dag / pi) - gamma) / gamma
    e4 = abs(((erpa + pmpa) / pi) - delta) / delta
    e5 = abs((cdpdag / pi) - epsilon) / epsilon
    return e1 + e2 + e3 + e4 + e5


def save_para(all_enz, error):
    data = {}
    for value in all_enz.values():
        data[value.name] = {
            "v": round(value.v, 4),
            "k": round(value.k, 4),
            "kinetics": value.kinetics
        }
    save_values = {
        "Error": round(error, 4),
        "Enzymes": data
    }
    OUTPUT.info(json.dumps(save_values, sort_keys=True))


def get_parameters():
    all_para = []
    with open("input.txt", "r") as f:
        for line in f:
            parameters = {}
            enzymes = extract_enz_from_log(line)
            for key in enzymes.keys():
                parameters[key] = Enzyme.make_with_values(key, enzymes[key])

            all_para.append(parameters)

    return all_para


def extract(system: str, kinetics: str):
    log_data = {
        "UID": CURRENT_JOB,
        "System": system,
        "Kinetics": kinetics,
        "TotalLipid": total_lipid_concentration,
        "Analysis": "Converting Old Parameters",
        "version": "3.0"}
    LOG.info(json.dumps(log_data, sort_keys=True))
    progress_counter = 0
    all_para = get_parameters()
    for para in all_para:
        initial_con = get_random_concentrations(total_lipid_concentration,
                                                system)
        sanity_counter = 0
        update_progress(progress_counter / len(all_para),
                        "Extracting Old Parameters")
        while sanity_counter < 10000:
            para[E_SOURCE].k = np.random.uniform(min_k, max_k)
            output = get_concentration_profile(system, initial_con, para,
                                               ode_end_time, ode_slices)
            e = calculate_wild_type_error(output[-1])
            if e < save_cutoff:
                save_para(para, e)
                break
            sanity_counter += 1
        progress_counter += 1
