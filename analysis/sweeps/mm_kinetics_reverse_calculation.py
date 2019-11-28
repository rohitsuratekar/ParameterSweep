"""
Generate parameter set for OPEN CYCLE 2 with MM kinetics.
"""

import math

from analysis.helper import *
from utils.functions import update_progress
from utils.log import LOG, CURRENT_JOB


def save_parameters(k_values, v_values, error):
    k_pitp, k_pi4k, k_pip5k, k_plc, k_dagk, k_sink, k_laza, k_patp, k_source, k_cds, k_pis = k_values
    v_pitp, v_pi4k, v_pip5k, v_plc, v_dagk, v_sink, v_laza, v_patp, v_source, v_cds, v_pis = v_values
    kinetics = MICHAELIS_MENTEN

    pitp = Enzyme(E_PITP, kinetics=kinetics, k=k_pitp, v=v_pitp)
    pi4k = Enzyme(E_PI4K, kinetics=kinetics, k=k_pi4k, v=v_pi4k)
    pip5k = Enzyme(E_PIP5K, kinetics=kinetics, k=k_pip5k, v=v_pip5k)
    plc = Enzyme(E_PLC, kinetics=kinetics, k=k_plc, v=v_plc)
    dagk = Enzyme(E_DAGK, kinetics=kinetics, k=k_dagk, v=v_dagk)
    laza = Enzyme(E_LAZA, kinetics=kinetics, k=k_laza, v=v_laza)
    patp = Enzyme(E_PATP, kinetics=kinetics, k=k_patp, v=v_patp)
    cds = Enzyme(E_CDS, kinetics=kinetics, k=k_cds, v=v_cds)
    pis = Enzyme(E_PIS, kinetics=kinetics, k=k_pis, v=v_pis)
    sink = Enzyme(E_SINK, kinetics=kinetics, k=k_sink, v=v_sink)
    source = Enzyme(E_SOURCE, kinetics=kinetics, k=k_source, v=v_source)
    p4tase = Enzyme(E_P4TASE, kinetics=kinetics, k=default_k)
    p5tase = Enzyme(E_P5TASE, kinetics=kinetics, k=default_k)
    ip3tase = Enzyme(E_IP3_PTASE, kinetics=kinetics, k=default_k)

    all_enz = {x.name: x for x in
               [pitp, pi4k, pip5k, plc, dagk, laza, patp, cds, pis, sink,
                source, p4tase, p5tase, ip3tase]}

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


def calculate_wild_type_error(parameters):
    pmpi, pi4p, pip2, dag, pmpa, erpa, cdpdag, erpi = parameters
    pi = pmpi + erpi
    e1 = abs((pip2 / pi) - alpha) / alpha
    e2 = abs((pi4p / pi) - beta) / beta
    e3 = abs((dag / pi) - gamma) / gamma
    e4 = abs(((erpa + pmpa) / pi) - delta) / delta
    e5 = abs((cdpdag / pi) - epsilon) / epsilon
    return e1 + e2 + e3 + e4 + e5


def get_pmpa_related(dag):
    pmpa, v_patp, k_patp, v_laza, k_laza, v_dagk, k_dagk = range(7)
    sanity_counter = 0
    while pmpa == 0:
        v_patp, v_laza, v_dagk = np.random.uniform(min_v, max_v, 3)
        k_patp, k_laza, k_dagk = np.random.uniform(min_k, max_k, 3)

        m = v_dagk * dag / (k_dagk + dag)
        a = (v_patp + v_laza - m)
        b = (v_patp * k_laza + v_laza * k_patp - (k_patp + k_laza) * m)
        c = -1 * m * k_patp * k_laza
        try:
            pmpa = (-b + math.sqrt((b * b) - (4 * a * c))) / (2 * a)
            if pmpa < 0:
                pmpa = 0
        except ValueError:
            pass
        sanity_counter += 1
        if sanity_counter > 10000:
            break
    return pmpa, v_patp, k_patp, v_laza, k_laza, v_dagk, k_dagk


def calculate(system: str, kinetics: str):
    if system != S_OPEN_2 and kinetics != MICHAELIS_MENTEN:
        raise Exception(
            "This analysis is only for OPEN Cycle 2 Michaelis-Menten "
            "reactions")

    log_data = {
        "UID": CURRENT_JOB,
        "System": system,
        "Kinetics": kinetics,
        "TotalLipid": "N/A",
        "Analysis": "Michealis-Menten Reverse Calculations",
        "version": "3.0"}
    LOG.info(json.dumps(log_data, sort_keys=True))

    smallest_error = 10000
    for i in range(outer_iterations):
        v_sink = np.random.uniform(min_v, max_k)
        k_sink = np.random.uniform(min_k, max_k)

        k_source = np.random.uniform(min_v, v_sink)

        dag = (k_sink * k_source) / (v_sink - k_source)

        pmpa, v_patp, k_patp, v_laza, k_laza, v_dagk, k_dagk = get_pmpa_related(
            dag)

        if pmpa > 0:
            m2 = v_patp * pmpa / (k_patp + pmpa)

            v_cds = np.random.uniform(m2, m2 + 10)
            k_cds = np.random.uniform(min_k, max_k)

            erpa = k_cds * m2 / (v_cds - m2)

            m = v_cds * erpa / (k_cds + erpa)

            k_pis, k_pitp, k_pi4k, k_pip5k, k_plc = np.random.uniform(min_k,
                                                                      max_k, 5)
            v_pis, v_pitp, v_pi4k, v_pip5k, v_plc = np.random.uniform(m,
                                                                      m + 10,
                                                                      5)

            cdpdag = k_pis * m / (v_pis - m)
            erpi = k_pitp * m / (v_pitp - m)
            pmpi = k_pi4k * m / (v_pi4k - m)
            pi4p = k_pip5k * m / (v_pip5k - m)
            pip2 = k_plc * m / (v_plc - m)

            e = calculate_wild_type_error(
                [pmpi, pi4p, pip2, dag, pmpa, erpa, cdpdag, erpi])

            if e < smallest_error:
                smallest_error = e
                if smallest_error < save_cutoff:
                    all_k = k_pitp, k_pi4k, k_pip5k, k_plc, k_dagk, k_sink, k_laza, k_patp, k_source, k_cds, k_pis
                    all_v = v_pitp, v_pi4k, v_pip5k, v_plc, v_dagk, v_sink, v_laza, v_patp, 1, v_cds, v_pis
                    save_parameters(all_k, all_v, smallest_error)

        update_progress(i / outer_iterations,
                        "lowest_error : %s" % smallest_error)
