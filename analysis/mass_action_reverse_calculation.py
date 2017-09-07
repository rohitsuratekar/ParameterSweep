"""
Calculates parameters for Mass Action kind of reactions from analytical solution
"""

from utils.functions import update_progress
from utils.log import LOG, CURRENT_JOB
from .helper import *


def calculate_steady_states(parameters):
    pitp, pi4k, pip5k, plc, dagk, sink, laza, patp, source, cds, pis = parameters
    dag = source / sink
    pip2 = dag * ((dagk * patp / (laza + patp)) + sink) / plc
    pmpa = pip2 * plc / (patp + ((laza + patp) * sink / dagk))
    erpa = pip2 * plc / cds
    cdpdag = pip2 * plc / pis
    erpi = pip2 * plc / pitp
    pmpi = pip2 * plc / pi4k
    pi4p = pip2 * plc / pip5k
    return pmpi, pi4p, pip2, dag, pmpa, erpa, cdpdag, erpi


def calculate_wild_type_error(parameters):
    pmpi, pi4p, pip2, dag, pmpa, erpa, cdpdag, erpi = parameters
    pi = pmpi + erpi
    e1 = abs((pip2 / pi) - alpha) / alpha
    e2 = abs((pi4p / pi) - beta) / beta
    e3 = abs((dag / pi) - gamma) / gamma
    e4 = abs(((erpa + pmpa) / pi) - delta) / delta
    e5 = abs((cdpdag / pi) - epsilon) / epsilon
    return e1 + e2 + e3 + e4 + e5


def get_pa(para):
    return (para[4] + para[5]) / (para[0] + para[7])


def get_dag(para):
    return para[3] / (para[0] + para[7])


def laza_pa_level_correction(error):
    if error < 2:
        return 5
    else:
        return 0


def error_with_tolerance(mt, wt, expected, percentage_tol=5):
    # We will give 5% tolerance
    e1 = abs((mt / wt) - expected)
    e2 = abs((mt / wt) - (expected + (expected * percentage_tol / 100)))
    e3 = abs((mt / wt) - (expected - (expected * percentage_tol / 100)))
    return min(e1, e2, e3)


def calculate_rdga3_error(wt, rdga3):
    wt_pa_ratio = get_pa(wt)
    rdga3_pa_ratio = get_pa(rdga3)

    wt_dag_ratio = get_dag(wt)
    rdga3_dag_ratio = get_dag(rdga3)

    return error_with_tolerance(rdga3_pa_ratio, wt_pa_ratio, 1) + error_with_tolerance(rdga3_dag_ratio, wt_dag_ratio, 1)


def calculate_laza22_error(wt, laza22):
    wt_pa_ratio = get_pa(wt)
    laza22_pa_ratio = get_pa(laza22)

    wt_dag_ratio = get_dag(wt)
    laza22_dag_ratio = get_dag(laza22)

    return error_with_tolerance(laza22_pa_ratio, wt_pa_ratio, 2.4, 10) + error_with_tolerance(laza22_dag_ratio,
                                                                                              wt_dag_ratio, 1)


def calculate_total_error(wt, rdga3, laza22):
    e = calculate_rdga3_error(wt, rdga3) + calculate_laza22_error(wt, laza22)
    return e


def save_parameters(parameters, error):
    epitp, epi4k, epip5k, eplc, edagk, esink, elaza, epatp, esource, ecds, epis = parameters
    kinetics = MASS_ACTION

    pitp = Enzyme(E_PITP, kinetics=kinetics, k=epitp)
    pi4k = Enzyme(E_PI4K, kinetics=kinetics, k=epi4k)
    pip5k = Enzyme(E_PIP5K, kinetics=kinetics, k=epip5k)
    plc = Enzyme(E_PLC, kinetics=kinetics, k=eplc)
    dagk = Enzyme(E_DAGK, kinetics=kinetics, k=edagk)
    laza = Enzyme(E_LAZA, kinetics=kinetics, k=elaza)
    patp = Enzyme(E_PATP, kinetics=kinetics, k=epatp)
    cds = Enzyme(E_CDS, kinetics=kinetics, k=ecds)
    pis = Enzyme(E_PIS, kinetics=kinetics, k=epis)
    sink = Enzyme(E_SINK, kinetics=kinetics, k=esink)
    source = Enzyme(E_SOURCE, kinetics=kinetics, k=esource)
    p4tase = Enzyme(E_P4TASE, kinetics=kinetics, k=default_k)
    p5tase = Enzyme(E_P5TASE, kinetics=kinetics, k=default_k)
    ip3tase = Enzyme(E_IP3_PTASE, kinetics=kinetics, k=default_k)

    all_enz = {x.name: x for x in
               [pitp, pi4k, pip5k, plc, dagk, laza, patp, cds, pis, sink, source, p4tase, p5tase, ip3tase]}

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


def calculate(system: str, kinetics: str):
    if system != S_OPEN_2 and kinetics != MASS_ACTION:
        raise Exception(" This analysis is only for OPEN Cycle 2 Mass action reactions")

    log_data = {
        "UID": CURRENT_JOB,
        "System": system,
        "Kinetics": kinetics,
        "TotalLipid": "N/A",
        "Analysis": "Mass Action Reverse Calculations",
        "version": "3.0"}
    LOG.info(json.dumps(log_data, sort_keys=True))

    lowest_error = 100000
    progress_counter = 0
    update_progress(progress_counter / outer_iterations, "lowest_error : %s" % lowest_error)
    for i in range(outer_iterations):
        B, E = np.random.uniform(0.01, 20, 2)
        A = np.random.uniform(alpha, 20)
        D = np.random.uniform(0.01, alpha / gamma)
        C = np.random.uniform((alpha - (gamma * D)) / delta, 20)
        plc = 7.7
        pitp = A * plc
        laza = B * plc
        patp = C * plc
        sink = D * plc
        source = E * plc
        pi4k = (alpha / (1 - (alpha / A))) * plc
        pip5k = alpha * plc / beta
        dagk = ((alpha / gamma) - D) * (1 + (B / C)) * plc
        pis = alpha * plc / epsilon
        cds = C * alpha * plc / ((C * delta) - alpha + (D * gamma))

        wt_ss = calculate_steady_states([pitp, pi4k, pip5k, plc, dagk, sink, laza, patp, source, cds, pis])
        dagk *= 0.1
        rdga3_ss = calculate_steady_states([pitp, pi4k, pip5k, plc, dagk, sink, laza, patp, source, cds, pis])
        dagk *= 10
        laza *= 0.1
        laza_ss = calculate_steady_states([pitp, pi4k, pip5k, plc, dagk, sink, laza, patp, source, cds, pis])

        current_error = calculate_total_error(wt_ss, rdga3_ss, laza_ss)

        if current_error < lowest_error:
            lowest_error = current_error
            if current_error < save_cutoff:
                laza *= 10
                save_parameters([pitp, pi4k, pip5k, plc, dagk, sink, laza, patp, source, cds, pis], current_error)

        progress_counter += 1
        update_progress(progress_counter / outer_iterations, "lowest_error : %s" % lowest_error)
