"""
Based on Nelder-Mead Simplex method of optimization.
Best-Good-Worst sets are altered differently than that of original method.
"""

from utils.functions import update_progress
from utils.log import LOG, CURRENT_JOB
from analysis.helper import *


# Sort parameter set in order of their error
def sort_best_set(initial_condition, system: str, s1, s2, s3, total_lipid) -> list:
    o1 = get_concentration_profile(system, initial_condition, s1, ode_end_time, ode_slices)
    o2 = get_concentration_profile(system, initial_condition, s2, ode_end_time, ode_slices)
    o3 = get_concentration_profile(system, initial_condition, s3, ode_end_time, ode_slices)
    e1 = Error(total_lipid, o1[-1], s1)
    e2 = Error(total_lipid, o2[-1], s2)
    e3 = Error(total_lipid, o3[-1], s3)
    final_list = [e1, e2, e3]
    final_list.sort(key=lambda x: x.total_error)
    return final_list


# Create initial random parameter set
def make_initial_set(initial_condition, system: str, kinetics: str, total_lipid) -> list:
    s1 = get_random_enzymes(kinetics)
    s2 = get_random_enzymes(kinetics)
    s3 = get_random_enzymes(kinetics)
    return sort_best_set(initial_condition, system, s1, s2, s3, total_lipid)


def mid_point(best, good) -> dict:
    best_enzyme = best.enzymes
    good_enzyme = good.enzymes
    mid_enzyme = {}
    key_rand = np.random.choice(list(best_enzyme.keys()))
    for key in best_enzyme.keys():
        current_enz1 = best_enzyme[key]  # type: Enzyme
        current_enz2 = good_enzyme[key]  # type: Enzyme
        if key == key_rand:
            vmax = abs((current_enz1.v - current_enz2.v) / 2)
            km = abs((current_enz1.k - current_enz2.k) / 2)
        else:
            vmax = current_enz1.v
            km = current_enz1.k

        mid_enzyme[key] = Enzyme(key, v=vmax, k=km, kinetics=current_enz1.kinetics)

    return mid_enzyme


def get_reflection(mid, worst):
    mid_enzyme = mid.enzymes
    worst_enzyme = worst.enzymes
    reflected_enzyme = {}
    key_rand = np.random.choice(list(mid_enzyme.keys()))
    for key in mid_enzyme.keys():
        current_enz1 = mid_enzyme[key]  # type: Enzyme
        current_enz2 = worst_enzyme[key]  # type: Enzyme
        if key == key_rand:
            vmax = abs(2 * current_enz1.v - current_enz2.v)
            km = abs(2 * current_enz1.k - current_enz2.k)
        else:
            vmax = current_enz1.v
            km = current_enz1.k
        reflected_enzyme[key] = Enzyme(key, v=vmax, k=km, kinetics=current_enz1.kinetics)

    return reflected_enzyme


def get_expansion(reflected, mid):
    return get_reflection(reflected, mid)


def do_sweep(system: str, kinetics: str, total_lipid=total_lipid_concentration):
    # Initial setup to start sweep
    log_data = {
        "UID": CURRENT_JOB,
        "System": system,
        "Kinetics": kinetics,
        "TotalLipid": total_lipid,
        "Analysis": "Nelder-Mead Simplex Sweep",
        "version": "3.0"}
    LOG.info(json.dumps(log_data, sort_keys=True))

    lowest_error = 1000000
    progress_counter = 0
    update_progress(progress_counter / outer_iterations, "lowest_error : %s" % lowest_error)
    for i in range(outer_iterations):
        initial_condition = get_random_concentrations(total_lipid, system)
        best, good, worst = make_initial_set(initial_condition, system, kinetics, total_lipid)
        para_skip = 0
        while para_skip < para_skip_threshold:
            mid = mid_point(best, good)
            mid_out = get_concentration_profile(system, initial_condition, mid, ode_end_time, ode_slices)
            mid_error = Error(total_lipid, mid_out[-1], mid)
            reflection = get_reflection(mid_error, worst)
            reflection_out = get_concentration_profile(system, initial_condition, reflection, ode_end_time, ode_slices)
            reflection_error = Error(total_lipid, reflection_out[-1], reflection)

            analysis_list = [best, good, worst, mid_error, reflection_error]
            before_sort = best
            analysis_list.sort(key=lambda x: x.total_error)
            best, good, worst = analysis_list[0], analysis_list[1], analysis_list[2]

            if before_sort == best:
                para_skip += 1
            if best.total_error < lowest_error:
                lowest_error = best.total_error
                best.record()

        progress_counter += 1
        update_progress(progress_counter / outer_iterations, "lowest_error : %s" % lowest_error)
