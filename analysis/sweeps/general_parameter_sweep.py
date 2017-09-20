"""General parameter sweep based on randomize Monte-Carlo type estimation.
It is also similar to "Simulated Annealing" but Instead of temperature,
we used constant number of iterations. """

from analysis.helper import *
from utils.functions import update_progress
from utils.log import LOG, CURRENT_JOB


def update_enzymes(enzyme_list: dict) -> None:
    """
    Updates enzyme sweeps. Use this if you want to accept changes to enzyme
    properties after randomization :param enzyme_list: dict of enzyme to be
    udated :return: updated enzymes
    """
    for e in enzyme_list.values():
        e.use_current()


def randomize(enzyme_list: dict):
    enzyme = enzyme_list[np.random.choice(list(enzyme_list.keys()))]
    if enzyme.kinetics == MASS_ACTION:
        enzyme.randomize_k()
    else:
        if np.random.choice([True, False]):
            enzyme.randomize_v()
        else:
            enzyme.randomize_k()


def do_sweep(system: str, kinetics: str,
             total_lipid=total_lipid_concentration):
    # Initial setup to start sweep
    log_data = {
        "UID": CURRENT_JOB,
        "System": system,
        "Kinetics": kinetics,
        "TotalLipid": total_lipid_concentration,
        "Analysis": "General Parameter Sweep",
        "version": "3.0"}
    LOG.info(json.dumps(log_data, sort_keys=True))

    # Start sweep
    progress_counter = 0
    lowest_error = 10000000
    update_progress(progress_counter / outer_iterations,
                    "lowest_error : %s" % lowest_error)
    for i in range(outer_iterations):
        current_error = 10000000
        initial_conditions = get_random_concentrations(total_lipid, system)
        enzymes = get_random_enzymes(kinetics)
        para_skip = 0
        for j in range(inner_iterations):
            output = get_concentration_profile(system, initial_conditions,
                                               enzymes, 20000, 2000)
            error = Error(total_lipid, list(output[-1]), enzymes)
            if error.total_error < current_error:
                error.record()
                current_error = error.total_error
                update_enzymes(enzymes)
            else:
                para_skip += 1
                if para_skip > para_skip_threshold:
                    break
                else:
                    for e in enzymes.values():
                        e.reset()
                    randomize(enzymes)

        if current_error < lowest_error:
            lowest_error = current_error
        progress_counter += 1
        update_progress(progress_counter / outer_iterations,
                        "lowest_error : %s" % lowest_error)
