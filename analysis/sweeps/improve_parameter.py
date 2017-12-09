"""It improves given parameter values by randoming around single parameter.
This will also check mutant levels"""

from analysis.best_parameters import get_open2
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


def randomize(enzyme_list: dict, at_a_time: int = 1):
    for i in range(at_a_time):
        enzyme = enzyme_list[np.random.choice(list(enzyme_list.keys()))]
        if enzyme.kinetics == MASS_ACTION:
            enzyme.randomize_k()
        else:
            if np.random.choice([True, False]):
                enzyme.randomize_v()
            else:
                enzyme.randomize_k()


class MutantError:
    def __init__(self, system: str, initial_condition: list, enzymes: dict,
                 time: int, mutant_level=0.1, cut_off=save_cutoff):
        self.system = system
        self.initial_condition = initial_condition
        self.time = time
        self.enzymes = enzymes
        self.slices = ode_slices
        self.cut_off = cut_off
        if cut_off > 1:
            self.cut_off = save_cutoff
        self.wt = \
            get_concentration_profile(self.system, self.initial_condition,
                                      self.enzymes, self.time, self.slices)[-1]

        self.enzymes.get(E_DAGK).mutate(mutant_level)

        self.rdga = \
            get_concentration_profile(self.system, self.initial_condition,
                                      self.enzymes, self.time, self.slices)[
                -1]
        self.enzymes.get(E_DAGK).mutate(1 / mutant_level)
        self.enzymes.get(E_LAZA).mutate(mutant_level)
        self.laza = get_concentration_profile(self.system,
                                              self.initial_condition,
                                              self.enzymes, self.time,
                                              self.slices)[
            -1]
        self.enzymes.get(E_LAZA).mutate(1 / mutant_level)

    @property
    def rdga_error(self):
        return abs(self.get_pa_ratio(self.rdga) - 1) + abs(
            self.get_dag_ratio(self.rdga) - 1)

    @property
    def laza_error(self):
        return abs(self.get_pa_ratio(self.laza) - 2.5) / 2.5 + abs(
            self.get_dag_ratio(self.laza) - 1)

    @property
    def wt_error(self):
        return Error(self.initial_condition, self.wt, self.enzymes).total_error

    @property
    def total_error(self):
        return self.rdga_error + self.laza_error + self.wt_error

    def get_pa_ratio(self, mt):
        return ((mt[4] + mt[5]) / (mt[0] + mt[7])) / (
            (self.wt[4] + self.wt[5]) / (self.wt[0] + self.wt[7]))

    def get_dag_ratio(self, mt):
        return (mt[3] / (mt[0] + mt[7])) / (
            self.wt[3] / (self.wt[0] + self.wt[7]))

    def record(self) -> None:
        """
        Records results if satisfying condition is true
        """
        if self.total_error < self.cut_off:
            data = {}
            for value in self.enzymes.values():
                data[value.name] = {
                    "v": round(value.v, 4),
                    "k": round(value.k, 4),
                    "kinetics": value.kinetics
                }
            save_values = {
                "Error": round(self.total_error, 4),
                "Enzymes": data
            }
            OUTPUT.info(json.dumps(save_values, sort_keys=True))


def improve(system: str, kinetics: str, filename: str,
            mutant_expression: float, total_lipid=total_lipid_concentration):
    # Initial setup to start sweep
    log_data = {
        "UID": CURRENT_JOB,
        "System": system,
        "Kinetics": kinetics,
        "TotalLipid": total_lipid_concentration,
        "Analysis": "Improve Parameter",
        "ExpressionLevel": mutant_expression,
        "version": "3.0"}
    LOG.info(json.dumps(log_data, sort_keys=True))

    # Start sweep
    progress_counter = 0
    lowest_error = 10000000
    update_progress(progress_counter / outer_iterations,
                    "lowest_error : %s" % lowest_error)

    time_end = ode_end_time
    if len(get_parameter_set(filename)) > 1:
        LOG.info("Multiple parameters found. Only first parameter will "
                 "be used")
    counter_error = 10000
    for i in range(outer_iterations):
        current_error = 10000000
        initial_conditions = get_random_concentrations(total_lipid, system)
        # enzymes = convert_to_enzyme(get_parameter_set(filename)[0])
        # enzymes = get_random_enzymes(kinetics)
        enzymes = get_open2(MASS_ACTION)
        para_skip = 0
        for j in range(inner_iterations):
            error = MutantError(system, initial_conditions, enzymes,
                                time_end, mutant_level=mutant_expression,
                                cut_off=0.3)
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
                    randomize(enzymes, 4)

        if current_error < lowest_error:
            lowest_error = current_error
        if current_error < counter_error:
            counter_error = current_error
        progress_counter += 1
        update_progress(progress_counter / outer_iterations,
                        "lowest_error : %s" % lowest_error)
