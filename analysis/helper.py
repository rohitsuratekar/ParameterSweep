"""
Helper methods used in all analysis
"""
import json

import numpy as np
from scipy.integrate import odeint

from biology.component import *
from biology.systems.closed_irreversible import get_equations as irreversible
from biology.systems.closed_reversible import get_equations as reversible
from biology.systems.open_cycle_1 import get_equations as open1
from biology.systems.open_cycle_2 import get_equations as open2
from constants.namespace import *
from utils.log import OUTPUT


def get_equations(system: str):
    """
    Returns equation of specific system
    :param system: topology or known model
    :return: set of equation function
    """
    if system == S_OPEN_2:
        return open2
    elif system == S_OPEN_1:
        return open1
    elif system == S_CLOSED_REVERSIBLE:
        return reversible
    elif system == S_ONLY_FORWARD:
        return irreversible
    else:
        raise Exception("No such system found :%s" % system)


def extract_enz_from_log(log_text: str):
    return json.loads(log_text.split(":", 1)[1])["Enzymes"]


def get_concentration_profile(system: str, initial_condition, parameters: dict,
                              ode_time: int, slices: int):
    """
    Solves ODE and returns output array
    :param system: topology or known model
    :param initial_condition: array of initial condition concentrations
    :param parameters: dict of parameter set
    :param ode_time: end time of ODE, initial time will be always 0
    :param slices: number of points to integrate
    :return: Output array
    """
    time = np.linspace(0, ode_time, slices)
    output = odeint(get_equations(system), initial_condition, time,
                    args=(parameters, 0))
    return output


def get_random_concentrations(total: float, system: str) -> list:
    """
    Creates random lipid distribution who's total is equal to given value
    :param system: Type of Cycle
    :param total: total amount of lipid
    :return: randomly distributing lipids
    """

    no_of_lipids = 8
    if system == S_CLASSICAL_REVERSIBLE:  # We need extra concentration for IP3
        no_of_lipids = 9

    all_ratios = np.random.uniform(0, 1, no_of_lipids)
    all_concentration = []
    for i in range(no_of_lipids):
        all_concentration.append(all_ratios[i] * total / sum(all_ratios))

    return all_concentration


def get_random_enzymes(kinetics) -> dict:
    """
    Makes new random enzyme based on kinetics

    :param kinetics: type of kinetics
    :return: enzyme list
    """
    pitp = RandomEnzyme(E_PITP, kinetics=kinetics)
    pi4k = RandomEnzyme(E_PI4K, kinetics=kinetics)
    pip5k = RandomEnzyme(E_PIP5K, kinetics=kinetics)
    plc = RandomEnzyme(E_PLC, kinetics=kinetics)
    dagk = RandomEnzyme(E_DAGK, kinetics=kinetics)
    laza = RandomEnzyme(E_LAZA, kinetics=kinetics)
    patp = RandomEnzyme(E_PATP, kinetics=kinetics)
    cds = RandomEnzyme(E_CDS, kinetics=kinetics)
    pis = RandomEnzyme(E_PIS, kinetics=kinetics)
    sink = RandomEnzyme(E_SINK, kinetics=kinetics)
    source = RandomEnzyme(E_SOURCE, kinetics=kinetics)
    p4tase = RandomEnzyme(E_P4TASE, kinetics=kinetics)
    p5tase = RandomEnzyme(E_P5TASE, kinetics=kinetics)
    ip3tase = RandomEnzyme(E_IP3_PTASE, kinetics=kinetics)
    return {x.name: x for x in
            [pitp, pi4k, pip5k, plc, dagk, laza, patp, cds, pis, sink, source,
             p4tase, p5tase, ip3tase]}


class Error:
    """
    Simple error class to handle result saving and keeping track of plots saved
    """
    total_files_saved = 0

    def __init__(self, initial_total: float, all_concentrations: list,
                 enzymes: dict,
                 cut_off: float = save_cutoff):
        self.all_concentrations = all_concentrations
        self.total_concentration = sum(all_concentrations)
        self.pmpi, self.pi4p, self.pip2, self.dag, self.pmpa, self.erpa, self.cdpdag, self.erpi = all_concentrations[
                                                                                                  :8]
        self.pa = self.erpa + self.pmpa
        self.pi = self.erpi + self.pmpi
        self.cut_off = cut_off
        self.enzymes = enzymes
        self.cdpdag_weight = 0.5  # Weight for CDPDAG error
        self.initial_total = initial_total

    @property
    def pip2_error(self) -> float:
        return abs((self.pip2 / self.pi) - alpha) / alpha

    @property
    def pi4p_error(self) -> float:
        return abs((self.pi4p / self.pi) - beta) / beta

    @property
    def dag_error(self) -> float:
        return abs((self.dag / self.pi) - gamma) / gamma

    @property
    def cdpdag_error(self) -> float:
        return abs(
            (self.cdpdag / self.pi) - epsilon) * self.cdpdag_weight / epsilon

    @property
    def pa_error(self) -> float:
        return abs((self.pa / self.pi) - delta) / delta

    @property
    def total_error(self) -> float:
        return self.pip2_error + self.pi4p_error + self.dag_error + self.cdpdag_error + self.pa_error

    def print_errors(self):
        print(self.total_error, self.pip2_error, self.pi4p_error,
              self.pa_error, self.dag_error, self.cdpdag_error)

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
