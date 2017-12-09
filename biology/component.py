"""
General objects used in all other simulations.
"""

from numpy.random import uniform

from constants.namespace import MASS_ACTION, MICHAELIS_MENTEN, E_SOURCE
from utils.functions import clamp
from .default_values import *


def randomize_within(number, around_range, min_allowed, max_allowed):
    """
    Randomize parameter within given range
    :param number: original number
    :param around_range: range around which it should be randomize
    :param min_allowed: minimum allowed value
    :param max_allowed: maximum allowed value
    :return: new value which is randomized according to above rules
    """
    min_number = clamp(number - around_range, min_allowed, max_allowed)
    max_number = clamp(number + around_range, min_allowed, max_allowed)
    return uniform(min_number, max_number)


class Enzyme:
    """
    General Enzyme class for with all essential features for analysis and
    parameter sweep
    """

    def __init__(self, name, **kwargs) -> None:
        super().__init__()
        self.name = name
        self.v = kwargs.get("v") or default_v
        self.k = kwargs.get("k") or default_k
        self.kinetics = kwargs.get("kinetics") or MASS_ACTION
        self.original_v = self.v
        self.original_k = self.k
        self.feedback_factor = 1
        self.feedback_amount = 1
        self.feedback_substrate = None

    @classmethod
    def make_with_values(cls, name: str, values: dict):
        enz = cls(name, v=values["v"], k=values["k"],
                  kinetics=values["kinetics"])
        return enz

    @property
    def properties(self) -> dict:
        return {"v": self.v, "k": self.k, "kinetics": self.kinetics,
                "name": self.name}

    def copy_properties(self, properties: dict) -> None:
        self.v = properties.get("v")
        self.k = properties.get("k")
        self.kinetics = properties.get("kinetics")

    def react_with(self, substrate: float) -> float:
        if substrate is not None:
            if self.kinetics == MASS_ACTION:
                return self.k * substrate
            elif self.kinetics == MICHAELIS_MENTEN:
                return (self.v * substrate) / (self.k + substrate)
        else:
            return self.k  # For source

    def react_with_feedback(self, substrate: float,
                            feedback_substrate: float) -> float:

        f_n = 1 + self.feedback_factor * (
            feedback_substrate / self.feedback_amount)

        f_d = 1 + (feedback_substrate / self.feedback_amount)

        feedback_factor = f_n / f_d

        if substrate is not None:
            if self.kinetics == MASS_ACTION:
                return self.k * substrate * feedback_factor
            elif self.kinetics == MICHAELIS_MENTEN:
                return (self.v * feedback_factor * substrate) / (
                    self.k + substrate)
        else:
            return self.k * feedback_factor  # For source

    def randomize_v(self):
        self.v = randomize_within(self.v, rand_around_v, min_v, max_v)

    def randomize_k(self):
        self.k = randomize_within(self.k, rand_around_k, min_k, max_k)

    def randomize(self):
        self.randomize_v()
        self.randomize_k()

    def use_current(self):
        self.original_v = self.v
        self.original_k = self.k

    def reset(self):
        self.v = self.original_v
        self.k = self.original_k
        self.feedback_amount = 1
        self.feedback_substrate = None
        self.feedback_factor = 1

    def mutate(self, factor):
        if self.kinetics == MASS_ACTION or self.name == E_SOURCE:
            self.k *= factor
        else:
            self.v *= factor


class RandomEnzyme(Enzyme):
    """
       Creates enzyme with random V as well as Km value.
       IMPORTANT: any V or K value set from constructor will be overridden
       """

    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)
        self.k = uniform(min_k, max_k)
        self.v = uniform(min_v, max_v)
