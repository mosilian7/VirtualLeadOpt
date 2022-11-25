from collections.abc import Callable
from typing import Dict

import pandas as pd
import math
import core.util as util


class Constraint(dict):
    """
    Constraint on the searching target.
    The object itself is a dictionary which contains key-value pairs, where keys are the name of considered property,
    and values are (desirability function, weight, default value) tuples. PREDEFINED_CONSTRAINT shows an example.
    """

    def __init__(self, evaluator: Callable, extra_args: dict = None):
        super(Constraint, self).__init__()
        self.evaluator = evaluator
        if extra_args is None:
            self.extra_args = {}
        else:
            self.extra_args = extra_args

    def update_weight(self, weight: Dict[str, float]):
        for k in weight:
            self[k] = (self[k][0], weight[k], self[k][2])

    def eval_property(self, prop: pd.Series) -> float:
        """
        Runs distance estimation.
        :param prop: the pd.DataFrame returned by Predictor.run()
        :return: a list of floats
        """
        return self.evaluator(self, prop, **self.extra_args)


class NanHandler:
    # Deprecated. Something goes wrong when pickle.dump is called.
    def __init__(self, func):
        self.func = func

    def __call__(self, val, **kwargs):
        if math.isnan(val):
            return 1e-8
        return self.func(val, **kwargs)


def geometric_mean_evaluator(constraint: Constraint, prop: pd.Series, **kwargs) -> float:
    desirabilities = []
    weight = []
    for item in constraint:
        desirabilities.append(constraint[item][0](prop.get(item, default=constraint[item][2]), **kwargs))
        weight.append(constraint[item][1])
    return util.geometric_mean(desirabilities, weight)


def check_nan(func):
    # Deprecated. This function cannot be dumped.
    def wrap(val, **kwargs):
        if math.isnan(val):
            return 1e-8
        return func(val, **kwargs)
    return wrap


def sigmoid_delta_dock_score(dock_score: float, **kwargs) -> float:
    if math.isnan(dock_score):
        return 1e-8
    a = 2  # a = 0.75 when vina is used
    return util.sigmoid(-(dock_score - kwargs["orig_dock_score"]) * a)


def eval_herg_log_ic50(herg_log_ic50: float, **kwargs) -> float:
    if math.isnan(herg_log_ic50):
        return 1e-8
    if herg_log_ic50 < -6 + 1e-8:
        return 1e-8
    elif herg_log_ic50 < -5:
        return herg_log_ic50 + 6
    else:
        return 1


def eval_qed(qed: float, **kwargs) -> float:
    if math.isnan(qed):
        return 1e-8
    return qed


def eval_sa_score(sa_score: float, **kwargs) -> float:
    if math.isnan(sa_score):
        return 1e-8
    return util.sigmoid(-(sa_score - 3))


def eval_pcaco(pcaco: float, **kwargs) -> float:
    if math.isnan(pcaco):
        return 1e-8
    if pcaco < 25:
        return pcaco / (100/3)
    elif pcaco < 500:
        return 0.75 + (pcaco - 25) / 1900
    else:
        return 1


def eval_num_rot_bond(num_rot_bond, **kwargs) -> float:
    if math.isnan(num_rot_bond):
        return 1e-8
    if num_rot_bond <= 10:
        return 1
    else:
        return 1e-8


PREDEFINED_CONSTRAINT = Constraint(geometric_mean_evaluator)
PREDEFINED_CONSTRAINT["dock_score"] = (sigmoid_delta_dock_score, 1, -float('inf'))
PREDEFINED_CONSTRAINT["QPlogHERG"] = (eval_herg_log_ic50, 0.3, float('nan'))
PREDEFINED_CONSTRAINT["qed"] = (eval_qed, 0.2, float('nan'))
PREDEFINED_CONSTRAINT["sa_score"] = (eval_sa_score, 0.3, float('nan'))
PREDEFINED_CONSTRAINT["QPPCaco"] = (eval_pcaco, 0.3, float('nan'))
PREDEFINED_CONSTRAINT["num_rot_bond"] = (eval_num_rot_bond, 0, float('nan'))

