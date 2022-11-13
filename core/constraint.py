from collections.abc import Callable
import pandas as pd
import math
import common.util as util


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

    def eval_property(self, prop: pd.Series) -> float:
        """
        Runs distance estimation.
        :param prop: the pd.DataFrame returned by Predictor.run()
        :return: a list of floats
        """
        return self.evaluator(self, prop, **self.extra_args)


def geometric_mean_evaluator(constraint: Constraint, prop: pd.Series, **kwargs) -> float:
    desirabilities = []
    weight = []
    for item in constraint:
        desirabilities.append(constraint[item][0](prop.get(item, default=constraint[item][2]), **kwargs))
        weight.append(constraint[item][1])
    return util.geometric_mean(desirabilities, weight)


def sigmoid_delta_dock_score(dock_score: float, **kwargs) -> float:
    if math.isnan(dock_score):
        return 1e-8
    return 1 / (1 + math.exp((dock_score - kwargs["orig_dock_score"]) * 0.75))


def eval_herg_log_ic50(herg_log_ic50: float, **kwargs) -> float:
    if math.isnan(herg_log_ic50):
        return 1e-8

    if herg_log_ic50 < -5.99:
        return 0.01
    elif herg_log_ic50 < -5:
        return herg_log_ic50 + 6
    else:
        return 1


def eval_qed(qed: float, **kwargs) -> float:
    if math.isnan(qed):
        return 1e-8

    return qed


PREDEFINED_CONSTRAINT = Constraint(geometric_mean_evaluator)
PREDEFINED_CONSTRAINT["dock_score"] = (sigmoid_delta_dock_score, 1, -float('inf'))
PREDEFINED_CONSTRAINT["QPlogHERG"] = (eval_herg_log_ic50, 0.5, float('nan'))
PREDEFINED_CONSTRAINT["qed"] = (eval_qed, 0.6, float('nan'))

