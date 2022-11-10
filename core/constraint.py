from collections.abc import Callable
import pandas as pd
import math
import common.util as util


class Constraint(dict):
    """
    Constraint on the searching target.
    The object itself is a dictionary which contains key-value pairs, where keys are the name of considered property,
    and values are tuples (the desirability function of the corresponding property, weight). PREDEFINED_CONSTRAINT shows
    an example.

    """

    def __init__(self, dis_estimator: Callable, extra_args: dict = None):
        super(Constraint, self).__init__()
        self.dis_estimator = dis_estimator
        if extra_args is None:
            self.extra_args = {}
        else:
            self.extra_args = extra_args

    def calculate_dis(self, prop: pd.Series) -> float:
        """
        Runs distance estimation.
        :param prop: the pd.DataFrame returned by Predictor.run()
        :return: a list of floats
        """
        return self.dis_estimator(self, prop, **self.extra_args)


def geometric_mean_estimator(constraint: Constraint, prop: pd.Series, **kwargs) -> float:
    desirabilities = []
    weight = []
    for item in constraint:
        desirabilities.append(constraint[item][0](prop[item], **kwargs))
        weight.append(constraint[item][1])
    return util.geometric_mean(desirabilities, weight)


def sigmoid_delta_dock_score(dock_score: float, **kwargs) -> float:
    return 1 / (1 + math.exp(dock_score - kwargs["orig_dock_score"]))


def eval_herg_log_ic50(herg_log_ic50: float, **kwargs) -> float:
    if herg_log_ic50 < -5.99:
        return 0.01
    elif herg_log_ic50 < -5:
        return herg_log_ic50 + 6
    else:
        return 1


def eval_qed(qed: float, **kwargs) -> float:
    return qed


PREDEFINED_CONSTRAINT = Constraint(geometric_mean_estimator)
PREDEFINED_CONSTRAINT["dock_score"] = (sigmoid_delta_dock_score, 1)
PREDEFINED_CONSTRAINT["QPlogHERG"] = (eval_herg_log_ic50, 0.5)
PREDEFINED_CONSTRAINT["qed"] = (eval_qed, 0.8)

