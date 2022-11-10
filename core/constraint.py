from collections.abc import Callable
import pandas as pd
import math
from abc import ABCMeta, abstractmethod
import common.util as util


class Constraint(dict):
    """
    constraint on the searching target

    """

    def __init__(self, dis_estimator: Callable):
        self.dis_estimator = dis_estimator

    def calculate_dis(self, prop: pd.Series) -> float:
        """
        Runs distance estimation.
        :param prop: the pd.DataFrame returned by Predictor.run()
        :return: a list of floats
        """
        return self.dis_estimator(self, prop)


def sigmoid_delta_dock_score(delta_dock_score: float) -> float:
    return 1 / (1 + math.exp(delta_dock_score))


def eval_herg_log_ic50(herg_log_ic50: float) -> float:
    if herg_log_ic50 < -6:
        return 0.01
    elif herg_log_ic50 < -5:
        return herg_log_ic50 + 6
    else:
        return 1


def eval_qed(qed: float) -> float:
    return qed


PREDEFINED_CONSTRAINT = Constraint(None)
PREDEFINED_CONSTRAINT["delta_dock_score"] = sigmoid_delta_dock_score
PREDEFINED_CONSTRAINT["QPlogHERG"] = eval_herg_log_ic50
PREDEFINED_CONSTRAINT["qed"] = eval_qed


class DisEstimator:
    """
    distance estimator
    """
    __metaclass__ = ABCMeta

    @staticmethod
    @abstractmethod
    def estimate(constraint: Constraint, prop: pd.Series) -> float:
        pass


class GeometricMean(DisEstimator):
    @staticmethod
    def estimate(constraint: Constraint, prop: pd.Series) -> float:
        desirabilities = []
        for item in constraint:
            desirabilities.append(constraint[item](prop[item]))
        return util.geometric_mean(desirabilities)




