from collections.abc import Callable
import pandas as pd


class Constraint(dict):
    """
    constraint type
    """

    def __init__(self, dis_estimator: Callable):
        self.dis_estimator = dis_estimator

    def calculate_dis(self, prop: pd.Series) -> float:
        """
        Runs distance estimation.
        :param prop: the pd.DataFrame returned by Predictor.run()
        :return: a list of floats
        """
        # TODO: complete this


class DisEstimator:
    """
    distance estimator
    """

    @staticmethod
    def euler(constraint: Constraint, prop: dict) -> float:
        """
        euler distance
        :param constraint: a set of constrains
        :param prop: dictionary of properties
        :return: distance estimation
        """
        # TODO: complete this
