

class Constraint(dict):
    """
    constraint type
    """

    def __init__(self):
        pass


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
