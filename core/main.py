from core.constraint import Constraint
import common.util


class Solver:
    """
    solver who finds the path
    """

    def __init__(self, start: str, constraint: Constraint, protein: str = None):
        """
        Creates a Solver object.
        :param start: a SMILES string of 1st generation lead compound
        :param constraint: a set of constraints
        :param protein: the structure of protein
        """
        self.start = start
        self.constraint = constraint
        self.protein = protein

    def dis_estimate(self, mol: str) -> float:
        """
        Estimates the distance from the given molecule to the target area.
        :param mol: SMILES string of the valued ligand molecule
        :return: the distance
        """
        # TODO: complete this
