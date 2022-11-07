import heapq
import os

import common.util as util
import pandas as pd
from rdkit import Chem
from core.predictor import Predictor, PredictorWrapper
from core.constraint import Constraint

DEFAULT_MMPDB = "data/ChEMBL_mmp/random_500000.mmpdb"


class Node:
    """
    node in searching algorithm
    """

    def __init__(self, mol: str, prop: pd.Series = None, neighbors: list = None):
        self.mol = Chem.CanonSmiles(mol)
        self.prop = prop
        self.neighbors = neighbors
        self.dis_to_target = float('inf')

    def __str__(self):
        return self.mol

    def __lt__(self, other):
        return self.dis_to_target < other.dis_to_target


class Graph:
    """
    a set of nodes
    """

    def __init__(self, pred_wrapper: PredictorWrapper, constraint: Constraint, mols: list = [],
                 mmpdb: str = DEFAULT_MMPDB, max_variable_size: int = 4, prediction_workers: int = 1):
        self.pred_wrapper = pred_wrapper
        self.constraint = constraint
        self.pq = mols
        self.mmpdb = mmpdb
        self.max_variable_size = max_variable_size
        self.prediction_workers = prediction_workers

    @staticmethod
    def _parse_mmpdb_trans_out(mmpdb_trans_out: str) -> list:
        lines = mmpdb_trans_out.split('\n')
        result = []
        for l in lines:
            record = l.split('\t')
            if len(record) > 1 and record[1] != "SMILES":
                result.append(record[1])
        return result

    def find_neighbors(self, node: Node):
        mmpdb_trans_out = util.run_args(["python", util.MMPDB, "transform",
                                         "--smiles", node.mol,
                                         self.mmpdb,
                                         "--max-variable-size", str(self.max_variable_size)])
        # run python mmpdb transform --help to get help

        mols = self._parse_mmpdb_trans_out(mmpdb_trans_out)
        node.neighbors = [Node(m) for m in mols]

    def run_prediction_on(self, nodes: list) -> None:
        os.chdir("./scratch")
        nodes_split = util.split_list(nodes, self.prediction_workers)
        predictors = [None for i in range(self.prediction_workers)]
        for i in range(self.prediction_workers):
            predictors[i] = self.pred_wrapper.predictor_instance([str(n) for n in nodes_split[i]])
            predictors[i].start()

        for i in range(self.prediction_workers):
            predictors[i].join()
            for j in range(len(nodes_split[i])):
                nodes_split[i][j].prop = predictors[i].predictions.loc[j]

        os.chdir("..")

    def estimate_dis_on(self, nodes: list):
        self.run_prediction(nodes)
        for n in nodes:
            n.dis_to_target = self.constraint.calculate_dis(n.prop)


class BeamSearchSolver(Graph):
    """
    search the chemical space by beam search
    """

    def __init__(self, pred_wrapper: PredictorWrapper, constraint: Constraint, iter_num: int,
                 mols: list = [], mmpdb: str = DEFAULT_MMPDB, max_variable_size: int = 4):
        super(BeamSearchSolver, self).__init__(pred_wrapper, constraint, mols, mmpdb, max_variable_size)
        self.iter_num = iter_num
        self.discard = set()
        heapq.heapify(self.pq)

    def run(self):
        # TODO: thread safety
        pass



