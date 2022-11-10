import heapq
import os
import threading

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
        """
        Each node correspond to a molecule in the searching algorithm. It contains some information about the molecule.
        :param mol: the SMILES string of the molecule
        :param prop: the properties of the molecule predicted by Predictor
        :param neighbors: "similar" molecules discovered by mmp transform
        """
        self.mol = Chem.CanonSmiles(mol)
        self.prop = prop
        self.neighbors = neighbors
        self.dis_to_target = float('inf')

    def __str__(self):
        return self.mol

    def __lt__(self, other):
        return self.dis_to_target < other.dis_to_target

    def __eq__(self, other):
        return self.mol == other.mol

    def __hash__(self):
        return hash(self.mol)


class MultiThreadPredictor:
    """
    use multiple predictors to perform prediction on a set of nodes
    """
    def __init__(self, nodes: list, predictor_wrapper: PredictorWrapper, prediction_workers: int):
        """
        Initializer of a MultiThreadPredictor
        :param nodes: corresponds to a set of molecules
        :param predictor_wrapper: a PredictorWrapper object
        :param prediction_workers: number of threads
        """
        self.nodes_split = util.split_list(nodes, prediction_workers)
        self.predictor_wrapper = predictor_wrapper
        self.prediction_workers = prediction_workers
        self.predictors = []
        for ns in self.nodes_split:
            self.predictors.append(self.predictor_wrapper.predictor_instance([str(n) for n in ns]))

    def run_predictor_method(self, method: str, update_node: bool = True) -> None:
        """
        runs predictor method on multiple threads
        :param method: a Predictor method
        :param update_node: whether update prop attribute of nodes in self.nodes_split
        TODO: is it possible to implement this without os.chdir?
        """
        os.chdir("./scratch")

        threads = [None for i in range(self.prediction_workers)]
        for i in range(self.prediction_workers):
            threads[i] = threading.Thread(target=getattr(self.predictors[i], method))
            threads[i].start()

        for i in range(self.prediction_workers):
            threads[i].join()
            if update_node:
                for j in range(len(self.nodes_split[i])):
                    self.nodes_split[i][j].prop = self.predictors[i].predictions.loc[j]

        os.chdir("..")


class Graph:
    """
    a set of nodes
    """

    def __init__(self, predictor_wrapper: PredictorWrapper, constraint: Constraint, mols: list = [],
                 mmpdb: str = DEFAULT_MMPDB, max_variable_size: int = 4, prediction_workers: int = 1):
        """
        A graph contains some information that is going to be used when implementing searching algorithm.
        :param predictor_wrapper: a PredictorWrapper object
        :param constraint: a Constraint object
        :param mols:
        :param mmpdb: mmp database used in find_neighbors
        :param max_variable_size: parameter in mmp transform
        :param prediction_workers: number of threads
        """
        self.predictor_wrapper = predictor_wrapper
        self.constraint = constraint
        self.pq = mols
        self.mmpdb = mmpdb
        self.max_variable_size = max_variable_size

        self.prediction_workers = prediction_workers
        self.predictors = [None for i in range(self.prediction_workers)]
        self.nodes_split = None

    def multi_thread_predictor(self, nodes: list) -> MultiThreadPredictor:
        """
        Returns a MultiThreadPredictor that performs predictions on a set of nodes
        :param nodes: a list of Node
        """
        return MultiThreadPredictor(nodes, self.predictor_wrapper, self.prediction_workers)

    @staticmethod
    def _parse_mmpdb_trans_out(mmpdb_trans_out: str) -> list:
        lines = mmpdb_trans_out.split('\n')
        result = []
        for l in lines:
            record = l.split('\t')
            if len(record) > 1 and record[1] != "SMILES":
                result.append(record[1])
        return result

    def find_neighbors(self, node: Node) -> None:
        """
        Finds neighbors of a node.
        :param node: a Node
        """
        mmpdb_trans_out = util.run_args(["python", util.MMPDB, "transform",
                                         "--smiles", node.mol,
                                         self.mmpdb,
                                         "--max-variable-size", str(self.max_variable_size)])
        # run python mmpdb transform --help to get help

        mols = self._parse_mmpdb_trans_out(mmpdb_trans_out)
        node.neighbors = [Node(m) for m in mols]

    def run_prediction_on(self, nodes: list) -> None:
        """
        Runs prediction on a list of nodes.
        :param nodes: a list of Node
        """
        os.chdir("./scratch")
        nodes_split = util.split_list(nodes, self.prediction_workers)
        for i in range(self.prediction_workers):
            self.predictors[i] = self.predictor_wrapper.predictor_instance([str(n) for n in nodes_split[i]])
            self.predictors[i].start()

        for i in range(self.prediction_workers):
            self.predictors[i].join()
            for j in range(len(nodes_split[i])):
                nodes_split[i][j].prop = self.predictors[i].predictions.loc[j]

        os.chdir("..")

    def estimate_dis_on(self, nodes: list) -> None:
        for n in nodes:
            n.dis_to_target = self.constraint.calculate_dis(n.prop)


class BeamSearchSolver(Graph):
    """
    search the chemical space by beam search
    """

    def __init__(self, predictor_wrapper: PredictorWrapper, constraint: Constraint, iter_num: int,
                 mols: list = [], mmpdb: str = DEFAULT_MMPDB, max_variable_size: int = 4):
        super(BeamSearchSolver, self).__init__(predictor_wrapper, constraint, mols, mmpdb, max_variable_size)
        self.iter_num = iter_num
        self.discard = set()
        heapq.heapify(self.pq)

    def run(self):
        # TODO: thread safety
        pass



