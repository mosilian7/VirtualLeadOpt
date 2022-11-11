import os
import threading
import random
import pickle

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
        self.prop = pd.Series([], dtype=float) if prop is None else prop
        self.neighbors = neighbors
        self.evaluation = 0

    def __str__(self):
        return self.mol

    def __lt__(self, other):
        return self.evaluation < other.evaluation

    def __eq__(self, other):
        return self.mol == other.mol

    def __hash__(self):
        return hash(self.mol)

    def qed_done(self) -> bool:
        return "qed" in self.prop

    def qikprop_done(self) -> bool:
        return "QPlogS" in self.prop

    def docking_done(self) -> bool:
        return "dock_score" in self.prop

    def prediction_done(self) -> bool:
        return self.qed_done() and self.qikprop_done() and self.docking_done()


class MultiThreadPredictor:
    """
    use multiple predictors to perform prediction on a set of nodes
    """
    def __init__(self, nodes: list, predictor_wrapper: PredictorWrapper, prediction_workers: int, mask: dict = None):
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
        for i in range(self.prediction_workers):
            self.predictors.append(self.predictor_wrapper
                                   .predictor_instance([str(n) for n in self.nodes_split[i]]))
        mask = {} if mask is None else mask
        self.update_mask(mask)

    def update_mask(self, mask: dict) -> None:
        """
        Update masking of each predictor in self.predictors.
        :param mask: masking directory
        """
        mask_split = {}
        for k in mask:
            mask_split[k] = util.split_list(mask[k], self.prediction_workers)
        for i in range(self.prediction_workers):
            mask_dict = {}
            for k in mask:
                mask_dict[k] = mask_split[k][i]
            self.predictors[i].mask = mask_dict

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
        :param mols: source molecules
        :param mmpdb: mmp database used in find_neighbors
        :param max_variable_size: parameter in mmp transform
        :param prediction_workers: number of threads
        """
        self.predictor_wrapper = predictor_wrapper
        self.constraint = constraint
        self.fringe = mols
        self.mmpdb = mmpdb
        self.max_variable_size = max_variable_size

        self.prediction_workers = prediction_workers
        self.predictors = [None for i in range(self.prediction_workers)]
        self.nodes_split = None

    def multi_thread_predictor(self, nodes: list, mask: dict = None) -> MultiThreadPredictor:
        """
        Returns a MultiThreadPredictor that performs predictions on a set of nodes
        :param nodes: a list of Node
        :param mask: masking dictionary
        """
        if mask is None:
            mask = {"docking": [False] * len(nodes)}
        return MultiThreadPredictor(nodes, self.predictor_wrapper, self.prediction_workers, mask=mask)

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
        Runs prediction on a list of nodes. Deprecated?
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

    def evaluate(self, nodes: list) -> None:
        for n in nodes:
            n.evaluation = self.constraint.eval_property(n.prop)


class BeamSearchSolver(Graph):
    """
    search the chemical space by beam search
    """

    def __init__(self, predictor_wrapper: PredictorWrapper, constraint: Constraint,
                 iter_num: int, beam_width: int, exhaustiveness: float, epsilon: float, source_mol: Node,
                 pass_line: float = 0.9, checkpoint: str = "bss_checkpoint",
                 mmpdb: str = DEFAULT_MMPDB, max_variable_size: int = 4, prediction_workers: int = 1):
        super(BeamSearchSolver, self).__init__(predictor_wrapper, constraint, [source_mol],
                                               mmpdb, max_variable_size, prediction_workers)
        self.source_mol = source_mol
        self.iter_num = iter_num
        self.beam_width = beam_width
        self.exhaustiveness = exhaustiveness
        self.retry_chance = int(beam_width * epsilon)
        self.discard = set()
        self.pass_line = pass_line
        self.result = set()
        self.error = set()
        self.checkpoint = checkpoint
        self.history = []

    def run(self):
        """
            Run the algorithm.
        """
        # TODO: thread safety
        self.start_up()
        self.history = [(self.fringe.copy(), self.result.copy())]
        for i in range(self.iter_num):
            self.one_step()
            print(f"\033[32mIteration {i} finished\033[0m")
            self.history.append((self.fringe.copy(), self.result.copy()))
            self._save_checkpoint()
            self.give_chance_to_discarded()

    def _save_checkpoint(self):
        # Serialize this.
        with open(self.checkpoint, "wb") as f:
            pickle.dump(self, f)
        print("checkpoint saved")

    def _to_be_predicted(self) -> list:
        # Returns molecules that will undergo prediction next.
        curr = set()
        run_prediction = []

        # Find neighbors for each node in self.fringe.
        for n in self.fringe:
            if n.neighbors is None:
                self.find_neighbors(n)
            curr.update(set(n.neighbors))

        for n in curr:
            if n.prediction_done() or n in self.discard or n in self.result:
                continue
            run_prediction.append(n)

        random.shuffle(run_prediction)

        # Because of limited resource, prediction will only run on a part of molecules.
        run_prediction = run_prediction[:int(self.exhaustiveness * len(run_prediction))]

        return run_prediction

    def _clipping(self, run_prediction) -> dict:
        # Clip the out-of-bound molecules, so that they don't have to dock.
        mask = {"docking": [False] * len(run_prediction)}
        bound = min(self.fringe).evaluation
        for i in range(len(run_prediction)):
            if run_prediction[i].evaluation < bound:
                self.discard.add(run_prediction[i])
                mask["docking"][i] = True
        return mask

    def _dedup_and_sort(self):
        # Deduplicate and sort the fringe.
        no_dup = set(self.fringe)
        self.fringe = list(no_dup)
        self.fringe = sorted(self.fringe, reverse=True)

    def _update_result(self):
        # If a molecule has evaluation > pass line, then adds it to self.result.
        i = 0
        while self.fringe[i].evaluation > self.pass_line:
            self.result.add(self.fringe[i])
            i += 1

    def one_step(self) -> None:
        """
            One iteration step in the algorithm.
        """

        run_prediction = self._to_be_predicted()

        mtp = self.multi_thread_predictor(run_prediction)
        mtp.run_predictor_method("prepare_sdf", update_node=False)
        mtp.run_predictor_method("qed")
        mtp.run_predictor_method("qikprop")

        # this is not the true distance!
        self.evaluate(run_prediction)

        mask = self._clipping(run_prediction)
        mtp.update_mask(mask)

        mtp.run_predictor_method("dock_score")
        mtp.run_predictor_method("delete_scratch", update_node=False)

        self.evaluate(run_prediction)
        self.fringe.extend(run_prediction)
        self._dedup_and_sort()
        self._update_result()

        for n in self.fringe[self.beam_width:]:
            self.discard.add(n)

        # update fringe
        self.fringe = self.fringe[:self.beam_width]

    def give_chance_to_discarded(self):
        """
            Adds some discarded molecules to the fringe.
        """
        if len(self.discard) > 0:
            s = set()
            l = list(self.discard)
            for i in range(self.retry_chance):
                s.add(random.choice(l))
            no_dup = set(self.fringe)
            no_dup.update(s)
            self.fringe = list(no_dup)

    def start_up(self) -> None:
        """
            Starts up from source molecule.
        """

        # there is only one source molecule, so more dock cpu is used
        orig_dock_cpu = self.predictor_wrapper.dock_cpu
        self.predictor_wrapper.dock_cpu = 12

        mtp = MultiThreadPredictor([self.source_mol], self.predictor_wrapper, prediction_workers=1)
        mtp.run_predictor_method("prepare_sdf", update_node=False)
        mtp.run_predictor_method("qed")
        mtp.run_predictor_method("qikprop")
        mtp.run_predictor_method("dock_score")
        mtp.run_predictor_method("delete_scratch", update_node=False)
        self.predictor_wrapper.dock_cpu = orig_dock_cpu
        self.constraint.extra_args["orig_dock_score"] = self.source_mol.prop["dock_score"]
        self.evaluate([self.source_mol])




