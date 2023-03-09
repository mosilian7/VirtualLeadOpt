import math
import os
import threading
import random
import pickle
import time
from functools import reduce
from typing import Dict, List

import core.util as util
import pandas as pd
from rdkit import Chem
from core.predictor import PredictorWrapper
from core.constraint import Constraint

DEFAULT_MMPDB = "data/ChEMBL_mmp/random_500000.mmpdb"


class Node:
    """
        node in searching algorithm
    """

    def __init__(self, mol: str, prop: pd.Series = None, neighbors: List["Node"] = None):
        """
        Each node correspond to a molecule in the searching algorithm. It contains some information about the molecule.
        :param mol: the SMILES string of the molecule
        :param prop: the properties of the molecule predicted by Predictor
        :param neighbors: "similar" molecules discovered by mmp transform
        """
        self.mol = Chem.CanonSmiles(mol)
        self.prop = pd.Series([], dtype=float) if prop is None else prop
        self.neighbors = neighbors if neighbors is not None else []
        self.evaluation = 0
        self.evaluation_done = False

    def __str__(self) -> str:
        return self.mol

    def __lt__(self, other) -> bool:
        return self.evaluation < other.evaluation

    def __eq__(self, other) -> bool:
        return self.mol == other.mol

    def __hash__(self) -> int:
        return hash(self.mol)

    def qed_done(self) -> bool:
        return "qed" in self.prop

    def qikprop_done(self) -> bool:
        return "QPlogS" in self.prop

    def sa_score_done(self) -> bool:
        return "sa_score" in self.prop

    def docking_done(self) -> bool:
        return "dock_score" in self.prop

    def prediction_done(self) -> bool:
        return self.sa_score_done() and self.docking_done()


class MultiThreadPredictor:
    """
    use multiple predictors to perform prediction on a set of nodes
    """

    def __init__(self, nodes: List[Node], predictor_wrapper: PredictorWrapper, prediction_workers: int,
                 mask: Dict[str, List[bool]] = None):
        """
        Initializer of a MultiThreadPredictor
        :param nodes: corresponds to a set of molecules
        :param predictor_wrapper: a PredictorWrapper object
        :param prediction_workers: number of threads
        """
        self.predictor_wrapper = predictor_wrapper
        self.prediction_workers = len(nodes) if prediction_workers > len(nodes) else prediction_workers
        self.nodes_split = util.split_list(nodes, self.prediction_workers)
        self.predictors = []
        for i in range(self.prediction_workers):
            self.predictors.append(self.predictor_wrapper
                                   .predictor_instance([str(n) for n in self.nodes_split[i]]))
        mask = {} if mask is None else mask
        self.update_mask(mask)

    def update_mask(self, mask: Dict[str, List[bool]]) -> None:
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
                 mmpdb: str = DEFAULT_MMPDB, substructure: str = None, max_variable_size: int = 4,
                 g2g_translation_model: str = util.G2G_TRANSLATION_MODEL, g2g_num_decode: int = 0,
                 prediction_workers: int = 1):
        """
        A graph contains some information that is going to be used when implementing searching algorithm
        :param predictor_wrapper: a PredictorWrapper object
        :param constraint: a Constraint object
        :param mols: source molecules
        :param mmpdb: mmp database used in find_neighbors
        :param substructure: SMARTS pattern, which describe the substructure that is wished to be preserved
        :param max_variable_size: parameter in mmp transform
        :param g2g_translation_model: path of graph-to-graph translation model
        :param g2g_num_decode: number of molecules generated by g2g translation
        :param prediction_workers: number of threads
        """
        self.predictor_wrapper = predictor_wrapper
        self.constraint = constraint
        self.fringe = mols
        self.mmpdb = mmpdb
        self.substructure = substructure
        self.max_variable_size = max_variable_size
        self.g2g_translation_model = g2g_translation_model
        self.g2g_num_decode = g2g_num_decode

        self.prediction_workers = prediction_workers
        self.predictors = [None for i in range(self.prediction_workers)]
        self.nodes_split = None

    @classmethod
    def _debug_helper(cls):
        return cls(None, None)

    def multi_thread_predictor(self, nodes: list, mask: dict = None) -> MultiThreadPredictor:
        """
        Returns a MultiThreadPredictor that performs predictions on a set of nodes
        :param nodes: a list of Node
        :param mask: masking dictionary
        """
        if mask is None:
            mask = {"dock_score": [False] * len(nodes)}
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
        substructure_param = [] if self.substructure is None else ["--substructure", self.substructure]
        mmpdb_trans_out = util.run_args(["python", util.MMPDB, "transform",
                                         "--smiles", node.mol,
                                         self.mmpdb,
                                         "--max-variable-size", str(self.max_variable_size),
                                         *substructure_param])
        # run python mmpdb transform --help to get help

        mols = self._parse_mmpdb_trans_out(mmpdb_trans_out)
        node.neighbors = [Node(m) for m in mols]

    @staticmethod
    def g2g_translate(mol: str, g2g_translation_model, num_decode: int = 20) -> List[str]:
        os.chdir("./scratch")
        mol_file = os.path.abspath(f"{hash(mol)}.txt")
        with open(mol_file, 'w') as f:
            f.write(f"{mol}\n\n")

        response = util.run_args(["python", "translate.py",
                                  "--test", mol_file,
                                  "--vocab", util.G2G_TRANSLATION_VOCAB,
                                  "--model", g2g_translation_model,
                                  "--num_decode", str(num_decode)],
                                 remote_endpoint="http://localhost:8112")

        lines = response.strip().split('\n')
        results = [line.split()[1] for line in lines]

        os.remove(mol_file)
        os.chdir("..")
        return results

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
                 iter_num: int, beam_width: int, exhaustiveness: float, retry_chance: float, source_mol: Node,
                 simulated_annealing_temperature: float = 1/40,
                 pass_line: float = 0.9, checkpoint: str = "bss_checkpoint",
                 mmpdb: str = DEFAULT_MMPDB, substructure: str = None, max_variable_size: int = 4,
                 g2g_translation_model: str = util.G2G_TRANSLATION_MODEL, g2g_num_decode: int = 0,
                 prediction_workers: int = 1, debug_info: bool = False):
        super(BeamSearchSolver, self).__init__(predictor_wrapper, constraint, [source_mol],
                                               mmpdb, substructure, max_variable_size,
                                               g2g_translation_model, g2g_num_decode, prediction_workers)
        self.source_mol = source_mol
        self.iter_num = iter_num
        self.beam_width = beam_width
        self.exhaustiveness = exhaustiveness
        self.retry_chance = int(beam_width * retry_chance) if retry_chance < 1 else int(retry_chance)
        self.simulated_annealing_temperature = simulated_annealing_temperature
        self.discard = set()
        self.discard = set()
        self.pass_line = pass_line
        self.result = set()
        self.error = set()
        self.checkpoint = checkpoint
        self.history = []
        self.last_elapse_time = None
        self.clip_bound = 0
        if debug_info:
            self.debug = []

    @staticmethod
    def load_from(checkpoint: str):
        with open(checkpoint, "rb") as f:
            return pickle.load(f)

    def run(self, start_from_checkpoint: bool = False, iter_num: int = -1):
        """
            Run the algorithm.
        """

        tic = time.time()

        if not start_from_checkpoint:
            self.start_up()
            self.history = [(self.fringe.copy(), self.result.copy())]
            iter_num = self.iter_num

        for i in range(iter_num):
            self.give_chance_to_discarded()
            self.one_step()
            self.history.append((self.fringe.copy(), self.result.copy()))
            self._save_checkpoint()
            print(f"\033[32mIteration {i} finished\033[0m")

        toc = time.time()
        self.last_elapse_time = toc - tic

    def _save_checkpoint(self):
        # Serialize this.

        if hasattr(self, "debug"):
            # predictor cannot be serialized
            debug_info = self.debug
            self.debug = []

        with open(self.checkpoint, "wb") as f:
            pickle.dump(self, f)
        print("checkpoint saved")

        if hasattr(self, "debug"):
            self.debug = debug_info

    def _update_discard(self, n: Node) -> None:
        self.discard.add(n)
        n.neighbors = []

    def _to_be_predicted(self) -> list:
        # Returns molecules that will undergo prediction next.
        curr = set()
        run_prediction = []

        # Find neighbors for each node in self.fringe.
        for n in self.fringe:
            if len(n.neighbors) == 0:
                self.find_neighbors(n)
            curr.update(set(n.neighbors))

        for n in curr:
            if n.prediction_done() or n in self.discard or n in self.result:
                continue
            run_prediction.append(n)

        random.shuffle(run_prediction)

        if self.exhaustiveness < 1:
            len_run_prediction = int(self.exhaustiveness * len(run_prediction))
        else:
            len_run_prediction = int(self.exhaustiveness)
        # Because of limited resource, prediction will only run on a part of molecules.
        run_prediction = run_prediction[:len_run_prediction]

        if self.g2g_num_decode != 0:
            for n in self.fringe:
                g2g_hints = self.g2g_translate(str(n), self.g2g_translation_model, self.g2g_num_decode)
                run_prediction.extend([Node(n) for n in g2g_hints])

        no_dup = set(run_prediction)
        run_prediction = list(no_dup)
        random.shuffle(run_prediction)

        if hasattr(self, "debug"):
            self.debug.append(run_prediction)
        return run_prediction

    def _clipping(self, run_prediction) -> dict:
        # Clip the out-of-bound molecules, so that they don't have to dock.
        # idea of branch and bound
        mask = {"dock_score": [False] * len(run_prediction)}
        for i in range(len(run_prediction)):
            if run_prediction[i].evaluation < self.clip_bound:
                self._update_discard(run_prediction[i])
                mask["dock_score"][i] = True
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

    def _update_fringe_and_discard(self):
        new_fringe = []
        i = 0
        while len(new_fringe) < self.beam_width and i < len(self.fringe):
            if self.fringe[i].evaluation_done:
                new_fringe.append(self.fringe[i])
            i += 1

        for n in self.fringe:
            if n not in new_fringe:
                self._update_discard(n)

        self.fringe = new_fringe
        self.clip_bound = min(self.fringe).evaluation

    def one_step(self) -> None:
        """
            One iteration step in the algorithm.
        """

        run_prediction = self._to_be_predicted()

        mtp = self.multi_thread_predictor(run_prediction)
        if hasattr(self, "debug"):
            self.debug.append(mtp)
        mtp.run_predictor_method("prepare_sdf", update_node=False)
        mtp.run_predictor_method("before_dock")

        # this is not the true distance!
        self.evaluate(run_prediction)

        # branch and bound
        mask = self._clipping(run_prediction)
        # mask the molecules that don't have to dock
        mtp.update_mask(mask)

        mtp.run_predictor_method("dock_score")

        # mark that evaluation is done (for the molecules that haven't been masked)
        for i in range(len(run_prediction)):
            # if a node is masked on any property, then its evaluation is not done
            assert "dock_score" in run_prediction[i].prop, "Hint: check if the server is running"
            run_prediction[i].evaluation_done = not reduce(lambda x, y: x or y, [mask[k][i] for k in mask])

        mtp.run_predictor_method("delete_scratch", update_node=False)

        self.evaluate(run_prediction)
        self.fringe.extend(run_prediction)
        self._dedup_and_sort()
        self._update_result()
        self._update_fringe_and_discard()

    def give_chance_to_discarded(self, max_num_sample: int = 100000):
        """
            Adds some discarded molecules to the fringe.
        """
        if len(self.discard) > 0:
            s = set()
            l = list(self.discard)
            cnt = 0
            while len(s) < self.retry_chance and cnt < max_num_sample:
                n = random.choice(l)
                r = random.random()
                if r < math.exp((n.evaluation - self.fringe[-1].evaluation) / self.simulated_annealing_temperature):
                    s.add(n)
                cnt += 1

            no_dup = set(self.fringe)
            no_dup.update(s)
            self.fringe = list(no_dup)

    def start_up(self, update_constraint: bool = True) -> None:
        """
            Starts up from source molecule.
        """

        # there is only one source molecule, so more dock cpu is used
        orig_dock_cpu = self.predictor_wrapper.dock_cpu
        self.predictor_wrapper.dock_cpu = 12

        mtp = MultiThreadPredictor([self.source_mol], self.predictor_wrapper, prediction_workers=1)
        mtp.run_predictor_method("prepare_sdf", update_node=False)
        mtp.run_predictor_method("before_dock")
        mtp.run_predictor_method("dock_score")
        mtp.run_predictor_method("delete_scratch", update_node=False)
        self.predictor_wrapper.dock_cpu = orig_dock_cpu
        if update_constraint:
            self.constraint.extra_args["orig_dock_score"] = self.source_mol.prop["dock_score"]
        self.evaluate([self.source_mol])
        self.clip_bound = self.source_mol.evaluation

    def summary(self, nodes: list) -> str:
        result = ""
        schema = self.constraint.keys()
        result += "evaluation "
        for k in schema:
            result += f"| {k} "
        result += "\n"

        for n in nodes:
            line = ""
            line += f"{n.evaluation:.4f} "
            for k in schema:
                line += f"| {n.prop.get(k, default=float('nan')):.4f} "
            line += "\n"
            result += line

        return result
