import core.graph as graph
import core.predictor as predictor
import core.constraint as constraint
import core.util as util

import pickle


class LeadOptimizer:
    """
    Wrap everything up.
    """

    def __init__(self, lead_candidate: str, receptor: str, dock_config: str, opt_target: constraint.Constraint,
                 iter_num: int, beam_width: int, exhaustiveness: float, epsilon: float = 1.0,
                 pass_line: float = 1, checkpoint: str = "optimizer_checkpoint", mmpdb: str = graph.DEFAULT_MMPDB,
                 max_variable_size: int = 4, prediction_workers: int = 1, dock_method: str = "vina", dock_cpu: int = 1):
        """
        Initialize a LeadOptimizer
        :param lead_candidate: A SMILES string of the lead candidate
        :param receptor: A pdb or pdbqt file that contains the structure of the receptor protein
        :param dock_config: A txt file that contains docking configuration
        :param opt_target: A set of constraints that defines the target of the optimization
        :param iter_num: Number of iterations
        :param beam_width: Beam width in the searching algorithm
        :param exhaustiveness: Defines the exhaustiveness of the searching algorithm (between 0 and 1)
        :param epsilon: Defines the preference to give retry chances to the suboptimal solutions
        :param pass_line: Deprecated
        :param checkpoint: Checkpoint filename
        :param mmpdb: The mmpdb file used by mmpdb transform.
                    See more on: https://github.com/rdkit/mmpdb
        :param max_variable_size: Argument used in mmpdb transform
        :param prediction_workers: Number of threads making prediction on molecular properties
        :param dock_method: Docking method
        :param dock_cpu: Number of cpu cores used in docking procedure.
        """
        pw = predictor.PredictorWrapper(receptor, dock_config, dock_method, dock_cpu)
        source_mol = graph.Node(lead_candidate)
        self.beam_search_solver = graph.BeamSearchSolver(pw, opt_target, iter_num, beam_width, exhaustiveness,
                                                         epsilon, source_mol, pass_line, checkpoint, mmpdb,
                                                         max_variable_size, prediction_workers)

    def run(self):
        self.beam_search_solver.run()

    def summary(self):
        result = "Lead candidate:\n"
        result += self.beam_search_solver.summary([self.beam_search_solver.source_mol])
        result += "Optimization result:\n"
        result += self.beam_search_solver.summary(self.beam_search_solver.fringe)
        return result

    def visualize(self):
        if self.beam_search_solver.last_elapse_time is None:
            return util.visualize([self.beam_search_solver.source_mol])
        else:
            return util.visualize(self.beam_search_solver.fringe)

    def load_checkpoint(self, checkpoint: str):
        with open(checkpoint, "rb") as f:
            self.beam_search_solver = pickle.load(f)


