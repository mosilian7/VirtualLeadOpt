import common.util as util
import pandas as pd


class Node:
    """
    node in searching algorithm
    """

    def __init__(self, mol: str, prop: dict = pd.DataFrame, neighbors: list = None):
        self.mol = mol
        self.prop = prop
        self.neighbors = neighbors

    def __str__(self):
        return self.mol


class Graph:
    """
    a set of nodes
    """

    def __init__(self, mols: list = [], mmpdb: str = "data/PubChem_mmp/random_200000.mmpdb", max_variable_size: int = 4):
        self.mols = mols
        self.mmpdb = mmpdb
        self.max_variable_size = max_variable_size

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
        self.mols.extend(node.neighbors)
