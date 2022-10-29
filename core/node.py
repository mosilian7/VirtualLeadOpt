class Node:
    """
    node in searching algorithm
    """

    def __init__(self, mol: str, prop: dict = None, neighbors: list = None):
        self.mol = mol
        self.prop = prop
        self.neighbors = neighbors

