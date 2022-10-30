from schrodinger.pipeline.stages.qikprop import QikPropStage as qp
from rdkit import Chem
import schrodinger.pipeline.pipeio as pipeio
from schrodinger import structure
import os
import re
import pandas as pd

VINA = "D:\\Program Files (x86)\\The Scripps Research Institute\\Vina"
OBABEL = "D:\\Program Files\\OpenBabel-3.1.1"


def qikprop(mols: list) -> pd.DataFrame:
    """
    Makes qikprop prediction.
    :param mols: SMILES strings of molecules
    :return: qikprop result
    """
    # list of 3d structures
    m_list = [structure.SmilesStructure(m).get3dStructure(require_stereo=False) for m in mols]
    # TODO: fix this for cases with chirality
    id = abs(hash(''.join(mols)))  # code will crash if id is negative. why? ask Inc. Schrodinger!
    sdf_filename = f"{id}.sdf"

    os.chdir("./scratch")
    with open(sdf_filename, 'w') as file:
        pass  # create empty file
    with structure.StructureWriter(sdf_filename) as writer:
        for m in m_list:
            writer.append(m)

    qp_solver = qp(f"{id}")
    ligandsobj = pipeio.Structures([sdf_filename])
    # magic number 1? ask Inc. Schrodinger!
    qp_solver.setInput(1, 'INPUT', ligandsobj)
    qp_solver.setOutputName(1, 'OUTPUT')
    qp_solver.run()

    result = pd.read_csv(f"{id}-001.CSV")  # why -001? ask Inc. Schrodinger!
    # _delete_scratch(str(id))
    os.chdir("..")

    return result


def _delete_scratch(id: str):
    file_list = os.listdir(".")
    for f in file_list:
        if re.match(id, f) is not None:
            os.remove(f)


def dock_score(mol: str, protein, method: str) -> float:
    """
    Returns docking score.
    :param mol: the ligand
    :param protein: the protein
    :param method: docking method
    :return: docking score
    """
    _prepare_docking_file(mol, protein, method)
    # TODO: complete this


def _prepare_docking_file(mol: str, protein, method: str) -> int:
    id = abs(hash(mol))

    pass


def visualize(mol: str):
    pass

