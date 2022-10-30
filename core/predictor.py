from schrodinger.pipeline.stages.qikprop import QikPropStage as qp
from rdkit import Chem
import schrodinger.pipeline.pipeio as pipeio
from schrodinger import structure
import os
import re
import pandas as pd

VINA = "D:\\Program Files (x86)\\The Scripps Research Institute\\Vina"
OBABEL = "D:\\Program Files\\OpenBabel-3.1.1"


class Predictor:
    """ predict the properties of a set of molecules """

    def __init__(self, mols: list, protein: str, dock_config: str, dock_method: str = "vina"):
        """
        initializer. make sure that open babel is installed and added to PATH
        :param mols: a list of SMILES string
        :param protein: a .pdb or .pdbqt filename
        :param dock_config: a .txt configure filename
        :param dock_method: docking method, vina by default.
                            supported methods: vina
        """
        self.mols = mols
        self.id = abs(hash(''.join(self.mols)))  # code will crash if id is negative. why? ask Inc. Schrodinger!
        self.protein = protein
        self.dock_config = dock_config
        self.dock_method = dock_method
        self.predictions = None

    def qikprop(self):
        """
        Makes qikprop prediction.
        """
        # list of 3d structures
        m_list = [structure.SmilesStructure(m).get3dStructure(require_stereo=False) for m in self.mols]
        # TODO: fix this for cases with chirality
        sdf_filename = f"{self.id}.sdf"

        os.chdir("./scratch")
        with open(sdf_filename, 'w') as file:
            pass  # create empty file under scratch folder
        with structure.StructureWriter(sdf_filename) as writer:
            for m in m_list:
                writer.append(m)
        qp_solver = qp(f"{self.id}")
        ligandsobj = pipeio.Structures([sdf_filename])
        qp_solver.setInput(1, 'INPUT', ligandsobj) # magic number 1? ask Inc. Schrodinger!
        qp_solver.setOutputName(1, 'OUTPUT')
        qp_solver.run()
        result = pd.read_csv(f"{self.id}-001.CSV")  # why -001? ask Inc. Schrodinger!
        # _delete_scratch(str(id))
        os.chdir("..")

        self.predictions = result

    def _delete_scratch(self):
        file_list = os.listdir(".")
        for f in file_list:
            if re.match(self.id, f) is not None:
                os.remove(f)

    def dock_score(self):
        """
        Returns docking score.
        """

        self._prepare_docking_file()

        # vina usage: see https://vina.scripps.edu/manual/#usage
        # TODO: complete this

    def _prepare_docking_file(self):

        # convert .pdb to .pdbqt
        if re.search(".pdbqt", self.protein) is None:
            os.system(f"obabel -ipdb {self.protein} -opdbqt -O {self.protein}qt -p -xr")
            # for -p: see https://openbabel.org/docs/dev/Command-line_tools/babel.html
            # for -xr: see http://openbabel.org/docs/current/FileFormats/Overview.html
            self.protein = f"{self.protein}qt"

        # convert .sdf ligands to .pdbqt
        os.system(f"obabel -isdf ./scratch/{self.id}.sdf -opdbqt -O {self.id}_.pdbqt -m")
