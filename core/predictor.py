import os
import re
import sys
import subprocess
import threading

import pandas as pd
import schrodinger.pipeline.pipeio as pipeio
from schrodinger import structure
from schrodinger.pipeline.stages.qikprop import QikPropStage as qp

VINA = "D:\\Program Files (x86)\\The Scripps Research Institute\\Vina"
OBABEL = "D:\\Program Files\\OpenBabel-3.1.1"


class Predictor:
    """
    predict the properties of a set of molecules
    methods should run in order: prepare_sdf, qikprop, dock_score
    """

    def __init__(self, mols: list, protein: str, dock_config: str, dock_method: str = "vina", log: str = "stdout"):
        """
        initializer. make sure that open babel is installed and added to PATH
        :param mols: a list of SMILES string
        :param protein: a .pdb or .pdbqt filename
        :param dock_config: a .txt configure filename. contains information about the docking box etc.
        :param dock_method: docking method, vina by default.
                            supported methods: vina
        :param log: log file. will be set to stdout if "stdout" is passed
        """
        self.mols = mols
        self.id = abs(hash(''.join(self.mols)))  # code will crash if id is negative. why? ask Inc. Schrodinger!
        self.protein = protein
        self.dock_config = dock_config
        self.dock_method = dock_method
        self.predictions = pd.DataFrame({}, index=[i for i in range(len(mols))])
        self.predictions_lock = threading.Lock()

        if log == "stdout":
            self.log = sys.stdout
        else:
            try:
                self.log = open(log, "w")
                sys.stdout = self.log
            except IOError:
                self.log = sys.stdout

    def run(self) -> pd.DataFrame:
        """
        Makes predictions.
        :return: a pd.DataFrame contains the predictions
        """
        self.prepare_sdf()

        os.chdir("./scratch")
        t_qikprop = threading.Thread(target=self.qikprop)
        t_dock = threading.Thread(target=self.dock_score)
        t_qikprop.start()
        t_dock.start()
        t_qikprop.join()
        t_dock.join()

        self.delete_scratch()
        os.chdir("..")

        return self.predictions

    def prepare_sdf(self):
        """
        Prepares sdf file of the molecules
        """
        # list of 3d structures
        m_list = [structure.SmilesStructure(m).get3dStructure(require_stereo=False) for m in self.mols]
        # TODO: fix this for cases with chirality

        with open(f"./scratch/{self.id}.sdf", 'w') as file:
            pass  # create empty file under scratch folder
        with structure.StructureWriter(f"./scratch/{self.id}.sdf") as writer:
            for m in m_list:
                writer.append(m)

    def qikprop(self):
        """
        Makes qikprop prediction. Requires sdf files of the molecules.
        TODO: use self._run_args(["qikprop",...]) may make abstraction clearer
        """
        print("-qikprop STDOUT/STDERR:")

        qp_solver = qp(f"{self.id}")
        ligandsobj = pipeio.Structures([f"{self.id}.sdf"])
        qp_solver.setInput(1, 'INPUT', ligandsobj)  # magic number 1? ask Inc. Schrodinger!
        qp_solver.setOutputName(1, 'OUTPUT')
        qp_solver.run(verbosity="quiet")
        qp_result = pd.read_csv(f"{self.id}-001.CSV")  # why -001? ask Inc. Schrodinger!

        self.predictions_lock.acquire()
        self.predictions = pd.concat([self.predictions, qp_result], axis=1, join="inner")
        self.predictions_lock.release()

        print("-qikprop DONE")

    def delete_scratch(self):
        # assume we are under ./scratch directory now
        file_list = os.listdir(".")
        for f in file_list:
            if re.match(str(self.id), f) is not None:
                os.remove(f)

    def _run_args(self, args):
        cp = subprocess.run(args, shell=True, capture_output=True, encoding="utf-8", errors="ignore")
        self.log.write(f"-{args[0]} STDOUT:\n" + cp.stdout)
        self.log.write(f"-{args[0]} STDERR:\n" + cp.stderr)

    def dock_score(self):
        """
        Computes docking score. Requires sdf file for the molecules in current directory.
        NOTE: everything is done under ./scratch directory
        TODO: how to specify the docking box in the config file?
        """

        self._prepare_docking_file()
        results = []

        for i in range(len(self.mols)):
            out_file = f"{self.id}_{i+1}_out.pdbqt"
            self._run_args(["vina",
                            "--receptor", f"../{self.protein}",
                            "--ligand", f"{self.id}_{i+1}.pdbqt",
                            "--out", out_file,
                            "--config", f"../{self.dock_config}"])
            # vina usage: see https://vina.scripps.edu/manual/#usage
            # i+1: ligands pdbqt files split by openbabel start at index 1
            # ../self.protein: ugliest code in this class. schrodinger api doesn't allow to change output filename
            # TODO: about schrodinger api

            result = self._parse_vina_out(out_file)
            results.append(result)

        self.predictions_lock.acquire()
        self.predictions['dock_score'] = results
        self.predictions_lock.release()

    @staticmethod
    def _parse_vina_out(out_file) -> float:
        # parse output .pdbqt file
        f = open(out_file)
        lines = f.readlines()
        f.close()
        try:
            line = lines[1]
            result = float(line.split(':')[1].split()[0])
        except IndexError:
            result = float('NaN')
        return result

    def _prepare_docking_file(self):

        # convert .pdb to .pdbqt
        if re.search(".pdbqt", self.protein) is None:
            protein_file = f"../{self.protein}"
            self._run_args(["obabel",
                            "-ipdb", protein_file,
                            "-opdbqt", "-O", f"{protein_file}qt",
                            "-p", "-xr"])
            # for -p: see https://openbabel.org/docs/dev/Command-line_tools/babel.html
            # for -xr: see http://openbabel.org/docs/current/FileFormats/Overview.html
            # -xr eliminates flexibility on side chains
            # TODO: fix this for flexible docking
            self.protein = f"{self.protein}qt"

        # convert .sdf ligands to .pdbqt
        self._run_args(["obabel",
                        "-isdf", f"{self.id}.sdf",
                        "-opdbqt", "-O", f"{self.id}_.pdbqt",
                        "-m"])
