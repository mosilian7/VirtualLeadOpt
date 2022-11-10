import os
import re
import sys
import subprocess
import threading

import pandas as pd
import schrodinger.pipeline.pipeio as pipeio
import common.util as util
from schrodinger import structure
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import QED
from schrodinger.pipeline.stages.qikprop import QikPropStage as qp


class Predictor(threading.Thread):
    """
        predict the properties of a set of molecules
        NOTE: public methods should run under ./scratch directory
    """

    def __init__(self, mols: list, protein: str, dock_config: str,
                 dock_method: str = "vina", dock_cpu: int = 4, log: str = "stdout"):
        """
        initializer. make sure that open babel is installed and added to PATH
        :param mols: a list of SMILES string
        :param protein: a .pdb or .pdbqt filename (relative path)
        :param dock_config: a .txt configure filename. contains information about the docking box etc. (relative path)
        :param dock_method: docking method, vina by default.
                            supported methods: vina
        :param log: log file. will be set to stdout if "stdout" is passed
        """
        super(Predictor, self).__init__()
        self.mols = mols
        self.id = abs(hash(''.join(self.mols)))  # code will crash if id is negative. why? ask Inc. Schrodinger!
        self.protein = protein
        self.dock_config = dock_config
        self.dock_method = dock_method
        self.dock_cpu = dock_cpu
        self.predictions = pd.DataFrame({}, index=[i for i in range(len(mols))])
        self.predictions_lock = threading.Lock()
        self.mol_list = [Chem.MolFromSmiles(smiles) for smiles in self.mols]

        if log == "stdout":
            self.log = sys.stdout
        else:
            try:
                self.log = open(log, "w")
                sys.stdout = self.log
            except IOError:
                self.log = sys.stdout

    def run(self) -> None:
        """
            Makes predictions.
            TODO: split this
            NOTE: Should be nested in
                os.chdir("./scratch")
                os.chdir("..")
        """
        self.prepare_sdf()

        t_qikprop = threading.Thread(target=self.qikprop)
        t_dock = threading.Thread(target=self.dock_score)
        t_qikprop.start()
        t_dock.start()
        t_qikprop.join()
        t_dock.join()

        self.delete_scratch()

    def prepare_sdf(self) -> None:
        """
            Prepares sdf file of the molecules
        """
        m_list = [Chem.MolFromSmiles(smi) for smi in self.mols]
        hm_list = [Chem.AddHs(m) for m in m_list]
        writer = Chem.SDWriter(f'{self.id}.sdf')
        for mol in hm_list:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol, 1000)
            writer.write(mol)
        writer.close()

    def qed(self) -> None:
        qed_list = [QED.qed(m) for m in self.mol_list]
        self.predictions_lock.acquire()
        self.predictions['qed'] = qed_list
        self.predictions_lock.release()

    def qikprop(self) -> None:
        """
            Makes qikprop prediction. Requires sdf files of the molecules.
        """
        util.run_args([util.QIKPROP,
                       "-fast", "-nosa", "-WAIT",
                       f"{self.id}.sdf"], log=self.log)
        qp_result = pd.read_csv(f"{self.id}.CSV")

        self.predictions_lock.acquire()
        self.predictions = pd.concat([self.predictions, qp_result], axis=1, join="inner")
        self.predictions_lock.release()

    def num_arom_ring(self) -> None:
        """
            Counts the number of aromatic rings.
        """
        num_arom_ring_list = [Chem.rdMolDescriptors.CalcNumAromaticRings(m) for m in self.mol_list]

        self.predictions_lock.acquire()
        self.predictions['num_arom_ring'] = num_arom_ring_list
        self.predictions_lock.release()

    def delete_scratch(self) -> None:
        """
            Deletes scratch files.
        """

        # assume we are under ./scratch directory now
        file_list = os.listdir(".")
        for f in file_list:
            if re.match(str(self.id), f) is not None:
                os.remove(f)

    def _run_args(self, args):
        # deprecated. can be replaced by util.ren_args(args, log=self.log)
        cp = subprocess.run(args, shell=True, capture_output=True, encoding="utf-8", errors="ignore")
        self.log.write(f"-{args[0]} STDOUT:\n" + cp.stdout)
        self.log.write(f"-{args[0]} STDERR:\n" + cp.stderr)

    def dock_score(self) -> None:
        """
        Computes docking score. Requires sdf file for the molecules in current directory.
        NOTE: everything is done under ./scratch directory
        TODO: how to specify the docking box in the config file?
        """

        self._prepare_docking_file()
        results = []

        for i in range(len(self.mols)):
            out_file = f"{self.id}_{i+1}_out.pdbqt"
            util.run_args(["vina",
                            "--receptor", f"../{self.protein}",
                            "--ligand", f"{self.id}_{i+1}.pdbqt",
                            "--out", out_file,
                            "--cpu", str(self.dock_cpu),
                            "--config", f"../{self.dock_config}"], log=self.log)
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

    def _prepare_docking_file(self) -> None:
        # prepare .pdbqt files

        # convert .pdb to .pdbqt
        if re.search(".pdbqt", self.protein) is None:
            protein_file = f"../{self.protein}"
            util.run_args(["obabel",
                            "-ipdb", protein_file,
                            "-opdbqt", "-O", f"{protein_file}qt",
                            "-p", "-xr"], log=self.log)
            # for -p: see https://openbabel.org/docs/dev/Command-line_tools/babel.html
            # for -xr: see http://openbabel.org/docs/current/FileFormats/Overview.html
            # -xr eliminates flexibility on side chains
            # TODO: fix this for flexible docking
            self.protein = f"{self.protein}qt"

        # convert .sdf ligands to .pdbqt
        util.run_args(["obabel",
                        "-isdf", f"{self.id}.sdf",
                        "-opdbqt", "-O", f"{self.id}_.pdbqt",
                        "-m"], log=self.log)


class PredictorWrapper:
    """
    Wrapper of predictor class.
    """

    def __init__(self, protein: str, dock_config: str, dock_method: str = "vina",
                 dock_cpu: int = 1, log: str = "stdout"):
        self.protein = protein
        self.dock_config = dock_config
        self.dock_method = dock_method
        self.dock_cpu = dock_cpu
        self.log = log

    def predictor_instance(self, mols: list) -> Predictor:
        """
        Creates a Predictor instance which makes prediction on properties of a list of molecules.
        :param mols: a list of SMILES strings
        :return: a Predictor instance
        """
        return Predictor(mols, self.protein, self.dock_config, self.dock_method, self.dock_cpu, self.log)

