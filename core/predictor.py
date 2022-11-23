import math
import os
import re
import subprocess
import sys
import threading

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import QED
from rdkit.Chem import RDConfig

import core.util as util

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


class Predictor(threading.Thread):
    """
        predict the properties of a set of molecules
        NOTE: public methods should run under ./scratch directory
    """

    def __init__(self, mols: list, protein: str, dock_config: str, mask: dict = None,
                 dock_method: str = "vina", dock_cpu: int = 4, log: str = "stdout"):
        """
        initializer. make sure that open babel is installed and added to PATH
        :param mols: a list of SMILES string. Should never be empty!
        :param protein: a .pdb or .pdbqt filename (relative path)
        :param dock_config: a .txt configure filename. contains information about the docking box etc. (relative path)
        :param mask: masking dictionary
        :param dock_method: docking method, vina by default.
                            supported methods: vina
        :param log: log file. will be set to stdout if "stdout" is passed
        """
        super(Predictor, self).__init__()
        self.mols = mols
        self.id = abs(hash(''.join(self.mols)))  # code will crash if id is negative.
        self.protein = protein
        self.mask = {} if mask is None else mask
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
            Makes predictions. Deprecated?
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

    def before_dock(self) -> None:
        self.qed()
        self.qikprop()
        self.sa_score()

    def qed(self) -> None:
        """
            Calculates QED score.
            see https://www.nature.com/articles/nchem.1243 Quantifying the chemical beauty of drugs
        """
        qed_list = [QED.qed(m) for m in self.mol_list]

        assert (len(qed_list) == len(self.mols))
        self.predictions_lock.acquire()
        self.predictions['qed'] = qed_list
        self.predictions_lock.release()

    @staticmethod
    def _fix_bad_csv(csv_file: str):
        with open(csv_file, "r") as f:
            ls = f.readlines()
        for i in range(len(ls)):
            ls[i] = ls[i].replace(',' * 500, ',' * 51)
        with open(csv_file, "w") as f:
            f.writelines(ls)

    def qikprop(self) -> None:
        """
            Makes qikprop prediction. Requires sdf files of the molecules.
        """
        util.run_args([util.QIKPROP,
                       "-fast", "-nosa", "-WAIT",
                       f"{self.id}.sdf"], log=self.log)
        qp_result = pd.read_csv(f"{self.id}.CSV", on_bad_lines=lambda x: x, engine='python')

        if math.isnan(qp_result["QPlogHERG"][0]):
            # fix and retry
            self._fix_bad_csv(f"{self.id}.CSV")
            qp_result = pd.read_csv(f"{self.id}.CSV", on_bad_lines=lambda x: x, engine='python')

        assert (qp_result.shape[0] == len(self.mols))
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

    def sa_score(self) -> None:
        sa_scores = [sascorer.calculateScore(m) for m in self.mol_list]

        self.predictions_lock.acquire()
        self.predictions['sa_score'] = sa_scores
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

        ligands = []
        for i in range(len(self.mols)):
            if "dock_score" in self.mask and self.mask["dock_score"][i]:
                continue
            ligands.append(f"{self.id}_{i + 1}.pdbqt")

        util.run_args([util.VINA,
                       "--receptor", f"../{self.protein}",
                       "--batch", *ligands,
                       "--dir", ".",
                       "--cpu", str(self.dock_cpu),
                       "--config", f"../{self.dock_config}"], log=self.log)
        # vina usage: see https://vina.scripps.edu/manual/#usage
        # i+1: ligands pdbqt files split by openbabel start at index 1
        # ../self.protein: ugliest code in this class. schrodinger api doesn't allow to change output filename

        results = []

        for i in range(len(self.mols)):
            if "dock_score" in self.mask and self.mask["dock_score"][i]:
                results.append(float('nan'))
                continue
            result = self._parse_vina_out(f"{self.id}_{i + 1}_out.pdbqt")
            results.append(result)

        assert (len(results) == len(self.mols))
        self.predictions_lock.acquire()
        self.predictions['dock_score'] = results
        self.predictions_lock.release()

    def dock_score_old(self) -> None:
        """
        Computes docking score. Requires sdf file for the molecules in current directory.
        NOTE: everything is done under ./scratch directory
        TODO: how to specify the docking box in the config file?
        """

        self._prepare_docking_file()
        results = []

        for i in range(len(self.mols)):
            if "dock_score" in self.mask and self.mask["dock_score"][i]:
                results.append(float('nan'))
                continue
            out_file = f"{self.id}_{i + 1}_out.pdbqt"
            util.run_args([util.VINA,
                           "--receptor", f"../{self.protein}",
                           "--ligand", f"{self.id}_{i + 1}.pdbqt",
                           "--out", out_file,
                           "--cpu", str(self.dock_cpu),
                           "--config", f"../{self.dock_config}"], log=self.log)
            # vina usage: see https://vina.scripps.edu/manual/#usage
            # i+1: ligands pdbqt files split by openbabel start at index 1
            # ../self.protein: ugliest code in this class. schrodinger api doesn't allow to change output filename

            result = self._parse_vina_out(out_file)
            results.append(result)

        self.predictions_lock.acquire()
        self.predictions['dock_score'] = results
        self.predictions_lock.release()

    @staticmethod
    def _parse_vina_out(out_file) -> float:
        # parse output .pdbqt file
        try:
            f = open(out_file)
            lines = f.readlines()
            f.close()
            try:
                line = lines[1]
                result = float(line.split(':')[1].split()[0])
            except IndexError:
                result = float('NaN')
            return result
        except FileNotFoundError:
            # when vina cannot process the input
            return float('nan')

    def _prepare_docking_file(self) -> None:
        # prepare .pdbqt files

        # convert .pdb to .pdbqt
        if re.search(".pdbqt", self.protein) is None:
            protein_file = f"../{self.protein}"
            util.run_args([util.OBABEL,
                           "-ipdb", protein_file,
                           "-opdbqt", "-O", f"{protein_file}qt",
                           "-p", "-xr"], log=self.log)
            # for -p: see https://openbabel.org/docs/dev/Command-line_tools/babel.html
            # for -xr: see http://openbabel.org/docs/current/FileFormats/Overview.html
            # -xr eliminates flexibility on side chains
            # TODO: fix this for flexible docking
            self.protein = f"{self.protein}qt"

        # convert .sdf ligands to .pdbqt
        util.run_args([util.OBABEL,
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

    def predictor_instance(self, mols: list, mask: dict = None) -> Predictor:
        """
        Creates a Predictor instance which makes prediction on properties of a list of molecules.
        :param mols: a list of SMILES strings
        :param mask: masking dictionary
        :return: a Predictor instance
        """
        return Predictor(mols, self.protein, self.dock_config, mask, self.dock_method, self.dock_cpu, self.log)
