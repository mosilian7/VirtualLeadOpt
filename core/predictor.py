import math
import os
import re
import subprocess
import sys
import threading
from typing import List, Dict

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

    def __init__(self, mols: List[str], protein: str, dock_config: str, mask: dict = None,
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
        self.protein = os.path.abspath(protein)
        self.mask = {} if mask is None else mask
        self.dock_config = os.path.abspath(dock_config)
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

    def _run_args(self, args: List[str]):
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

        if self.dock_method == "vina":
            ligands = []
            for i in range(len(self.mols)):
                if "dock_score" in self.mask and self.mask["dock_score"][i]:
                    continue
                ligands.append(f"{self.id}_{i + 1}.pdbqt")

            util.run_args([util.VINA,
                           "--receptor", self.protein,
                           "--batch", *ligands,
                           "--dir", ".",
                           "--cpu", str(self.dock_cpu),
                           "--config", self.dock_config], log=self.log)
            # vina usage: see https://vina.scripps.edu/manual/#usage
            # i+1: ligands pdbqt files split by openbabel start at index 1

            for i in range(len(self.mols)):
                if "dock_score" in self.mask and self.mask["dock_score"][i]:
                    results.append(float('nan'))
                    continue
                result = self._parse_vina_out(f"{self.id}_{i + 1}_out.pdbqt")
                results.append(result)

        elif self.dock_method == "gnina":
            # gnina must run on Linux system

            # file paths on remote machine.
            # TODO: modify this for generalization
            ligands_path = util.wsl_path_of(f"{self.id}_gnina_ligands.sdf")
            protein_path = util.wsl_path_of(self.protein)
            config_path = util.wsl_path_of(self.dock_config)
            output_path = util.wsl_path_of(f"{self.id}_gnina_out.sdf")

            wsl_ip = util.run_args(["wsl", "hostname", "-I"], logging=False).strip()
            util.run_args(["gnina",
                           "-r", protein_path,
                           "-l", ligands_path,
                           "--config", config_path,
                           "--cpu", str(self.dock_cpu),
                           "-o", output_path],
                          log=self.log, remote_endpoint=f"http://{wsl_ip}:8111")
            results = self._parse_gnina_out(f"{self.id}_gnina_out.sdf")

            for i in range(len(self.mols)):
                if "dock_score" in self.mask and self.mask["dock_score"][i]:
                    results.insert(i, float('nan'))

        assert(len(results) == len(self.mols))
        self.predictions_lock.acquire()
        self.predictions['dock_score'] = results
        self.predictions_lock.release()

    def dock_score_old(self) -> None:
        """
        Computes docking score. Requires sdf file for the molecules in current directory.
        NOTE: everything is done under ./scratch directory
        """

        self._prepare_docking_file()
        results = []

        for i in range(len(self.mols)):
            if "dock_score" in self.mask and self.mask["dock_score"][i]:
                results.append(float('nan'))
                continue
            out_file = f"{self.id}_{i + 1}_out.pdbqt"
            util.run_args([util.VINA,
                           "--receptor", self.protein,
                           "--ligand", f"{self.id}_{i + 1}.pdbqt",
                           "--out", out_file,
                           "--cpu", str(self.dock_cpu),
                           "--config", self.dock_config], log=self.log)
            # vina usage: see https://vina.scripps.edu/manual/#usage
            # i+1: ligands pdbqt files split by openbabel start at index 1

            result = self._parse_vina_out(out_file)
            results.append(result)

        self.predictions_lock.acquire()
        self.predictions['dock_score'] = results
        self.predictions_lock.release()

    @staticmethod
    def _parse_gnina_out(out_file: str) -> List[float]:
        # parser is a state machine
        flag = 0
        results = []
        last_smiles = ""

        with open(out_file, "r") as f:
            line = f.readline()
            while line != '':
                if flag == 0:
                    if re.search("<smiles>", line) is not None:
                        flag = 1
                elif flag == 1:
                    smiles = str(line.strip())
                    if smiles != last_smiles:
                        flag = 2
                    else:
                        flag = 0
                    last_smiles = smiles
                elif flag == 2:
                    if re.search("<CNNaffinity>", line) is not None:
                        flag = 3
                elif flag == 3:
                    cnn_affinity = float(line.strip())
                    results.append(-cnn_affinity)
                    flag = 0

                line = f.readline()

        return results



    @staticmethod
    def _parse_vina_out(out_file: str) -> float:
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
            util.run_args([util.OBABEL,
                           "-ipdb", self.protein,
                           "-opdbqt", "-O", f"{self.protein}qt",
                           "-p", "-xr"], log=self.log)
            # for -p: see https://openbabel.org/docs/dev/Command-line_tools/babel.html
            # for -xr: see http://openbabel.org/docs/current/FileFormats/Overview.html
            # -xr eliminates flexibility on side chains
            # TODO: fix this for flexible docking
            self.protein = os.path.abspath(f"{self.protein}qt")

        if self.dock_method == "vina":
            # convert .sdf ligands to .pdbqt
            util.run_args([util.OBABEL,
                           "-isdf", f"{self.id}.sdf",
                           "-opdbqt", "-O", f"{self.id}_.pdbqt",
                           "-m"], log=self.log)

        elif self.dock_method == "gnina":
            # prepare a .sdf file for all unmasked molecules
            writer = Chem.SDWriter(f'{self.id}_gnina_ligands.sdf')
            for i in range(len(self.mols)):
                if "dock_score" in self.mask and self.mask["dock_score"][i]:
                    continue
                mol = Chem.MolFromSmiles(self.mols[i])
                hmol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(hmol, AllChem.ETKDG())
                AllChem.UFFOptimizeMolecule(hmol, 1000)
                hmol.SetProp("smiles", self.mols[i])
                writer.write(hmol)
            writer.close()


class PredictorWrapper:
    """
    Wrapper of predictor class.
    """

    def __init__(self, protein: str, dock_config: str, dock_method: str = "vina",
                 dock_cpu: int = 1, log: str = "stdout"):
        self.protein = os.path.abspath(protein)
        self.dock_config = os.path.abspath(dock_config)
        self.dock_method = dock_method
        self.dock_cpu = dock_cpu
        self.log = log

    def predictor_instance(self, mols: List[str], mask: Dict[str, List[bool]] = None) -> Predictor:
        """
        Creates a Predictor instance which makes prediction on properties of a list of molecules.
        :param mols: a list of SMILES strings
        :param mask: masking dictionary
        :return: a Predictor instance
        """
        return Predictor(mols, self.protein, self.dock_config, mask, self.dock_method, self.dock_cpu, self.log)
