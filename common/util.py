import subprocess
import sys
from rdkit import Chem
from rdkit.Chem import Draw


VINA = "D:\\Program Files (x86)\\The Scripps Research Institute\\Vina\\vina"
OBABEL = "D:\\Program Files\\OpenBabel-3.1.1\\obabel"
MMPDB = "D:\\py_projects\\virtual_lead_opt\\dependency\\mmpdb\\mmpdb"


def run_args(args, logging=True, log=sys.stdout):
    """
    Run args with command line.
    :param args: list of arguments
    :param logging: if logging, stdout and stderr will be writen to log
    :param log: opened log file
    :return: stdout of the process
    """
    cp = subprocess.run(args, shell=True, capture_output=True, encoding="utf-8", errors="ignore")

    if logging:
        log.write(f"-{args[0]} STDOUT:\n" + cp.stdout)
        log.write(f"-{args[0]} STDERR:\n" + cp.stderr)
    return cp.stdout


def visualize(mols: list):
    mlist = [Chem.MolFromSmiles(str(m)) for m in mols]
    return Draw.MolsToImage(mlist, subImgSize=(300, 300))

