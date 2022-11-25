import subprocess
import sys
from rdkit import Chem
from rdkit.Chem import Draw
import time
import math

# path of dependencies
VINA = "D:\\py_projects\\VirtualLeadOpt\\dependency\\vina"
OBABEL = "D:\\Program Files\\OpenBabel-3.1.1\\obabel"
MMPDB = "D:\\py_projects\\VirtualLeadOpt\\dependency\\mmpdb\\mmpdb"
QIKPROP = "D:\\Program Files\\Schrodinger\\qikprop"


def run_args(args, logging=True, log=sys.stdout, simplify_bound: int = 40):
    """
    Run args with command line.
    :param args: list of arguments
    :param logging: if logging, stdout and stderr will be writen to log
    :param log: opened log file
    :param simplify_bound: if the number of lines in stdout is bigger than simplify_bound, then the output will be
                            simplified.
    :return: stdout of the process
    """
    cp = subprocess.run(args, shell=True, capture_output=True, encoding="utf-8", errors="ignore")

    if logging:
        curr_time = time.strftime('%H:%M:%S', time.localtime(time.time()))
        if len(cp.stdout) > 0:
            log.write(f"\033[32m[{curr_time}] {args[0]} STDOUT:\033[0m\n")
            so = str(cp.stdout).split('\n')
            if len(so) < simplify_bound:
                log.write(cp.stdout)
            else:
                for ln in so[:int(simplify_bound/2)]:
                    log.write(ln + '\n')
                for i in range(3):
                    log.write('...\n')
                for ln in so[-int(simplify_bound/2):]:
                    log.write(ln + '\n')
        if len(cp.stderr) > 0:
            log.write(f"\033[31m[{curr_time}] {args[0]} STDERR:\033[0m\n" + cp.stderr)
    return cp.stdout


def visualize(mols: list):
    mlist = [Chem.MolFromSmiles(str(m)) for m in mols]
    return Draw.MolsToImage(mlist, subImgSize=(300, 300))


def split_list(l: list, share: int) -> list:
    if share > len(l):
        share = len(l)
    result = []
    len_per_share = int(len(l)/share)
    for i in range(share):
        result.append(l[i * len_per_share: (i + 1) * len_per_share])
    for i in range(len_per_share * share, len(l)):
        result[i - len_per_share * share].append(l[i])
    return result


def geometric_mean(l: list, weight: list = None) -> float:
    if weight is None:
        weight = [1] * len(l)
    numerator = sum([weight[i] * math.log(l[i]) for i in range(len(l))])
    denominator = sum(weight)
    return math.exp(numerator / denominator)


def sigmoid(x: float) -> float:
    return 1 / (1 + math.exp(-x))
