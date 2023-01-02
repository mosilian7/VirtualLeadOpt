import os.path
import subprocess
import sys
from typing import List

from rdkit import Chem
from rdkit.Chem import Draw
import time
import math
import xmlrpc.client

# path of dependencies
VINA = "D:\\py_projects\\VirtualLeadOpt\\dependency\\vina"
OBABEL = "D:\\Program Files\\OpenBabel-3.1.1\\obabel"
MMPDB = "D:\\py_projects\\VirtualLeadOpt\\dependency\\mmpdb\\mmpdb"
QIKPROP = "D:\\Program Files\\Schrodinger\\qikprop"
G2G_TRANSLATION_MODEL = "D:\\py_projects\\VirtualLeadOpt\\dependency\\hgraph2graph\\ckpt\\translation\\model.8"
G2G_TRANSLATION_VOCAB = "D:\\py_projects\\VirtualLeadOpt\\dependency\\hgraph2graph\\data\\qed\\vocab.txt"


def run_args(args, logging=True, log=sys.stdout, simplify_bound: int = 40, remote_endpoint: str = None):
    """
    Run args with command line.
    :param args: list of arguments
    :param logging: if logging, stdout and stderr will be writen to log
    :param log: opened log file
    :param simplify_bound: if the number of lines in stdout is bigger than simplify_bound, then the output will be
                            simplified.
    :param remote_endpoint: endpoint of remote host. If remote_endpoint is None, then run_args will run locally.
    :return: stdout of the process
    """
    if remote_endpoint is None:
        completed_process = subprocess.run(args, shell=True, capture_output=True, encoding="utf-8", errors="ignore")
        out = completed_process.stdout
        err = completed_process.stderr
    else:
        server = xmlrpc.client.ServerProxy(remote_endpoint)
        response = server.run_args(args)
        out = response["stdout"]
        err = response["stderr"]

    if logging:
        host_info = "localhost" if remote_endpoint is None else remote_endpoint
        curr_time = time.strftime('%H:%M:%S', time.localtime(time.time()))
        if len(out) > 0:
            log.write(f"\033[32m[{curr_time}][{host_info}] {args[0]} STDOUT:\033[0m\n")
            so = str(out).split('\n')
            if len(so) < simplify_bound:
                log.write(out)
            else:
                for ln in so[:int(simplify_bound / 2)]:
                    log.write(ln + '\n')
                for i in range(3):
                    log.write('...\n')
                for ln in so[-int(simplify_bound / 2):]:
                    log.write(ln + '\n')
        if len(err) > 0:
            log.write(f"\033[31m[{curr_time}][{host_info}] {args[0]} STDERR:\033[0m\n" + err)

    return out


def visualize(mols: List[str]):
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


def geometric_mean(l: List[float], weight: List[float] = None) -> float:
    if weight is None:
        weight = [1] * len(l)
    numerator = sum([weight[i] * math.log(l[i]) for i in range(len(l))])
    denominator = sum(weight)
    return math.exp(numerator / denominator)


def sigmoid(x: float) -> float:
    return 1 / (1 + math.exp(-x))


def wsl_path_of(windows_path: str) -> str:
    windows_abs_path = os.path.abspath(windows_path)
    try:
        drive_letter, abs_path = windows_abs_path.split(":\\")
    except ValueError:
        raise ValueError("invalid windows path. maybe the path is not absolute")

    drive_letter = drive_letter.lower()
    mounted_path = '/'.join(abs_path.split('\\'))
    return '/'.join(['/mnt', drive_letter, mounted_path])


