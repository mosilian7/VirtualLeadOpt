import os
import re
import sys
import subprocess
import threading
import time

import pandas as pd
import schrodinger.pipeline.pipeio as pipeio
import common.util as util
from schrodinger import structure
from schrodinger.pipeline.stages.qikprop import QikPropStage as qp
from core.predictor import Predictor


p = Predictor(['Cn1c(CN2CCN(CC2)c3ccc(Cl)cc3)nc4ccccc14', 'COc1cc(OC)c(cc1NC(=O)CSCC(=O)O)S(=O)(=O)N2C(C)CCc3ccccc23'],
             "test/pdbdemo.pdb","test/configtest.txt")


class foo(threading.Thread):
    def __init__(self, size):
        super(foo, self).__init__()
        self.size = size

    def run(self):
        time.sleep(1)
        self.size += 1

    def bar(self):
        print(self.size)

threads = [None for i in range(10)]

for i in range(10):
    threads[i] = foo(i)
    threads[i].start()


for i in range(10):

    threads[i].join()
    threads[i].bar()
