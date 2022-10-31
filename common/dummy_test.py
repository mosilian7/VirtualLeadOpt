
from core.predictor import Predictor

p = Predictor(['Cn1c(CN2CCN(CC2)c3ccc(Cl)cc3)nc4ccccc14', 'COc1cc(OC)c(cc1NC(=O)CSCC(=O)O)S(=O)(=O)N2C(C)CCc3ccccc23'],
             "test/pdbdemo.pdb","test/configtest.txt")
p.run()



