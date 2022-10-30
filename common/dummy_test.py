from openbabel import openbabel
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("smi", "mdl")

mol = openbabel.OBMol()
if obConversion.ReadFile(mol, "smiDemo.smi"):
    print("hello")

print('Should not print 0 (atoms)')
print(mol.NumAtoms())

mol.AddHydrogens()
print('Should not print 0 (atoms) after adding hydrogens')
print(mol.NumAtoms())


