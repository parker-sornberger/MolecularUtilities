# MolecularUtilities


**BestMol**
* Utility wrapper for RDKit

**BestMol Examples**

```python
from molutils import Molecule

molecule = Molecule("C1=CC=C(C2=C1C5=C(C3=C2C=CC(=C3)C4=CC=CC=C4)C=CC(=C5)C6=CC=CC=C6)C7=CC=CC=C7", name = "Your cool molecule")

lgfr_molecule = molecule.get_LGFR_BM()

print(F"Molecule SMILES: {molecule.smiles}\nLGFR backbone SMILES: {lgfr_molecule.smiles}")

other_molecule = Molecule("C1=CC=CC2=C1C=CC=C2")

print(F"Molecule has this substructure {molecule.has(other_molecule.mol)}")

print(F"This molecule, {molecule['name']}, has {len(molecule)} atoms and {molecule['all_rings']} rings")

```




**MolInference**
* A wrapper for Molecule where a few properties are prebaked for convenience



