import copy
import pathlib

import pytest
from context import CONFIG, RESOURCES, SCR
from rdkit import Chem
from rdkit.Chem import AllChem

from ppqm import chembridge, orca, tasks, units
from molcalc_lib import orca_calculations

ORCA_OPTIONS = {
    "scr": SCR,
    "cmd": CONFIG["orca"]["cmd"],
}


TEST_SMILES = ["C"]

TEST_SMILES_COORD = [
    ("CCC", -23.62341),
]
TEST_ENERGIES_PM3 = [
    ("O", -11.935809225486 * units.hartree_to_kcalmol),
    ("CC", -12.124869353328 * units.hartree_to_kcalmol),
    ("[NH4+]", -7.972867788142 * units.hartree_to_kcalmol),
]

TEST_ERROR_SDF = ["wrong_methane.sdf", "wrong_benzene.sdf"]


TEST_SMILES_SOLVATION = [
    "C",
    "CCCBr",
    "C[NH3+]",
]


def _prepare_molobj(smiles):
    """
    Helper function for getting 3D coordinates from SMILES
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    _ = AllChem.EmbedMolecule(mol)
    _ = AllChem.UFFOptimizeMolecule(mol)
    return mol


@pytest.mark.parametrize("smiles, test_energy", TEST_ENERGIES_PM3)
def test_optimize_coordinates(smiles, test_energy):

    # molobj = _prepare_molobj(smiles)
    # The two lines below replace the line above
    molobj = chembridge.smiles_to_molobj(smiles)
    molobj = tasks.generate_conformers(molobj, n_conformers=1)

    assert molobj is not None
    properties = orca_calculations.optimize_coordinates(
        molobj, ORCA_OPTIONS
    )

    print('\n' + 80*'*')
    print(properties)
    print(80*'*')
    assert properties[orca.COLUMN_SCF_ENERGY] == pytest.approx(
        test_energy,
        1.0e-7 #2.0e-3
    )


if __name__ == "__main__":

    tmpdir = pathlib.Path(".tmptest")
    tmpdir.mkdir(parents=True, exist_ok=True)

    # test_optimize_coordinates("C", 5.0)
    # test_calculate_all_properties("C")
    # test_error_smiles(tmpdir, TEST_ERROR_SDF[0])
