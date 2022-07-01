import copy
import pathlib

import pytest
from context import CONFIG, RESOURCES, SCR
from rdkit import Chem
from rdkit.Chem import AllChem

import ppqm
from molcalc_lib import orca_calculations
from ppqm import chembridge

ORCA_OPTIONS = {
    "scr": SCR,
    "cmd": CONFIG["orca"]["cmd"],
}


TEST_SMILES = ["C"]

TEST_SMILES_COORD = [
    ("CCC", -23.62341),
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


@pytest.mark.parametrize("smiles, test_energy", TEST_SMILES_COORD)
def test_optimize_coordinates(smiles, test_energy):

    molobj = _prepare_molobj(smiles)
    properties = orca_calculations.optimize_coordinates(
        molobj, ORCA_OPTIONS
    )

    assert properties[ppqm.constants.COLUMN_ENERGY] == pytest.approx(
        test_energy
    )


@pytest.mark.parametrize("smiles", TEST_SMILES_SOLVATION)
def test_calculate_solvation(smiles):
    pass


@pytest.mark.parametrize("smiles", TEST_SMILES)
def test_calculate_all_properties(smiles):
    pass


@pytest.mark.parametrize("filename", TEST_ERROR_SDF)
def test_error_smiles(tmpdir, filename):
    pass


if __name__ == "__main__":

    tmpdir = pathlib.Path(".tmptest")
    tmpdir.mkdir(parents=True, exist_ok=True)

    # test_optimize_coordinates("C", 5.0)
    # test_calculate_all_properties("C")
    # test_error_smiles(tmpdir, TEST_ERROR_SDF[0])
