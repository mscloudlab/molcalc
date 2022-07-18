import copy
import datetime
import logging
import pathlib

import pprint  # pretty print dictionaries

import models

import ppqm
from molcalc_lib import gamess_calculations, qchem_calculations
from ppqm import chembridge, misc
from ppqm.constants import COLUMN_COORDINATES

_logger = logging.getLogger("molcalc:pipe")


def calculation_pipeline(molinfo, settings):
    """

    Assumed that rdkit understands the molecule

    args:
        molinfo - dict
        settings -

    """

    # Read input
    molobj = molinfo["molobj"]
    sdfstr = molinfo["sdfstr"]
    hashkey = molinfo["hashkey"]

    scratch_dir = settings["scr.scr"]
    scratch_dir = pathlib.Path(scratch_dir)

    # TODO Get molecule names

    # Get that smile on your face
    try:
        smiles = chembridge.molobj_to_smiles(molobj, remove_hs=True)
    except Exception:
        smiles = chembridge.molobj_to_smiles(molobj)

    # Start respond message
    msg = {"smiles": smiles, "hashkey": hashkey}

    atoms = chembridge.get_atoms(molobj)
    _logger.info(f"{hashkey} '{smiles}' {atoms}")

    # Create new calculation
    calculation = models.GamessCalculation()

    # Switch to scrdir / hashkey
    hashdir = scratch_dir / hashkey
    hashdir.mkdir(parents=True, exist_ok=True)

    gamess_options = {
        "cmd": settings["gamess.rungms"],
        "gamess_scr": settings["gamess.scr"],
        "gamess_userscr": settings["gamess.userscr"],
        "scr": hashdir,
        "filename": hashkey,
    }

    qchem_options = copy.deepcopy(gamess_options)
    qchem_options['orca_cmd'] = settings['orca_cmd']

    # TODO Add error messages when gamess fails
    # TODO add timeouts for all gamess calls

    # Optimize molecule
    try:
        # software_to_use = 'orca'
        # match software_to_use:
        #     case 'orca':
        #         print('Hello, trying to use Orca!')
        #         properties = orca_calculations.optimize_coordinates(
        #             molobj, orca_options
        #         )
        #     case _:
        #         properties = gamess_calculations.optimize_coordinates(
        #             molobj, gamess_options
        #         )
        # orca_properties = orca_calculations.optimize_coordinates(
        #     molobj, orca_options
        # )
        properties_opt = qchem_calculations.optimize_coordinates(
            molobj, copy.deepcopy(qchem_options)
        )
        gamess_props_opt, orca_props_opt = properties_opt
        # print('\n'+80*'*')
        # print('orca_properties:\n\t')
        # print(orca_properties)
        # print('\n'+80*'*')
        # print('properties:\n\t')
        # print(properties)
        # print(80*'*')
    except Exception:
        # TODO Logger + rich should store these exceptions somewhere. One file
        # per exception for easy debugging.
        # TODO Should store SDF of the molecule if exception
        sdfstr = chembridge.molobj_to_sdfstr(molobj)
        _logger.error(f"{hashkey} OptimizationError", exc_info=True)
        _logger.error(sdfstr)
        gamess_props_opt = None

    if gamess_props_opt is None:
        return {
            "error": "Error g-80 - gamess optimization error",
            "message": "Error. Unable to optimize molecule",
        }, None

    if "error" in gamess_props_opt:
        return {
            "error": "Error g-93 - gamess optimization error known",
            "message": gamess_props_opt["error"],
        }, None

    if (
        COLUMN_COORDINATES not in gamess_props_opt
        or gamess_props_opt[COLUMN_COORDINATES] is None
    ):
        return {
            "error": "Error g-104 - gamess optimization error",
            "message": "Error. Unable to optimize molecule",
        }, None

    _logger.info(f"{hashkey} OptimizationSuccess")

    # Save and set coordinates
    coord = gamess_props_opt[ppqm.constants.COLUMN_COORDINATES]
    calculation.coordinates = misc.save_array(coord)
    calculation.enthalpy = gamess_props_opt[ppqm.constants.COLUMN_ENERGY]
    chembridge.molobj_set_coordinates(molobj, coord)

    # Optimization is finished, do other calculation async-like

    (
        properties_vib,
        properties_orb,
        properties_sol,
    ) = qchem_calculations.calculate_all_properties(molobj, qchem_options)

    gamess_props_vib, orca_props_vib = properties_vib
    gamess_props_orb, orca_props_orb = properties_orb
    print('\n' + 80*'*')
    print('properties_vib_orca:\n', pprint.pprint(orca_props_vib))
    print()
    print('properties_vib_gamess:\n', pprint.pprint(gamess_props_vib))
    # print(80*'*')
    # print('properties_orb:', properties_orb)
    # print(80*'*')
    # print('properties_sol:', properties_sol)
    print(80*'*' + '\n')

    # Check results

    if gamess_props_vib is None or "error" in gamess_props_vib:
        return {
            "error": "Error g-104 - gamess vibration error",
            "message": "Error. Unable to vibrate molecule",
        }, None

    _logger.info(f"{hashkey} VibrationSuccess")

    # TODO Make a custom reader and move this out of ppqm
    calculation.islinear = gamess_props_vib["linear"]
    calculation.vibjsmol = gamess_props_vib["jsmol"]
    calculation.vibfreq = misc.save_array(gamess_props_vib["freq"])
    calculation.vibintens = misc.save_array(gamess_props_vib["intens"])
    calculation.thermo = misc.save_array(gamess_props_vib["thermo"])

    if gamess_props_orb is None or "error" in gamess_props_orb:
        return {
            "error": "Error g-128 - gamess orbital error",
            "message": "Error. Unable to calculate molecular orbitals",
        }, None

    _logger.info(f"{hashkey} OrbitalsSuccess")
    calculation.orbitals = misc.save_array(gamess_props_orb["orbitals"])
    calculation.orbitalstxt = gamess_props_orb["stdout"]

    if properties_sol is None or "error" in properties_sol:

        # Is okay solvation didn't converge, just warn.
        _logger.warning(f"{hashkey} SolvationError")

    else:
        # 'charges', 'solvation_total', 'solvation_polar',
        # 'solvation_nonpolar', 'surface', 'total_charge', 'dipole',
        # 'dipole_total'
        _logger.info(f"{hashkey} SolvationSuccess")

        charges = properties_sol["charges"]
        calculation.charges = misc.save_array(charges)
        calculation.soltotal = properties_sol["solvation_total"]
        calculation.solpolar = properties_sol["solvation_polar"]
        calculation.solnonpolar = properties_sol["solvation_nonpolar"]
        calculation.solsurface = properties_sol["surface"]
        calculation.soldipole = misc.save_array(properties_sol["dipole"])
        calculation.soldipoletotal = properties_sol["dipole_total"]

        # Save mol2 fmt
        mol2 = chembridge.molobj_to_mol2(molobj, charges=charges)
        calculation.mol2 = mol2

    # Saveable sdf and reset title
    sdfstr = chembridge.molobj_to_sdfstr(molobj)
    sdfstr = chembridge.clean_sdf_header(sdfstr)

    # Get a 2D Picture
    svgstr = chembridge.molobj_to_svgstr(molobj, removeHs=True, use_2d=True)

    # Success, store results database
    calculation.smiles = smiles
    calculation.hashkey = hashkey
    calculation.sdf = sdfstr
    calculation.svg = svgstr
    calculation.created = datetime.datetime.now()

    return msg, calculation


def update_smiles_counter(request, smiles):

    # Add smiles to counter
    countobj = (
        request.dbsession.query(models.Counter)
        .filter_by(smiles=smiles)
        .first()
    )

    if countobj is None:
        counter = models.Counter()
        counter.smiles = smiles
        counter.count = 1
        request.dbsession.add(counter)
    else:
        countobj.count += 1

    return
