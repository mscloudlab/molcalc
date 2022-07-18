from context import CONFIG, RESOURCES, SCR

from molcalc import pipelines
from ppqm import chembridge

GAMESS_OPTIONS = {
    "scr": SCR,
    "cmd": CONFIG["gamess"].get("rungms"),
    "gamess_scr": CONFIG["gamess"].get("scr"),
    "gamess_userscr": CONFIG["gamess"].get("userscr"),
    "debug": True,
}

ORCA_OPTIONS = {
    "scr": SCR,
    "cmd": CONFIG["orca"]["cmd"],
}


def test_pipelines():

    settings = dict()
    settings["scr.scr"] = SCR
    settings["gamess.rungms"] = GAMESS_OPTIONS["cmd"]
    settings["gamess.scr"] = GAMESS_OPTIONS["gamess_scr"]
    settings["gamess.userscr"] = GAMESS_OPTIONS["gamess_userscr"]
    settings["orca_cmd"] = ORCA_OPTIONS["cmd"]

    print(settings)

    sdf = """


  2  1  0  0  0  0  0  0  0  0999 V2000
    0.7500    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  3  0
M  END """
#
#     sdf = """
#
#
#   5  4  0  0  0  0  0  0  0  0999 V2000
#     0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#     0.0000   -0.8900   -0.6293 H   0  0  0  0  0  0  0  0  0  0  0  0
#     0.0000    0.8900   -0.6293 H   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.8900   -0.0000    0.6293 H   0  0  0  0  0  0  0  0  0  0  0  0
#     0.8900   -0.0000    0.6293 H   0  0  0  0  0  0  0  0  0  0  0  0
#   1  2  1  0  0  0  0
#   1  3  1  0  0  0  0
#   1  4  1  0  0  0  0
#   1  5  1  0  0  0  0
# M  END
# $$$$
#     """

    # smi = "N#N"
    # molobj = Chem.MolFromSmiles(smi)
    # AllChem.Compute2DCoords(molobj)
    molobj = chembridge.sdfstr_to_molobj(sdf)
    sdf = chembridge.molobj_to_sdfstr(molobj)
    print(sdf)

    # NOTE: The string value corresponding to "hashkey" in the molecular_info
    # dictionary will serve as the file/dir prefix used by pipelines.py
    molecule_info = {"sdfstr": sdf, "molobj": molobj, "hashkey": "TEST"}

    results = pipelines.calculation_pipeline(molecule_info, settings)

    print(results)

    # Necessary to ensure that pytest doesn't pass when calculation_pipeline()
    #   returns with an error message
    # msg, calculation = results
    # assert 'error' not in msg #and calculation is not None

    return


# def test_partial_pipeline():
#     """
#
#     Mg will fail solvation calculation
#
#     """
#
#     filename = RESOURCES / "mg.sdf"
#     with open(filename, "r") as f:
#         sdfstr = f.read()
#
#     settings = dict()
#     settings["scr.scr"] = SCR
#     settings["gamess.rungms"] = GAMESS_OPTIONS["cmd"]
#     settings["gamess.scr"] = GAMESS_OPTIONS["gamess_scr"]
#     settings["gamess.userscr"] = GAMESS_OPTIONS["gamess_userscr"]
#
#     molobj = chembridge.sdfstr_to_molobj(sdfstr)
#     sdf = chembridge.molobj_to_sdfstr(molobj)
#
#     molecule_info = {"sdfstr": sdf, "molobj": molobj, "hashkey": "TEST"}
#
#     results = pipelines.calculation_pipeline(molecule_info, settings)
#
#     # TODO results is valid
#     # TODO results does not contain solvation
#
#     print(results)
#
#     return


if __name__ == "__main__":
    test_pipelines()
