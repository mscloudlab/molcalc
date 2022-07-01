import copy
import logging
from multiprocessing import Pipe, Process

import ppqm

_logger = logging.getLogger("molcalc:calc")

MAX_TIME = 10  # seconds


def optimize_coordinates(molobj, orca_options):

    calculation_options = {"pm3": None}

    # FIX
    orca_options.get("filename", None)

    calc_obj = ppqm.orca.OrcaCalculator(**orca_options)
    results = calc_obj.calculate(molobj, calculation_options)

    properties = results[0]

    return properties

