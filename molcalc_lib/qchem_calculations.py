import copy
import logging
from collections import ChainMap
from multiprocessing import Pipe, Process

from molcalc_lib import gamess_calculations, orca_calculations
import ppqm

_logger = logging.getLogger("molcalc:calc")

MAX_TIME = 10  # seconds


def _get_options(existing_options, tmp_path):
    orca_options = {"scr": tmp_path, "n_cores": 1, "memory": 2}
    options_prime = ChainMap(orca_options, existing_options)
    options_prime = dict(options_prime)
    return options_prime


def optimize_coordinates(molobj, qchem_options, engine='gamess'):
    options_prime = copy.deepcopy(qchem_options)

    if engine == 'gamess':
        options_prime['cmd'] = options_prime['gamess.cmd']
        properties = gamess_calculations.optimize_coordinates(molobj, options_prime)
    elif engine == 'orca':
        options_prime['cmd'] = options_prime['orca.cmd']
        properties = orca_calculations.optimize_coordinates(molobj, options_prime)
    else:
        raise ValueError(f'Error: keyword argument engine = "{engine}" unknown.')
    return properties


def calculate_vibrations(molobj, qchem_options, engine='gamess'):
    options_prime = copy.deepcopy(qchem_options)
    engine = options_prime.pop('engine', None)

    if engine == 'gamess':
        options_prime['cmd'] = options_prime['gamess.cmd']
        properties = gamess_calculations.calculate_vibrations(molobj, options_prime)
    elif engine == 'orca':
        options_prime['cmd'] = options_prime['orca.cmd']
        properties = orca_calculations.calculate_vibrations(molobj, options_prime)
    else:
        raise ValueError(f'Error: keyword argument engine = "{engine}" unknown.')
    return properties


def calculate_orbitals(molobj, qchem_options, engine='gamess'):
    options_prime = copy.deepcopy(qchem_options)
    engine = options_prime.pop('engine', None)

    if engine == 'gamess':
        options_prime['cmd'] = options_prime['gamess.cmd']
        properties = gamess_calculations.calculate_orbitals(molobj, options_prime)
    else:
        raise ValueError(f'Error: keyword argument engine = "{engine}" unknown.')

    return properties


def calculate_solvation(molobj, qchem_options, engine='gamess'):
    options_prime = copy.deepcopy(qchem_options)
    engine = options_prime.pop('engine', None)

    if engine == 'gamess':
        options_prime['cmd'] = options_prime['gamess.cmd']
        properties = gamess_calculations.calculate_solvation(molobj, options_prime)
    else:
        raise ValueError(f'Error: keyword argument engine = "{engine}" unknown.')

    return properties


def calculate_all_properties(molobj, qchem_options):
    engines = {
        'vibrations': 'orca',
        'orbitals': 'gamess',
        'solvation': 'gamess'
    }
    funcs = {
        'vibrations': calculate_vibrations,
        'orbitals': calculate_orbitals,
        'solvation': calculate_solvation,
    }

    def procfunc(conn, func, *args, **kwargs):
        properties = func(*args, **kwargs)
        conn.send(properties)
        conn.close()

    procs = []
    conns = []

    for calc in engines.keys():

        func = funcs[calc]
        engine = engines[calc]

        # Change scr
        filename = qchem_options.get("filename", "qchem_calc")
        qchem_options = copy.deepcopy(qchem_options)
        qchem_options["filename"] = filename + "_" + func.__name__
        qchem_options["engine"] = engine

        parent_conn, child_conn = Pipe()
        p = Process(
            target=procfunc,
            args=(child_conn, func, molobj, qchem_options),
        )
        p.start()

        procs.append(p)
        conns.append(parent_conn)

    for proc in procs:
        proc.join(timeout=MAX_TIME)

    properties_vib = conns[0].recv()
    properties_orb = conns[1].recv()
    properties_sol = conns[2].recv()

    return properties_vib, properties_orb, properties_sol
