
MoleCalc - The Molecule Calculator
=================================

MoleCalc is now publically available at molecalc.cloud_!

MoleCalc is based on the MolCalc project (`https://github.com/jensengroup/molcalc`_), a web-based
chemistry teaching tool available at molcalc.org_ that was originally built by the Jensen Group at
the University of Copenhagen. MoleCalc is a separate project currently being developed by Dr. Sean L.
Seyler and led by Prof. Jeffery Yarger in the
Yarger Research Group at the School of Molecular Sciences at Arizona State University. MoleCalc allows
chemists to build small molecules and estimate key molecular properties using research-grade quantum
chemistry software. Currently, MoleCalc can provide, in a matter of minutes or seconds, a decent sense
of the physical and chemical properties of a chosen molecule, such as

* energy-minimized molecular structure
* electronic structure and molecular orbitals
* thermochemical properties of the molecular gas form
* vibrational modes and their frequencies


|screenshot|

**Important**: MoleCalc,  is heavily based on the original molcalc.org_ project
, is *under active development*! This
ReadMe is also not guaranteed to reflect the most recent state of the project,
so please be patient!

MoleCalc is a small web-based interface for doing small-scale
quantum chemistry calculation with the intent of giving chemical and physical intuition to
students, from high-school to university.

.. _molecalc.cloud: https://molecalc.cloud

.. _molcalc.org: http://molcalc.org

.. _`https://github.com/jensengroup/molcalc`: https://github.com/jensengroup/molcalc

.. |screenshot| image:: https://raw.githubusercontent.com/mscloudlab/molecalc/chm343-beta/molecalc_v1.jpg


Installation
------------

MoleCalc is a Python based web-service, so dependencies include
python packages, javascript modules and at least one backend quantum chemistry program
(currently GAMESS_, with plans to incorporate new capabilities using Orca_).

To setup the Python environment please use Anaconda_, because we use RDKit in the background.

.. _GAMESS: https://www.msg.chem.iastate.edu/gamess/

.. _Orca: https://www.faccts.de/orca/

.. _Anaconda: https://www.anaconda.com/download

.. code-block:: bash

   # Install anaconda (only needed if you don't already have a Python enviroment with conda)
   wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
   bash miniconda3.sh -b -p /opt/miniconda3

with the Python environment we can setup MoleCalc. Note that most of the steps are inserted into the `Makefile`.

1. Clone down the repository

.. code-block:: bash

    git clone https://github.com/mscloudlab/molecalc --depth 1
    cd molecalc


2. Create the Python environment using `conda` and `pip` and the provided
   `environment.yml` and `requirements.txt` files.

.. code-block:: bash

    # make env chemhelp
    conda env create -f environment.yml -p env
    pip install -r requirements.txt

3. Install `molecalc` in `./env` and `ppqm` locally in the molecalc directory:

.. code-block:: bash

    pip install -e .
    git clone https://github.com/ppqm/ppqm ppqm.git --depth 1
    ln -s ppqm.git/ppqm ppqm

4. Download the JavaScript and frontend libraries, using the scripts.
   You need `unzip` and `wget` installed.
   All JavaScript libraries will be installed in the `molcalc/static` folder.

.. code-block:: bash

    # make setup_assets
    bash scripts/setup_chemdoodle.sh
    bash scripts/setup_jsmol.sh
    bash scripts/setup_fontawesome.sh
    bash scripts/setup_jquery.sh
    bash scripts/setup_rdkit.sh

5. Set up GAMESS. You need to download_ and `compile GAMESS`__.


.. _download: https://www.msg.chem.iastate.edu/gamess/download.html
.. __: http://computerandchemistry.blogspot.com/2014/02/compiling-and-setting-up-gamess.html

6. Setup configuration by copying the example and edit.
   Especially note to edit the GAMESS section to reflect the setup of your setup.

.. code-block:: bash

    cp example.development.ini development.ini
    # edit development.ini


7. Test. Use the unittest to check that the configuration for GAMESS is setup correctly

.. code-block:: bash

    # make test
    python -m pytest tests


8. You are ready. Serve the server by

.. code-block:: bash

    # make serve
    env/bin/pserve development.ini --reload

MoleCalc should now be available on ``localhost:6543``, based on the settings of development.ini.


Dependencies
------------

rdkit,
pyramid,
fontawesome,
jquery,
chemdoodle,
jsmol,
gamess


TODO
----

Remove connections from javascript libs

.. code-block::

    Failed to load resource: net::ERR_INTERNET_DISCONNECTED
    ichemlabs.cloud.chemdoodle.com/icl_cdc_v070001/WebHQ


TODO computation
----------------

Extend the computations for molcalc to include

* spectrum
** H/C-NMR
** mass spectrum
** vibrational

* open shell systems


TODO Better texts
-----------------

Tutorials and assignment examples (with answers)

Better FAQ interface


Known Problems
--------

If rdkit has problems finding `libxrender.so` then you need to install

.. code-block:: bash

    apt install -y libxrender-dev

or

.. code-block:: bash

    ./env/bin/conda install nox
    ./env/bin/conda install cairo

