
MoleCalc - The Molecule Calculator
=================================

|screenshot|

Important! Molecalc, which is based on molcalc.org_, is *under development*!

The molecule calculator is a small web-based interface for doing small-scale
quantum chemistry calculation with the intent of giving chemical intuition to
students, from high-school to university.
Hosted at molecalc.cloud_.

.. _molecalc.cloud: https://_molecalc.cloud

.. _molcalc.org: http://molcalc.org

.. |screenshot| image:: https://raw.githubusercontent.com/jensengroup/molcalc/master/screenshot.jpg

Installation
------------

MoleCalc is a Python based web-service, so dependencies includes
python-packages, javascript-modules and a backend quantum chemistry program (for now it will be GAMESS).

To setup the Python environment please use Anaconda_, because we use RDKit in the background.

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


2. Create the Python environment using `conda` and `pip`.

.. code-block:: bash

    # make env chemhelp
    conda env create -f environment.yml -p env
    pip install -r requirements.txt
    git clone https://github.com/ppqm/ppqm ppqm.git --depth 1
    ln -s ppqm.git/ppqm ppqm

3. Download the JavaScript and frontend libraries, using the scripts.
   You need `unzip` and `wget` installed.
   All JavaScript libraries will be installed in the `molcalc/static` folder.

.. code-block:: bash

    # make setup_assets
    bash scripts/setup_chemdoodle.sh
    bash scripts/setup_jsmol.sh
    bash scripts/setup_fontawesome.sh
    bash scripts/setup_jquery.sh
    bash scripts/setup_rdkit.sh

4. Setup GAMESS_. You need to download and `compile GAMESS`__.


.. _GAMESS: https://www.msg.chem.iastate.edu/gamess/download.html
.. __: http://computerandchemistry.blogspot.com/2014/02/compiling-and-setting-up-gamess.html

5. Setup configuration by copying the example and edit.
   Especially note to edit the GAMESS section to reflect the setup of your setup.

.. code-block:: bash

    cp example.development.ini development.ini
    # edit development.ini


6. Test. Use the unittest to check that the configuration for GAMESS is setup correctly

.. code-block:: bash

    # make test
    python -m pytest tests


7. You are ready. Serve the server by

.. code-block:: bash

    # make serve
    env/bin/pserve development.ini --reload


And molecalc should now be available on ``localhost:6543``, based on the settings of development.ini.


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

