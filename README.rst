CxDB web-apps
=============

Web-apps for:

* old CMR-projects
* C2DB
* BiDB
* CRYSP
* QPOD
* OQMD12345


Important links
---------------

* `Bottle <https://bottlepy.org/docs/dev/index.html>`__
* `Plotly <https://plotly.com/python/>`__
* `Bootstrap
  <https://getbootstrap.com/docs/5.3/getting-started/introduction/>`__
* `ASE <https://wiki.fysik.dtu.dk/ase/index.html>`__
* `CMR-projects <https://cmrdb.fysik.dtu.dk/>`__
* `CMR-repo <https://gitlab.com/camd/cmr>`__
* `ASR-repo <https://gitlab.com/asr-dev/asr>`__


Test URLs
---------

* `CXDB-test <https://fysik-cmr02.fysik.dtu.dk:8081/>`__
* `C2DB-test <https://c2db-test.fysik.dtu.dk/>`__
* `CMR-test <https://cmrdb-test.fysik.dtu.dk/>`__


Installation
------------

CxDB-web needs Python_ version 3.9 or later.

.. code:: bash

    $ python -m venv venv
    $ source venv/bin/activate
    $ # ASE master:
    $ pip install --no-deps git+https://gitlab.com/ase/ase.git
    $ # A minimal ASR installation:
    $ pip install --no-deps git+https://gitlab.com/asr-dev/asr.git
    $ pip install click
    $ And CXDB:
    $ git clone git@gitlab.com:jensj/cxdb-web
    $ pip install -e cxdb-web


.. _Python: https://python.org/

Command-line interface
----------------------

usage: cxdb [-h] filename [filename ...]

positional arguments:
  filename    Filename of atomic structure file.

options:
  -h, --help  show this help message and exit


CMR-app for old projects
------------------------

::

    $ python -m cxdb.cmr.app *.db


C2DB-app
--------

::

    $ python -m cxdb.c2db.copy_files ~cmr/C2DB-ASR "tree/*/*/*/" ...
    $ python -m cxdb.c2db.app

