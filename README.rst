CxDB web-apps
=============

Web-apps for:

* old CMR-projects
* C2DB
* BiDB
* CRYSP
* QPOD


Important links
---------------

* `Bottle <https://bottlepy.org/docs/dev/index.html>`__
* `Plotly <https://plotly.com/python/>`__
* `ASE <https://wiki.fysik.dtu.dk/ase/index.html>`__
* `CMR-projects <https://cmrdb.fysik.dtu.dk/>`__
* `CMR-repo <https://gitlab.com/camd/cmr>`__
* `ASR-repo <https://gitlab.com/asr-dev/asr>`__
* `CXDB-test <http://fysik-cmr02.fysik.dtu.dk:8081/>`__

Installation
------------

CxDB-web needs Python_ version 3.9 or later.

::

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


CMR-app for old projects
------------------------

::

    $ python -m cxdb.cmr.app *.db


C2DB-app
--------

::

    $ python -m cxdb.c2db.copy_files ~cmr/C2DB-ASR "tree/*/*/*/" ...
    $ python -m cxdb.c2db.app

