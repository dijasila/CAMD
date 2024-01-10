CxDB web-apps
=============

Web-apps for:

* old CMR-projects
* C2DB
* BiDB
* CRYSP
* QPOD


Installation
------------

CxDB-web needs Python_ version 3.9 or later.

::

    $ python -m venv venv
    $ source venv/bin/activate
    $ pip install --no-deps git+https://gitlab.com/ase/ase.git
    $ pip install --no-deps git+https://gitlab.com/asr-dev/asr.git
    $ pip install click
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

