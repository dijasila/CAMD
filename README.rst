=============
CAMd web-apps
=============

Web-apps for:

* C2DB
* old CMR-projects
* QPOD
* BiDB
* CRYSP
* OQMD12345 (CrystalBank)

.. contents::


Important links
===============

* `Bottle <https://bottlepy.org/docs/dev/index.html>`__
* `uWSGI <https://uwsgi-docs.readthedocs.io/en/latest/index.html>`__
* `Plotly <https://plotly.com/python/>`__
* `Bootstrap
  <https://getbootstrap.com/docs/5.3/getting-started/introduction/>`__
* `WSGI <https://peps.python.org/pep-3333/>`_
* `ASE <https://wiki.fysik.dtu.dk/ase/index.html>`__
* `CMR-projects <https://cmrdb.fysik.dtu.dk/>`__
* `CMR-repo <https://gitlab.com/camd/cmr>`__
* `ASR-repo <https://gitlab.com/asr-dev/asr>`__
* `Scaling to many users <https://workchronicles.com/white-lies/>`__


Test URLs
=========

Forwarding to ``https://fysik-cmr02.fysik.dtu.dk:<port>``:

======================================================  ====
Link                                                    port
======================================================  ====
`C2DB-test <https://c2db-test.fysik.dtu.dk/>`__         8081
`CMR-test <https://cmrdb-test.fysik.dtu.dk/>`__         8082
`QPOD <https://qpod.fysik.dtu.dk/>`__                   8083
`BiDB <https://bidb.fysik.dtu.dk/>`__                   8084
`CRYSP <https://crysp.fysik.dtu.dk/>`__                 8085
`OQMD12345 <https://oqmd12345.fysik.dtu.dk/>`__         8086
`OPTIMADE <https://c2db-test.fysik.dtu.dk/optimade>`__
======================================================  ====


Installation
============

CAMd-web needs Python_ version 3.9 or later.

.. code:: bash

    $ python -m venv venv
    $ source venv/bin/activate
    $ git clone git@gitlab.com:ase/ase
    $ pip install -e ase
    $ git clone git@gitlab.com:asr-dev/asr
    $ pip install -e asr
    $ git clone git@gitlab.com:camd/camd-web
    $ pip install -e camd-web[test,optimade]


.. _Python: https://python.org/


Command-line interface
======================

usage: camd-web [-h] filename [filename ...]

positional arguments:
  filename    Filename of atomic structure file.

options:
  -h, --help  show this help message and exit


CMR-app for old projects
========================

::

    $ python -m camdweb.cmr.app *.db


C2DB-app
========

::

    $ python -m camdweb.c2db.copy ~cmr/C2DB-ASR "tree/*/*/*/" ...
    $ python -m camdweb.c2db.app A*/

Folder structure for UIDs ``1MoS2-1`` and ``1MoS2-2``::

  AB2/1MoS2/1/data.json
              structure.xyz
              results-asr.<property1>.json
              results-asr.<property2>.json
              ...
              ...
            2/data.json
              structure.xyz
              ...
  oqmd123.json.gz
  convex-hulls/MoS.json
               Mo.json
               S.json


Testing the C2DB-app
--------------------

For development work, just copy one or a few meterial folders from Niflheim
to your local machine::

    $ mkdir C2DB-test
    $ ssh sylg
    $ cd /home/niflheim2/cmr/C2DB-ASR/tree/AB2/MoS2
    $ scp -r MoS2-b3b4685fb6e1 <your-machine>:C2DB-test/
    $ ^D

Then you can play with those files like this::

    $ cd C2DB-test
    $ python -m camdweb.c2db.copy . "MoS2*/"
    $ python -m camdweb.c2db.app AB2


Development
===========

Please run the following checks on your code::

    $ cd <root-of-repo>
    $ mypy
    $ flake8 camdweb
    $ camd-web-coverage

If 100 % coverage is not possible then you can make CI pass by adding
``# pragma: no cover`` or ``# pragma: no branch`` comments.


Deployment
==========

On the ``fysik-cmr02`` server run uWSGI like this::

    $ uwsgi -w "camdweb.c2db.app:create_app()" --http :8081 --master --threads=2 --enable-threads --daemonize=c2db.log
    $ uwsgi -w "camdweb.cmr.app:create_app()" --http :8082 --master --threads=2 --enable-threads --daemonize=cmr.log
    $ uwsgi -w "camdweb.bidb.app:create_app()" --http :8084 --master --threads=2 --enable-threads --daemonize=bidb.log
    $ uwsgi -w "camdweb.oqmd12345.app:create_app()" --http :8086 --master --threads=2 --enable-threads --daemonize=oqmd12345.log


How it works
============

In the picture below, ``camd.app`` is the WSGI_ app::

  camd
    |
    v
  +---------+
  | CAMDApp |   app   +------------+
  |         |-------->| bottle.App |
  |   and   |         +------------+
  |   sub-  |
  | classes |   sessions   +----------+
  |         |------------->| Sessions |
  |         |              +----------+
  +---------+
     |
     |materials
     |
     v
  +------------+  index   +-------+
  | Materials  |--------->| Index |
  |            |          +-------+
  |  --------  |
  | |Material| |
  |  --------  |   panels   +--------------+
  | |Material| |----------->| list[Panel]  |
  |  --------  |            |              |
  |     :      |            |  ----------  |
  |     :      |            | |AtomsPanel| |
  +------------+            |  ----------  |
                            | |OtherPanel| |
                            |  ----------  |
                            |      :       |
                            |      :       |
                            +--------------+


Objects
-------

:bottle.App:
    WSGI_ app.  Defines the end-points ``/`` and ``/material/<uid>/``.

:Material:
    Attributes: ``uid: str``, ``atoms: Atoms``, ``folder: Path``
    ``columns: dict[str, bool | int | float | str]`` and
    ``count: dict[str, int]``.
    The ``columns`` dictionary stores key-value pairs for displaying
    in the landing-page table.  The ``count`` dictionary stores the
    number of each species present.

:Panel:
    Has a ``get_data(material)`` method that can produce a ``PanelData``
    object (html, title, pop-up text, javascript and possibly sub-panels)
    to be assembled in the ``/material/<uid>/`` end-point.

:Index:
    Handles efficient filtering of materials using the values in
    ``Material.columns`` and ``Material.count``.

:Materials:
    Keeps track of all the ``Panel`` and ``Material`` objects
    (``materials[uid]`` will give you the material with the given ``uid``
    and ``for material in materials:`` will loop over them all).
    Also handles two dictionaries that it shares with the panels:

    * ``column_descriptions: dict[str, str]`` for longer descriptions of the
      short column names.
    * ``html_column_formatters: dict[str, Callable[..., str]]`` for converting
      bool, int, float and str values to HTML strings.

:Sessions:
    Handles ``Session`` objects for clients (one for each browser-tab).

:Session:
    Remembers selected columns, sorting information, ...
