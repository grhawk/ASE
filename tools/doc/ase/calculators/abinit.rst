.. module:: abinit

======
ABINIT
======

Introduction
============

ABINIT_ is a density-functional theory code based on pseudopotentials
and a planewave basis.


.. _ABINIT: http://www.abinit.org



Environment variables
=====================

.. highlight:: bash

The environment variable :envvar:`ASE_ABINIT_COMMAND` must be set to
something like this::

  abinis < PREFIX.files > PREFIX.log

where ``abinis`` is the executable.

A directory containing the pseudopotential files (at least of
:file:`.fhi` type) is also needed, and it is to be put in the
environment variable :envvar:`ABINIT_PP_PATH`. Abinit does not provide
tarballs of pseudopotentials so the easiest way is to download and
unpack
http://wiki.fysik.dtu.dk/abinit-files/abinit-pseudopotentials-2.tar.gz

Set the environment variables in your in your shell configuration file::

  export ASE_ABINIT_COMMAND="abinis < PREFIX.files > PREFIX.log"
  PP=${HOME}/abinit-pseudopotentials-2
  export ABINIT_PP_PATH=$PP/LDA_FHI
  export ABINIT_PP_PATH=$PP/GGA_FHI:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/LDA_HGH:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/LDA_PAW:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/LDA_TM:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/GGA_FHI:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/GGA_HGHK:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/GGA_PAW:$ABINIT_PP_PATH

.. highlight:: python


ABINIT Calculator
=================

Abinit does not specify a default value for the plane-wave cutoff
energy.  You need to set them as in the example at the bottom of the
page, otherwise calculations will fail.  Calculations wihout k-points
are not parallelized by default and will fail! To enable band
paralellization specify ``Number of BanDs in a BLOCK`` (``nbdblock``)
--- see
`<http://www.abinit.org/Infos_v5.2/tutorial/lesson_parallelism.html>`_.


Pseudopotentials
================

Pseudopotentials in the ABINIT format are available on the
`pseudopotentials`_ website.  A database of user contributed
pseudopotentials is also available there.

.. _pseudopotentials: http://www.abinit.org/Psps/?text=psps


Example 1
=========

Here is an example of how to calculate the total energy for bulk Silicon:

.. literalinclude:: Si_abinit.py
