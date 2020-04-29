.. SensorPlacement documentation master file, created by
   sphinx-quickstart on Fri May 12 07:10:21 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Sensor placement for the analysis of seismic surface waves
==========================================================

Sensor placement is in general a really hard problem. It is difficult to obtain optimal array geometries. However it is possible to find very good array geometries. Here is a Python implementation of the algorithm presented in

* Stefano Maranò, Donat Fäh, and Yue M. Lu, **"Sensor Placement for the Analysis of Seismic Surface Waves: Sources of Error, Design Criterion and Array Design Algorithms"**,  Geophys. J. Int. (2014) 197 (3): 1566--1581. `Free access online <http://academic.oup.com//gji/article/197/3/1566/652229/Sensor-placement-for-the-analysis-of-seismic?guestAccessKey=fbe9d23c-d168-4b3a-8946-118c73b482ad>`__, `doi:10.1093/gji/ggt489 <https://doi.org/10.1093/gji/ggt489>`__.

In the proposed work, the sensor placement problem is addressed as a `mixed integer program  <https://en.wikipedia.org/wiki/Mixed_integer_program>`__. The proposed array design technique was developed for usage in seismology. Two-dimensional arrays for other applications, including acoustic and radar, may also be designed.

The algorithm was developed at the `Swiss Seismological Service <http://www.seismo.ethz.ch/>`__ of `ETH Zurich <https://www.ethz.ch>`__.

.. note::  The Python code is working but I am aware that some features are not polished. Please feel to contact me if there is anything you wish to discuss or improve :math:`\textrm{wavedec@gmail.com}`.


Contents
========

.. toctree::
   :maxdepth: 3

   scientific.rst
   installation.rst
   userguide.rst
   examples.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

