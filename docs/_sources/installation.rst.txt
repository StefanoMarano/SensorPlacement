=========================
Download and Installation
=========================

This page describes the steps necessary to install software required for the sensor placement algorithm.

Download
########

The code can be cloned from the `repository <https://github.com/StefanoMarano/SensorPlacement>`__ ::

  git clone https://github.com/StefanoMarano/SensorPlacement.git

or downloaded as `ZIP archive <https://github.com/StefanoMarano/SensorPlacement/archive/master.zip>`__.


Install Python and packages
###########################

* Install Python 3 and the Python package installer by running::

    sudo apt-get install python3 python3-pip

* The following Python libraries are needed: ``glob, os, errno, sys, yaml, time, csv, logging, argparse, numpy, scipy, matplotlib, pulp``

  On Ubuntu many packages are installed with the default Python 3 installation. You may need to run install them manualy with::

    sudo pip install <package_name>

Install a solver
################

The computationally intense job is performed by an external, thirt-party, solver. Different solvers can be used including `Gurobi <http://www.gurobi.com/>`__, `CPLEX <https://www.ibm.com/analytics/cplex-optimizer>`__, `GNU Linear Programming Toolkit <https://www.gnu.org/software/glpk/>`__. Any solver supported by `PuLP <https://coin-or.github.io/pulp>`__ can be used.

To install the GNU Linear Programming Toolkit type::

  sudo apt-get install glpk-utils

For windows check `GLPK for Windows <http://winglpk.sourceforge.net/>`__.

For using a different solver, check the documentation of the specific solver.

Legacy MATLAB code with Gurobi
##############################

In the folder ``src_matlab`` MATLAB scripts for sensor placement are found.

 * Install MATLAB
 * Install and configure the solver `Gurobi <http://www.gurobi.com/>`__. You need to get a license to use Gurobi, the license is free for academic use.
 * The code was tested with Gurobi 6.5.1 and Matlab 2013b.

To run the code note the following:

 * The folder ``./utils/`` needs to be in the path searchable by MATLAB
 * From the MATLAB prompt, configure Gurobi with the script ``gurobi_setup.m``
 * Run the matlab script ``SensorPlacement.m``. Changing some variables in the code allows to modify design parameters including number of sensors and spatial bandwidth. See comments in the file itself.
 * Output is saved in the folder ``./output/``. Each output MAT file can be loaded into matlab. The variable `pos_sol` contains the array layout found by the algorithm. The array layout can be plotted with `plotArray_C(pos_sol)`
