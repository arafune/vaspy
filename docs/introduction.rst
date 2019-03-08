VASPy: Introduction
=====================

What is this?
---------------
    This is VASP post-processing scripts for python.
    Translated from scRipt4VASP written by Ryuichi Arafune.

Requirement
------------
  *   Python3
  *   NumPy (any version)
  *   matplotlib (only when using outcar and locpot script)
  *   Scipy

Installation
-------------
  *   In this directory type following command:

    - All users (require admin)


    .. code-block:: bash

      $ python setup.py install


    - Personal

    .. code-block:: bash

      $ python setup.py install --user

Use scripts
-----------
    Add directory which contains scripts to your PATH.
    The default directory is:

    - Linux/Mac

      :All user: /usr/local/bin
      :Personal: ~/.local/bin

    - Windows

      :All user: [PythonDirectory]\Scripts

      :Personal(Vista and later): C:\Users\[YourName]\AppData\Roaming\Python\Scripts

      :Personal(XP): C:\Documents and Settings\[YourName]\Application Data\Python\[PythonXY]\Scripts

    [PythonDirectory] is the directory of python.exe; Default is C:\[PythonXY]

    [PythonXY] varies with python version, e.g. Python27 for ver. 2.7

Uninstall
-----------
    W/ pip:
    If you have installed pip package to your python,

    .. code-block:: bash

      $ pip uninstall vaspy

    W/o pip:
    Use rmvaspy.py and log file (install-vaspy.txt by default) in this directory.

    .. code-block:: bash

      $ python rmvaspy.py install-vaspy.txt
