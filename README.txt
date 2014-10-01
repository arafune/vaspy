VASPy

What is this?
    This is VASP post-processing scripts for python.
    Translated from scRipt4VASP written by Ryuichi Arafune.

Requirement
    Python (2.6 or newer)
    NumPy (any version)
    argparse (for Python 2.6, 3.0, 3.1)
    matplotlib (only when using outcar script)

Installation
    In this directory type following command:
    - All users (require admin)
    $ python setup.py install
    
    - Personal
    $ python setup.py install --user

Use scripts
    *** SCRIPTS OF THIS VERSION DOES NOT WORK CORRECTLY (yet). ***
    Add directory which contains scripts to your PATH.
    The default directory is:
    - Linux
    All user: /usr/local/bin
    Personal: ~/.local/bin
    
    - Windows
    All user: [PythonDirectory]\Scripts
    Personal(Vista and later): C:\Users\[YourName]\AppData\Roaming\Python\Scripts
    Personal(XP): C:\Documents and Settings\[YourName]\Application Data\Python\[PythonXY]\Scripts
    
    [PythonDirectory] is the directory of python.exe; Default is C:\[PythonXY]
    [PythonXY] varies with python version, e.g. Python27 for ver. 2.7

Uninstall
    With pip:
    If you have installed pip package to your python,
    $ pip uninstall vaspy
    
    Without pip:
    Use rmvaspy.py and log file (install-vaspy.txt by default) in this directory.
    $ python rmvaspy.py install-vaspy.txt
    
