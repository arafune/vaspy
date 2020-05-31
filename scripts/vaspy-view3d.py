#! /usr/bin/env python3

"""Script to generate 3D view with mlab (Mayavi).

In generarl, VESTA is pretty nice. However, due to GUI, it is
difficult to make a series of figures of the model structure.
"""


import argparse
from logging import DEBUG, INFO, Formatter, StreamHandler, getLogger

from vaspy import tools
from vaspy.poscar import POSCAR
import vaspy

LOGLEVEL = DEBUG
logger = getLogger(__name__)
fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
formatter = Formatter(fmt)
handler = StreamHandler()
handler.setLevel(LOGLEVEL)
logger.setLevel(LOGLEVEL)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False


def view3d(
    vaspy_obj,
    repeat=(1, 1, 1),
    output="output.tiff",
    figsize=(800, 800),
    bgcolor=(255, 255, 255),
    line_width=0.05,
    line_color=(0, 0, 0),
    is_parallel=True,
    theta=90,
    phi=0,
):
    """Drow 3D view by using mayavi.
    
    Parameters
    -------------
    vaspy_obj: vaspy.poscar.POSCAR or vaspy.CHGCAR.CHGCAR, vaspy.CHGCAR.LOCPOT
        VASPY object that include POSCAR object
    repeat: tuple of three int. default: (1, 1, 1)
        number of the unit cell repeat. For drawing the supercell structure
        ##FUTURE
    output: str
        output file name
    figsize: tuple of two float. default: (800, 800)
        image size
    bgcolor: tuple of three int.
        background color
    line_thickness: float
        thickness of the line of the cell box.
    line_color: (0, 0, 0)
        color of the line of the cell box.
    is_parallel: Boolean 
        if True, the view is parallel.
    theta: float
        the polar angle of the camera.
    phi: float
        the azimuth angle of the camera.
    """
