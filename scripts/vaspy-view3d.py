#! /usr/bin/env python3

"""Script to generate 3D view with mlab (Mayavi).

In generarl, VESTA is pretty nice. However, due to GUI, it is
difficult to make a series of figures of the model structure.
"""


import argparse
from logging import DEBUG, INFO, Formatter, StreamHandler, getLogger
import numpy as np
from mayavi import mlab

from vaspy import tools, const
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
    vaspy_poscar,
    repeat=(1, 1, 1),
    output="output.tiff",
    figsize=(800, 800),
    bgcolor=(255, 255, 255),
    line_thickness=0.05,
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
    line_color: tuple
        color of the line of the cell box. default: (0, 0, 0)
    is_parallel: Boolean
        if True, the view is parallel.
    theta: float
        the polar angle of the camera.
    phi: float
        the azimuth angle of the camera.
    """

    poscar = vaspy_poscar
    poscar.to_cartesian()
    poscar.tune_scaling_factor(1.0)
    unit_cell = poscar.cell_vecs
    ## poscar.atomnums = [144, 288, 19, 18, 2, 3]
    ## poscar.atomtypes = ['Mo', 'S', 'C', 'H', 'N', 'O']
    atom_symbols = tools.atomtypes_atomnums_to_atoms(poscar.atomtypes, poscar.atomnums)
    uniq_atom_symbols = list(set(atom_symbols))
    site_indexes = {}
    for atom in uniq_atom_symbols:
        site_indexes[atom] = [i for i, x in enumerate(atom_symbols) if x == atom]
    cell_box = draw_cell_box(unit_cell, line_thickness, line_color)
    #
    for atom in uniq_atom_symbols:
        atom_positions = np.array(poscar.positions)[site_indexes[atom]]
        mlab.points3d(
            atom_positions[:, 0],
            atom_positions[:, 1],
            atom_positions[:, 2],
            np.ones(len(site_indexes[atom])) * const.radii[atom],
            color=const.colors[atom],
            resolution=60,  #  tunable ?
            scale_factor=1.0,  # tunable ?
            name="Atom_{}".format(atom),
        )


def draw_cell_box(unit_cell, line_thickness, line_color):
    """Return mlab.plot3d object of cell box
    
    unit_cell: three vector
        Represent the unit cell
    line_thickness: float
        thickness of the box line
    line_color: tuple of int
        color of the box line
    """
    corner_indexes = np.array(np.meshgrid(range(2), range(2), range(2), indexing="ij"))
    corner_coordinates = np.array(
        np.tensordot(unit_cell, corner_indexes, axes=(0, 0))
    ).reshape((3, -1))
    corner_indexes = corner_indexes.reshape((3, -1))
    connections = []
    for i in range(corner_indexes.shape[1]):
        for j in range(i):
            L = corner_indexes[:, i] - corner_indexes[:, j]
            if list(L).count(0) == 2:
                connections.append((i, j))
    cell_box = mlab.plot3d(
        corner_coordinates[0],
        corner_coordinates[1],
        corner_coordinates[2],
        tube_radius=line_thickness,
        color=line_color,
        name="CellBox",
    )
    cell_box.mlab_source.dataset.lines = np.array(connections)
    return cell_box
