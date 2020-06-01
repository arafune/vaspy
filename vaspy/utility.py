#! /usr/bin/env python3

"""Script to generate 3D view with mlab (Mayavi).

In generarl, VESTA is pretty nice. However, due to GUI, it is
difficult to make a series of figures of the model structure.
"""


import argparse
from logging import DEBUG, INFO, Formatter, StreamHandler, getLogger
import itertools
import numpy as np
from mayavi import mlab

from vaspy import tools, const
from vaspy.poscar import POSCAR
import vaspy

LOGLEVEL = INFO
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
    bgcolor=(1, 1, 1),
    line_thickness=0.05,
    line_color=(0, 0, 0),
    is_parallel=True,
    scale=20.0,
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
    # Read POSCAR
    poscar = vaspy_poscar
    poscar.repack_in_cell()
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
    fig = mlab.figure(1, bgcolor=bgcolor, size=figsize)
    mlab.clf()
    logger.debug("The type of 'fig' is {}".format(type(fig)))
    # Draw cell box
    draw_cell_box(unit_cell, line_thickness, line_color)
    #
    # Draw atoms
    #
    for atom in uniq_atom_symbols:
        atom_positions = np.array(poscar.positions)[site_indexes[atom]]
        mlab.points3d(
            atom_positions[:, 0],
            atom_positions[:, 1],
            atom_positions[:, 2],
            np.ones(len(site_indexes[atom])) * const.radii[atom],
            color=const.colors[atom],
            resolution=60,  # tunable ?
            scale_factor=1.0,  # tunable ?
            name="Atom_{}".format(atom),
        )
    # Drow bonds
    # step 1: find the bond to draw
    active_bonds = {}
    for i in itertools.product(uniq_atom_symbols, uniq_atom_symbols):
        if i in const.max_bond_length:
            active_bonds[i] = const.max_bond_length[i]
    # step 2: find the connections
    logger.debug(active_bonds.keys())
    for bonds, length in active_bonds.items():
        atom_a, atom_b = bonds
        connect_sites = []
        if atom_a == atom_b:
            for i in itertools.combinations(site_indexes[atom_a], 2):
                if (
                    np.linalg.norm(poscar.positions[i[0]] - poscar.positions[i[1]])
                    < length
                ):
                    connect_sites.append(i)
        else:
            for i in itertools.product(site_indexes[atom_a], site_indexes[atom_b]):
                if (
                    np.linalg.norm(poscar.positions[i[0]] - poscar.positions[i[1]])
                    < length
                ):
                    connect_sites.append(i)
        if not connect_sites:
            continue
        positions_a = np.array(poscar.positions)[[i[0] for i in connect_sites]]
        positions_b = np.array(poscar.positions)[[i[1] for i in connect_sites]]
        positions_middle = (positions_a + positions_b) / 2.0
        positions_bonding = np.zeros((len(connect_sites) * 2, 3))
        positions_bonding[1::2, :] = positions_middle
        bond_connectivity = np.vstack(
            [range(0, 2 * len(connect_sites), 2), range(1, 2 * len(connect_sites), 2)]
        ).T
        # plot the first half of the bond (from atom_a to the half.)
        positions_bonding[0::2, :] = positions_a
        bond_a = mlab.plot3d(
            positions_bonding[:, 0],
            positions_bonding[:, 1],
            positions_bonding[:, 2],
            tube_radius=0.1,  # should be tunable ?
            color=const.colors[atom_a],
            name="Bonds_{}-{}".format(atom_a, atom_b),
        )
        bond_a.mlab_source.dataset.lines = bond_connectivity
        # plot the last half of the bond (from atom_b to the half.)
        positions_bonding[0::2, :] = positions_b
        bond_b = mlab.plot3d(
            positions_bonding[:, 0],
            positions_bonding[:, 1],
            positions_bonding[:, 2],
            tube_radius=0.1,
            color=const.colors[atom_b],
            name="Bonds_{}-{}".format(atom_a, atom_b),
        )
        bond_b.mlab_source.dataset.lines = bond_connectivity
    logger.debug("ALL setting should be done.")
    if is_parallel:
        fig.scene.parallel_projection = True
        fig.scene.camera.parallel_scale = scale
        logger.debug(" fig.scene.camera.parallel_scale is {}".format(scale))
    else:
        fig.scene.parallel_projection = False
    mlab.orientation_axes()
    mlab.view(azimuth=phi, elevation=theta, distance=10)
    logger.debug("mlab.view ... done")
    mlab.savefig(output, magnification=5)
    logger.debug("mlab.save ... done. output file name is {}".format(output))


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
    return mlab.gcf()
