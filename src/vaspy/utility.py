#! /usr/bin/env python3

"""Script to generate 3D view with mlab (Mayavi).

In general, VESTA is pretty nice. However, due to GUI, it is
difficult to make a series of figures of the model structure.
"""
from __future__ import annotations

import argparse
import itertools
from logging import DEBUG, INFO, Formatter, StreamHandler, getLogger
from typing import TYPE_CHECKING

import numpy as np
from mayavi import mlab

import vaspy
from vaspy import const, tools

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from vaspy.poscar import POSCAR

LOGLEVELS = (DEBUG, INFO)
LOGLEVEL = LOGLEVELS[1]
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
    vaspy_poscar: POSCAR,
    repeat: tuple[int, int, int] = (1, 1, 1),
    output: str | None = None,
    figsize: tuple[float, float] = (800.0, 800.0),
    fgcolor: tuple[int, int, int] = (0, 0, 0),
    bgcolor: tuple[int, int, int] = (1, 1, 1),
    line_thickness: float = 0.05,
    line_color: tuple[int, int, int] = (0, 0, 0),
    is_parallel: bool = True,
    scale: float = 20.0,
    theta: float = 90,
    phi: float = 0,
):
    """Drow 3D view by using mayavi.

    Parameters
    ----------
    vaspy_obj: vaspy.poscar.POSCAR, vaspy.chgcar.CHGCAR or vaspy.locpot.LOCPOT
        VASPY object that includes POSCAR object
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
    atom_symbols = tools.atom_types_atomnums_to_atoms(
        poscar.atom_types,
        poscar.atomnums,
    )
    uniq_atom_symbols = list(set(atom_symbols))
    site_indexes = {}
    for atom in uniq_atom_symbols:
        site_indexes[atom] = [i for i, x in enumerate(atom_symbols) if x == atom]
    fig = mlab.figure(1, bgcolor=bgcolor, fgcolor=fgcolor, size=figsize)
    mlab.clf()
    logger.debug(f"The type of 'fig' is {type(fig)}")
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
            name=f"Atom_{atom}",
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
            [range(0, 2 * len(connect_sites), 2), range(1, 2 * len(connect_sites), 2)],
        ).T
        # plot the first half of the bond (from atom_a to the half.)
        positions_bonding[0::2, :] = positions_a
        bond_a = mlab.plot3d(
            positions_bonding[:, 0],
            positions_bonding[:, 1],
            positions_bonding[:, 2],
            tube_radius=0.1,  # should be tunable ?
            color=const.colors[atom_a],
            name=f"Bonds_{atom_a}-{atom_b}",
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
            name=f"Bonds_{atom_a}-{atom_b}",
        )
        bond_b.mlab_source.dataset.lines = bond_connectivity
    logger.debug("ALL setting should be done.")
    if is_parallel:
        fig.scene.parallel_projection = True
        fig.scene.camera.parallel_scale = scale
        logger.debug(f" fig.scene.camera.parallel_scale is {scale}")
    else:
        fig.scene.parallel_projection = False
    """
    orientation_axes = mlab.orientation_axes()
    ## for coloring the text of the axis, (X, Y, Z)
    ## https://github.com/enthought/mayavi/issues/750
    orientation_axes.axes.x_axis_caption_actor2d.caption_text_property.color = fgcolor
    orientation_axes.axes.y_axis_caption_actor2d.caption_text_property.color = fgcolor
    orientation_axes.axes.z_axis_caption_actor2d.caption_text_property.color = fgcolor

    if output:
        mlab.view(azimuth=phi, elevation=theta, distance=10)
        mlab.savefig(output, magnification=5)
    logger.debug("mlab.save ... done. output file name is {}".format(output))
    """
    return mlab.gcf()


def view_atom_with_surface(
    vaspy_chgcar: vaspy.chgcar.CHGCAR,
    repeat: tuple[int, int, int] = (1, 1, 1),
    output: str | None = None,
    figsize: tuple[float, float] = (800, 800),
    fgcolor: tuple[int, int, int] = (0, 0, 0),
    bgcolor: tuple[int, int, int] = (1, 1, 1),
    line_thickness: float = 0.05,
    line_color: tuple[int, int, int] = (0, 0, 0),
    is_parallel: bool = True,
    volume_iso_color: tuple[int, int, int] = (1, 1, 1),
    scale: float = 20.0,
    theta: float = 90.0,
    phi: float = 0,
):
    """Draw contour surface from CHGCAR, LOCPOT, etc.

    Important:
    ---------
    The unit cell must be rectangle.  This limitation comes from mayavi.

    """
    poscar = vaspy_chgcar.poscar
    poscar.repack_in_cell()
    poscar.to_cartesian()
    poscar.tune_scaling_factor(1.0)
    unit_cell = poscar.cell_vecs
    grid_size = vaspy_chgcar.grid.shape
    lab_grid_size = grid_nums(grid_size, poscar.cell_vecs)
    if grid_size == lab_grid_size:
        volume_data = vaspy_chgcar.grid.frame(0).data.reshape(grid_size, order="F")
    else:
        volume_data = reallocate_to_labframe(
            grid_size,
            unit_cell,
            vaspy_chgcar.grid.frame(0).data,
        )
    #
    view3d(
        poscar,
        output=None,
        figsize=figsize,
        fgcolor=fgcolor,
        bgcolor=bgcolor,
        line_thickness=line_thickness,
        line_color=line_color,
        is_parallel=is_parallel,
        scale=scale,
        theta=theta,
        phi=phi,
    )

    cuboid = tools.cuboid(unit_cell)
    positions = np.mgrid[
        cuboid[0][0] : cuboid[0][1] : lab_grid_size[0] * 1j,
        cuboid[1][0] : cuboid[1][1] : lab_grid_size[1] * 1j,
        cuboid[2][0] : cuboid[2][1] : lab_grid_size[2] * 1j,
    ]
    #
    mlab.contour3d(
        positions[0],
        positions[1],
        positions[2],
        volume_data,
        color=volume_iso_color,
        transparent=True,
        name="VolumeData",
    )
    """"    mlab.view(azimuth=phi, elevation=theta, distance=10)
    if output:
        mlab.savefig(output, magnification=5)"""
    return mlab.gcf()


def draw_cell_box(unit_cell, line_thickness, line_color) -> mlab.plot3d:
    """Return mlab.plot3d object of cell box.

    unit_cell: three vector
        Represent the unit cell
    line_thickness: float
        thickness of the box line
    line_color: tuple of int
        color of the box line
    """
    corner_indexes = np.array(np.meshgrid(range(2), range(2), range(2), indexing="ij"))
    corner_coordinates = np.array(
        np.tensordot(unit_cell, corner_indexes, axes=(0, 0)),
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


def _grid_nums(n_grids, crystal_axes, lab_axes):
    """Return tuple of the Grid determined from the original grids.

    Parameters
    ----------
    original_grids:
    tetrahedron_vectors:
    cuboid_vectors

    """
    new_grids = []
    for lab_axis in lab_axes:
        cosins = [
            np.dot(crystal_axis, lab_axis)
            / np.linalg.norm(crystal_axis)
            / np.linalg.norm(lab_axis)
            for crystal_axis in crystal_axes
        ]
        max_index = np.argmax(cosins)
        new_grids.append(
            int(
                round(
                    n_grids[max_index]
                    / (
                        np.linalg.norm(crystal_axes[max_index])
                        * cosins[max_index]
                        / np.linalg.norm(lab_axis)
                    ),
                ),
            ),
        )
    return tuple(new_grids)


def grid_nums(n_grids, crystal_axes) -> tuple[int, ...]:
    """Return the number of gris of the Laboratory frame.

    Parameters
    ----------
    n_grids: tuple
        number of grid points. (Written in CHGCAR)
    crystal_axes: array_like
        Three vectors for the laboratory frame.

    Return: tuple
        number of grid points of the laboratory frame.

    """
    lab_box = tools.cuboid(crystal_axes)
    lab_axes = np.array(
        (
            (next(p[1] - p[0] for p in lab_box), 0, 0),
            (0, [p[1] - p[0] for p in lab_box][1], 0),
            (0, 0, [p[1] - p[0] for p in lab_box][2]),
        ),
    )
    return tuple(_grid_nums(n_grids, crystal_axes, lab_axes))


def _find_diagonal_indexes(n_grids, crystal_axes):
    """Return the tuple of the index corresponds to origin and its diagonal point.

    Parameters
    ----------
    crystal_axes: array-like
        Three vectors that represents a-, b- and c-axis.
    n_grids

    Returns
    -------

    """
    cuboid = tools.cuboid(crystal_axes)
    n_labgrids = grid_nums(n_grids, crystal_axes)
    o = (0, 0, 0)
    diag = crystal_axes[0] + crystal_axes[1] + crystal_axes[2]
    logger.debug(f"diag: {diag}")
    index_o = []
    index_diag = []
    for i in range(3):
        tmp = np.abs(np.mgrid[cuboid[i][0] : cuboid[i][1] : n_labgrids[i] * 1j] - o[i])
        index_o.append(np.where(tmp == tmp.min())[0][0])
        tmp = np.abs(
            np.mgrid[cuboid[i][0] : cuboid[i][1] : n_labgrids[i] * 1j] - diag[i],
        )
        index_diag.append(np.where(tmp == tmp.min())[0][0])
    return tuple(index_o), tuple(index_diag)


def reallocate_to_labframe(
    mesh_in_direct_coor: tuple[int, ...],
    crystal_axes: NDArray[np.float64],
    volume_data: NDArray[np.float64],
    no_roll: bool = False,
):
    """Return the volume mesh data in lab frame (Cartesian coordinate).

    Parameters
    ----------
    mesh_in_direct_coor: tuple of int
        Number of mesh size
    crystal_axes: array-like
        Vector of crystal axes.
    volume_data: array-like
        Volume data

    Returns
    -------
    np.array

    """
    if mesh_in_direct_coor != volume_data.shape:
        volume_data = volume_data.reshape(mesh_in_direct_coor, order="F")

    lab_grid = grid_nums(mesh_in_direct_coor, crystal_axes)
    cuboid = tools.cuboid(crystal_axes)
    len_cuboid = np.array(
        (
            cuboid[0][1] - cuboid[0][0],
            cuboid[1][1] - cuboid[1][0],
            cuboid[2][1] - cuboid[2][0],
        ),
    )
    logger.debug(f"cuboid: {cuboid}, len_cuboid: {len_cuboid}")
    len_cyrstal_axes = np.linalg.norm(crystal_axes, axis=1)
    logger.debug(f"len_crystal_axes: {len_cyrstal_axes}")
    lab_frame = np.empty(lab_grid, dtype=float)
    lab_frame[:, :, :] = np.nan
    logger.debug(f"Shape of lab_frame: {lab_frame.shape}")
    nx = np.linspace(0, 1, mesh_in_direct_coor[0], endpoint=False)
    ny = np.linspace(0, 1, mesh_in_direct_coor[1], endpoint=False)
    nz = np.linspace(0, 1, mesh_in_direct_coor[2], endpoint=False)
    index_o = []
    for i in range(3):
        tmp = np.abs(np.linspace(cuboid[i][0], cuboid[i][1], lab_grid[i]))
        idx = np.where(tmp == tmp.min())[0][0]
        index_o.append(idx)
    logger.debug(f"index_o: {index_o}")
    logger.debug(f"mesh_in_direct_coor is {mesh_in_direct_coor}")
    for i_x, i_y, i_z in itertools.product(
        range(mesh_in_direct_coor[0]),
        range(mesh_in_direct_coor[1]),
        range(mesh_in_direct_coor[2]),
    ):
        lab_coordinate = crystal_axes.transpose().dot(
            np.array((nx[i_x], ny[i_y], nz[i_z])),
        )
        relative_lab_coordinate = lab_coordinate / len_cuboid
        lab_index = (np.array(lab_grid) * relative_lab_coordinate).round().astype(int)
        logger.debug(f"lab_index at {i_x}, {i_y}, {i_z} is {lab_index}")
        lab_frame[lab_index[0], lab_index[1], lab_index[2]] = volume_data[i_x, i_y, i_z]
    if no_roll:
        return lab_frame
    return np.roll(lab_frame, index_o, (0, 1, 2))


if __name__ == "__main__":
    arg = argparse.ArgumentParser()
    arg.add_argument(
        "-o",
        action="store",
        dest="outImg",
        type=str,
        default="output.tiff",
        help="Output image name",
    )
    arg.add_argument(
        "-size",
        action="store",
        dest="outImgSize",
        nargs=2,
        type=int,
        default=(400, 400),
        help="Output image size",
    )
