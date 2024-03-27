"""Module for WAVECAR class."""

from __future__ import annotations

from pathlib import Path
from typing import IO, TYPE_CHECKING

import numpy as np
from scipy.fftpack import ifftn

from vaspy import mesh3d
from vaspy.poscar import POSCAR

if TYPE_CHECKING:
    from numpy.typing import DTypeLike, NDArray

    from vaspy.mesh3d import VASPGrid

Ry_in_eV: float = 13.605826
Au_in_AA: float = 0.529177249

# If parallel version vasp, set True
# See OUTCAR file
# If you find:
#   use parallel FFT for wavefunctions z direction half grid
# ('z' direction is important),
# set True.
#   use serial FFT for wavefunctions x direction half grid
# set False

PARALLEL: bool = False


class WAVECAR:
    """Class for storing the data of WAVECAR file.

    Parameters
    ----------
    WAEVCAR_file: str
        File name of 'WAVECAR'


    Attributes
    ----------
    recl: numpy.int
        Record length
    n_spin: numpy.int
        Number of spins
    rtag: int
        Tag for precision in WAVECAR
    nplwvs: int
        Number of plane waves.
    num_k: int
        Number of k-points
    n_bands: int
        Number of bands
    encut: float
        Cut off energy in eV unit.
    realcell: NDArray
        Vectors for the unit cell in real space
    fermi: float
        fermi level
    rcpcell: numpy.array
        Vectors for the unit cell in reciprocal space
    k_vectors: NDArray
        kvector
    bands: NDArray
        Energy
    occs: NDArray
        Occupation
    gamma: boolean
        True if Wavecar by vasp with -DwNGZHalf

    """

    def __init__(self, filename: str | Path = "WAVECAR") -> None:
        """Initialize WAVECAR class.

        Parameters
        ----------
        filename: str | Path (Default: "WAVECAR")
            File name of the 'WAVECAR'

        """
        self.wfc: IO[bytes] = Path(filename).open("rb")
        self.gamma: bool = False
        #
        self.header()
        self.band()
        if self.num_k == 1:
            self.check_DwNGZHalf()

    def header(self) -> None:
        """Read the information of the system.

        Information of the system is stored in the first two record
        in WAVECAR file

        rec1: recl, n_spin, rtag
        rec2: num_k, n_bands ,encut, ((cell(i, j) i=1, 3), j=1, 3), efermi

        """
        self.wfc.seek(0)
        self.recl: int
        self.n_spin: int
        self.rtag: int
        self.recl, self.n_spin, self.rtag = np.array(
            np.fromfile(self.wfc, dtype=float, count=3),
            dtype=int,
        )
        self.wfc.seek(self.recl)
        #
        dump: NDArray[np.float64] = np.fromfile(self.wfc, dtype=float, count=13)
        #
        self.num_k: int = int(dump[0])
        self.n_bands: int = int(dump[1])
        self.encut: float = dump[2]
        self.realcell: NDArray[np.float64] = dump[3:12].reshape((3, 3))
        self.efermi: float = dump[12]
        self.volume: float = np.linalg.det(self.realcell)
        self.rcpcell: NDArray[np.float64] = np.linalg.inv(self.realcell).T
        unit_cell_vector_magnitude: NDArray[np.float64] = np.linalg.norm(
            self.realcell,
            axis=1,
        )
        cutoff: NDArray[np.float64] = np.ceil(
            np.sqrt(self.encut / Ry_in_eV)
            / (2 * np.pi / (unit_cell_vector_magnitude / Au_in_AA)),
        )
        # FFT Minimum grid size. Always odd!!
        self.ngrid: NDArray[np.int64] = np.array(2 * cutoff + 1, dtype=int)

    def check_DwNGZHalf(self) -> bool:  # noqa: N802
        r"""Check DwNGZHalf.

        self.gamma = True if self gvectors(0)[0] :math:`\neq` nplwvs[0] and
        about half of the number of gvectors equals number of plane waves.
        """
        if self.gamma:
            return True
        if (
            self.gvectors(0).shape[0] // 2 == self.nplwvs[0]
            or self.gvectors(0).shape[0] // 2 + 1 == self.nplwvs[0]
            or self.gvectors(0).shape[0] // 2 - 1 == self.nplwvs[0]
        ):
            self.gamma = True
            return True
        return False

    @property
    def prec(self) -> DTypeLike:
        """Return precision determined from self.rtag."""
        if self.rtag == 45200:  # noqa: PLR2004
            return np.complex64
        if self.rtag == 45210:  # noqa: PLR2004
            return np.complex128
        msg = f"Invalid TAG value: {self.rtag}"
        raise ValueError(msg)

    def band(self) -> None:
        """Read the information about the band from WAVECAR file.

        The information obtained by this method is as follows:

        * Number of plane waves (nplwv)
        * A integer set for k-vectors
        * energy of the band (as a function of spin-, k-, and band index)
        * occupation  (as a function of spin-, k-, and band index)

        """
        self.k_vectors: NDArray[np.float64] = np.zeros((self.num_k, 3), dtype=float)
        self.bands: NDArray[np.float64] = np.zeros(
            (self.n_spin, self.num_k, self.n_bands),
            dtype=float,
        )
        self.nplwvs: NDArray[np.float64] = np.zeros(self.num_k, dtype=int)
        self.occs: NDArray[np.float64] = np.zeros(
            (self.n_spin, self.num_k, self.n_bands),
            dtype=float,
        )
        for spin_i in range(self.n_spin):
            for k_i in range(self.num_k):
                pos = 2 + spin_i * self.num_k * (self.n_bands + 1)
                pos += k_i * (self.n_bands + 1)
                self.wfc.seek(pos * self.recl)
                dump = np.fromfile(self.wfc, dtype=float, count=4 + 3 * self.n_bands)
                if spin_i == 0:
                    self.nplwvs[k_i] = int(dump[0])
                    self.k_vectors[k_i] = dump[1:4]
                dump = dump[4:].reshape((-1, 3))
                self.bands[spin_i, k_i, :] = dump[:, 0]
                self.occs[spin_i, k_i, :] = dump[:, 2]
        if self.num_k == 1:
            self.kpath = None
        else:
            self.kpath = np.concatenate(
                (
                    [
                        0,
                    ],
                    np.cumsum(
                        np.linalg.norm(
                            np.dot(np.diff(self.k_vectors, axis=0), self.rcpcell),
                            axis=1,
                        ),
                    ),
                ),
            )

    def gvectors(self, k_i: float = 0) -> NDArray[np.int64]:
        r"""Return G vector.

        G-vectors :math:`G` is determined by the following condition:
        :math:`\frac{(G+k)^2}{ 2} < E_{cut}`

        note: hbar2over2m is :math:`\frac{\hbar^2}{2m}`

        * 0.529177249 is au unit in AA
        * 13.605826 is Ry unit in eV

        Parameters
        ----------
        k_i: int, optional
            k index :math:`k_i` (the default value is 0).

        Returns
        -------
        NDArray
            G vectors

        """
        k_vector = self.k_vectors[k_i]
        k_grid: NDArray[np.float64] = make_k_grid(
            self.ngrid,
            gamma=self.gamma,
            para=PARALLEL,
        )
        hbar2over2m = 13.605826 * 0.529177249 * 0.529177249
        energy_k = (
            hbar2over2m
            * np.linalg.norm(
                np.dot(k_grid + k_vector[np.newaxis, :], 2 * np.pi * self.rcpcell),
                axis=1,
            )
            ** 2
        )
        g_vec = k_grid[np.where(energy_k < self.encut)[0]]
        return np.array(g_vec, dtype=int)

    def bandcoeff(
        self,
        spin_i: int = 0,
        k_i: int = 0,
        band_i: int = 0,
        *,
        norm: bool = False,
    ) -> NDArray[np.complex128]:
        """Read the coefficient of the planewave of the KS states.

        The KS states is specified by the `spin_i`, `k_i` and `band_i`.

        Parameters
        ----------
        spin_i: int
            spin index :math:`s_i` (0 or 1) (default value is 0)
        k_i: int
            k index :math:`k_i`. Starts with 0 (default value is 0)
        band_i: int
            band index :math:`b_i`. starts with 0 (default value is 0)
        norm: bool
            If true the Band coefficients are normliazed (default is false)

        """
        irec = 3 + spin_i * self.num_k * (self.n_bands + 1)
        irec += k_i * (self.n_bands + 1) + band_i
        self.wfc.seek(irec * self.recl)
        nplw = self.nplwvs[k_i]
        dump = np.fromfile(self.wfc, dtype=self.prec, count=nplw)
        cg: NDArray[np.complex128] = np.array(dump, dtype=np.complex128)
        if norm:
            cg /= np.linalg.norm(cg)
        return cg

    def realspace_wfc(
        self,
        spin_i: int = 0,
        k_i: int = 0,
        band_i: int = 0,
        gvector: NDArray[np.int64] | None = None,
        ngrid: NDArray[np.int64] | None = None,
        *,
        norm: bool = False,
        poscar: POSCAR | None = None,
    ) -> NDArray[np.complex128] | tuple[NDArray[np.float64]] | VASPGrid:
        r"""Return the pseudo-wavefunction in real space.

        Calculate the pseudo-wavefunction of the KS states in
        the real space by using FFT transformation of the reciprocal
        space planewave coefficients.

        The 3D FE grid size is determined by ngrid, which defaults
        to self.ngrid if it is not provided.  GVectors of the KS
        states is used to put 1D plane wave coefficient back to 3D
        grid.

        Parameters
        ----------
        spin_i: int
            spin index (0 or 1). default is 0
        k_i: int
            k index :math:`k_i`. Starts with 0. default is 0
        band_i: int
            band index :math:`b_i`. starts with 0. default is 0.
        norm: bool
            If true the Band coefficients are normalized
        gvector: numpy.array, optional
            G-vector for calculation. (default is self.gvectors(k_i))
        ngrid: numpy.array, optional
            Ngrid for calculation. (default is self.ngrid).
        poscar: vaspy.poscar.POSCAR, optional
            POSCAR object (default is blank POSCAR object)

        Returns
        -------
        numpy.array
            If poscar is not specified, for Collinear-wavecar file.
            data for the wavefunction in the real space.
            .T is due to the original data is made by fortran program
            (i.e. 'VASP', of course).

        tuple
            If poscar is not specified, for SOI-wavecar file.
            the first item is for 'up' wavefunction in real space.
            the second item is for 'down' wavefunction in real space.

        vaspy.mesh3d.VASPGrid
            Returns VASPGrid object, if poscar is specified.  The former frame
            represents the real part of the wavefunction at :math:`k_i` and
            :math:`b_i` in the real space, the latter frame the imaginary
            part. For the SOI wavecar, 4 frames.
            The first and second are for the "up" wavefunction (Real, Imaginary),
            and the third and fourth are "down" wavefunction(Real, Imaginary).
            (Judging SOI by gvectors(k_i).shape[0] :math:`\neq` bandcoeff(k_i).size)

        """
        poscar = POSCAR() if poscar is None else poscar
        ngrid = self.ngrid.copy() if ngrid is None else np.array(ngrid, dtype=int)
        gvec = self.gvectors(k_i) if gvector is None else gvector
        gvec %= ngrid[np.newaxis, :]
        if self.gamma:
            phi_k = (
                np.zeros(
                    (ngrid[0], ngrid[1], ngrid[2] // 2 + 1),
                    dtype=np.complex128,
                )
                if PARALLEL
                else np.zeros(
                    (ngrid[0] // 2 + 1, ngrid[1], ngrid[2]),
                    dtype=np.complex128,
                )
            )
        else:
            phi_k = np.zeros(ngrid, dtype=np.complex128)
        try:  # Collininear
            phi_k[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = self.bandcoeff(
                spin_i,
                k_i,
                band_i,
                norm=norm,
            )
        except ValueError:  # SOI:
            bandcoeff = self.bandcoeff(spin_i, k_i, band_i, norm=norm)
            phi_k = np.zeros((2, ngrid[0], ngrid[1], ngrid[2]), dtype=np.complex128)
            phi_k[0][gvec[:, 0], gvec[:, 1], gvec[:, 2]] = bandcoeff[
                : bandcoeff.size // 2
            ]
            phi_k[1][gvec[:, 0], gvec[:, 1], gvec[:, 2]] = bandcoeff[
                bandcoeff.size // 2 :
            ]
        #
        if self.gamma:
            if PARALLEL:
                for ix in range(ngrid[0]):
                    for iy in range(ngrid[1]):
                        fx = ix if ix < ngrid[0] // 2 + 1 else ix - ngrid[0]
                        fy = iy if iy < ngrid[1] // 2 + 1 else iy - ngrid[1]
                        if (fy > 0) or (fy == 0 and fx >= 0):
                            continue
                        phi_k[ix, iy, 0] = phi_k[-ix, -iy, 0].conjugate()
            else:
                for iz in range(ngrid[2]):
                    for iy in range(ngrid[1]):
                        fz = iz if iz < ngrid[2] // 2 + 1 else iz - ngrid[2]
                        fy = iy if iy < ngrid[1] // 2 + 1 else iy - ngrid[1]
                        if (fy > 0) or (fy == 0 and fz >= 0):
                            continue
                        phi_k[0, iy, iz] = phi_k[0, -iy, -iz].conjugate()
            phi_k /= np.sqrt(2.0)
            phi_k[0, 0, 0] *= np.sqrt(2.0)
            phi_k = restore_gamma_grid(phi_k)
        #
        self.phi_k: NDArray[np.complex128] = phi_k  # For debug
        phi_r: NDArray[np.complex128] = ifftn(phi_k)
        if poscar.scaling_factor == 0:  # poscar is not given.
            if phi_r.ndim == 3:  # noqa: PLR2004
                return phi_r.T
            return (phi_r[0] + phi_r[1]).T, (phi_r[0] - phi_r[1]).T
        vaspgrid = mesh3d.VASPGrid()
        vaspgrid.poscar, vaspgrid.grid.shape = poscar, ngrid
        # checking consistency between POSCAR and WAVECAR
        np.testing.assert_array_almost_equal(
            poscar.scaling_factor * poscar.cell_vecs,
            self.realcell,
        )
        re, im = np.real(phi_r), np.imag(phi_r)
        if phi_r.ndim == 3:  # noqa: PLR2004
            vaspgrid.grid.data = np.concatenate((re.flatten("F"), im.flatten("F")))
        else:  # SOI
            vaspgrid.grid.data = np.concatenate(
                (
                    (re[0] + re[1]).flatten("F"),
                    (im[0] + im[1]).flatten("F"),
                    (re[0] - re[1]).flatten("F"),
                    (im[0] - im[1]).flatten("F"),
                ),
            )
        return vaspgrid

    def __str__(self) -> str:
        """Print out the system parameters."""
        the1stline = "record length  =       {0}  "
        the1stline += "spins =           {1}  "
        the1stline += "prec flag        {2}"
        string = the1stline.format(self.recl, self.n_spin, self.rtag)
        string += f"\nno. k points =          {self.num_k}"
        string += f"\nno. bands =          {self.n_bands}"
        string += "\nreal space lattice vectors:"
        for i in range(3):
            string += "\na" + str(i + 1)
            string += f" = {self.realcell[i][0]}"
            string += f"    {self.realcell[i][1]}"
            string += f"    {self.realcell[i][2]}"
        string += "\n"
        string += f"\nvolume unit cell =   {self.volume}"
        string += "\nReciprocal lattice vectors:"
        for i in range(3):
            string += "\nb" + str(i + 1)
            string += f" = {self.rcpcell[i][0]}    "
            string += f"{self.rcpcell[i][1]}    "
            string += f"{self.rcpcell[i][2]}"
        return string


def make_k_grid(
    ngrid: tuple[int, ...] | NDArray[np.int64] = (),
    *,
    gamma: bool = False,
    para: bool = PARALLEL,
) -> NDArray[np.float64]:
    """Return k_grid.

    Parameters
    ----------
    ngrid: tuple or NDArray
        Grid size
    gamma: bool, default, false
        Set true if only gamma calculations (use vasp with -DwNGZHalf)
    para: bool, optional
        default is global variable `PARALLEL`

    Returns
    -------
    numpy.array

    """
    fx = [
        ii if ii < ngrid[0] // 2 + 1 else ii - ngrid[0]  # <<< // or / (?)
        for ii in range(ngrid[0])
    ]
    fy = [ii if ii < ngrid[1] // 2 + 1 else ii - ngrid[1] for ii in range(ngrid[1])]
    fz = [ii if ii < ngrid[2] // 2 + 1 else ii - ngrid[2] for ii in range(ngrid[2])]
    if gamma and para:
        k_grid: NDArray[np.float64] = np.array(
            [
                (fx[ix], fy[iy], fz[iz])
                for iz in range(ngrid[2])
                for iy in range(ngrid[1])
                for ix in range(ngrid[0])
                if (
                    (fz[iz] > 0)
                    or (fz[iz] == 0 and fy[iy] > 0)
                    or (fz[iz] == 0 and fy[iy] == 0 and fx[ix] >= 0)
                )
            ],
            dtype=float,
        )
    elif gamma and not para:
        k_grid = np.array(
            [
                (fx[ix], fy[iy], fz[iz])
                for iz in range(ngrid[2])
                for iy in range(ngrid[1])
                for ix in range(ngrid[0])
                if (
                    (fz[ix] > 0)
                    or (fz[ix] == 0 and fy[iy] > 0)
                    or (fz[ix] == 0 and fy[iy] == 0 and fx[iz] >= 0)
                )
            ],
            dtype=float,
        )

    else:
        k_grid = np.array(
            [
                (fx[ix], fy[iy], fz[iz])
                for iz in range(ngrid[2])
                for iy in range(ngrid[1])
                for ix in range(ngrid[0])
            ],
            dtype=float,
        )
    return k_grid


def check_symmetry(grid3d: NDArray[np.float64]) -> bool:
    """Return True if grid3d(G) == np.conjugate(grid3d(-G)) for all G.

    Parameters
    ----------
    grid3d: numpy.array
        3D grid data

    Returns
    -------
    Boolean

    """
    assert grid3d.ndim == 3, "Must be 3D Grid"  # noqa: PLR2004
    grid = grid3d.shape
    k_grid = make_k_grid(grid)
    for k in k_grid:
        ix, iy, iz = int(k[0]), int(k[1]), int(k[2])
        if (
            ix >= 0
            and iy >= 0
            and iz >= 0
            and (grid3d[ix][iy][iz] != np.conjugate(grid3d[-ix][-iy][-iz]))
        ):
            print(f"[{ix} {iy} {iz}] is {grid3d[ix][iy][iz]}\n")
            print(
                f"[{-ix} {-iy} {-iz}] is {grid3d[-ix][-iy][-iz]}\n",
            )
            print("check the value\n")
            return False
    return True


def restore_gamma_grid(
    grid3d: NDArray[np.complex128],
    *,
    para: bool = PARALLEL,
) -> NDArray[np.complex128]:
    """Return Grid from the size-reduced matrix for gammareal Wavecar.

    Parameters
    ----------
    grid3d: numpy.array
        3D grid data created with gamma-only version VASP
    para: bool, optional
        default is global variable `PARALLEL`

    """
    assert grid3d.ndim == 3, "Must be 3D Grid"  # noqa: PLR2004
    if para:
        toconj = np.copy(grid3d[:, :, 1:])
        # x=0 slice
        x0slice = toconj[0]
        x0slice = x0slice[:, ::-1]
        x0slice[1:] = x0slice[1:][::-1]
        toconj[0] = x0slice
        # y=0 slice
        toconj[1:, 0, :] = toconj[1:, 0, :][::-1, ::-1]
        # block part
        toconj[1:, 1:, :] = toconj[1:, 1:, :][::-1, ::-1, ::-1]
        return np.concatenate((grid3d, np.conjugate(toconj)), axis=-1)
    toconj = np.copy(grid3d[1:, :, :])
    # z=0 slice
    z0slice = toconj[:, :, 0]
    z0slice = z0slice[::-1, :]
    z0slice[::, 1:] = z0slice[::, 1:][:, ::-1]
    toconj[:, :, 0] = z0slice
    # y = 0 slice
    toconj[:, 0, 1:] = toconj[:, 0, 1:][::-1, ::-1]
    # block part
    toconj[:, 1:, 1:] = toconj[:, 1:, 1:][::-1, ::-1, ::-1]
    return np.concatenate((grid3d, np.conjugate(toconj)), axis=0)
