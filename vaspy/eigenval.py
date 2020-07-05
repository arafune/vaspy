"""This module provides EIGENVAL."""

from __future__ import division, print_function

import csv
import sys
from logging import DEBUG, INFO, Formatter, StreamHandler, getLogger

import numpy as np
from nptyping import NDArray
from typing import Optional, Tuple

from vaspy.tools import open_by_suffix

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.stderr.write("Install matplotlib, or you cannot use methods relating to draw\n")

# logger
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


class EnergyBand(object):
    """Simple band structure object for analyzing by using ipython.

    Class for band structure

    Attributes
    ----------
    kvecs: numpy.ndarray
        kvectors
    kdistances: numpy.ndarray
        kdisance
    numk: int
        number of kpoints
    nbands: int
        number of bands
    energies: numpy.ndarray
        energies[spin_i, k_i, band_i], where spin_i, k_i, and band_i are spin-,
        k- and band-index, respectively.
    label: dict
        used as a label (data 'title' such as '#k', 'Energy') in str format


    Parameters
    ----------
    kvecs: numpy.ndarray
            1D array data of k-vectors.
    energies: numpy.ndarray
            1D array data of energies
    nspin: int
            number of spin: '1' means No-spin.  '2' means collinear spin,
            '4' means noncollinear spin.
            In this class does not distinguish non-collinear spin
            and No-spin.  (default is 1)

    """

    def __init__(self, kvecs=(), energies=(), nspin=1):
        """Initialize."""
        self.kvecs: NDArray = np.array(kvecs)
        self.numk = len(self.kvecs)
        self.label = {}
        try:
            self.nbands = len(energies) // len(kvecs)
        except ZeroDivisionError:
            self.nbands = 0
        self.energies = energies
        self.nspin = nspin
        if self.nspin == 1:  # standard
            self.label["spin"] = [""]
            self.label["energy"] = ["Energy"]
        elif self.nspin == 2:  # spin-polarized
            self.label["energy"] = ["Energy_up", "Energy_down"]
            self.label["spin"] = ["_up", "_down"]
        elif self.nspin == 4:  # non-collinear
            self.label["energy"] = ["Energy"]
            self.label["spin"] = ["_mT", "_mX", "_mY", "_mZ"]
        self.label["k"] = ["#k"]

    @property
    def kdistances(self) -> NDArray:
        """Return kdistances."""
        return np.cumsum(
            np.linalg.norm(
                np.concatenate((np.array([[0, 0, 0]]), np.diff(self.kvecs, axis=0))),
                axis=1,
            )
        )

    def fermi_correction(self, fermi: float) -> None:
        """Correct the Fermi level.

        Parameters
        ----------
        fermi: float
            value of the Fermi level.

        """
        self.energies -= fermi

    def make_label(self, *keys) -> List[str]:
        """Return array the used for label for CSV-like data.

        Parameters
        ----------
        keys: tuple
            key tuple used for label

        """
        label_list = []
        for key in keys:
            for tmp in self.label[key]:
                label_list.append(tmp)
        return label_list

    def to_3dlist(self, **kwargs):
        """Return 3D mentional list.

        list[band_i, [k_i, energy, (energy_down)]]

        This list format would be useful for str output

        """
        bandstructure = []
        for energies in self.energies.T.tolist():
            band = []
            for k, energy in zip(self.kdistances[:, np.newaxis].tolist(), energies):
                k.extend(energy)
                band.append(k)
            bandstructure.append(band)
        return bandstructure

    def to_csv(self, csv_file: str, blankline: bool = True) -> None:
        """Write data to csv file.

        Parameters
        ------------
        csv_file: str
            filename for output
        label_str: str
            string for label (put it on the first line)
        blankline: boolean
            It True (default), the blank line is inserted between band data

        """
        label_str: str = "\t".join(self.make_label("k", "energy")) + "\n"
        with open(csv_file, "w") as fhandle:
            fhandle.writelines(label_str)
            writer = csv.writer(fhandle, delimiter="\t")
            for band_i in self.to_3dlist():
                writer.writerows(band_i)
                if blankline:
                    fhandle.writelines("\n")

    def __str__(self) -> str:
        """Return the str object.

        Returns
        --------
        str
            a string represntation of EnergyBand.
            **Useful for gnuplot and Igor**.

        """
        labels = self.make_label("k", "energy")
        output = labels[0]
        for label in labels[1:]:
            output += "\t" + label
        output += "\n"
        list3d = self.to_3dlist()
        for band_i in list3d:
            for line in band_i:
                output += "{0:.8e}".format(line[0])
                for energy in line[1:]:
                    output += "\t{0:.8e}".format(energy)
                output += "\n"
            output += "\n"
        return output

    def figure(self, color: str = "blue", spin_i: int = 0):
        """Return Axes object of the energy band.

        Parameters
        -----------
        color: str, optional (default is 'blue')
            color of the band line

        spin_i: spin_index
            default is 0

        Returns
        ---------
        matplotlib.pyplot.Axes

        Example
        --------
        Here is a typical code::

            fig = plt.figure()
            ax = band.figure(color='blue')
            ax.set_ylabel('Energy  ( eV )')
            ax.set_ylim(-5, 5)
            ax.set_xlim(0, 4)
            plt.show()

        """
        [
            plt.plot(self.kdistances, self.energies[spin_i, :, band_i], color=color)
            for band_i in range(self.energies.shape[2])
        ]
        return plt.gca()

    def show(
        self, yrange: Optional[Tuple[float, float]] = None, spin_i: int = 0
    ) -> None:  # How to set default value?
        """Draw band structure by using maptlotlib.

        For 'just seeing' use.

        Parameters
        ----------
        yrange: tuple, optional  (default: all range)
            Minimum and maximum value of the y-axis.
            If not specified, use the matplotlib default value.

        spin_i: int  (default is 0 for no spin or 'up' spin)
            Spin index. For spin-polarized collinear band

        """
        for band_i in range(self.energies.shape[2]):
            plt.plot(self.kdistances, self.energies[spin_i, :, band_i], color="blue")
        if yrange is not None:
            plt.ylim([yrange[0], yrange[1]])
        plt.xlim([self.kdistances[0], self.kdistances[-1]])
        plt.ylabel(self.label["energy"][spin_i] + " (eV)")
        plt.show()

    def to_physical_kvector(
        self, recvec=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    ) -> None:
        """Change kvec unit to inverse AA.

        Parameters
        -----------
        recvec: array, numpy.ndarray, optional (default is the unit vector)
            reciprocal vector

        Notes
        -----
            Don't forget that the reciprocal vector used
            in VASP needs 2PI to match  the conventional
            unit of the wavevector.

        """
        logger.debug("recvec: {}".format(recvec))
        logger.debug("self.kvecs: {}".format(self.kvecs))
        recvec = np.array(recvec)
        self.kvecs = np.array([recvec.dot(kvecs) for kvecs in self.kvecs])


class EIGENVAL(EnergyBand):
    """Class for storing the data of EIGENVAL file.

    Parameters
    -----------
    filename: str
        File name of 'EIGENVAL'

    Attributes
    ----------
    natom: int
        Number of atoms

    """

    def __init__(self, filename: Optional[str] = None) -> None:
        """Initialize."""
        super(EIGENVAL, self).__init__()
        self.natom = 0
        #
        if filename:
            self.load_file(open_by_suffix(filename))

    def load_file(self, thefile):
        """Parse EIGENVAL."""
        self.natom, _, _, self.nspin = [int(i) for i in next(thefile).split()]
        if self.nspin == 2:
            self.label["energy"] = ["Energy_up", "Energy_down"]
        else:
            self.label["energy"] = ["Energy"]
        next(thefile)
        next(thefile)
        next(thefile)
        next(thefile)
        _, self.numk, self.nbands = [int(i) for i in next(thefile).split()]
        self.kvecs = []
        self.energies = []
        for _ in range(self.numk):
            # the first line in the sigleset begins with the blank
            next(thefile)
            self.kvecs.append([float(i) for i in next(thefile).split()[0:3]])
            for _ in range(self.nbands):
                self.energies.append(
                    [float(i) for i in next(thefile).split()[1 : self.nspin + 1]]
                )
        self.kvecs = np.array(self.kvecs)
        self.energies = np.array(self.energies).T.reshape(
            self.nspin, self.numk, self.nbands
        )
        thefile.close()
