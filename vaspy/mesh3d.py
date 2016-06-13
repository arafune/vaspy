''' .. py:module:: mesh3D
mesh3D Module to provide class FFTGRID for FFT-grid NG(X,Y,Z)F 
used in CHGCAR, ELFCAR, LOCPOT and the others...
'''

import numpy as np

class FFTGRID(object):
    '''.. py:class:: FFTGRID(meshsize, meshdata)
    Class for NG(X,Y,Z)F in VASP

    This class is used chg_array in CHGCAR, Potential in LOCPOT,
    electron localization function (ELF) in ELFCAR

    Parameters
    ----------

    meshsize: tuple
        mesh_x, mesh_y, mesh_z
    meshdata: list or numpy.array
        1D-list or 1D-numpy array.  The length of grid is mesh_x * mesh_y * mesh_z
    '''
    def __init__(self, meshsize, meshdata):
        if len(meshdata) == meshsize(0) * meshsize(1) * meshsize(2):
            self.mesh_x = meshsize(0)
            self.mesy_y = meshsize(1)
            self.mesh_z = meshsize(2)
            self.meshdata = np.array(meshdata)
        else:
            raise RuntimeError

    
    def slice(self, axis, postition):
        '''.. py:method:: slice(axis, position)

        Parameters
        ----------

        axis: str
            'x', 'y', or 'z'.  Case insensitive.
        position: int
            position for slice

        Return
        ------

        numpy.array
            2D numpy array that sliced from 3D mesh data.
        '''
        pass

    
    def integrate(self, axis, from_coor, to_coor):
        '''.. py:method:: integrate(axis, from_coor, to_coor)
        Return 2D data integrated occupacy along the 'axis' from_coor to to_coor.

        Parameters
        ----------

        axis: str
            'x', 'y', or 'z'.  Case insensitive

        from_coor: int
            'from' value of range of interval integration

        to_coor: int
            'to' value of range interval integration

        Return
        ------
        
        numpy.array
            2D numpy array that integrated from 3D mesh data
        
        '''
        pass

    def __str__(self):
        ''' .. py:method:: __str__()

        Returns
        -------

        str
            a string representation of occupancy.
        '''
        pass

