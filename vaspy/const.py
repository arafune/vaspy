# -*- coding: utf-8 -*-
import numpy as np
'''
Physical constants used in VASP
'''

#  Some important Parameters, to convert to a.u.
#  - AUTOA     =  1. a.u. in Angstroem
#  - RYTOEV    =  1 Ry in Ev
#  - EVTOJ     =  1 eV in Joule
#  - AMTOKG    =  1 atomic mass unit ("proton mass") in kg
#  - BOLKEV    =  Boltzmanns constant in eV/K
#  - BOLK      =  Boltzmanns constant in Joule/K

au_to_A    = 0.529177249
RytoeV   = 13.605826
CLIGHT   =  137.037          # speed of light in a.u.
EVTOJ    = 1.60217733E-19
AMTOKG   = 1.6605402E-27
BOLKEV   = 8.6173857E-5
BOLK     = BOLKEV * EVTOJ
EVTOKCAL = 23.06

# FELECT    =  (the electronic charge)/(4*pi*the permittivity of free space)
#         in atomic units this is just e^2
# EDEPS    =  electron charge divided by the permittivity of free space
#         in atomic units this is just 4 pi e^2
# HSQDTM    =  (plancks CONSTANT/(2*PI))**2/(2*ELECTRON MASS)
#

TPI    = 2 * np.pi
CITPI  = 1j * np.pi
FELECT = 2 * au_to_A * RytoeV
EDEPS  = 4 * np.pi * 2 * RytoeV * au_to_A
HSQDTM = RytoeV * au_to_A * au_to_A

# vector field A times momentum times e/ (2 m_e c) is an energy
# magnetic moments are supplied in Bohr magnetons
# e / (2 m_e c) A(r) p(r)    =  energy
# e / (2 m_e c) m_s x ( r - r_s) / (r-r_s)^3 hbar nabla    = 
# e^2 hbar^2 / (2 m_e^2 c^2) 1/ lenght^3    =  energy
# conversion factor from magnetic moment to energy
# checked independently in SI by Gilles de Wijs

MAGMOMTOENERGY  = 1 / CLIGHT**2 * au_to_A**3 * RytoeV

# dimensionless number connecting input and output magnetic moments
# au_to_A e^2 (2 m_e c^2)
MOMTOMOM   = au_to_A / CLIGHT / CLIGHT / 2
au_to_A2   = au_to_A * au_to_A
au_to_A3   = au_to_A2 * au_to_A
au_to_A4   = au_to_A2 * au_to_A2
au_to_A5   = au_to_A3 * au_to_A2

# dipole moment in atomic units to Debye
AUTDEBYE = 2.541746 
