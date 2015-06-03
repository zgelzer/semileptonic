# Created by Zechariah Gelzer (University of Iowa) on 2015-06-03.
# Copyright (C) 2015 Zechariah Gelzer.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or any later version (see
# <http://www.gnu.org/licenses/>).
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


"""
--------------------------------------------------------------------------------
Defines functions for use in relating particle energies.
--------------------------------------------------------------------------------
Definitions
-----------
__all__ : list of strs
    Functions to be imported during 'from calculators.fcns import *'.
E_out : function
    Calculates outgoing meson energy from lepton momentum transfer.
q2 : function
    Calculates lepton momentum transfer from outgoing meson energy.
--------------------------------------------------------------------------------
"""


__all__ = ['E_out', 'q2']


from settings.constants import hbarc, r1
from settings.fit import decayname, formfactor


def E_out(q2):
    """
    ----------------------------------------------------------------------------
    Calculates outgoing meson energy from lepton momentum transfer q2 in
    (GeV)^2, using conservation of four-momentum in rest frame of incoming
    meson.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    q2 : float
        Lepton momentum transfer q^2 in (GeV)^2.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Outgoing meson energy.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    decayname : str, from settings.fit
    hbarc : float, from settings.constants
    r1 : float, from settings.constants
    ----------------------------------------------------------------------------
    """
    if decayname == 'B2K':
        from settings.constants import mB as M_in, mK as M_out
    elif decayname == 'B2pi':
        from settings.constants import mB as M_in, mpi as M_out
    return (M_in ** 2 + M_out ** 2 - q2 * 1e06 * (r1 / hbarc) ** 2) / (2 * M_in)


def q2(E_out):
    """
    ----------------------------------------------------------------------------
    Calculates lepton momentum transfer from outgoing meson energy E_out in MeV,
    using conservation of four-momentum in rest frame of incoming meson.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    E_out : float
        Outgoing meson energy in MeV.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Lepton momentum transfer.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    decayname : str, from settings.fit
    hbarc : float, from settings.constants
    r1 : float, from settings.constants
    ----------------------------------------------------------------------------
    """
    if decayname == 'B2K':
        from settings.constants import mB as M_in, mK as M_out
    elif decayname == 'B2pi':
        from settings.constants import mB as M_in, mpi as M_out
    return M_in ** 2 + M_out ** 2 - 2 * M_in * E_out * r1 / hbarc

