# Created by Zechariah Gelzer (University of Iowa) on 2015-03-30.
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
Defines functions for use in calculating terms related to chiral logs.
--------------------------------------------------------------------------------
Definitions
-----------
__all__ : list of strs
    Functions to be imported during 'from calculators.chilogs.fcns import *'.
Deltabar : function
    Calculates root-mean-square taste splitting for lattice spacing.
Deltas : function
    Returns taste splittings for lattice spacing.
F : function
    Calculates arctan(h)s for use in functions I2, J1sub.
I1 : function
    Calculates chiral loop integral.
I2 : function
    Calculates chiral loop integral with pole.
J1sub : function
    Calculates singularity-subtracted chiral loop integral with pole.
R31 : function
    Calculates Euclidean space residue R_j^{[3, 1]}.
a_fermi : function
    Calculates lattice spacing in fm.
dD : function
    Calculates corrections to pole for self-energy contribution to f_perp.
deltaps : function
    Returns hairpin parameters for lattice spacing.
deltaps_calc : function
    Calculates hairpin parameters for lattice spacing.
df_para : function
    Calculates loop corrections (incl. wavefunction renormalizations) to f_para.
df_perp : function
    Calculates loop corrections (incl. wavefunction renormalizations) to f_perp.
mu : function
    Calculates mass slope constant for lattice spacing.
--------------------------------------------------------------------------------
"""


__all__ = ['Deltabar', 'Deltas', 'F', 'I1', 'I2', 'J1sub', 'R31', 'deltaps',
           'mu']


from math import atan, atanh, log, pi, sqrt
from settings.constants import Lambda, a_atol, nearzero, tastemults
from settings.fit import decayname, hardpiK


def Deltabar(a_fm):
    """
    ----------------------------------------------------------------------------
    Calculates root-mean-square taste splitting for lattice spacing a_fm in fm.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    a_fm : float
        Lattice spacing in fm.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Root-mean-square taste splitting scaled with square of lattice spacing a
        in r_1 units.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    Deltas : function
    numpy : module, as np
    sqrt : function, from math
    tastemults : dict, from settings.constants
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function Deltas.
    ----------------------------------------------------------------------------
    """
    Deltabar = 0.
    for Xi in tastemults.keys():
        Deltabar += tastemults[Xi] * Deltas(a_fm)[Xi] ** 2
    return sqrt(Deltabar / 16.)


def Deltas(a_fm):
    """
    ----------------------------------------------------------------------------
    Returns taste splittings for lattice spacing a_fm in fermi.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    a_fm : float
        Lattice spacing in fm.
    ----------------------------------------------------------------------------
    Returns
    -------
    dict of floats
        Taste splittings scaled with square of lattice spacing a in r_1 units
        [1]. See constants.py for description of taste names.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    a_atol : float, from settings.constants
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'invalid lattice spacing'
        Lattice spacing in fm must be one of: 0.12, 0.09, 0.06, 0.045, or 0.0.
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Bernard, et al. (Fermilab Lattice and MILC Collaborations),
        "mass^2 splittings in massind r1 units",
        <http://physics.wustl.edu/~cb/Fermilab-MILC/secure/>, accessed
        2015-03-30.
    ----------------------------------------------------------------------------
    """
    if abs(a_fm - 0.12) < a_atol:
        return {'P': 0.000000e+00,
                'A': 2.270460e-01,
                'T': 3.661620e-01,
                'V': 4.802591e-01,
                'I': 6.008212e-01}
    elif abs(a_fm - 0.09) < a_atol:
        return {'P': 0.000000e+00,
                'A': 7.469220e-02,
                'T': 1.237760e-01,
                'V': 1.593220e-01,
                'I': 2.206520e-01}
    elif abs(a_fm - 0.06) < a_atol:
        return {'P': 0.000000e+00,
                'A': 2.634800e-02,
                'T': 4.297780e-02,
                'V': 5.743780e-02,
                'I': 7.038790e-02}
    elif abs(a_fm - 0.045) < a_atol:
        return {'P': 0.000000e+00,
                'A': 1.040930e-02,
                'T': 1.697920e-02,
                'V': 2.269190e-02,
                'I': 2.780810e-02}
    elif abs(a_fm) < a_atol:
        return {'P': 0.000000e+00,
                'A': 0.000000e+00,
                'T': 0.000000e+00,
                'V': 0.000000e+00,
                'I': 0.000000e+00}
    else:
        raise ValueError('invalid lattice spacing')


def F(x):
    """
    ----------------------------------------------------------------------------
    Calculates arctan(h)s for use in functions I2, J1sub.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    x : float
        As implemented in I2, J1sub, (x = {meson mass} / {pole mass}).
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        [(-1 * ) -arg * arctan(h){arg}] for argument (arg = sqrt(diff{x^2, 1})),
        where diff is positive difference [1].
    ----------------------------------------------------------------------------
    Requirements
    ------------
    atan : function, from math
    atanh : function, from math
    sqrt : function, from math
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError
        x must be positive.
    ----------------------------------------------------------------------------
    Notes
    -----
    + See functions I2, J1sub.
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Aubin and C. Bernard, "Heavy-Light Semileptonic Decays in Staggered
        Chiral Perturbation Theory", Phys. Rev. D 76, 014002 (2007)
        [arXiv:0704.0795 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    if 0 <= x <= 1:
        arg = sqrt(1. - x ** 2)
        return arg * atanh(arg)
    elif x >= 1:
        arg = sqrt(x ** 2 - 1.)
        return -arg * atan(arg)
    else:
        raise ValueError('x must be positive')


def I1(m):
    """
    ----------------------------------------------------------------------------
    Calculates chiral loop integral for given mass m.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    m : float
        Meson mass.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Chiral loop integral (I1 = m^2 * ln(m^2 / Lambda^2)) [1].
    ----------------------------------------------------------------------------
    Requirements
    ------------
    Lambda : float, from settings.constants
    log : function, from math
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Aubin and C. Bernard, "Heavy-Light Semileptonic Decays in Staggered
        Chiral Perturbation Theory", Phys. Rev. D 76, 014002 (2007)
        [arXiv:0704.0795 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    m2 = m ** 2
    return m2 * log(m2 / Lambda ** 2)


def I2(m, Delta):
    """
    ----------------------------------------------------------------------------
    Calculates chiral loop integral for given mass m and pole Delta.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    x : float
        Meson mass.
    Delta : float
        Pole mass; as implemented in B2{K or pi}, energy of pion/kaon.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Chiral loop integral; simplifies to -I1(m) if using hard pion/kaon [1].
    ----------------------------------------------------------------------------
    Requirements
    ------------
    F : function
    I1 : function
    Lambda : float, from settings.constants
    hardpiK : bool, from settings.fit
    log : function, from math
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Aubin and C. Bernard, "Heavy-Light Semileptonic Decays in Staggered
        Chiral Perturbation Theory", Phys. Rev. D 76, 014002 (2007)
        [arXiv:0704.0795 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    if hardpiK:
        return -I1(m)
    else:
        Delta2 = Delta ** 2
        return (-2 * Delta2 * log(m ** 2 / Lambda ** 2) -
                4 * Delta2 * F(m / Delta) + 2 * Delta2)


def J1sub(m, Delta):
    """
    ----------------------------------------------------------------------------
    Calculates singularity-subtracted chiral loop integral for given mass m and
    pole Delta.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    x : float
        Meson mass.
    Delta : float
        Pole mass; as implemented in B2{K or pi}, energy of pion/kaon.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Singularity-subtracted chiral loop integral for cross-term; vanishes if
        using hard pion/kaon [1].
    ----------------------------------------------------------------------------
    Requirements
    ------------
    F : function
    Lambda : float, from settings.constants
    hardpiK : bool, from settings.fit
    log : function, from math
    pi : float, from math
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Aubin and C. Bernard, "Heavy-Light Semileptonic Decays in Staggered
        Chiral Perturbation Theory", Phys. Rev. D 76, 014002 (2007)
        [arXiv:0704.0795 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    if hardpiK:
        return 0
    else:
        m2 = m ** 2
        Delta2 = Delta ** 2
        return ((-m2 + (2 / 3.) * Delta2) * log(m2 / Lambda ** 2) +
                (4 / 3.) * (Delta2 - m2) * F(m / Delta) - (10 / 9.) * Delta2 +
                (4 / 3.) * m2 - (2 * pi * m ** 3) / (3 * Delta))


def R31(m_mu, m1, m2, m3, j):
    """
    ----------------------------------------------------------------------------
    Calculates Euclidean space residue R_j^{[3, 1]}({m1, m2, m3}; m_mu).
    ----------------------------------------------------------------------------
    Parameters
    ----------
    m_mu : float
        Mass of (numerator) meson.
    m1 : float
        Mass of (denominator) meson.
    m2 : float
        Mass of (denominator) meson.
    m3 : float
        Mass of (denominator) meson.
    j : int
        Index representing particular meson.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Euclidean space residue for sets of 3 denominator, 1 numerator meson
        masses [1].
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError
        (1 <= j <= 3) is required for R_j^{[3, 1]}.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    nearzero : float, from settings.constants
    ----------------------------------------------------------------------------
    Notes
    -----
    + Could lead to delicate cancellation between large terms.
    + SU(3) fits may result in vanishing denominators in R_j^{[3, 1]} during
      continuum extrapolation. This is avoided by instead returning
      (1 / {approx. zero}) via settings.constants.nearzero. Since R_j^{[3, 1]}
      is scaled by (a^2 * delta^{\\prime}s = deltaps(a_fm)), which are all zero
      in continuum, continuum SU(3) results from R31 are ultimately annihilated.
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Aubin and C. Bernard, "Heavy-Light Semileptonic Decays in Staggered
        Chiral Perturbation Theory", Phys. Rev. D 76, 014002 (2007)
        [arXiv:0704.0795 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    if j==1:
        mj2 = m1 ** 2
        try:
            return (m_mu ** 2 - mj2) / ((m2 ** 2 - mj2) * (m3 ** 2 - mj2))
        except ZeroDivisionError:
            return 1 / nearzero
    elif j==2:
        mj2 = m2 ** 2
        try:
            return (m_mu ** 2 - mj2) / ((m1 ** 2 - mj2) * (m3 ** 2 - mj2))
        except ZeroDivisionError:
            return 1 / nearzero
    elif j==3:
        mj2 = m3 ** 2
        try:
            return (m_mu ** 2 - mj2) / ((m1 ** 2 - mj2) * (m2 ** 2 - mj2))
        except ZeroDivisionError:
            return 1 / nearzero
    else:
        raise ValueError('1<=j<=3 is required for R31')


def a_fermi(a_str):
    """
    ----------------------------------------------------------------------------
    Calculates lattice spacing in fm corresponding to given string name.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    a_str : str
        Lattice spacing name.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Lattice spacing in fm.
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'invalid lattice spacing name'
        Lattice spacing name must be one of: 'coarse', 'fine', 'super_fine',
        'ultra_fine', or 'continuum'.
    ----------------------------------------------------------------------------
    """
    if a_str == 'coarse':
        return 0.12
    elif a_str == 'fine':
        return 0.09
    elif a_str == 'super_fine':
        return 0.06
    elif a_str == 'ultra_fine':
        return 0.045
    elif a_str == 'continuum':
        return 0
    else:
        raise ValueError('invalid lattice spacing name')


def dD(inputs, g_pi):
    """
    ----------------------------------------------------------------------------
    Calculates corrections to pole for self-energy contribution to f_perp of
    particular decay, given inputs and heavy-light coupling constant g_pi.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : dict of floats
        Input floats for particular form factor. See fileIOs.readers.inputs for
        complete list of inputs.
    g_pi : float
        Heavy-light coupling constant. See constants.py.
    ----------------------------------------------------------------------------
    Returns
    -------
    float, from function B2{K or pi}.dD
        Corrections to pole for self-energy contribution to f_perp.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    decayname : str, from settings.constants
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function B2{K or pi}.dD.
    ----------------------------------------------------------------------------
    """
    if decayname == 'B2K':
        from calculators.chilogs import B2K
        return B2K.dD(inputs, g_pi)
    elif decayname == 'B2pi':
        from calculators.chilogs import B2pi
        return B2pi.dD(inputs, g_pi)


def deltaps(a_fm):
    """
    ----------------------------------------------------------------------------
    Returns hairpin parameters for lattice spacing a_fm in fm.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    a_fm : float
        Lattice spacing in fm.
    ----------------------------------------------------------------------------
    Returns
    -------
    dict of floats
        Hairpin parameters of axial-vector taste 'A' and vector taste 'V'.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    a_atol : float, from settings.constants
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'invalid lattice spacing'
        Lattice spacing in fm must be one of: 0.12, 0.09, 0.06, 0.045, or 0.0.
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function deltaps_calc.
    ----------------------------------------------------------------------------
    """
    if abs(a_fm - 0.12) < a_atol:
        return {'A': -2.800000e-01,
                'V':  0.000000e+00}
    elif abs(a_fm - 0.09) < a_atol:
        return {'A': -9.506479e-02,
                'V':  0.000000e+00}
    elif abs(a_fm - 0.06) < a_atol:
        return {'A': -3.307382e-02,
                'V':  0.000000e+00}
    elif abs(a_fm - 0.045) < a_atol:
        return {'A': -1.306645e-02,
                'V':  0.000000e+00}
    elif abs(a_fm) < a_atol:
        return {'A':  0.000000e+00,
                'V':  0.000000e+00}
    else:
        raise ValueError('invalid lattice spacing')


def deltaps_calc(a_fm):
    """
    Calculates hairpin parameters for lattice spacing a_fm in fermi.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    a_fm : float
        Lattice spacing in fm.
    ----------------------------------------------------------------------------
    Returns
    -------
    deltaps : dict of floats
        Hairpin parameters of axial-vector taste 'A' and vector taste 'V'.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    Deltabar : function
    a_atol : float, from settings.constants
    a_fermi : function
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'invalid lattice spacing'
        Lattice spacing in fm must be one of: 0.12, 0.09, 0.06, 0.045, or 0.0.
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function deltaps.
    ----------------------------------------------------------------------------
    """
    deltaps = {'A': -0.28, 'V': 0.00}
    Deltabar_coarse = Deltabar(a_fermi('coarse'))
    if ((abs(a_fm - 0.12) < a_atol) or (abs(a_fm - 0.09) < a_atol) or
          (abs(a_fm - 0.06) < a_atol) or (abs(a_fm - 0.045) < a_atol) or
          (abs(a_fm) < a_atol)):
        ratio = Deltabar(a_fm) / Deltabar_coarse
    else:
        raise ValueError('invalid lattice spacing')
    for taste in deltaps.keys():
        deltaps[taste] *= ratio
    return deltaps


def df_para(inputs, g_pi):
    """
    ----------------------------------------------------------------------------
    Calculates loop corrections (incl. wavefunction renormalizations) to f_para
    of particular decay, given inputs and heavy-light coupling constant g_pi.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : dict of floats
        Input floats for particular form factor. See fileIOs.readers.inputs for
        complete list of inputs.
    g_pi : float
        Heavy-light coupling constant. See constants.py.
    ----------------------------------------------------------------------------
    Returns
    -------
    float, from function B2{K or pi}.df_para
        Loop corrections (incl. wavefunction renormalizations) to f_para.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    decayname : str, from settings.constants
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function B2{K or pi}.df_para.
    ----------------------------------------------------------------------------
    """
    if decayname == 'B2K':
        from calculators.chilogs import B2K
        return B2K.df_para(inputs, g_pi)
    elif decayname == 'B2pi':
        from calculators.chilogs import B2pi
        return B2pi.df_para(inputs, g_pi)


def df_perp(inputs, g_pi):
    """
    ----------------------------------------------------------------------------
    Calculates loop corrections (incl. wavefunction renormalizations) to f_perp
    of particular decay, given inputs and heavy-light coupling constant g_pi.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : dict of floats
        Input floats for particular form factor. See fileIOs.readers.inputs for
        complete list of inputs.
    g_pi : float
        Heavy-light coupling constant. See constants.py.
    ----------------------------------------------------------------------------
    Returns
    -------
    float, from function B2{K or pi}.df_perp
        Loop corrections (incl. wavefunction renormalizations) to f_perp.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    decayname : str, from settings.constants
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function B2{K or pi}.df_perp.
    ----------------------------------------------------------------------------
    """
    if decayname == 'B2K':
        from calculators.chilogs import B2K
        return B2K.df_perp(inputs, g_pi)
    elif decayname == 'B2pi':
        from calculators.chilogs import B2pi
        return B2pi.df_perp(inputs, g_pi)


def mu(a_fm):
    """
    ----------------------------------------------------------------------------
    Calculates mass slope constant for lattice spacing a_fm in fermi.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    a_fm : float
        Lattice spacing in fm.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Coefficient of (m_l + m_l) in (m_pi^2) in r_1 units [1].
    ----------------------------------------------------------------------------
    Requirements
    ------------
    a_atol : float, from settings.constants
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'invalid lattice spacing'
        Lattice spacing in fm must be one of: 0.12, 0.09, 0.06, 0.045, or 0.0.
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Bernard, et al. (Fermilab Lattice and MILC Collaborations),
        "slopes in mass-independent r1 units",
        <http://physics.wustl.edu/~cb/Fermilab-MILC/secure/>, accessed
        2015-03-30.
    ----------------------------------------------------------------------------
    """
    if abs(a_fm - 0.12) < a_atol:
        return 6.831904e+00
    elif abs(a_fm - 0.09) < a_atol:
        return 6.638563e+00
    elif abs(a_fm - 0.06) < a_atol:
        return 6.486649e+00
    elif abs(a_fm - 0.045) < a_atol:
        return 6.417427e+00
    elif abs(a_fm) < a_atol:
        return 6.015349e+00
    else:
        raise ValueError('invalid lattice spacing')

