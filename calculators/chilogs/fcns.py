"""Defines functions for use in calculating terms related to chiral logs."""


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


__all__ = ['a_fermi', 'Deltabar', 'deltas', 'Deltas', 'F', 'I1', 'I2', 'J1Sub',
           'mu', 'R31']


from settings.constants import *
from settings.fit import *
from math import atan, atanh, log, pi, sqrt


def a_fermi(a_str):
    """Calculates the lattice spacing in fermi corresponding to a given string
    name."""
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


def D(inputs, g_pi, decayname):
    """Calculates correction to pole for self-energy contribution to f_perp of
    a given decay, according to given inputs and g_pi."""
    if decayname == 'B2K':
        from calculators.chilogs import B2K
        return B2K.D(inputs, g_pi)
    elif decayname == 'B2pi':
        from calculators.chilogs import B2pi
        return B2pi.D(inputs, g_pi)


def Deltabar(a_fm):
    """Calculates root-mean-square taste splitting for some lattice spacing a_fm
    in fermi."""
    Deltabar = 0.
    for Xi in tastemults.keys():
        Deltabar += tastemults[Xi] * Deltas(a_fm)[Xi] ** 2
    return sqrt(Deltabar / 16.)


def deltas(a_fm):
    """Returns hairpin parameters for some lattice spacing a_fm in fermi."""
    if abs(a_fm - 0.12) < a_atol:
        return {'A': -2.800000e-01,
                'V':  1.000000e-10}
    elif abs(a_fm - 0.09) < a_atol:
        return {'A': -9.506479e-02,
                'V':  3.395171e-11}
    elif abs(a_fm - 0.06) < a_atol:
        return {'A': -3.307382e-02,
                'V':  1.181208e-11}
    elif abs(a_fm - 0.045) < a_atol:
        return {'A': -1.306645e-02,
                'V':  4.666589e-12}
    elif abs(a_fm) < a_atol:
        return {'A':  0.000000e+00,
                'V':  0.000000e+00}
    else:
        raise ValueError('invalid lattice spacing in inputs')


def Deltas(a_fm):
    """Returns taste splittings for some lattice spacing a_fm in fermi."""
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
        raise ValueError('invalid lattice spacing in inputs')


def deltas_calc(a_fm):
    """Calculates hairpin parameters for some lattice spacing a_fm in fermi."""
    deltas = {'V': nearzero, 'A': -0.28}
    Deltabar_coarse = Deltabar(a_fermi('coarse'))
    if ((abs(a_fm - 0.12) < a_atol) or (abs(a_fm - 0.09) < a_atol) or
          (abs(a_fm - 0.06) < a_atol) or (abs(a_fm - 0.045) < a_atol) or
          (abs(a_fm) < a_atol)):
        ratio = Deltabar(a_fm) / Deltabar_coarse
    else:
        raise ValueError('invalid lattice spacing in inputs')
    for key in deltas.keys():
        deltas[key] *= ratio
    return deltas


def df_para(inputs, g_pi, decayname):
    """Calculates loop corrections to f_para of a given decay, according to
    given inputs and g_pi."""
    if decayname == 'B2K':
        from calculators.chilogs import B2K
        return B2K.df_para(inputs, g_pi)
    elif decayname == 'B2pi':
        from calculators.chilogs import B2pi
        return B2pi.df_para(inputs, g_pi)


def df_perp(inputs, g_pi, decayname):
    """Calculates loop corrections to f_perp of a given decay, according to
    given inputs and g_pi."""
    if decayname == 'B2K':
        from calculators.chilogs import B2K
        return B2K.df_perp(inputs, g_pi)
    elif decayname == 'B2pi':
        from calculators.chilogs import B2pi
        return B2pi.df_perp(inputs, g_pi)


def F(x):
    """Calculates arctan(h)'s for use in calculators.chilogs.fcns.{I1, I2}."""
    if 0 <= x <= 1:
        arg = sqrt(1. - x ** 2)
        return arg * atanh(arg)
    elif x >= 1:
        arg = sqrt(x ** 2 - 1.)
        return -arg * atan(arg)
    else:
        raise ValueError('x must be positive')


def I1(m):
    """Calculates chiral loop integral term for given mass m."""
    m2 = m ** 2
    return m2 * log(m2 / Lambda ** 2)


def I2(m, Delta):
    """Calculates chiral loop integral term for given mass m and pole Delta."""
    if hardPiK:
        return -I1(m)
    else:
        Delta2 = Delta ** 2
        return (-2 * Delta2 * log(m ** 2 / Lambda ** 2) -
                4 * Delta2 * F(m / Delta) + 2 * Delta2)


def J1Sub(m, Delta):
    """Calculates singularity-subtracted chiral loop integral cross-term for
    given mass m and pole Delta."""
    if hardPiK:
        return 0
    else:
        m2 = m ** 2
        Delta2 = Delta ** 2
        return ((-m2 + (2 / 3.) * Delta2) * log(m2 / Lambda**2) +
                (4 / 3.) * (Delta2 - m2) * F(m / Delta) - (10 / 9.) * Delta2 +
                (4 / 3.) * m2 - (2 * pi * m**3) / (3 * Delta))


def mu(a_fm):
    """Calculates mass slope constant for some lattice spacing a_fm in fermi."""
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
        raise ValueError('invalid lattice spacing in inputs')


def R31(m_mu, m1, m2, m3, j):
    """Calculates Euclidean space residue R_j^[3,1]({m1, m2, m3}; m_mu)."""
    if j==1:
        mj2 = m1 ** 2
        return (m_mu ** 2 - mj2) / ((m2 ** 2 - mj2) * (m3 ** 2 - mj2))
    elif j==2:
        mj2 = m2 ** 2
        return (m_mu ** 2 - mj2) / ((m1 ** 2 - mj2) * (m3 ** 2 - mj2))
    elif j==3:
        mj2 = m3 ** 2
        return (m_mu ** 2 - mj2) / ((m1 ** 2 - mj2) * (m2 ** 2 - mj2))
    else:
        raise ValueError('1<=j<=3 is required for R31')

