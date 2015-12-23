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
Calculates pole and loop contributions to B-->pi form factors.
--------------------------------------------------------------------------------
Definitions
-----------
D : function
    Calculates corrections to pole for self-energy contribution to f_perp.
df_para : function
    Calculates loop corrections (incl. wavefunction renormalizations) to f_para.
df_perp : function
    Calculates loop corrections (incl. wavefunction renormalizations) to f_perp.
--------------------------------------------------------------------------------
References
----------
[1] C. Aubin, C. Bernard, "Heavy-Light Semileptonic Decays in Staggered Chiral
    Perturbation Theory", Phys. Rev. D 76, 014002 (2007) [arXiv:0704.0795
    [hep-lat]].
--------------------------------------------------------------------------------
"""


from calculators.chilogs.fcns import *
from math import pi, sqrt
from settings.constants import fpi, tastemults
from settings.fit import SU3


def D(inputs, g_pi):
    """
    ----------------------------------------------------------------------------
    Calculates corrections to pole for self-energy contribution to f_perp, given
    inputs and heavy-light coupling constant g_pi.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : dict of floats
        Input floats for particular experiment. See fileIOs.readers.inputs for
        complete list of inputs.
    g_pi : float
        Heavy-light coupling constant. See constants.py.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Corrections to pole for self-energy contribution to f_perp [1]; has
        fewer terms if not using SU(3) fit.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    Deltas : function, from calculators.chilogs.fcns
    J1sub : function, from calculators.chilogs.fcns
    R31 : function, from calculators.chilogs.fcns
    SU3: bool, from settings.fit
    deltaps : function, from calculators.chilogs.fcns
    fpi : float, from settings.constants
    mu : function, from calculators.chilogs.fcns
    pi : float, from math
    sqrt : function, from math
    tastemults : dict, from settings.constants
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Aubin, C. Bernard, "Heavy-Light Semileptonic Decays in Staggered
        Chiral Perturbation Theory", Phys. Rev. D 76, 014002 (2007)
        [arXiv:0704.0795 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    E = inputs['E']
    m_l = inputs['ml_sea']
    m_h = inputs['mh_sea']
    g_pi2 = g_pi ** 2
    Xisum = 0.
    for Xi in tastemults.keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_K = sqrt(mu(inputs['a_fm']) * (m_l + m_h) +
                   Deltas(inputs['a_fm'])[Xi])
        Xisum += ((2 * J1sub(m_pi, E) + SU3 * J1sub(m_K, E)) * tastemults[Xi])
    jsum = 0.
    for Xi in deltaps(inputs['a_fm']).keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])[Xi])
        m_pi2 = m_pi ** 2
        m_S2 = m_S ** 2
        Z = sqrt((m_S2 - m_pi2) ** 2 -
                 (deltaps(inputs['a_fm'])[Xi] / 2) * (m_S2 - m_pi2) +
                 (9. / 16) * deltaps(inputs['a_fm'])[Xi] ** 2)
        m_eta = sqrt((m_pi2 + m_S2 + (3. / 4) *
                      deltaps(inputs['a_fm'])[Xi] - Z) / 2)
        m_etap = sqrt((m_pi2 + m_S2 + (3. / 4) *
                       deltaps(inputs['a_fm'])[Xi] + Z) / 2)
        if SU3:
            for j, m_j in [[1, m_pi], [2, m_eta], [3, m_etap]]:
                jsum -= (deltaps(inputs['a_fm'])[Xi] *
                         R31(m_S, m_pi, m_eta, m_etap, j) * J1sub(m_j, E))
        else:
            jsum += 2 * (J1sub(m_eta, E) - J1sub(m_pi, E))
    m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])['I'])
    m_eta = sqrt(mu(inputs['a_fm']) * ((2 * m_l + 4 * m_h) / 3) +
                 Deltas(inputs['a_fm'])['I'])
    D = ((Xisum / 16) + jsum - (J1sub(m_pi, E) / 2) +
         SU3 * (J1sub(m_eta, E) / 6))
    return (-3 * g_pi2 * E * D) / (4 * pi * fpi) ** 2


def df_para(inputs, g_pi):
    """
    ----------------------------------------------------------------------------
    Calculates loop corrections (incl. wavefunction renormalizations) to f_para
    of particular decay, given inputs and heavy-light coupling constant g_pi.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : dict of floats
        Input floats for particular experiment. See fileIOs.readers.inputs for
        complete list of inputs.
    g_pi : float
        Heavy-light coupling constant. See constants.py.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Loop corrections (incl. wavefunction renormalizations) to f_para [1];
        has fewer terms if not using SU(3) fit.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    Deltas : function, from calculators.chilogs.fcns
    I1 : function, from calculators.chilogs.fcns
    I2 : function, from calculators.chilogs.fcns
    R31 : function, from calculators.chilogs.fcns
    SU3: bool, from settings.fit
    deltaps : function, from calculators.chilogs.fcns
    fpi : float, from settings.constants
    mu : function, from calculators.chilogs.fcns
    pi : float, from math
    sqrt : function, from math
    tastemults : dict, from settings.constants
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Aubin, C. Bernard, "Heavy-Light Semileptonic Decays in Staggered
        Chiral Perturbation Theory", Phys. Rev. D 76, 014002 (2007)
        [arXiv:0704.0795 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    E = inputs['E']
    m_l = inputs['ml_sea']
    m_h = inputs['mh_sea']
    g_pi2 = g_pi ** 2
    Xisum = 0.
    for Xi in tastemults.keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_K = sqrt(mu(inputs['a_fm']) * (m_l + m_h) +
                   Deltas(inputs['a_fm'])[Xi])
        Xisum += ((((1 - 3 * g_pi2) / 2) * (2 * I1(m_pi) + SU3 * I1(m_K)) +
                   2 * I2(m_pi, E) + SU3 * I2(m_K, E)) * tastemults[Xi])
    jsum = 0.
    for Xi in deltaps(inputs['a_fm']).keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])[Xi])
        m_pi2 = m_pi ** 2
        m_S2 = m_S ** 2
        Z = sqrt((m_S2 - m_pi2) ** 2 -
                 (deltaps(inputs['a_fm'])[Xi] / 2) * (m_S2 - m_pi2) +
                 (9. / 16) * deltaps(inputs['a_fm'])[Xi] ** 2)
        m_eta = sqrt((m_pi2 + m_S2 + (3. / 4) *
                      deltaps(inputs['a_fm'])[Xi] - Z) / 2)
        m_etap = sqrt((m_pi2 + m_S2 + (3. / 4) *
                       deltaps(inputs['a_fm'])[Xi] + Z) / 2)
        if SU3:
            for j, m_j in [[1, m_pi], [2, m_eta], [3, m_etap]]:
                jsum += (deltaps(inputs['a_fm'])[Xi] *
                         R31(m_S, m_pi, m_eta, m_etap, j) *
                         ((3. / 2) * (g_pi2 - 1) * I1(m_j) - 2 * I2(m_j, E)))
        else:
            jsum += (3 * (g_pi2 - 1) * (I1(m_pi) - I1(m_eta)) +
                     4 * (I2(m_eta, E) - I2(m_pi, E)))
    m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])['I'])
    m_eta = sqrt(mu(inputs['a_fm']) * ((2 * m_l + 4 * m_h) / 3) +
                 Deltas(inputs['a_fm'])['I'])
    df = ((Xisum / 16) + jsum +
          ((1 + 3 * g_pi2) / 4) * (I1(m_pi) - SU3 * (I1(m_eta) / 3)))
    return df / (4 * pi * fpi) ** 2


def df_perp(inputs, g_pi):
    """
    ----------------------------------------------------------------------------
    Calculates loop corrections (incl. wavefunction renormalizations) to f_perp
    of particular decay, given inputs and heavy-light coupling constant g_pi.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : dict of floats
        Input floats for particular experiment. See fileIOs.readers.inputs for
        complete list of inputs.
    g_pi : float
        Heavy-light coupling constant. See constants.py.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Loop corrections (incl. wavefunction renormalizations) to f_perp [1];
        has fewer terms if not using SU(3) fit.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    Deltas : function, from calculators.chilogs.fcns
    I1 : function, from calculators.chilogs.fcns
    J1sub : function, from calculators.chilogs.fcns
    R31 : function, from calculators.chilogs.fcns
    SU3: bool, from settings.fit
    deltaps : function, from calculators.chilogs.fcns
    fpi : float, from settings.constants
    mu : function, from calculators.chilogs.fcns
    pi : float, from math
    sqrt : function, from math
    tastemults : dict, from settings.constants
    ----------------------------------------------------------------------------
    References
    ----------
    [1] C. Aubin, C. Bernard, "Heavy-Light Semileptonic Decays in Staggered
        Chiral Perturbation Theory", Phys. Rev. D 76, 014002 (2007)
        [arXiv:0704.0795 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    E = inputs['E']
    m_l = inputs['ml_sea']
    m_h = inputs['mh_sea']
    g_pi2 = g_pi ** 2
    Xisum = 0.
    for Xi in tastemults.keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_K = sqrt(mu(inputs['a_fm']) * (m_l + m_h) +
                   Deltas(inputs['a_fm'])[Xi])
        Xisum -= (2 * I1(m_pi) + SU3 * I1(m_K)) * tastemults[Xi]
    Xisum *= (1 + 3 * g_pi2) / 2
    jsum = 0.
    for Xi in deltaps(inputs['a_fm']).keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])[Xi])
        m_pi2 = m_pi ** 2
        m_S2 = m_S ** 2
        Z = sqrt((m_S2 - m_pi2) ** 2 -
                 (deltaps(inputs['a_fm'])[Xi] / 2) * (m_S2 - m_pi2) +
                 (9. / 16) * deltaps(inputs['a_fm'])[Xi] ** 2)
        m_eta = sqrt((m_pi2 + m_S2 + (3. / 4) *
                      deltaps(inputs['a_fm'])[Xi] - Z) / 2)
        m_etap = sqrt((m_pi2 + m_S2 + (3. / 4) *
                       deltaps(inputs['a_fm'])[Xi] + Z) / 2)
        if SU3:
            for j, m_j in [[1, m_pi], [2, m_eta], [3, m_etap]]:
                jsum += (deltaps(inputs['a_fm'])[Xi] *
                         R31(m_S, m_pi, m_eta, m_etap, j) *
                         (g_pi2 * J1sub(m_j, E) +
                          ((1 + 3 * g_pi2) / 2) * I1(m_j)))
        else:
            jsum += (2 * g_pi2 * (J1sub(m_pi, E) - J1sub(m_eta, E)) +
                     (1 + 3 * g_pi2) * (I1(m_pi) - I1(m_eta)))
    m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])['I'])
    m_eta = sqrt(mu(inputs['a_fm']) * ((2 * m_l + 4 * m_h) / 3) +
                 Deltas(inputs['a_fm'])['I'])
    df = ((Xisum / 16) + jsum -
          (g_pi2 / 2) * J1sub(m_pi, E) +
          SU3 * (g_pi2 / 6) * J1sub(m_eta, E) +
          ((1 + 3 * g_pi2) / 12) * (3 * I1(m_pi) - SU3 * I1(m_eta)))
    return df / (4 * pi * fpi) ** 2

