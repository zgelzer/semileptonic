"""Calculates pole and loop contributions to B-->K form factors."""


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


from calculators.chilogs.fcns import *
from settings.constants import *
from settings.fit import *
from math import pi, sqrt


def D(inputs, g_pi):
    """Calculates correction to pole for self-energy contribution to f_perp,
    according to given inputs, g_pi (the heavy-light coupling constant), and
    global constants."""
    E = inputs['E']
    m_l = inputs['ml_val']
    m_h = inputs['mh_val']
    g_pi2 = g_pi ** 2
    if not SU3:
        return 0
    else:
        XiSum = 0.
        for Xi in tastemults.keys():
            m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h +
                       Deltas(inputs['a_fm'])[Xi])
            m_K = sqrt(mu(inputs['a_fm']) * (m_l + m_h) +
                       Deltas(inputs['a_fm'])[Xi])
            XiSum += (2 * J1Sub(m_K, E) + J1Sub(m_S, E)) * tastemults[Xi]
        jSum = 0.
        for Xi in deltas(inputs['a_fm']).keys():
            m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l +
                        Deltas(inputs['a_fm'])[Xi])
            m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h +
                       Deltas(inputs['a_fm'])[Xi])
            m_pi2 = m_pi ** 2
            m_S2 = m_S ** 2
            Z = sqrt((m_S2 - m_pi2) ** 2 -
                     (deltas(inputs['a_fm'])[Xi] / 2) * (m_S2 - m_pi2) +
                     (9. / 16) * deltas(inputs['a_fm'])[Xi] ** 2)
            m_eta = sqrt((m_pi2 + m_S2 + (3. / 4) *
                          deltas(inputs['a_fm'])[Xi] - Z) / 2)
            m_etap = sqrt((m_pi2 + m_S2 + (3. / 4) *
                           deltas(inputs['a_fm'])[Xi] + Z) / 2)
            for j, m_j in [[1, m_S], [2, m_eta], [3, m_etap]]:
                jSum -= (deltas(inputs['a_fm'])[Xi] *
                         R31(m_pi, m_S, m_eta, m_etap, j) * J1Sub(m_j, E))
        m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])['I'])
        m_eta = sqrt(mu(inputs['a_fm']) * ((2 * m_l + 4 * m_h) / 3) +
                     Deltas(inputs['a_fm'])['I'])
        D = (XiSum / 16) + jSum + (2. / 3) * J1Sub(m_eta, E) - J1Sub(m_S, E)
        return (-3 * g_pi2 * E * D) / (4 * pi * fpi) ** 2


def df_para(inputs, g_pi):
    """Calculates loop corrections to f_para (including contributions from
    wavefunction renormalizations), according to given inputs, g_pi (the
    heavy-light coupling constant), and global constants."""
    E = inputs['E']
    m_l = inputs['ml_val']
    m_h = inputs['mh_val']
    g_pi2 = g_pi ** 2
    XiSum = 0.
    for Xi in tastemults.keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])[Xi])
        m_K = sqrt(mu(inputs['a_fm']) * (m_l + m_h) +
                   Deltas(inputs['a_fm'])[Xi])
        XiSum -= (3 * g_pi2 * I1(m_pi))* tastemults[Xi]
        if SU3:
            XiSum += ((((2 - 3 * g_pi2) / 2) * I1(m_K) + (I1(m_S) / 2) +
                       2 * I2(m_K, E) + I2(m_S, E)) * tastemults[Xi])
    jSum = 0.
    for Xi in deltas(inputs['a_fm']).keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])[Xi])
        m_pi2 = m_pi ** 2
        m_S2 = m_S ** 2
        Z = sqrt((m_S2 - m_pi2) ** 2 -
                 (deltas(inputs['a_fm'])[Xi] / 2) * (m_S2 - m_pi2) +
                 (9. / 16) * deltas(inputs['a_fm'])[Xi] ** 2)
        m_eta = sqrt((m_pi2 + m_S2 + (3. / 4) *
                      deltas(inputs['a_fm'])[Xi] - Z) / 2)
        m_etap = sqrt((m_pi2 + m_S2 + (3. / 4) *
                       deltas(inputs['a_fm'])[Xi] + Z) / 2)
        if SU3:
            for j, m_j in [[1, m_pi], [2, m_eta], [3, m_etap]]:
                jSum += (deltas(inputs['a_fm'])[Xi] *
                         R31(m_S, m_pi, m_eta, m_etap, j) *
                         (3. / 2) * g_pi2 * I1(m_j))
            for j, m_j in [[1, m_S], [2, m_eta], [3, m_etap]]:
                jSum -= (deltas(inputs['a_fm'])[Xi] *
                         R31(m_pi, m_S, m_eta, m_etap, j) *
                         ((I1(m_j) / 2) + I2(m_j, E)))
            jSum += ((deltas(inputs['a_fm'])[Xi] / (m_etap ** 2 - m_eta ** 2)) *
                     (I1(m_etap) - I1(m_eta) + I2(m_etap, E) - I2(m_eta, E)))
        else:
            jSum += 3 * g_pi2 * (I1(m_pi) - I1(m_eta))
    m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])['I'])
    m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])['I'])
    m_eta = sqrt(mu(inputs['a_fm']) * ((2 * m_l + 4 * m_h) / 3) +
                 Deltas(inputs['a_fm'])['I'])
    df = (XiSum / 16) + jSum + (3. / 4) * g_pi2 * I1(m_pi)
    if SU3:
        df += (((8 - 3 * g_pi2) / 12) * I1(m_eta) - (I1(m_S) / 2) +
               I2(m_eta, E) - I2(m_S, E))
    return df / (4 * pi * fpi) ** 2


def df_perp(inputs, g_pi):
    """Calculates loop corrections to f_perp (including contributions from
    wavefunction renormalizations), according to given inputs, g_pi (the
    heavy-light coupling constant), and global constants."""
    E = inputs['E']
    m_l = inputs['ml_val']
    m_h = inputs['mh_val']
    g_pi2 = g_pi ** 2
    XiSum = 0.
    for Xi in tastemults.keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])[Xi])
        m_K = sqrt(mu(inputs['a_fm']) * (m_l + m_h) +
                   Deltas(inputs['a_fm'])[Xi])
        XiSum -= ((3 * g_pi2 * I1(m_pi) +
                   SU3 * ((I1(m_S) / 2) + ((2 + 3 * g_pi2) / 2) * I1(m_K))) *
                  tastemults[Xi])
    jSum = 0.
    for Xi in deltas(inputs['a_fm']).keys():
        m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])[Xi])
        m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])[Xi])
        m_pi2 = m_pi ** 2
        m_S2 = m_S ** 2
        Z = sqrt((m_S2 - m_pi2) ** 2 -
                 (deltas(inputs['a_fm'])[Xi] / 2) * (m_S2 - m_pi2) +
                 (9. / 16) * deltas(inputs['a_fm'])[Xi] ** 2)
        m_eta = sqrt((m_pi2 + m_S2 + (3. / 4) *
                      deltas(inputs['a_fm'])[Xi] - Z) / 2)
        m_etap = sqrt((m_pi2 + m_S2 + (3. / 4) *
                       deltas(inputs['a_fm'])[Xi] + Z) / 2)
        if SU3:
            jSum += (deltas(inputs['a_fm'])[Xi] * (g_pi2 / (m_etap ** 2 - m_eta ** 2)) *
                     (J1Sub(m_eta, E) - J1Sub(m_etap, E)))
            for j, m_j in [[1, m_pi], [2, m_eta], [3, m_etap]]:
                jSum += (deltas(inputs['a_fm'])[Xi] *
                         R31(m_S, m_pi, m_eta, m_etap, j) * I1(m_j) * 3 *
                         (g_pi2 / 2))
            for j, m_j in [[1, m_S], [2, m_eta], [3, m_etap]]:
                jSum += (deltas(inputs['a_fm'])[Xi] *
                         R31(m_pi, m_S, m_eta, m_etap, j) * (I1(m_j) / 2))
        else:
            jSum += 3 * g_pi2 * (I1(m_pi) - I1(m_eta))
    m_pi = sqrt(mu(inputs['a_fm']) * 2 * m_l + Deltas(inputs['a_fm'])['I'])
    m_S = sqrt(mu(inputs['a_fm']) * 2 * m_h + Deltas(inputs['a_fm'])['I'])
    m_eta = sqrt(mu(inputs['a_fm']) * ((2 * m_l + 4 * m_h) / 3) +
                 Deltas(inputs['a_fm'])['I'])
    df = (XiSum / 16) + jSum + (3. / 4) * g_pi2 * I1(m_pi)
    if SU3:
        df += ((I1(m_S) / 2) - ((4 + 3 * g_pi2) / 12) * I1(m_eta) -
               (g_pi2 / 3) * J1Sub(m_eta, E))
    return df / (4 * pi * fpi) ** 2

