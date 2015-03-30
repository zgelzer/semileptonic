"""Defines chiral fit function for B decay form factors."""


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


from calculators.chilogs.fcns import D, df_para, df_perp, mu, Deltabar
from settings.constants import *
from settings.fit import *
from math import pi, sqrt
import numpy as np


def fitfcn(inputs, params, nouts=nexperiments):
    """Returns nouts chiral estimate(s) of pertinent form factor."""
    fpis = fpi * np.ones(nouts)
    Es = np.array([inputs[i]['E'] for i in range(nouts)])
    chi_l = np.array([(2 * mu(inputs[i]['a_fm']) * inputs[i]['ml_val']) /
                      (8 * pi ** 2 * fpi ** 2) for i in range(nouts)])
    chi_h = np.array([(2 * mu(inputs[i]['a_fm']) * inputs[i]['mh_val']) /
                      (8 * pi**2 * fpi ** 2) for i in range(nouts)])
    chi_a2 = np.array([Deltabar(inputs[i]['a_fm']) /
                       (8 * pi**2 * fpi ** 2) for i in range(nouts)])
    chi_E = (Es * sqrt(2)) / (4 * pi * fpi)
    if 'C0' in params.keys():
        C0 = params['C0']
    else:
        C0 = 0
    if 'Cl' in params.keys():
        Cl = params['Cl']
    else:
        Cl = 0
    if 'Cl2' in params.keys():
        Cl2 = params['Cl2']
    else:
        Cl2 = 0
    if 'CE' in params.keys():
        CE = params['CE']
    else:
        CE = 0
    if 'CE2' in params.keys():
        CE2 = params['CE2']
    else:
        CE2 = 0
    if 'Ca2' in params.keys():
        Ca2 = params['Ca2']
    else:
        Ca2 = 0
    if 'Cla2' in params.keys():
        Cla2 = params['Cla2']
    else:
        Cla2 = 0
    if 'ClE' in params.keys():
        ClE = params['ClE']
    else:
        ClE = 0
    if 'ClE2' in params.keys():
        ClE2 = params['ClE2']
    else:
        ClE2 = 0
    if 'Ca2E' in params.keys():
        Ca2E = params['Ca2E']
    else:
        Ca2E = 0
    if 'Ca2E2' in params.keys():
        Ca2E2 = params['Ca2E2']
    else:
        Ca2E2 = 0
    if 'CE3' in params.keys():
        CE3 = params['CE3']
    else:
        CE3 = 0
    if 'CE4' in params.keys():
        CE4 = params['CE4']
    else:
        CE4 = 0
    if 'Ca4' in params.keys():
        Ca4 = params['Ca4']
    else:
        Ca4 = 0
    if 'Ch' in params.keys():
        Ch = params['Ch']
    else:
        Ch = 0
    if 'gpi' in params.keys():
        g_pi = params['gpi']
    else:
        g_pi = gpi
    if formfactor == 'para':
        dfs = np.array([df_para(inputs[i], g_pi, decayname) for i in
                        range(nouts)])
        return ((C0 * (1. + dfs) +
                 Cl * chi_l +
                 Cl2 * chi_l ** 2 +
                 CE * chi_E +
                 CE2 * chi_E ** 2 +
                 Ca2 * chi_a2 +
                 Cla2 * chi_l * chi_a2 +
                 ClE * chi_l * chi_E +
                 ClE2 * chi_l * chi_E ** 2 +
                 Ca2E * chi_a2 * chi_E +
                 Ca2E2 * chi_a2 * chi_E ** 2 +
                 CE3 * chi_E ** 3 +
                 CE4 * chi_E ** 4 +
                 Ca4 * chi_a2 ** 2 +
                 Ch * chi_h) /
                fpis)
    elif formfactor == 'perp':
        dfs = np.array([df_perp(inputs[i], g_pi, decayname) for i in
                        range(nouts)])
        Ds = np.array([D(inputs[i], g_pi, decayname) for i in range(nouts)])
        Delta_Bs = Delta_B * np.ones(nouts)
        return ((C0 * (1. + dfs) +
                 Cl * chi_l +
                 Cl2 * chi_l ** 2 +
                 CE * chi_E +
                 CE2 * chi_E ** 2 +
                 Ca2 * chi_a2 +
                 Cla2 * chi_l * chi_a2 +
                 ClE * chi_l * chi_E +
                 ClE2 * chi_l * chi_E ** 2 +
                 Ca2E * chi_a2 * chi_E +
                 Ca2E2 * chi_a2 * chi_E ** 2 +
                 CE3 * chi_E ** 3 +
                 CE4 * chi_E ** 4 +
                 Ca4 * chi_a2 ** 2 +
                 Ch * chi_h) *
                (gpi / (fpis * (Es + Delta_Bs + Ds))))

