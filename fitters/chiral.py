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
Defines chiral fit function for decay form factors of B mesons.
--------------------------------------------------------------------------------
Definitions
-----------
fitfcn : function
    Returns chiral estimates of pertinent form factor for each experiment.
--------------------------------------------------------------------------------
"""


from calculators.chilogs.fcns import D, df_para, df_perp, mu, Deltabar
from settings.constants import Delta_B, fpi, gpi
from settings.fit import decayname, formfactor, nexperiments
from math import pi, sqrt
import numpy as np


def fitfcn(inputs, params, nouts=nexperiments):
    """
    ----------------------------------------------------------------------------
    Returns chiral estimates f(X, p) of pertinent form factor for all
    experiments in (X = inputs) given fit parameters (p = params).
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Sets any fit parameter (other than g_pi) to zero if not supplied in params;
    sets g_pi to settings.constants.gpi if g_pi not supplied in params. All
    (next-to-)next-to-leading-order (N)NLO fit parameters C^(i) (with (i) > 0)
    enter into chiral fit function as products of dimensionless expansion terms
    chi_(i). chi_E is proportional to E_(K/pi); chi_a2 is proportional to
    (a^2 * Deltabar = calculators.chilogs.fcns.Deltabar); chi_h is proportional
    to (mu * m_h), where m_h is heavy quark mass and mu is mass slope constant
    (see calculators.chilogs.fcns.mu); chi_l is proportional to (mu * m_l),
    where m_l is light quark mass. All chi_(i) have appropriate number of
    factors of f_pi in their denominators in order to remain dimensionless. NNLO
    expansion terms involve products of multiple chi_(i); e.g., ClE2
    parameterizes (chi_l * (chi_E)^2). NLO chiral fit function, similar in form
    to that in [1], is as follows:
        f_para = (1 / f_pi) * [C^(0) * (1 + df) + Sum_(i){C^(i) * chi_(i)}]
        f_perp = [g_pi / (f_pi * (E_(K/pi) + Delta_B + D))] *
                 [C^(0) * (1 + df) + Sum_(i){C^(i) * chi_(i)}]
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of {nouts} total input dictionaries, where each dictionary stores
        input floats for particular experiment. See fileIOs.readers.inputs for
        complete list of inputs.
    params : dict of floats or gvar.BufferDict
        Current fit parameters in dictionary-like container.
    nouts : int (optional; default is {number of experiments})
        Number of floats to return in output array; must match {size of inputs}.
    ----------------------------------------------------------------------------
    Returns
    -------
    numpy.ndarray of floats or of gvar.GVars
        Array of nouts total form factors.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    D : function, from calculators.chilogs.fcns
    Delta_B : float, from settings.constants
    Deltabar : function, from calculators.chilogs.fcns
    decayname : str, from settings.fit
    df_para : function, from calculators.chilogs.fcns
    df_perp : function, from calculators.chilogs.fcns
    formfactor : str, from settings.fit
    fpi : float, from settings.constants
    gpi : float, from settings.constants
    mu : function, from calculators.chilogs.fcns
    nexperiments : int, from settings.fit
    numpy : module, as np
    pi : function, from math
    sqrt : function, from math
    ----------------------------------------------------------------------------
    References
    ----------
    [1] J. Bailey et al. (Fermilab Lattice and MILC Collaborations), "The B -->
        pi l nu semileptonic form factor from three-flavor lattice QCD: A model-
        independent determination of |V(ub)|", Phys. Rev. D 79, 054507 (2009)
        [arXiv:0811.3640 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    Es = np.array([inputs[i]['E'] for i in range(nouts)])
    fpis = fpi * np.ones(nouts)
    chi_E = (Es * sqrt(2)) / (4 * pi * fpi)
    chi_a2 = np.array([Deltabar(inputs[i]['a_fm']) /
                       (8 * pi**2 * fpi ** 2) for i in range(nouts)])
    chi_h = np.array([(2 * mu(inputs[i]['a_fm']) * inputs[i]['mh_val']) /
                      (8 * pi**2 * fpi ** 2) for i in range(nouts)])
    chi_l = np.array([(2 * mu(inputs[i]['a_fm']) * inputs[i]['ml_val']) /
                      (8 * pi ** 2 * fpi ** 2) for i in range(nouts)])
    if 'C0' in params.keys():
        C0 = params['C0']
    else:
        C0 = 0
    if 'CE' in params.keys():
        CE = params['CE']
    else:
        CE = 0
    if 'CE2' in params.keys():
        CE2 = params['CE2']
    else:
        CE2 = 0
    if 'CE3' in params.keys():
        CE3 = params['CE3']
    else:
        CE3 = 0
    if 'CE4' in params.keys():
        CE4 = params['CE4']
    else:
        CE4 = 0
    if 'Ca2' in params.keys():
        Ca2 = params['Ca2']
    else:
        Ca2 = 0
    if 'Ca2E' in params.keys():
        Ca2E = params['Ca2E']
    else:
        Ca2E = 0
    if 'Ca2E2' in params.keys():
        Ca2E2 = params['Ca2E2']
    else:
        Ca2E2 = 0
    if 'Ca4' in params.keys():
        Ca4 = params['Ca4']
    else:
        Ca4 = 0
    if 'Ch' in params.keys():
        Ch = params['Ch']
    else:
        Ch = 0
    if 'Cl' in params.keys():
        Cl = params['Cl']
    else:
        Cl = 0
    if 'Cl2' in params.keys():
        Cl2 = params['Cl2']
    else:
        Cl2 = 0
    if 'ClE' in params.keys():
        ClE = params['ClE']
    else:
        ClE = 0
    if 'ClE2' in params.keys():
        ClE2 = params['ClE2']
    else:
        ClE2 = 0
    if 'Cla2' in params.keys():
        Cla2 = params['Cla2']
    else:
        Cla2 = 0
    if 'gpi' in params.keys():
        g_pi = params['gpi']
    else:
        g_pi = gpi
    if formfactor == 'para':
        dfs = np.array([df_para(inputs[i], g_pi, decayname) for i in
                        range(nouts)])
        return ((C0 * (1 + dfs) +
                 CE * chi_E +
                 CE2 * chi_E ** 2 +
                 CE3 * chi_E ** 3 +
                 CE4 * chi_E ** 4 +
                 Ca2 * chi_a2 +
                 Ca2E * chi_a2 * chi_E +
                 Ca2E2 * chi_a2 * chi_E ** 2 +
                 Ca4 * chi_a2 ** 2 +
                 Ch * chi_h +
                 Cl * chi_l +
                 Cl2 * chi_l ** 2 +
                 ClE * chi_l * chi_E +
                 ClE2 * chi_l * chi_E ** 2 +
                 Cla2 * chi_l * chi_a2) /
                fpis)
    elif formfactor == 'perp':
        Delta_Bs = Delta_B * np.ones(nouts)
        Ds = np.array([D(inputs[i], g_pi, decayname) for i in range(nouts)])
        dfs = np.array([df_perp(inputs[i], g_pi, decayname) for i in
                        range(nouts)])
        return ((C0 * (1 + dfs) +
                 CE * chi_E +
                 CE2 * chi_E ** 2 +
                 CE3 * chi_E ** 3 +
                 CE4 * chi_E ** 4 +
                 Ca2 * chi_a2 +
                 Ca2E * chi_a2 * chi_E +
                 Ca2E2 * chi_a2 * chi_E ** 2 +
                 Ca4 * chi_a2 ** 2 +
                 Ch * chi_h +
                 Cl * chi_l +
                 Cl2 * chi_l ** 2 +
                 ClE * chi_l * chi_E +
                 ClE2 * chi_l * chi_E ** 2 +
                 Cla2 * chi_l * chi_a2) *
                (gpi / (fpis * (Es + Delta_Bs + Ds))))

