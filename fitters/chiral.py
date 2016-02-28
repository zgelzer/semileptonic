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
Defines chiral fit functions for decay form factors of B mesons.
--------------------------------------------------------------------------------
Definitions
-----------
f_para : function
    Returns chiral estimates of parallel form factor for each experiment.
f_perp : function
    Returns chiral estimates of perpendicular form factor for each experiment.
f_scalar : function
    Returns chiral estimates of scalar form factor for each experiment.
f_vector : function
    Returns chiral estimates of vector form factor for each experiment.
fitfcn : function
    Returns chiral estimates of pertinent form factor for each experiment.
fitparse : function
    Sets absent fit parameters to zero (or to constant value, where applicable).
--------------------------------------------------------------------------------
"""


from calculators.chilogs.fcns import D, df_para, df_perp, mu, Deltabar
from math import pi, sqrt
from settings.constants import DeltaB_para, DeltaB_perp, fpi, gpi
from settings.fit import decayname, formfactor
import numpy as np


def f_para(inputs, params):
    """
    ----------------------------------------------------------------------------
    Returns chiral estimates f(X, p) of parallel form factor for all experiments
    in (X = inputs) given fit parameters (p = params).
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Sets absent fit parameters to zero via fitparse(params) (see fitparse). All
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
        f_para = norms * [C^(0) * (1 + df) + Sum_(i){C^(i) * chi_(i)}], where
        norms = (1 / f_pi) if DeltaB_para is 0, or
        norms = gpi / (fpi * (Es + DeltaB_para)) otherwise
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of input dictionaries, where each dictionary stores input floats
        for particular experiment. See fileIOs.readers.inputs for complete list
        of inputs.
    params : dict of floats or gvar.BufferDict
        Current fit parameters in dictionary-like container.
    ----------------------------------------------------------------------------
    Returns
    -------
    numpy.ndarray of floats or of gvar.GVars
        Array of form factors with size {size of inputs}.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    DeltaB_para : float, from settings.constants
    Deltabar : function, from calculators.chilogs.fcns
    df_para : function, from calculators.chilogs.fcns
    fitparse : function
    fpi : float, from settings.constants
    mu : function, from calculators.chilogs.fcns
    numpy : module, as np
    pi : float, from math
    sqrt : function, from math
    ----------------------------------------------------------------------------
    References
    ----------
    [1] J. Bailey, et al. (Fermilab Lattice and MILC Collaborations), "The B -->
        pi l nu semileptonic form factor from three-flavor lattice QCD: A model-
        independent determination of |V(ub)|", Phys. Rev. D 79, 054507 (2009)
        [arXiv:0811.3640 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    fitparams = fitparse(params)
    nouts = len(inputs)
    Es = np.array([inputs[i]['E'] for i in range(nouts)])
    chi_Es = (Es * sqrt(2)) / (4 * pi * fpi)
    chi_a2s = np.array([Deltabar(inputs[i]['a_fm']) /
                        (8 * pi ** 2 * fpi ** 2) for i in range(nouts)])
    chi_hs = np.array([(2 * mu(inputs[i]['a_fm']) * inputs[i]['mh_val']) /
                       (8 * pi ** 2 * fpi ** 2) for i in range(nouts)])
    chi_ls = np.array([(2 * mu(inputs[i]['a_fm']) * inputs[i]['ml_val']) /
                       (8 * pi ** 2 * fpi ** 2) for i in range(nouts)])
    dfs = np.array([df_para(inputs[i], fitparams['gpi']) for i in range(nouts)])
    if DeltaB_para == 0:
        norms = 1. / fpi
    else:
        norms = fitparams['gpi'] / (fpi * (Es + DeltaB_para))
    return ((fitparams['C0'] * (1 + dfs) +
             fitparams['CE'] * chi_Es +
             fitparams['CE2'] * chi_Es ** 2 +
             fitparams['CE3'] * chi_Es ** 3 +
             fitparams['CE4'] * chi_Es ** 4 +
             fitparams['Ca2'] * chi_a2s +
             fitparams['Ca2E'] * chi_a2s * chi_Es +
             fitparams['Ca2E2'] * chi_a2s * chi_Es ** 2 +
             fitparams['Ca4'] * chi_a2s ** 2 +
             fitparams['Ch'] * chi_hs +
             fitparams['Cl'] * chi_ls +
             fitparams['Cl2'] * chi_ls ** 2 +
             fitparams['ClE'] * chi_ls * chi_Es +
             fitparams['ClE2'] * chi_ls * chi_Es ** 2 +
             fitparams['Cla2'] * chi_ls * chi_a2s) *
            norms)


def f_perp(inputs, params):
    """
    ----------------------------------------------------------------------------
    Returns chiral estimates f(X, p) of perpendicular (or tensor) form factor
    for all experiments in (X = inputs) given fit parameters (p = params).
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Sets absent fit parameters to zero via fitparse(params) (see fitparse). All
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
        f_perp = [g_pi / (f_pi * (E_(K/pi) + DeltaB_perp + D))] *
                 [C^(0) * (1 + df) + Sum_(i){C^(i) * chi_(i)}]
        f_tensor = (same form as f_perp)
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of input dictionaries, where each dictionary stores input floats
        for particular experiment. See fileIOs.readers.inputs for complete list
        of inputs.
    params : dict of floats or gvar.BufferDict
        Current fit parameters in dictionary-like container.
    ----------------------------------------------------------------------------
    Returns
    -------
    numpy.ndarray of floats or of gvar.GVars
        Array of form factors with size {size of inputs}.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    D : function, from calculators.chilogs.fcns
    DeltaB_perp : float, from settings.constants
    Deltabar : function, from calculators.chilogs.fcns
    df_perp : function, from calculators.chilogs.fcns
    fitparse : function
    fpi : float, from settings.constants
    mu : function, from calculators.chilogs.fcns
    numpy : module, as np
    pi : float, from math
    sqrt : function, from math
    ----------------------------------------------------------------------------
    References
    ----------
    [1] J. Bailey, et al. (Fermilab Lattice and MILC Collaborations), "The B -->
        pi l nu semileptonic form factor from three-flavor lattice QCD: A model-
        independent determination of |V(ub)|", Phys. Rev. D 79, 054507 (2009)
        [arXiv:0811.3640 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    fitparams = fitparse(params)
    nouts = len(inputs)
    Es = np.array([inputs[i]['E'] for i in range(nouts)])
    chi_Es = (Es * sqrt(2)) / (4 * pi * fpi)
    chi_a2s = np.array([Deltabar(inputs[i]['a_fm']) /
                        (8 * pi ** 2 * fpi ** 2) for i in range(nouts)])
    chi_hs = np.array([(2 * mu(inputs[i]['a_fm']) * inputs[i]['mh_val']) /
                       (8 * pi ** 2 * fpi ** 2) for i in range(nouts)])
    chi_ls = np.array([(2 * mu(inputs[i]['a_fm']) * inputs[i]['ml_val']) /
                       (8 * pi ** 2 * fpi ** 2) for i in range(nouts)])
    Ds = np.array([D(inputs[i], fitparams['gpi']) for i in range(nouts)])
    dfs = np.array([df_perp(inputs[i], fitparams['gpi']) for i in range(nouts)])
    return ((fitparams['C0'] * (1 + dfs) +
             fitparams['CE'] * chi_Es +
             fitparams['CE2'] * chi_Es ** 2 +
             fitparams['CE3'] * chi_Es ** 3 +
             fitparams['CE4'] * chi_Es ** 4 +
             fitparams['Ca2'] * chi_a2s +
             fitparams['Ca2E'] * chi_a2s * chi_Es +
             fitparams['Ca2E2'] * chi_a2s * chi_Es ** 2 +
             fitparams['Ca4'] * chi_a2s ** 2 +
             fitparams['Ch'] * chi_hs +
             fitparams['Cl'] * chi_ls +
             fitparams['Cl2'] * chi_ls ** 2 +
             fitparams['ClE'] * chi_ls * chi_Es +
             fitparams['ClE2'] * chi_ls * chi_Es ** 2 +
             fitparams['Cla2'] * chi_ls * chi_a2s) *
            (fitparams['gpi'] / (fpi * (Es + DeltaB_perp + Ds))))


def f_scalar(inputs, params_para, params_perp):
    """
    ----------------------------------------------------------------------------
    Returns chiral estimates f(X, P) of scalar form factor for all experiments
    in (X = inputs) given fit parameters (P = params) for both parallel and
    perpendicular form factors.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Imports masses of pertinent incoming and outgoing mesons (M_in and M_out,
    resp.) for particular decay. Scalar form factor is as follows [1]:
        f_0(q^2) = [sqrt(2 * M_in) / (M_in^2 - M_out^2)] *
                   [(M_in - E_(K/pi)) * f_para(E_(K/pi)) +
                    (E_(K/pi)^2 - M_out^2) * f_perp(E_(K/pi))]
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of input dictionaries, where each dictionary stores input floats
        for particular experiment. See fileIOs.readers.inputs for complete list
        of inputs.
    params_para : dict of floats or gvar.BufferDict
        Current parallel fit parameters in dictionary-like container.
    params_perp : dict of floats or gvar.BufferDict
        Current perpendicular fit parameters in dictionary-like container.
    ----------------------------------------------------------------------------
    Returns
    -------
    numpy.ndarray of floats or of gvar.GVars
        Array of form factors with size {size of inputs}.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    decayname : str, from settings.fit
    f_para : function
    f_perp : function
    numpy : module, as np
    sqrt : function, from math
    ----------------------------------------------------------------------------
    Notes
    -----
    + See functions f_para, f_perp.
    ----------------------------------------------------------------------------
    References
    ----------
    [1] J. Bailey, et al. (Fermilab Lattice and MILC Collaborations), "The B -->
        pi l nu semileptonic form factor from three-flavor lattice QCD: A model-
        independent determination of |V(ub)|", Phys. Rev. D 79, 054507 (2009)
        [arXiv:0811.3640 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    if decayname == 'B2K':
        from settings.constants import mB as M_in, mK as M_out
    elif decayname == 'B2pi':
        from settings.constants import mB as M_in, mpi as M_out
    Es = np.array([inputs[i]['E'] for i in range(len(inputs))])
    return (((M_in - Es) * f_para(inputs, params_para) +
             (np.square(Es) - M_out ** 2) * f_perp(inputs, params_perp)) *
            sqrt(2 * M_in) / (M_in ** 2 - M_out ** 2))


def f_vector(inputs, params_para, params_perp):
    """
    ----------------------------------------------------------------------------
    Returns chiral estimates f(X, P) of vector form factor for all experiments
    in (X = inputs) given fit parameters (P = params) for both parallel and
    perpendicular form factors.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Imports masses of pertinent incoming and outgoing mesons (M_in and M_out,
    resp.) for particular decay. Vector form factor is as follows [1]:
        f_+(q^2) = [1 / sqrt(2 * M_in)] *
                   [f_para(E_(K/pi)) + (M_in - E_(K/pi)) * f_perp(E_(K/pi))]
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of input dictionaries, where each dictionary stores input floats
        for particular experiment. See fileIOs.readers.inputs for complete list
        of inputs.
    params_para : dict of floats or gvar.BufferDict
        Current parallel fit parameters in dictionary-like container.
    params_perp : dict of floats or gvar.BufferDict
        Current perpendicular fit parameters in dictionary-like container.
    ----------------------------------------------------------------------------
    Returns
    -------
    numpy.ndarray of floats or of gvar.GVars
        Array of form factors with size {size of inputs}.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    decayname : str, from settings.fit
    f_para : function
    f_perp : function
    numpy : module, as np
    sqrt : function, from math
    ----------------------------------------------------------------------------
    Notes
    -----
    + See functions f_para, f_perp.
    ----------------------------------------------------------------------------
    References
    ----------
    [1] J. Bailey, et al. (Fermilab Lattice and MILC Collaborations), "The B -->
        pi l nu semileptonic form factor from three-flavor lattice QCD: A model-
        independent determination of |V(ub)|", Phys. Rev. D 79, 054507 (2009)
        [arXiv:0811.3640 [hep-lat]].
    ----------------------------------------------------------------------------
    """
    if (decayname == 'B2K') or (decayname == 'B2pi'):
        from settings.constants import mB as M_in
    Es = np.array([inputs[i]['E'] for i in range(len(inputs))])
    return ((f_para(inputs, params_para) +
            (M_in - Es) * f_perp(inputs, params_perp)) /
            sqrt(2 * M_in))


def fitfcn(*args):
    """
    ----------------------------------------------------------------------------
    Returns chiral estimates f_{form factor}(X, P) of particular form factor for
    all experiments in (X = inputs) given one or more sets of fit parameters
    (P = params).
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Positional arguments *args are {inputs, params} for single (parallel,
    perpendicular, or tensor) form factors; *args are {inputs, params_para,
    params_perp} for combined (scalar or vector) form factors.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of input dictionaries, where each dictionary stores input floats
        for particular experiment. See fileIOs.readers.inputs for complete list
        of inputs.
    params : dict of floats or gvar.BufferDict
        Current fit parameters in dictionary-like container.
    params_para : dict of floats or gvar.BufferDict
        Current parallel fit parameters in dictionary-like container.
    params_perp : dict of floats or gvar.BufferDict
        Current perpendicular fit parameters in dictionary-like container.
    ----------------------------------------------------------------------------
    Returns
    -------
    numpy.ndarray of floats or of gvar.GVars
        Array of form factors with size {size of inputs}.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    f_para : function
    f_perp : function
    f_scalar : function
    f_vector : function
    formfactor : str, from settings.fit
    ----------------------------------------------------------------------------
    Notes
    -----
    + See functions f_para, f_perp, f_scalar, f_vector.
    ----------------------------------------------------------------------------
    """
    if formfactor == 'para':
        return f_para(*args)
    elif (formfactor == 'perp') or (formfactor == 'tensor'):
        return f_perp(*args)
    elif formfactor == 'scalar':
        return f_scalar(*args)
    elif formfactor == 'vector':
        return f_vector(*args)


def fitparse(params):
    """
    ----------------------------------------------------------------------------
    Sets any fit parameter (other than g_pi) to zero if not supplied in params;
    sets g_pi to settings.constants.gpi if not supplied in params.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    params : dict of floats or gvar.BufferDict
        Current fit parameters in dictionary-like container.
    ----------------------------------------------------------------------------
    Returns
    -------
    fitparams : dict of floats or gvar.BufferDict
        Parsed current fit parameters in dictionary-like container.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    gpi : float, from settings.constants
    ----------------------------------------------------------------------------
    """
    paramlist = ['C0', 'CE', 'CE2', 'CE3', 'CE4', 'Ca2', 'Ca2E', 'Ca2E2', 'Ca4',
                 'Ch', 'Cl', 'Cl2', 'ClE', 'ClE2', 'Cla2']
    fitparams = {}
    for param in paramlist:
        if param in params.keys():
            fitparams[param] = params[param]
        else:
            fitparams[param] = 0
    if 'gpi' in params.keys():
        fitparams['gpi'] = params['gpi']
    else:
        fitparams['gpi'] = gpi
    return fitparams

