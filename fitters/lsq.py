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
Performs least squares fit(s) to data using G. P. Lepage's lsqfit.nonlinear_fit.
--------------------------------------------------------------------------------
Definitions
-----------
one : function
    Performs one-time fit to inputs and data.
all : function
    Performs fits for all bootstrap/jackknife samples of inputs and data.
--------------------------------------------------------------------------------
"""


from fileIOs.readers import array2dict
from fitters.chiral import fitfcn
from gvar import gvar
from lsqfit import nonlinear_fit
from settings.fit import correlated, datatype, nenergies, nformfactors
import numpy as np


def one(inputs, data, fitfcn=fitfcn, inits=None, priors=None):
    """
    ----------------------------------------------------------------------------
    Performs one-time fit to ((X = inputs), (Y = data)), using given fit
    function fitfcn and one of inits or priors as fit parameters.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Passes (x, ymean, yinfo) as data for lsqfit.nonlinear_fit, with (x = X),
    (ymean = {central values of Y}) as determined by settings.fit.datatype, and
    (yinfo = {covariance matrix or errors of Y}) as determined by
    settings.fit.correlated. Correlated fits require covariance matrix with
    shape ({number of form factors}, {number of form factors}), where matrix
    elements involving form factors from different ensembles are set to zero.
    Uncorrelated fits require errors, as determined by settings.fit.datatype,
    with size {number of form factors}.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of {number of form factors} total dictionaries, where each
        dictionary stores input floats for particular form factor. See
        fileIOs.readers.inputs for complete list of inputs.
    data : numpy.ndarray of floats
        Array of particular form factor, with shape ({number of form factors},
        {number of bootstrap/jackknife samples}).
    fitfcn : function (optional; default is fitters.chiral.fitfcn)
        Fit function to be used by lsqfit.nonlinear_fit; must return values in
        same shape as that of nonlinear_fit's ymean. See fitters.chiral.fitfcn.
    priors : dict of gvar.GVars or NoneType (optional; default is None)
        A priori estimates of fit parameters, as Gaussian variables with widths,
        if specified; meant to be used with constrained fit.
    inits : dict of floats or NoneType (optional; default is None)
        Initial values of fit parameters, if specified; meant to be used with
        unconstrained fit.
    ----------------------------------------------------------------------------
    Returns
    -------
    fit.p : gvar.BufferDict
        Best-fit parameters from fit.
    fit.chi2 : float
        Minimum chi^2 from fit.
    fit.dof : float
        Degrees of freedom in fit; equates to size of data being fit if priors
        are specified, or to ({size of data being fit} - {number of fit
        parameters}) if no priors are specified.
    fit.Q : float
        p-value of fit, defined as probability that chi^2 from fit could have
        been larger assuming that best-fit model is correct.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    correlated : bool, from settings.fit
    datatype : str, from settings.fit
    nenergies : int, from settings.fit
    fitfcn : function, from fitters.chiral
    nformfactors : int, from settings.fit
    nonlinear_fit : class, from lsqfit
    numpy : module, as np
    ----------------------------------------------------------------------------
    Notes
    -----
    + (fit.chi2 / fit.dof) is of order unity for good fits; a value much greater
      than unity indicates a poor model, while a value much less than unity
      indicates overestimation of standard deviations in priors or data.
    + Good fits have p-values greater than about 0.1.
    ----------------------------------------------------------------------------
    """
    if datatype == 'bs':
        from calculators.stats.bootstrap import avg, err
    if correlated:
        cov = np.cov(data, ddof=1)
        data_info = np.zeros_like(cov)
        for i in range(0, nformfactors, nenergies):
            data_info[i : (i + nenergies), i : (i + nenergies)] = (
                  cov[i : (i + nenergies), i : (i + nenergies)])
    else:
        data_info = np.asarray([err(data[ff, :]) for ff in range(nformfactors)])
    data_cvals = np.asarray([avg(data[ff, :]) for ff in range(nformfactors)])
    fit = nonlinear_fit(data=(inputs, data_cvals, data_info), fcn=fitfcn,
                        prior=priors, p0=inits)
    return fit.p, fit.chi2, fit.dof, fit.Q


def all(inputs, data, fitfcn=fitfcn, inits=None, priors=None):
    """
    ----------------------------------------------------------------------------
    Performs {number of samples} total fits to ((X = inputs), (Y = data)), where
    ({number of samples} = {size of data} / {number of form factors}), using
    given fit function fitfcn and one of inits or priors as fit parameters.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    For each bootstrap/jackknife sample, passes (x, ymean, yinfo) as data for
    lsqfit.nonlinear_fit, with (x = X[sample]), (ymean = {central values of
    Y[sample]}) as determined by settings.fit.datatype, and (yinfo = {covariance
    matrix or errors of Y}) as determined by settings.fit.correlated. Correlated
    fits require full covariance matrix with shape ({number of form factors},
    {number of form factors}), where matrix elements involving form factors from
    different ensembles are set to zero. Uncorrelated fits require full errors,
    as determined by settings.fit.datatype, with size {number of form factors}.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of dictionaries with size {number of form factors}, where each
        dictionary stores input floats for particular form factor. See
        fileIOs.readers.inputs for complete list of inputs.
    data : numpy.ndarray
        Array of particular form factor, with shape ({number of form factors},
        {number of bootstrap/jackknife samples}).
    fitfcn : function (optional; default is fitters.chiral.fitfcn)
        Fit function to be used by lsqfit.nonlinear_fit; must return values in
        same shape as that of nonlinear_fit's ymean. See fitters.chiral.fitfcn.
    priors : dict of gvar.GVars or NoneType (optional; default is None)
        A priori estimates of fit parameters, as Gaussian variables with widths,
        if specified; meant to be used with constrained fit.
    inits : dict of floats or NoneType (optional; default is None)
        Initial values of fit parameters, if specified; meant to be used with
        unconstrained fit.
    ----------------------------------------------------------------------------
    Returns
    -------
    fit_cvals : numpy.ndarray of floats
        Central values of fit parameter results.
    fit_cov : numpy.ndarray of floats
        Covariance matrix of fit parameter results with shape ({number of
        parameters, number of parameters}).
    fit_gvars : dict of gvar.GVars
        Fit parameter results as Gaussian variables gvars with attached
        covariance matrix.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    array2dict : function, from fileIOs.readers
    correlated : bool, from settings.fit
    datatype : str, from settings.fit
    nenergies : int, from settings.fit
    fitfcn : function, from fitters.chiral
    gvar : class, from gvar
    nformfactors : int, from settings.fit
    nonlinear_fit : class, from lsqfit
    numpy : module, as np
    ----------------------------------------------------------------------------
    Notes
    -----
    + fit_cvals and fit_cov are returned so that they may be stored in pickled
      binary format for later exporting/importing of fit parameter results.
      Format of fit_gvars is not able to be pickled. See function
      fileIOs.writers.params.
    ----------------------------------------------------------------------------
    """
    params = sorted(priors.keys() if priors is not None else inits.keys())
    if datatype == 'bs':
        from calculators.stats.bootstrap import avg, err
    if correlated:
        cov = np.cov(data, ddof=1)
        data_info = np.zeros_like(cov)
        for i in range(0, nformfactors, nenergies):
            data_info[i : (i + nenergies), i : (i + nenergies)] = (
                  cov[i : (i + nenergies), i : (i + nenergies)])
    else:
        data_info = np.asarray([err(data[ff, :]) for ff in range(nformfactors)])
    nsamples = data.size / nformfactors
    fit_avgs = {}
    for param in params:
        fit_avgs[param] = np.empty(nsamples)
    for sample in range(nsamples):
        data_samples = np.array([data[ff, sample] for ff in
                                 range(nformfactors)])
        fit = nonlinear_fit(data=(inputs, data_samples, data_info), fcn=fitfcn,
                            prior=priors, p0=inits)
        for param in params:
            fit_avgs[param][sample] = fit.p[param].mean
    fit_cvals = np.asarray([avg(fit_avgs[param]) for param in params])
    fit_cov = np.cov(np.vstack([fit_avgs[param] for param in
                                params]), ddof=1)
    fit_gvars = array2dict(gvar(fit_cvals, fit_cov), params)
    return fit_gvars, fit_cvals, fit_cov

