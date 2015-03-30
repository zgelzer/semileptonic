"""Performs least squares fit(s) to data using Prof. G. Peter LePage's
lsqfit.nonlinear_fit."""


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


from fitters.chiral import fitfcn
from gvar import gvar
from lsqfit import nonlinear_fit
from fileIOs.readers import array2dict
from settings.fit import *
import numpy as np
import pickle


def one(inputs, data, fitfcn=fitfcn, priors=None, p0s=None):
    """Performs a one-time fit to X=inputs Y=data, using the supplied fit
    function (default: fitters.chiral.fitfcn) and one of priors or p0s."""
    if datatype == 'bs':
        from calculators.stats.bootstrap import avg, err
    elif datatype == 'raw':
        from calculators.stats.raw import avg, err
    if correlated:
        cov = np.cov(data, ddof=1)
        data_info = np.zeros_like(cov)
        for i in range(0, nexperiments, ensemblesize):
            data_info[i : (i + ensemblesize), i : (i + ensemblesize)] = (
                cov[i : (i + ensemblesize), i : (i + ensemblesize)])
    else:
        data_info = np.asarray([err(data[xpmt, :]) for xpmt in
                                range(nexperiments)])
    data_cvals = np.asarray([avg(data[xpmt, :]) for xpmt in
                             range(nexperiments)])
    fit = nonlinear_fit(data=(inputs, data_cvals, data_info), fcn=fitfcn,
                        prior=priors, p0=p0s)
    return fit.p, fit.chi2, fit.dof, fit.Q


def all(inputs, data, fitfcn=fitfcn, priors=None, p0s=None):
    """Performs nsamples total fits to X=inputs Y=data for a total, using the
    supplied fit function (default: fitters.chiral.fitfcn) and one of priors or
    p0s. Note that nsamples = data.size / nexperiments requires nexperiments to
    have been defined."""
    params = sorted(priors.keys() if priors is not None else p0s.keys())
    if datatype == 'bs':
        from calculators.stats.bootstrap import avg, err
    elif datatype == 'raw':
        from calculators.stats.raw import avg, err
    if correlated:
        cov = np.cov(data, ddof=1)
        data_info = np.zeros_like(cov)
        for i in range(0, nexperiments, ensemblesize):
            data_info[i : (i + ensemblesize), i : (i + ensemblesize)] = (
                cov[i : (i + ensemblesize), i : (i + ensemblesize)])
    else:
        data_info = np.asarray([err(data[xpmt, :]) for xpmt in
                                range(nexperiments)])
    nsamples = data.size / nexperiments
    fit_avgs = {}
    for param in params:
        fit_avgs[param] = np.empty(nsamples)
    for sample in range(nsamples):
        data_samples = np.array([data[xpmt, sample] for xpmt in
                                 range(nexperiments)])
        fit = nonlinear_fit(data=(inputs, data_samples, data_info), fcn=fitfcn,
                            prior=priors, p0=p0s)
        for param in params:
            fit_avgs[param][sample] = fit.p[param].mean
    fit_cvals = np.asarray([avg(fit_avgs[param]) for param in params])
    fit_cov = np.cov(np.vstack([fit_avgs[param] for param in
                                params]), ddof=1)
    return array2dict(gvar(fit_cvals, fit_cov), params), fit_cvals, fit_cov

