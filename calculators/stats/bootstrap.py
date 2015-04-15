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
Defines functions for use with bootstrapping raw data and reporting averages
and errors for bootstrapped data.
--------------------------------------------------------------------------------
Definitions
-----------
avg : function
    Returns central value of bootstrapped data.
err : function
    Returns standard estimate of error of bootstrapped data.
middle : function
    Returns half of range of middle percent of bootstrapped data.
middle_bounds : function
    Returns lower, upper error bars of middle percent of bootstrapped data.
resample : function
    Applies bootstrap resamplings on data.
--------------------------------------------------------------------------------
"""


import numpy as np


def avg(bsdata):
    """
    ----------------------------------------------------------------------------
    Returns central value of bootstrapped data bsdata.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    bsdata : numpy.ndarray or list of floats
        Bootstrap samples of data.
    ----------------------------------------------------------------------------
    Returns
    -------
    float, from result of function numpy.mean
        Central value of bootstrap samples of data.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    ----------------------------------------------------------------------------
    """
    return np.mean(bsdata)


def err(bsdata):
    """
    ----------------------------------------------------------------------------
    Returns standard estimate of error of bootstrapped data bsdata.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    bsdata : numpy.ndarray or list of floats
        Bootstrap samples of data.
    ----------------------------------------------------------------------------
    Returns
    -------
    float, from result of function middle
        Standard estimate of error of bootstrap samples of data.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    middle : function
    ----------------------------------------------------------------------------
    """
    #return np.std(bsdata, ddof=1)
    return middle(bsdata)


def middle(bsdata, percent=68):
    """
    ----------------------------------------------------------------------------
    Returns half of range of middle percent of bootstrapped data bsdata.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    bsdata : numpy.ndarray or list of floats
        Bootstrap samples of data.
    percent : int or float (optional; default is 68)
        Middle percentage of data to use for range.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Half of range of middle percent of bootstrap samples of data.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    middle_bounds : function
    ----------------------------------------------------------------------------
    """
    bounds = middle_bounds(bsdata, percent=percent)
    return (bounds[-1] - bounds[0]) / 2.


def middle_bounds(bsdata, percent=68):
    """
    ----------------------------------------------------------------------------
    Returns lower, upper error bars of middle percent of bootstrapped data.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    bsdata : numpy.ndarray or list of floats
        Bootstrap samples of data.
    percent : int or float (optional; default is 68)
        Middle percentage of data to use for range.
    ----------------------------------------------------------------------------
    Returns
    -------
    float
        Lower error bar (< 0) of middle percent of bootstrap samples of data.
    float
        Upper error bar (> 0) of middle percent of bootstrap samples of data.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    avg : function
    ----------------------------------------------------------------------------
    """
    skip = len(bsdata) * (1. - percent / 100.) / 2.
    middle = np.sort(bsdata)[skip: -skip]
    mean = avg(bsdata)
    return middle[0] - mean, middle[-1] - mean


def resample(data, niter=None, seed=None):
    """
    ----------------------------------------------------------------------------
    Applies niter bootstrap resamplings on data with randomness seed.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    data : numpy.ndarray or list of floats
        Raw data.
    niter : int (optional; default is None)
        Number of bootstrap resampling iterations to perform; defaults to {size
        of data} if not specified.
    seed : int (optional; default is None)
        Randomness seed to use (via function numpy.random.seed) in generating
        bootstrap samples; skipped if not specified.
    ----------------------------------------------------------------------------
    Returns
    -------
    numpy.ndarray of floats
        Bootstrap resamples of data.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    ----------------------------------------------------------------------------
    """
    if niter is None:
        niter = len(data)
    if seed is not None:
        np.random.seed(seed)
    resamples = np.floor(np.random.rand(niter * len(data)) *
                         len(data)).astype(int)
    return data[resamples]

