"""Defines functions for use with bootstrapping raw data and reporting averages
and errors for bootstrapped data."""


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


import numpy as np


def avg(bsdata):
    """Returns central value of bootstrapped data bsdata."""
    return np.mean(bsdata)


def err(bsdata):
    """Returns standard estimate of error of bootstrapped data bsdata."""
    #return middle(bsdata)
    return np.std(bsdata, ddof=1)


def middle(bsdata, percent=68):
    """Returns range of middle percent of bootstrapped data bsdata."""
    skip = len(bsdata) * (1 - percent / 100.) / 2
    middle = np.sort(bsdata)[skip: -skip]
    return (middle[-1] - middle[0]) / 2.


def resample(data, niter=None):
    """Applies niter bootstrap resamplings on data."""
    if niter is None:
        niter = len(data)
    resamples = np.floor(np.random.rand(niter * len(data)) *
                         len(data)).astype(int)
    return data[resamples]

