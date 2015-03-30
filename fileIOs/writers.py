"""Defines functions for writing and plotting data."""


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


from calculators.chilogs.fcns import a_fermi
from datetime import datetime as dt
from fitters import chiral
from math import copysign, floor, log10, sqrt
from matplotlib import pyplot as plt
import numpy as np
from os import getcwd
from settings import fit
from settings.fit import *
from sys import stdout
import pickle


def params(avgs, info, datetime=None):
    """Pickles parameter data, where avgs are their means and info is either
    their standard deviations or their covariance matrix, saving result with
    date/time appended if datetime is not False."""
    if datetime is None:
        dt = '_' + dt.now().strftime('%Y%m%d-%H%M')
    elif datetime is False:
        dt = ''
    else:
        dt = '_' + datetime
    pickle.dump((avgs, info), open('result' + dt + '.p', 'wb'))


def results(params, chi2, dof, Q, stdouts=True, fileouts=True,
            datetime=None):
    """Writes fit settings, qualities, and parameters to stdout and/or output
    file, saving the latter with date/time appended if datetime is not False."""
    if fileouts:
        stdout.write('\nWriting results to {}\n'.format(getcwd()))
    output = '\nFitSettings:\n'
    for setting in sorted(fit.__all__):
        output += '\t{0}\t{1}\n'.format(setting.rjust(12), eval(setting))
    if correlated:
        output += '\nCorrelated'
    else:
        output += '\nUncorrelated'
    output += ' Least Squares Fit:\n'
    output += '\tchi2/dof [dof] = {0} [{1}]\n'.format(sigfig(chi2 / dof, n=2),
                                                     int(np.round(dof)))
    output += '\tQ = {}\n'.format(sigfig(Q, n=2))
    output += '\nParameters:\n'
    for param in sorted(params.keys()):
        output += '\t{0}\t{1}\n'.format(param.rjust(5), params[param])
    output += '\n'
    if stdouts:
        stdout.write(output)
    if fileouts:
        if datetime is None:
            dt = '_' + dt.now().strftime('%Y%m%d-%H%M')
        elif datetime is False:
            dt = ''
        else:
            dt = '_' + datetime
        outfile = open('result' + dt + '.txt', 'w')
        outfile.write(output)
        outfile.close()


def plot_data(inputs, data, markershape='o', markerfill='none', colormap=None,
              axis=None):
    """Plots data vs. inputs for all experiments."""
    axis = axis if axis is not None else plt.gca()
    E2s = np.array([inputs[xpmt]['E']**2 for xpmt in range(nexperiments)])
    if datatype == 'bs':
        from calculators.stats.bootstrap import avg, err
    elif datatype == 'raw':
        from calculators.stats.raw import avg, err
    data_avgs = np.asarray([avg(data[xpmt, :]) for xpmt in range(nexperiments)])
    data_errs = np.asarray([err(data[xpmt, :]) for xpmt in range(nexperiments)])
    colormap = (colormap if colormap is not None else
                plt.get_cmap('gist_rainbow'))
    nensembles = int(nexperiments / ensemblesize)
    colors = [colormap(1. * ensemble / nensembles) for ensemble in
              range(nensembles)]
    for ensemble, xpmt in enumerate(range(0, nexperiments, ensemblesize)):
        label = (r'$a m_\ell / a m_h$ = ' +
                 str(sigfig(inputs[xpmt]['a'] * inputs[xpmt]['ml_val'])) +
                 ' / ' +
                 str(sigfig(inputs[xpmt]['a'] * inputs[xpmt]['mh_val'])))
        for ensxpmt in range(ensemblesize):
            plt.errorbar(E2s[xpmt + ensxpmt], data_avgs[xpmt + ensxpmt],
                         yerr=data_errs[xpmt + ensxpmt], label=label,
                         color=colors[ensemble], fmt=markershape,
                         fillstyle=markerfill)
            label = None


def plot_errfill(x, y, yerr, color=None, label=None, alphafill=0.3, axis=None):
    """Plots (y +- yerr) vs. x, with shading about y of the error yerr."""
    axis = axis if axis is not None else plt.gca()
    if color is None:
        color = axis._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = yerr
    else:
        raise ValueError('invalid shape for yerr')
    axis.plot(x, y, color=color, label=label)
    axis.fill_between(x, ymax, ymin, color=color, alpha=alphafill)


def plot_fit(lattspace, inputs, params, linelength, color=None, label=None,
             linewidth=1.0, colormap=None, axis=None):
    """Plots fit line with given length linelength and lattice spacing string
    lattspace, according to given inputs and fit parameter results params."""
    axis = axis if axis is not None else plt.gca()
    fit_a_fm = a_fermi(lattspace)
    if fit_a_fm != 0:
        fit_xpmt = np.where([inputs[xpmt]['a_fm'] == fit_a_fm for xpmt in
                             range(nexperiments)])[0][-1]
        fit_inputs = inputs[fit_xpmt]
    else:
        from settings.constants import ml_continuum, mh_continuum
        fit_inputs = np.array([{'extrapolation': 'continuum'}])
        for param in inputs[0].keys():
            fit_inputs[0][param] = inputs[0][param]
        fit_inputs[0]['a_fm'] = 0
        fit_inputs[0]['ml_sea'] = ml_continuum
        fit_inputs[0]['ml_val'] = ml_continuum
        fit_inputs[0]['mh_sea'] = mh_continuum
        fit_inputs[0]['mh_val'] = mh_continuum
    E2s = np.array([inputs[xpmt]['E'] ** 2 for xpmt in range(nexperiments)])
    if linelength == 'full':
        fit_E2s = np.linspace(min(E2s), max(E2s), num=50)
    elif (type(linelength) == list) and (len(linelength) == 3):
        fit_E2s = np.linspace(linelength[0], linelength[1], num=linelength[2])
    fit_ff_avgs = np.empty_like(fit_E2s)
    fit_ff_errs = np.empty_like(fit_E2s)
    fit_ff_gvars = []
    for E2 in fit_E2s:
        fit_inputs[0]['E'] = sqrt(E2)
        fit_ff_gvars.append(chiral.fitfcn(fit_inputs, params, nouts=1)[0])
    fit_ff_avgs = np.array([fit_ff_gvars[E2index].mean for E2index in
                            range(len(fit_E2s))])
    fit_ff_errs = np.array([fit_ff_gvars[E2index].sdev for E2index in
                            range(len(fit_E2s))])
    plot_errfill(fit_E2s, fit_ff_avgs, fit_ff_errs, color=color, label=label,
                 alphafill=None, axis=axis)


def plot_fitavgs(inputs, params, linelength, color=None, label=None,
                 linewidth=1.0, colormap=None, axis=None):
    """Plots fit lines with given length linelength for all experiments,
    according to given inputs and fit parameter results params."""
    axis = axis if axis is not None else plt.gca()
    colormap = (colormap if colormap is not None else
                plt.get_cmap('gist_rainbow'))
    nensembles = int(nexperiments / ensemblesize)
    colors = [colormap(1. * ensemble / nensembles) for ensemble in
              range(nensembles)]
    E2s = np.array([inputs[xpmt]['E'] ** 2 for xpmt in range(nexperiments)])
    if linelength == 'full':
        fit_E2s = np.linspace(min(E2s), max(E2s), num=50)
    elif (type(linelength) == list) and (len(linelength) == 3):
        fit_E2s = np.linspace(linelength[0], linelength[1], num=linelength[2])
    params_avg = {}
    for param in params.keys():
        params_avg[param] = params[param].mean
    for ensemble in range(nensembles):
        fit_ffs = np.empty_like(fit_E2s)
        fit_inputs = np.array([{'ensemble': ensemble}])
        for param in inputs[ensemble * ensemblesize].keys():
            fit_inputs[0][param] = inputs[ensemble * ensemblesize][param]
        for E2index, E2 in enumerate(fit_E2s):
            fit_inputs[0]['E'] = sqrt(E2)
            fit_ffs[E2index] = chiral.fitfcn(fit_inputs, params_avg, nouts=1)[0]
        axis.plot(fit_E2s, fit_ffs, color=colors[ensemble], label=label,
                  linewidth=linewidth)


def plot_labels(legendloc='upper right', legendsize='8'):
    """Adds axis labels and legend to existing plot."""
    if decayname == 'B2K':
        plt.xlabel(r'$\left( r_1 E_K \right) ^2$')
    elif decayname == 'B2pi':
        plt.xlabel(r'$\left( r_1 E_\pi \right) ^2$')
    if formfactor == 'para':
        plt.ylabel(r'$r_1^{1 / 2} f_\parallel$')
    elif formfactor == 'perp':
        plt.ylabel(r'$r_1^{-1 / 2} f_\perp$')
    plt.legend(loc=legendloc, prop={'size': legendsize}, numpoints=1)


def plot_save(datetime=None):
    """Saves plot with date/time appended if datetime is not False."""
    if datetime is None:
        dt = '_' + dt.now().strftime('%Y%m%d-%H%M')
    elif datetime is False:
        dt = ''
    else:
        dt = '_' + datetime
    plt.savefig('result_' + decayname + '_' + formfactor + dt + '.pdf')


def sigfig(x, n=6):
    """Rounds any real number x to n significant figures."""
    if x == 0:
        return 0
    else:
        return copysign(1, x) * round(abs(x), n - int(floor(log10(abs(x)))) - 1)

