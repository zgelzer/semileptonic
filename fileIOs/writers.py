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
Defines functions for writing and plotting data.
--------------------------------------------------------------------------------
Definitions
-----------
date_time : function
    Returns string of date/time to be appended to output file names.
params : function
    Pickles fit parameters.
plot_data : function
    Plots data vs. E_(K/pi) in r_1 units for all experiments.
plot_errfill : function
    Plots y vs. x, with shading about y of its error.
plot_fit : function
    Plots fit line vs. E_(K/pi) in r_1 units at some lattice spacing.
plot_fitavgs : function
    Plots fit lines vs. E_(K/pi) in r_1 units for all ensembles.
plot_labels : function
    Adds axis labels and legend to existing plot instance.
plot_save : function
    Saves plot instance.
results : function
    Writes results of fit settings, qualities, and parameters.
sigfig : function
    Rounds real numbers to some significant figure.
--------------------------------------------------------------------------------
"""


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


def date_time(use, sep='_'):
    """
    ----------------------------------------------------------------------------
    Returns string of date/time, to be appended with separation sep to output
    file names if use is not False.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    use : str or bool or NoneType
        Usage of date/time; may be string of date/time, or False, or None.
    sep : str (optional; default is '_')
        Separation to use before date/time.
    ----------------------------------------------------------------------------
    Returns
    -------
    str
        (sep + use) if use is string, {empty} if use is False, or (sep +
        {current date/time}) if use is None.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    datetime : class, from datetime, as dt
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError
        use must be one of: string, False, or None.
    ----------------------------------------------------------------------------
    """
    if use is None:
        return sep + dt.now().strftime('%Y%m%d-%H%M')
    elif use is False:
        return ''
    elif type(use) is str:
        return sep + use
    else:
        raise ValueError('use must be string, False, or None')


def params(avgs, cov, datetime=None):
    """
    ----------------------------------------------------------------------------
    Pickles fit parameters, where avgs are their central values and cov is their
    covariance matrix, saving result with date/time appended if datetime is not
    False.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    avgs : numpy.ndarray of floats
        Central values of results of fit parameters.
    cov : numpy.ndarray of floats
        Covariance matrix of results of fit parameters with shape ({number of
        parameters, number of parameters}).
    datetime : str or bool or NoneType (optional; default is None)
        Date/time to be appended to output file name.
    ----------------------------------------------------------------------------
    Results
    -------
    result.p : file
        Results of fit parameters in pickled binary format.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    date_time : function
    pickle : module
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function fitters.lsq.all.
    ----------------------------------------------------------------------------
    """
    pickle.dump((avgs, cov), open('result' + date_time(datetime) + '.p', 'wb'))


def plot_data(inputs, data, axis=None, colormap=None, datetime=None,
              markerfill='none', markershape='o'):
    """
    ----------------------------------------------------------------------------
    Plots (Y = data) vs. (X = inputs) for all experiments.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    All experiments of particular ensemble are given same color from colormap;
    these colors should match those in function plot_fitavgs. Labels are created
    only for first experiment of each ensemble, to avoid redundancy in legend
    (if created in future). Plot values and labels are ultimately saved via
    numpy.savetxt.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of {number of experiments} total dictionaries, where each
        dictionary stores input floats for particular experiment. See
        fileIOs.readers.inputs for complete list of inputs.
    data : numpy.ndarray of floats
        Array of particular form factor, with shape ({number of experiments},
        {number of bootstrap/jackknife samples}).
    axis : matplotlib.axes._subplots.AxesSubplot or NoneType (optional; default
           is None)
        Axis container for plot items.
    colormap : matplotlib.colors.LinearSegmentedColormap or NoneType (optional;
               default is None)
        Color map to use when creating successive plot line/marker colors.
    datetime : str or bool or NoneType (optional; default is None)
        Date/time to be appended to output file name.
    markerfill : str (optional; default is 'none')
        Type of fill to use for markers of data points.
    markershape : str (optional; default is 'o')
        Type of shape to use for markers of data points.
    ----------------------------------------------------------------------------
    Results
    -------
    result.dat : file
        Values from plot of input data with error bars.
    matplotlib.container.ErrorbarContainer
        Markers at (Y = data), (X = inputs) with error bars about Y.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    datatype : str, from settings.fit
    date_time : function
    ensemblesize : int, from settings.fit
    nexperiments : int, from settings.fit
    numpy : module, as np
    pyplot : module, from matplotlib, as plt
    sigfig : function
    ----------------------------------------------------------------------------
    """
    axis = axis if axis is not None else plt.gca()
    E2s = np.array([inputs[xpmt]['E']**2 for xpmt in range(nexperiments)])
    if datatype == 'bs':
        from calculators.stats.bootstrap import avg, middle_bounds as err
    data_avg = np.asarray([avg(data[xpmt, :]) for xpmt in range(nexperiments)])
    data_err = np.asarray([err(data[xpmt, :]) for xpmt in range(nexperiments)])
    data_errlower = np.asarray(-data_err[:, 0])
    data_errupper = np.asarray(data_err[:, -1])
    colormap = (colormap if colormap is not None else
                plt.get_cmap('gist_rainbow'))
    nensembles = int(nexperiments / ensemblesize)
    colors = [colormap(1. * ensemble / nensembles) for ensemble in
              range(nensembles)]
    aml_amhs = []
    for ensemble, xpmt in enumerate(range(0, nexperiments, ensemblesize)):
        aml = str(sigfig(inputs[xpmt]['a'] * inputs[xpmt]['ml_val'], n=3))
        amh = str(sigfig(inputs[xpmt]['a'] * inputs[xpmt]['mh_val'], n=3))
        label = r'$a m_\ell / a m_h$ = ' + aml + ' / ' + amh
        for ensxpmt in range(ensemblesize):
            plt.errorbar(E2s[xpmt + ensxpmt], data_avg[xpmt + ensxpmt],
                         yerr=([data_errlower[xpmt + ensxpmt]],
                               [data_errupper[xpmt + ensxpmt]]),
                         label=label, color=colors[ensemble], fmt=markershape,
                         fillstyle=markerfill)
            label = None
            aml_amhs.append(str(aml + ' / ' + amh).ljust(18))
    outputs = np.empty(nexperiments, dtype=[('aml_amh', 'a18'),
                                            ('E2', float),
                                            ('ff_avg', float),
                                            ('ff_err+', float),
                                            ('ff_err-', float)])
    outputs['aml_amh'] = aml_amhs
    outputs['E2'] = E2s
    outputs['ff_avg'] = data_avg
    outputs['ff_err+'] = data_errupper
    outputs['ff_err-'] = data_errlower
    np.savetxt('result' + date_time(datetime) + '.dat', outputs,
               fmt=('%18s', '%.6e', '%.6e', '%.6e', '%.6e'), delimiter='  ',
               header='  '.join(['a*ml / a*mh'.ljust(16), 'E2'.ljust(12),
                                 'ff_avg'.ljust(12), 'ff_err+'.ljust(12),
                                 'ff_err-'.ljust(12)]))


def plot_errfill(x, y, yerr, alphafill=0.3, axis=None, color=None, label=None):
    """
    ----------------------------------------------------------------------------
    Plots (y +- yerr) vs. x, with shading about y of the error yerr.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    x : numpy.ndarray
        Data for x values.
    y : numpy.ndarray
        Data for y values, with size {size of x}.
    yerr : numpy.ndarray or numpy.sctypeNA
        Data for y errors, as scalar, array with size 2, array with size {size
        of y}, or array with shape ({size of y}, 2).
    alphafill : float or NoneType (optional; default is 0.3)
        Alpha transparency of shading.
    axis : matplotlib.axes._subplots.AxesSubplot or NoneType (optional; default
           is None)
        Axis container for plot items.
    color : str or NoneType (optional; default is None)
        Color of plot line and shading.
    label : str or NoneType (optional; default is None)
        Label of plot line.
    ----------------------------------------------------------------------------
    Results
    -------
    matplotlib.lines.Line2D
        Plot line of y vs. x.
    matplotlib.collections.PolyCollection
        Plot shading of error bars yerr about y.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    pyplot : module, from matplotlib, as plt
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'invalid shape for yerr'
        yerr must be one of: scalar, array with size 2, array with size {size of
        y}, or array with shape ({size of y}, 2).
    ----------------------------------------------------------------------------
    """
    axis = axis if axis is not None else plt.gca()
    if color is None:
        color = axis._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin = y - yerr[0]
        ymax = y + yerr[1]
    else:
        raise ValueError('invalid shape for yerr')
    axis.plot(x, y, color=color, label=label)
    axis.fill_between(x, ymax, ymin, color=color, alpha=alphafill)


def plot_fit(lattspace, inputs, params, linelength, alphafill=0.3, axis=None,
             color=None, colormap=None, datetime=None, label=None):
    """
    ----------------------------------------------------------------------------
    Plots fit line with length linelength and lattice spacing name lattspace,
    given (X = inputs) and fit parameters params.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Creates fit inputs dictionary from last experiment of ensemble whose lattice
    spacing matches that of lattspace, or creates fit inputs dictionary from
    continuum constants if (lattspace = 0). Creates values for E_(K/pi) via
    numpy.linspace(linelength), where its bounds are made to match those of
    E_(K/pi) in inputs if (linelength = 'full'). Applies
    fitters.chiral.fitfcn(x, p) with (x = {fit inputs}) and (p = params), which
    returns an array of gvar.GVars with standard propagation of errors. Plot
    values and labels are ultimately saved via numpy.savetxt.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    lattspace : str
        Lattice spacing name.
    inputs : numpy.ndarray of dicts
        Array of {number of experiments} total dictionaries, where each
        dictionary stores input floats for particular experiment. See
        fileIOs.readers.inputs for complete list of inputs.
    params : dict of floats or gvar.BufferDict
        Fit parameters in dictionary-like container.
    linelength : list or str
        Length of fit line, as list of [min, max, numpoints] for use with
        numpy.linspace or as string 'full'. 'full' causes numpy.linspace bounds
        to match those of E_(K/pi) in inputs.
    alphafill : float or NoneType (optional; default is 0.3)
        Alpha transparency of shading.
    axis : matplotlib.axes._subplots.AxesSubplot or NoneType (optional; default
           is None)
        Axis container for plot items.
    color : str or NoneType (optional; default is None)
        Color of plot line.
    colormap : matplotlib.colors.LinearSegmentedColormap or NoneType (optional;
               default is None)
        Color map to use when creating successive plot line colors.
    datetime : str or bool or NoneType (optional; default is None)
        Date/time to be appended to output file name.
    label : str or NoneType (optional; default is None)
        Label of plot line.
    ----------------------------------------------------------------------------
    Results
    -------
    result_fit.dat : file
        Values from plot of fit line at some lattice spacing with errors.
    matplotlib.lines.Line2D
        Plot line of central values of fit.
    matplotlib.collections.PolyCollection
        Plot shading of error bars about fit line.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    a_fermi : function, from calculators.chilogs.fcns
    chiral : module, from fitters
    date_time : function
    nexperiments : int, from settings.fit
    numpy : module, as np
    plot_errfill : function
    pyplot : module, from matplotlib, as plt
    sigfig : function
    ----------------------------------------------------------------------------
    """
    axis = axis if axis is not None else plt.gca()
    fit_a_fm = a_fermi(lattspace)
    if fit_a_fm != 0:
        fit_xpmt = np.where([inputs[xpmt]['a_fm'] == fit_a_fm for xpmt in
                             range(nexperiments)])[0][-1]
        fit_inputs = np.array([{}])
        fit_inputs[0] = inputs[fit_xpmt]
    else:
        from settings.constants import r1_continuum, ml_continuum, mh_continuum
        fit_inputs = np.array([{'extrapolation': 'continuum'}])
        fit_inputs[0]['a_fm'] = 0
        fit_inputs[0]['a'] = 1 / r1_continuum
        fit_inputs[0]['ml_sea'] = ml_continuum
        fit_inputs[0]['ml_val'] = ml_continuum
        fit_inputs[0]['mh_sea'] = mh_continuum
        fit_inputs[0]['mh_val'] = mh_continuum
    E2s = np.array([inputs[xpmt]['E'] ** 2 for xpmt in range(nexperiments)])
    if linelength == 'full':
        fit_E2s = np.linspace(min(E2s), max(E2s), num=50)
    elif (type(linelength) == list) and (len(linelength) == 3):
        fit_E2s = np.linspace(linelength[0], linelength[1], num=linelength[2])
    fit_ffs_avg = np.empty_like(fit_E2s)
    fit_ffs_err = np.empty_like(fit_E2s)
    fit_ffs_gvar = []
    for E2 in fit_E2s:
        fit_inputs[0]['E'] = sqrt(E2)
        fit_ffs_gvar.append(chiral.fitfcn(fit_inputs, params, nouts=1)[0])
    fit_ffs_avg = np.array([fit_ffs_gvar[E2index].mean for E2index in
                            range(len(fit_E2s))])
    fit_ffs_err = np.array([fit_ffs_gvar[E2index].sdev for E2index in
                            range(len(fit_E2s))])
    plot_errfill(fit_E2s, fit_ffs_avg, fit_ffs_err, color=color, label=label,
                 alphafill=alphafill, axis=axis)
    aml = str(sigfig(fit_inputs[0]['a'] * fit_inputs[0]['ml_val'], n=3))
    amh = str(sigfig(fit_inputs[0]['a'] * fit_inputs[0]['mh_val'], n=3))
    np.savetxt('result_fit' + date_time(datetime) + '.dat',
               np.vstack((fit_E2s, fit_ffs_avg, fit_ffs_err)).T,
               fmt=('%.6e', '%.6e', '%.6e'), delimiter='  ',
               header='a*ml / a*mh = {}\n'.format(aml + ' / ' + amh) +
                      '  '.join(['E2'.ljust(10), 'ff_avg'.ljust(12),
                                 'ff_err'.ljust(12)]))


def plot_fitavgs(inputs, params, linelength, axis=None, colormap=None,
                 datetime=None, label=None, linewidth=1.0):
    """
    ----------------------------------------------------------------------------
    Plots fit lines with given length linelength for all ensembles, given
    (X = inputs) and fit parameters params.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Creates inputs dictionaries from last experiment of each ensemble. Creates
    values for E_(K/pi) via numpy.linspace(linelength), where its bounds are
    made to match those of E_(K/pi) in inputs if (linelength = 'full'). Applies
    fitters.chiral.fitfcn(x, p) with (x = {ensemble inputs}) and (p = {central
    values of params}) for each ensemble, which returns an array of floats. Plot
    values and labels are ultimately saved via numpy.savetxt.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    inputs : numpy.ndarray of dicts
        Array of {number of experiments} total dictionaries, where each
        dictionary stores input floats for particular experiment. See
        fileIOs.readers.inputs for complete list of inputs.
    params : dict of floats or gvar.BufferDict
        Fit parameters in dictionary-like container.
    linelength : list or str
        Length of fit line, as list of [min, max, numpoints] for use with
        numpy.linspace or as string 'full'. 'full' causes numpy.linspace bounds
        to match those of E_(K/pi) in inputs.
    axis : matplotlib.axes._subplots.AxesSubplot or NoneType (optional; default
           is None)
        Axis container for plot items.
    colormap : matplotlib.colors.LinearSegmentedColormap or NoneType (optional;
               default is None)
        Color map to use when creating successive plot line colors.
    datetime : str or bool or NoneType (optional; default is None)
        Date/time to be appended to output file name.
    label : str or NoneType (optional; default is None)
        Label of plot line.
    linewidth : float (optional; default is 1.0)
        Relative width of plot line.
    ----------------------------------------------------------------------------
    Results
    -------
    result_fits.dat : file
        Values from plot of ensemble fit averages.
    matplotlib.lines.Line2D
        Plot line of central values of fits.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    chiral : module, from fitters
    date_time : function
    ensemblesize : int, from settings.fit
    nexperiments : int, from settings.fit
    numpy : module, as np
    pyplot : module, from matplotlib, as plt
    ----------------------------------------------------------------------------
    """
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
    outputs = np.empty((nensembles, len(fit_E2s)))
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
        outputs[ensemble] = fit_ffs
    outputs = np.vstack((fit_E2s, outputs)).T
    np.savetxt('result_fits' + date_time(datetime) + '.dat', outputs,
               fmt='%.6e', delimiter='  ',
               header='  '.join(['E2'.ljust(10)] +
                                ['ff_avg(e={})'.format(ens).ljust(12)
                                 for ens in range(nensembles)]))


def plot_labels(legendloc='upper right', legendsize='8'):
    """
    ----------------------------------------------------------------------------
    Adds axis labels and legend to existing plot instance.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    legendloc : str or int or tuple of floats (optional; default is
                'upper right')
        Location of legend in plot figure.
    legendsize : str of int or of float (optional; default is '8')
        Size of fonts used in legend.
    ----------------------------------------------------------------------------
    Results
    -------
    matplotlib.legend.Legend
        Plot legend that will display all previously declared plot labels.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    decayname : str, from settings.fit
    formfactor : str, from settings.fit
    pyplot : module, from matplotlib, as plt
    ----------------------------------------------------------------------------
    """
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
    """
    ----------------------------------------------------------------------------
    Saves plot instance with date/time appended if datetime is not False.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    datetime : str or bool or NoneType (optional; default is None)
        Date/time to be appended to output file name.
    ----------------------------------------------------------------------------
    Results
    -------
    result.pdf : file
        Saves plot instance.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    date_time : function
    pyplot : module, from matplotlib, as plt
    ----------------------------------------------------------------------------
    Notes
    -----
    + plot_save should be called before plot instance is displayed via function
      matplotlib.pyplot.show.
    ----------------------------------------------------------------------------
    """
    plt.savefig('result' + date_time(datetime) + '.pdf')


def results(params, chi2, dof, Q, datetime=None, fileouts=True, stdouts=True):
    """
    ----------------------------------------------------------------------------
    Writes results of fit settings, qualities, and parameters to stdout and/or
    output file, saving latter with date/time appended if datetime is not False.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    params : dict of floats or gvar.BufferDict
        Fit parameters in dictionary-like container.
    chi2 : float
        Minimum chi^2 from fit.
    dof : float
        Degrees of freedom in fit.
    Q : float
        p-value of fit.
    datetime : str or bool or NoneType (optional; default is None)
        Date/time to be appended to output file name.
    fileouts : bool (optional; default is True)
        Determines if outputs are written to file.
    stdouts : bool (optional; default is True)
        Determines if outputs are written to stdout.
    ----------------------------------------------------------------------------
    Results
    -------
    result.txt : file (optional; created by default)
        Results of fit settings, qualities, and parameters.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    correlated : bool, from settings.fit
    date_time : function
    fit : module, from settings
    getcwd : function, from os
    numpy : module, as np
    sigfig : function
    stdout: file, from settings
    ----------------------------------------------------------------------------
    """
    if fileouts and stdouts:
        stdout.write('\nWriting results to {}\n'.format(getcwd()))
    output = '#\n# Settings:\n'
    for setting in sorted(fit.__all__):
        output += '#\t{0}\t{1}\n'.format(setting.rjust(12), eval(setting))
    if correlated:
        output += '#\n# Correlated'
    else:
        output += '#\n# Uncorrelated'
    output += ' Least Squares Fit:\n'
    output += '#\tchi2/dof [dof] = {0} [{1}]\n'.format(sigfig(chi2 / dof, n=2),
                                                     int(np.round(dof)))
    output += '#\tQ = {}\n'.format(sigfig(Q, n=2))
    output += '#\n# Parameters:\n'
    for param in sorted(params.keys()):
        output += '\t{0}\t{1}\n'.format(param.rjust(5), params[param])
    output += '\n'
    if stdouts:
        stdout.write(output)
    if fileouts:
        outfile = open('result' + date_time(datetime) + '.txt', 'w')
        outfile.write(output)
        outfile.close()


def sigfig(x, n=6):
    """
    ----------------------------------------------------------------------------
    Rounds any real number x to n significant figures.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    x : float or int
        Any real number.
    n : int (optional; default is 6)
        Number of significant figures to use.
    ----------------------------------------------------------------------------
    Returns
    -------
    float (or int if (x = 0))
    ----------------------------------------------------------------------------
    Requirements
    ------------
    copysign : function, from math
    floor : function, from math
    log10 : function, from math
    ----------------------------------------------------------------------------
    """
    if x == 0:
        return 0
    else:
        return copysign(1, x) * round(abs(x), n - int(floor(log10(abs(x)))) - 1)

