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
inputs_continuum : function
    Returns inputs for continuum extrapolations.
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
plot_fitcombo : function
    Plots continuum extrapolation for combined (scalar or vector) form factor.
plot_labels : function
    Adds axis labels and legend to existing plot instance.
results : function
    Writes results of fit settings, qualities, and parameters.
sigfig : function
    Rounds real numbers to some significant figure.
--------------------------------------------------------------------------------
"""


from calculators.chilogs.fcns import a_fermi
from calculators.fcns import E_out
from fitters import chiral
from math import copysign, floor, log10, sqrt
from matplotlib import pyplot as plt
from os import getcwd
from settings import fit
from settings.constants import (alphaV_continuum, mh_continuum, ml_continuum,
                                r1_a_continuum)
from settings.fit import *
from sys import stdout
import numpy as np
import pickle


def inputs_continuum():
    """
    ----------------------------------------------------------------------------
    Returns inputs for continuum extrapolations.
    ----------------------------------------------------------------------------
    Returns
    -------
    inputs : np.ndarray of dict
        Array of dictionary of input floats for continuum extrapolations.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    alphaV_continuum : float, from settings.constants
    mh_continuum : float, from settings.constants
    ml_continuum : float, from settings.constants
    numpy : module, as np
    ----------------------------------------------------------------------------
    Notes
    -----
    + See fileIOs.readers.inputs for complete list and description of inputs.
    ----------------------------------------------------------------------------
    """
    inputs = np.array([{'extrapolation': 'continuum'}])
    inputs[0]['a'] = 0
    inputs[0]['a_fm'] = 0
    inputs[0]['m0 * a'] = 0
    inputs[0]['alpha_V'] = alphaV_continuum
    inputs[0]['mh_sea'] = mh_continuum
    inputs[0]['mh_val'] = mh_continuum
    inputs[0]['ml_sea'] = ml_continuum
    inputs[0]['ml_val'] = ml_continuum
    return inputs


def params(avgs, cov, savename='result'):
    """
    ----------------------------------------------------------------------------
    Pickles fit parameters, where avgs are their central values and cov is their
    covariance matrix, saving result with name savename.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    avgs : numpy.ndarray of floats
        Central values of results of fit parameters.
    cov : numpy.ndarray of floats
        Covariance matrix of results of fit parameters with shape ({number of
        parameters, number of parameters}).
    savename : str (optional; default is 'result')
        Root name to use when saving output file.
    ----------------------------------------------------------------------------
    Results
    -------
    {savename}.p : file
        Results of fit parameters in pickled binary format.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    pickle : module
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function fitters.lsq.all.
    ----------------------------------------------------------------------------
    """
    pickle.dump((avgs, cov), open(savename + '.p', 'wb'))


def plot_data(inputs, data, axis=None, colormap=None, markerfill='none',
              markershape='o', savename='result'):
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
    markerfill : str (optional; default is 'none')
        Type of fill to use for markers of data points.
    markershape : str (optional; default is 'o')
        Type of shape to use for markers of data points.
    savename : str (optional; default is 'result')
        Root name to use when saving output file.
    ----------------------------------------------------------------------------
    Results
    -------
    {savename}.dat : file
        Values from plot of input data with error bars.
    matplotlib.container.ErrorbarContainer
        Markers at (Y = data), (X = inputs) with error bars about Y.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    datatype : str, from settings.fit
    ensemblesize : int, from settings.fit
    nexperiments : int, from settings.fit
    numpy : module, as np
    pyplot : module, from matplotlib, as plt
    sigfig : function
    ----------------------------------------------------------------------------
    """
    axis = axis if axis is not None else plt.gca()
    Es = np.array([inputs[xpmt]['E'] for xpmt in range(nexperiments)])
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
        label = '$a m_\\ell / a m_h =\\ {0} / {1}$'.format(aml, amh)
        for ensxpmt in range(ensemblesize):
            plt.errorbar(Es[xpmt + ensxpmt], data_avg[xpmt + ensxpmt],
                         yerr=([data_errlower[xpmt + ensxpmt]],
                               [data_errupper[xpmt + ensxpmt]]),
                         label=label, color=colors[ensemble], fmt=markershape,
                         fillstyle=markerfill)
            label = None
            aml_amhs.append(str(aml + ' / ' + amh).ljust(18))
    outputs = np.empty(nexperiments, dtype=[('aml_amh', 'a18'),
                                            ('E', float),
                                            ('ff_avg', float),
                                            ('ff_err+', float),
                                            ('ff_err-', float)])
    outputs['aml_amh'] = aml_amhs
    outputs['E'] = Es
    outputs['ff_avg'] = data_avg
    outputs['ff_err+'] = data_errupper
    outputs['ff_err-'] = data_errlower
    np.savetxt(savename + '.dat', outputs,
               fmt=('%18s', '%.6e', '%.6e', '%.6e', '%.6e'), delimiter='  ',
               header='  '.join(['a*ml / a*mh'.ljust(16), 'E'.ljust(12),
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
        color = axis._get_lines.get_next_color()
    if np.isscalar(yerr) or ((len(yerr) == len(y)) and (len(yerr.shape) == 1)):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin = y - yerr[0]
        ymax = y + yerr[1]
    elif len(yerr.shape) == 2:
        ymin = y - yerr[:, 0]
        ymax = y + yerr[:, 1]
    else:
        raise ValueError('invalid shape for yerr')
    axis.plot(x, y, color=color, label=label)
    axis.fill_between(x, ymax, ymin, color=color, alpha=alphafill)


def plot_fit(inputs, params, lattspace, alphafill=0.3, axis=None, color=None,
             label=None, linelength=None, savename='result'):
    """
    ----------------------------------------------------------------------------
    Plots fit line at some lattice spacing lattspace, given (X = inputs) and fit
    parameters params.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Creates fit inputs dictionary from last experiment of ensemble whose lattice
    spacing matches that of lattspace, or creates fit inputs dictionary from
    continuum constants if (lattspace = 0). Creates values for E_(K/pi) via
    numpy.linspace(linelength), where its bounds are made to match those of
    E_(K/pi) in inputs if (linelength = None). Applies
    fitters.chiral.fitfcn(x, p) with (x = {fit inputs}) and (p = params), which
    returns an array of gvar.GVars with standard propagation of errors. Plot
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
    lattspace : str
        Lattice spacing name.
    alphafill : float or NoneType (optional; default is 0.3)
        Alpha transparency of shading.
    axis : matplotlib.axes._subplots.AxesSubplot or NoneType (optional; default
           is None)
        Axis container for plot items.
    color : str or NoneType (optional; default is None)
        Color of plot line.
    label : str or NoneType (optional; default is None)
        Label of plot line.
    linelength : list or NoneType (optional; default is None)
        Length of fit line, as list of [min, max, numpoints] if specified, for
        use with numpy.linspace. If not specified, numpy.linspace bounds match
        those of E_(K/pi) in inputs.
    savename : str (optional; default is 'result')
        Root name to use when saving output file.
    ----------------------------------------------------------------------------
    Results
    -------
    {savename}_fit.dat : file
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
    nexperiments : int, from settings.fit
    numpy : module, as np
    plot_errfill : function
    pyplot : module, from matplotlib, as plt
    r1_a_continuum : float, from settings.constants
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
        aml = str(sigfig(fit_inputs[0]['a'] * fit_inputs[0]['ml_val'], n=3))
        amh = str(sigfig(fit_inputs[0]['a'] * fit_inputs[0]['mh_val'], n=3))
    else:
        fit_inputs = inputs_continuum()
        aml = str(sigfig(fit_inputs[0]['ml_val'] / r1_a_continuum, n=3))
        amh = str(sigfig(fit_inputs[0]['mh_val'] / r1_a_continuum, n=3))
    Es = np.array([inputs[xpmt]['E'] for xpmt in range(nexperiments)])
    if linelength is None:
        fit_Es = np.linspace(min(Es), max(Es), num=50)
    elif (type(linelength) == list) and (len(linelength) == 3):
        fit_Es = np.linspace(linelength[0], linelength[1], num=linelength[2])
    fit_ffs_avg = np.empty_like(fit_Es)
    fit_ffs_err = np.empty_like(fit_Es)
    fit_ffs_gvar = []
    for E in fit_Es:
        fit_inputs[0]['E'] = E
        fit_ffs_gvar.append(chiral.fitfcn(fit_inputs, params)[0])
    fit_ffs_avg = np.array([fit_ffs_gvar[Eindex].mean for Eindex in
                            range(len(fit_Es))])
    fit_ffs_err = np.array([fit_ffs_gvar[Eindex].sdev for Eindex in
                            range(len(fit_Es))])
    plot_errfill(fit_Es, fit_ffs_avg, fit_ffs_err, color=color, label=label,
                 alphafill=alphafill, axis=axis)
    np.savetxt(savename + '_fit.dat',
               np.vstack((fit_Es, fit_ffs_avg, fit_ffs_err)).T,
               fmt=('%.6e', '%.6e', '%.6e'), delimiter='  ',
               header='a*ml / a*mh = {}\n'.format(aml + ' / ' + amh) +
                      '  '.join(['E'.ljust(10), 'ff_avg'.ljust(12),
                                 'ff_err'.ljust(12)]))


def plot_fitavgs(inputs, params, axis=None, colormap=None, label=None,
                 linelength=None, linewidth=1.0, savename='result'):
    """
    ----------------------------------------------------------------------------
    Plots fit lines for all ensembles, given (X = inputs) and fit parameters
    params.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Creates inputs dictionaries from last experiment of each ensemble. Creates
    values for E_(K/pi) via numpy.linspace(linelength), where its bounds are
    made to match those of E_(K/pi) in inputs if (linelength = None). Applies
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
    axis : matplotlib.axes._subplots.AxesSubplot or NoneType (optional; default
           is None)
        Axis container for plot items.
    colormap : matplotlib.colors.LinearSegmentedColormap or NoneType (optional;
               default is None)
        Color map to use when creating successive plot line colors.
    label : str or NoneType (optional; default is None)
        Label of plot line.
    linelength : list or NoneType (optional; default is None)
        Length of fit line, as list of [min, max, numpoints] if specified, for
        use with numpy.linspace. If not specified, numpy.linspace bounds match
        those of E_(K/pi) in inputs.
    linewidth : float (optional; default is 1.0)
        Relative width of plot line.
    savename : str (optional; default is 'result')
        Root name to use when saving output file.
    ----------------------------------------------------------------------------
    Results
    -------
    {savename}_fits.dat : file
        Values from plot of ensemble fit averages.
    matplotlib.lines.Line2D
        Plot line of central values of fits.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    chiral : module, from fitters
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
    Es = np.array([inputs[xpmt]['E'] for xpmt in range(nexperiments)])
    if linelength is None:
        fit_Es = np.linspace(min(Es), max(Es), num=50)
    elif (type(linelength) == list) and (len(linelength) == 3):
        fit_Es = np.linspace(linelength[0], linelength[1], num=linelength[2])
    params_avg = {}
    for param in params.keys():
        params_avg[param] = params[param].mean
    outputs = np.empty((nensembles, len(fit_Es)))
    for ensemble in range(nensembles):
        fit_ffs = np.empty_like(fit_Es)
        fit_inputs = np.array([{'ensemble': ensemble}])
        for param in inputs[ensemble * ensemblesize].keys():
            fit_inputs[0][param] = inputs[ensemble * ensemblesize][param]
        for Eindex, E in enumerate(fit_Es):
            fit_inputs[0]['E'] = E
            fit_ffs[Eindex] = chiral.fitfcn(fit_inputs, params_avg)[0]
        axis.plot(fit_Es, fit_ffs, color=colors[ensemble], label=label,
                  linewidth=linewidth)
        outputs[ensemble] = fit_ffs
    outputs = np.vstack((fit_Es, outputs)).T
    np.savetxt(savename + '_fits.dat', outputs,
               fmt='%.6e', delimiter='  ',
               header='  '.join(['E'.ljust(10)] +
                                ['ff_avg(e={})'.format(ens).ljust(12)
                                 for ens in range(nensembles)]))


def plot_fitcombo(params_para, params_perp, alphafill=0.3, axis=None,
                  color=None, label=None, linelength=None, savename='result'):
    """
    ----------------------------------------------------------------------------
    Plots continuum extrapolation for combined (scalar or vector) form factor,
    given fit parameters for both parallel and perpendicular form factors.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Creates fit inputs dictionary from continuum constants. Creates values for
    q^2 via numpy.linspace(linelength), where its bounds are set to typical
    values if (linelength = None). Applies fitters.chiral.fitfcn(x, p1, p2) with
    (x = {fit inputs}), (p1 = params_para), and (p2 = params_perp), which
    returns an array of gvar.GVars with standard propagation of errors. Plot
    values and labels are ultimately saved via numpy.savetxt.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    params_para : dict of floats or gvar.BufferDict
        Fit parameters of parallel form factor in dictionary-like container.
    params_perp : dict of floats or gvar.BufferDict
        Fit parameters of perpendicular form factor in dictionary-like
        container.
    alphafill : float or NoneType (optional; default is 0.3)
        Alpha transparency of shading.
    axis : matplotlib.axes._subplots.AxesSubplot or NoneType (optional; default
           is None)
        Axis container for plot items.
    color : str or NoneType (optional; default is None)
        Color of plot line.
    label : str or NoneType (optional; default is None)
        Label of plot line.
    linelength : list or NoneType (optional; default is None)
        Length of fit line, as list of [min, max, numpoints] if specified, for
        use with numpy.linspace. If not specified, numpy.linspace bounds are set
        to (q^2 = {16, 23} (GeV)^2) to mirror bounds of typical datasets.
    savename : str (optional; default is 'result')
        Root name to use when saving output file.
    ----------------------------------------------------------------------------
    Results
    -------
    {savename}_fit.dat : file
        Values from plot of continuum extrapolation with errors.
    matplotlib.lines.Line2D
        Plot line of central values.
    matplotlib.collections.PolyCollection
        Plot shading of error bars about central line.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    E_out : function, from calculators.fcns
    chiral : module, from fitters
    numpy : module, as np
    plot_errfill : function
    pyplot : module, from matplotlib, as plt
    r1_a_continuum : float, from settings.constants
    sigfig : function
    ----------------------------------------------------------------------------
    """
    axis = axis if axis is not None else plt.gca()
    fit_inputs = inputs_continuum()
    if linelength is None:
        fit_q2s = np.linspace(16., 23., num=50)
    elif (type(linelength) == list) and (len(linelength) == 3):
        fit_q2s = np.linspace(linelength[0], linelength[1], num=linelength[2])
    fit_ffs_avg = np.empty_like(fit_q2s)
    fit_ffs_err = np.empty_like(fit_q2s)
    fit_ffs_gvar = []
    for q2 in fit_q2s:
        fit_inputs[0]['E'] = E_out(q2)
        fit_ffs_gvar.append(chiral.fitfcn(fit_inputs,
                                          params_para, params_perp)[0])
    fit_ffs_avg = np.array([fit_ffs_gvar[q2index].mean for q2index in
                            range(len(fit_q2s))])
    fit_ffs_err = np.array([fit_ffs_gvar[q2index].sdev for q2index in
                            range(len(fit_q2s))])
    plot_errfill(fit_q2s, fit_ffs_avg, fit_ffs_err, color=color, label=label,
                 alphafill=alphafill, axis=axis)
    aml = str(sigfig(fit_inputs[0]['ml_val'] / r1_a_continuum, n=3))
    amh = str(sigfig(fit_inputs[0]['mh_val'] / r1_a_continuum, n=3))
    np.savetxt(savename + '_fit.dat',
               np.vstack((fit_q2s, fit_ffs_avg, fit_ffs_err)).T,
               fmt=('%.6e', '%.6e', '%.6e'), delimiter='  ',
               header='a*ml / a*mh = {}\n'.format(aml + ' / ' + amh) +
                      '  '.join(['q2'.ljust(10), 'ff_avg'.ljust(12),
                                 'ff_err'.ljust(12)]))


def plot_labels(chi2=None, dof=None, legendloc='best', legendsize='12', p=None,
                xlims=None, ylims=None):
    """
    ----------------------------------------------------------------------------
    Adds axis labels and legend to existing plot instance, creates a title with
    chi^2, degrees of freedom, and p-value (if specified), and alters x- and/or
    y-axis limits (if specified).
    ----------------------------------------------------------------------------
    Parameters
    ----------
    chi2 : float
        Minimum chi^2 from fit.
    dof : int
        Degrees of freedom in fit.
    legendloc : str or int or tuple of floats or NoneType (optional; default is
                'best')
        Location of legend in plot figure.
    legendsize : str of int or float (optional; default is '12')
        Size of fonts used in legend.
    p : float
        p-value of fit.
    xlims : list or NoneType (optional; default is None)
        Custom bounds (min, max) to use with x-axis, if specified.
    ylims : list or NoneType (optional; default is None)
        Custom bounds (min, max) to use with y-axis, if specified.
    ----------------------------------------------------------------------------
    Results
    -------
    matplotlib.legend.Legend
        Plot legend that will display all previously declared plot labels.
    matplotlib.text.Text
        Plot text container(s) that comprise x- and y-axis labels, as well as
        title (if specified).
    ----------------------------------------------------------------------------
    Requirements
    ------------
    decayname : str, from settings.fit
    formfactor : str, from settings.fit
    numpy : module, as np
    pyplot : module, from matplotlib, as plt
    sigfig : function
    ----------------------------------------------------------------------------
    """
    if ((formfactor == 'para') or (formfactor == 'perp') or
        (formfactor == 'tensor')):
        if decayname == 'B2K':
            plt.xlabel('$r_1 E_K$')
        elif decayname == 'B2pi':
            plt.xlabel('$r_1 E_\\pi$')
    elif (formfactor == 'scalar') or (formfactor == 'vector'):
        plt.xlabel('$q^2\\! \\ (\\mathrm{{GeV}}^2\\!)$')
    if formfactor == 'para':
        plt.ylabel('$r_1^{1 / 2} f_\\parallel$')
    elif formfactor == 'perp':
        plt.ylabel('$r_1^{-1 / 2} f_\\perp$')
    elif formfactor == 'scalar':
        plt.ylabel('$f_0(q^2\\!)$')
    elif formfactor == 'tensor':
        plt.ylabel('$r_1^{-1 / 2} f_T$')
    elif formfactor == 'vector':
        plt.ylabel('$f_{\\!+}(q^2\\!)$')
    if xlims is not None:
        plt.xlim(xlims[0], xlims[1])
    if ylims is not None:
        plt.ylim(ylims[0], ylims[1])
    title = ''
    if (chi2 is not None) and (dof is not None):
        title += '$\\chi^2\\! / \\mathrm{{dof}} = {0} / {1}$'.format(
                                          sigfig(chi2, n=3), int(np.round(dof)))
        if p is not None:
            title += '$,\\ $'
    if p is not None:
        title += '$p = {}$'.format(sigfig(p, n=2))
    if title != '':
        plt.title(title)
    if legendloc is not None:
        plt.legend(loc=legendloc, prop={'size': legendsize}, numpoints=1)


def results(chi2=None, dof=None, fileouts=True, p=None, params=None,
            savename='result', stdouts=True):
    """
    ----------------------------------------------------------------------------
    Writes results of fit settings, qualities, and parameters to stdout and/or
    output file, saving latter with root name savename.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    All of chi2, dof, p, and params must be supplied if using single (parallel,
    perpendicular, or tensor) form factor. Only fit settings are written if
    using combined (scalar or vector) form factor.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    params : dict of floats or gvar.BufferDict (optional; default is None)
        Fit parameters in dictionary-like container.
    chi2 : float (optional; default is None)
        Minimum chi^2 from fit.
    dof : float (optional; default is None)
        Degrees of freedom in fit.
    p : float (optional; default is None)
        p-value of fit.
    fileouts : bool (optional; default is True)
        Determines if outputs are written to file.
    savename : str (optional; default is 'result')
        Root name to use when saving output file.
    stdouts : bool (optional; default is True)
        Determines if outputs are written to stdout.
    ----------------------------------------------------------------------------
    Results
    -------
    {savename}.txt : file (optional; created by default)
        Results of fit settings, qualities, and parameters.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    correlated : bool, from settings.fit
    fit : module, from settings
    getcwd : function, from os
    numpy : module, as np
    sigfig : function
    stdout: file, from settings
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'all of chi2, dof, p, and params must be supplied [...]'
        All of chi2, dof, p, and params must be supplied if using single
        (parallel, perpendicular, or tensor) form factor.
    ----------------------------------------------------------------------------
    """
    if fileouts and stdouts:
        stdout.write('\nWriting results to {}\n'.format(getcwd()))
    output = '#\n# Settings:\n'
    for setting in sorted(fit.__all__):
        output += '\t{0}\t{1}\n'.format(setting.rjust(12), eval(setting))
    if (fit.formfactor != 'vector') and (fit.formfactor != 'scalar'):
        if (chi2 is None) or (dof is None) or (p is None) or (params is None):
            raise ValueError('all of chi2, dof, p, and params must be ' +
                             'supplied if using f_para, f_perp, or f_T')
        if correlated:
            output += '#\n# Correlated'
        else:
            output += '#\n# Uncorrelated'
        output += ' Least Squares Fit:\n'
        output += '#\tchi^2 / dof = {0} / {1}\n'.format(sigfig(chi2, n=3),
                                                        int(np.round(dof)))
        output += '#\tp = {}\n'.format(sigfig(p, n=2))
        output += '#\n# Parameters:\n'
        for param in sorted(params.keys()):
            output += '\t{0}\t{1}\n'.format(param.rjust(5), params[param])
    output += '#\n'
    if stdouts:
        stdout.write(output)
    if fileouts:
        outfile = open(savename + '.txt', 'w')
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

