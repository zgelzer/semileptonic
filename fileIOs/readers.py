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
Defines functions for reading command-line arguments, inputs, data, initial fit
parameters, and previous results.
--------------------------------------------------------------------------------
Definitions
-----------
args : function
    Reads in command-line arguments of main script.
args_parse : function
    Parses command-line arguments of main script.
array2dict : function
    Converts array to dictionary.
data : function
    Reads in data from data source.
inputs : function
    Reads in inputs from inputs source.
params : function
    Reads a priori fit parameters.
results : function
    Reads results of fit parameters from pickled source.
--------------------------------------------------------------------------------
"""


from gvar import gvar
import argparse
import numpy as np
import os
import pickle


def args():
    """
    ----------------------------------------------------------------------------
    Reads in command-line arguments of main script.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Passes command-line arguments to function args_parse. Intended to be called
    before any other function in semileptonic module.
    ----------------------------------------------------------------------------
    Returns
    -------
    args : argparse.Namespace, from function args_parse
        Strings of command-line arguments along with their values.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    argparse : module
    args_parse : function
    ----------------------------------------------------------------------------
    """
    args = argparse.ArgumentParser(description='Perform chiral fits on ' +
                                               'semileptonic form factors.')
    args.add_argument('-c', '--constrained', dest='constrained',
                      action='store_true',
                      help='uses constrained fit')
    args.add_argument('-C', '--correlated', dest='correlated',
                      action='store_true',
                      help='uses correlated fit')
    args.add_argument('-d', '--decay', dest='decayname',
                      default=None,
                      help="name of decay (must be 'B2K' or 'B2pi')")
    args.add_argument('-D', '--dir', dest='datadir',
                      default='.',
                      help="data directory; default='.'")
    args.add_argument('-e', '--enssize', dest='ensemblesize',
                      default=3, type=int,
                      help='ensemble size; default=3')
    args.add_argument('-f', '--formfactor', dest='formfactor',
                      default=None,
                      help="form factor to be computed (must be 'para' or " +
                           "'perp')")
    args.add_argument('-H', '--hard', dest='hardpiK',
                      action='store_true',
                      help='uses hard pion/kaon')
    args.add_argument('-i', '--include', dest='include',
                      default=None,
                      help='includes only given experiment numbers (must be ' +
                           'comma-separated list)')
    args.add_argument('-l', '--length', dest='fitlength',
                      default=None,
                      help='sets curve length of averaged fit equation (must ' +
                           'be comma-separated list entered as min,max,' +
                           'numpoints)')
    args.add_argument('-L', '--load', dest='load',
                      default=None,
                      help='loads fit parameters from LOAD.p and fit ' +
                           'settings from LOAD.txt (except for ENSEMBLESIZE ' +
                           'and INCLUDE/EXCLUDE)')
    args.add_argument('-n', '--nsamples', dest='nsamples',
                      default=None,
                      help='number of samples to use from resampled data')
    args.add_argument('-o', '--outputdir', dest='outputdir',
                      default=None,
                      help="data directory; default='{decayname}" + os.sep +
                           "{formfactor}'")
    args.add_argument('-p', '--plot', dest='plot',
                      action='store_true',
                      help='displays plot')
    args.add_argument('-s', '--su3', dest='SU3',
                      action='store_true',
                      help='uses SU(3) gauge')
    args.add_argument('-t', '--type', dest='datatype',
                      default='bs',
                      help="type of data (must be 'bs'); default='bs'")
    args.add_argument('-x', '--exclude', dest='exclude',
                      default=None,
                      help='excludes given experiment numbers (must be comma-' +
                           'separated list)')
    args.add_argument('-X', '--inputs', dest='inputsource',
                      default='X.dat',
                      help="inputs source; default='X.dat'")
    args.add_argument('-Y', '--data', dest='datasource',
                      default='Y.dat',
                      help="data source; default='Y.dat'")
    args = args.parse_args()
    return args_parse(args)


def args_parse(args):
    """
    ----------------------------------------------------------------------------
    Parses command-line arguments args of main script.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Applies tests; alters formats if necessary; adds arguments relating to data
    and defaults; and saves important arguments to semileptonic/settings/fit.py.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    args : argparse.Namespace
        Strings of command-line arguments along with their values.
    ----------------------------------------------------------------------------
    Returns
    -------
    args : argparse.Namespace
        Strings of command-line arguments along with their values.
    ----------------------------------------------------------------------------
    Results
    -------
    fit.py : file
        Storage of important settings, located in script directory 'settings'.
    args : argparse.Namespace
        Alteration and addition of arguments, as follows:
        args.exclude or args.include is converted to list of integers if
        specified.
        args.fitlength is converted to list of floats if specified.
        args.nensembles is added.
        args.nexperiments is added.
        args.nexperiments_source is added.
        args.nsamples is converted to integer if specified.
        args.nsamples_source is added.
        args.outputdir is set to './{args.decayname}/{args.formfactor}' if not
        specified.
        args.workdir is added.
        args.xpmtlist is added.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    os : module
    results : function
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'must specify decay name'
        Must specify decay name if not loading previous results.
    ValueError : 'invalid decay name'
        Decay name must be one of: 'B2K' (for B-->K) or 'B2pi' (for B-->pi).
    ValueError : 'must specify form factor'
        Must specify form factor if not loading previous results.
    ValueError : 'invalid form factor'
        Form factor must be one of: 'para' (for parallel) or 'perp' (for
        perpendicular).
    ValueError : 'must specify data type'
        Must specify data type if not loading previous results.
    ValueError : 'invalid data type'
        Data type must be 'bs' (for bootstrap).
    ValueError : 'invalid fit length [...]'
        Fit length must be comma-separated list entered as: {min,max,numpoints}.
    ValueError : 'invalid experiment list [...]'
        Cannot simultaneously specify lists of experiments both to include and
        to exclude. Instead, summarize choices as single list of inclusions or
        of exclusions.
    ValueError : 'invalid number of experiments [...]'
        Number of experiments must be multiple of {ensemble size}.
    ValueError : 'invalid number of samples [...]'
        Number of bootstrap/jackknife samples must be greater than zero and less
        than {number of samples in data source}.
    ----------------------------------------------------------------------------
    Notes
    -----
    + Ignores previous values of args.ensemblesize and args.nexperiments if
      loading fit settings from results of previous run, so that user may have
      full control when importing data during each run.
    ----------------------------------------------------------------------------
    """
    if args.load is not None:
        fitsettings = np.loadtxt(args.load + '.txt', delimiter='\t', dtype=str,
                                 usecols=(1, 2))[:results.func_defaults[1]]
        names = [fitsetting.strip() for fitsetting in fitsettings[:, 0]]
        values = []
        for value in fitsettings[:, 1]:
            try:
                values.append(eval(value))
            except NameError:
                values.append(str(value))
        for name in names:
            if (name != 'ensemblesize') and (name != 'nexperiments'):
                value = values[np.where(np.asarray(names) == name)[0][0]]
                args.__setattr__(name, value)
    else:
        if args.decayname is None:
            raise ValueError('must specify decay name')
        elif args.decayname not in ['B2K', 'B2pi']:
            raise ValueError('invalid decay name')
        if args.formfactor is None:
            raise ValueError('must specify form factor')
        elif args.formfactor not in ['para', 'perp']:
            raise ValueError('invalid form factor')
        if args.datatype is None:
            raise ValueError('must specify data type')
        elif args.datatype not in ['bs']:
            raise ValueError('invalid data type')
    if args.fitlength is not None:
        args.fitlength = [float(n) for n in args.fitlength.split(',')]
        if len(args.fitlength) != 3:
            raise ValueError('invalid length of fit line (must be comma-' +
                             'separated list entered as min,max,numpoints)')
    if args.outputdir is None:
        args.outputdir = args.decayname + os.sep + args.formfactor
    if (args.include is not None) and (args.exclude is not None):
        raise ValueError('invalid experiment list (cannot specify both ' +
                         'include and exclude experiment lists)')
    args.workdir = os.getcwd()
    os.chdir(args.datadir)
    args.nexperiments_source = len(np.loadtxt(args.inputsource))
    args.nsamples_source = int(len(np.loadtxt(args.datasource)) /
                               args.nexperiments_source)
    os.chdir(args.workdir)
    if args.include is not None:
        args.include = [int(xpmt) for xpmt in args.include.split(',')]
        args.nexperiments = len(args.include)
        args.xpmtlist = args.include
    elif args.exclude is not None:
        args.exclude = [int(xpmt) for xpmt in args.exclude.split(',')]
        args.nexperiments = args.nexperiments_source - len(args.exclude)
        args.xpmtlist = [xpmt for xpmt in list(range(args.nexperiments_source))
                         if xpmt not in args.exclude]
    else:
        args.nexperiments = args.nexperiments_source
        args.xpmtlist = list(range(args.nexperiments_source))
    if args.nexperiments % args.ensemblesize == 0:
        args.nensembles = int(args.nexperiments / args.ensemblesize)
    else:
        raise ValueError('invalid number of experiments (must be multiple of ' +
                         'ensemble size)')
    if args.nsamples is not None:
        args.nsamples = int(args.nsamples)
        if (args.nsamples <= 0) or (args.nsamples > args.nsamples_source):
            raise ValueError('invalid number of samples (must be > 0 and <= ' +
                             'number of samples in data source)')
    else:
        args.nsamples = args.nsamples_source
    savefile = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                   '..', 'settings', 'fit.py'), 'w')
    savefile.write('# Fit Settings #\n')
    savefile.write('# -> generated by semileptonic/fileIOs/readers.py\n')
    savefile.write('# -> not for user input; all manual changes will be ' +
                   'overwritten\n')
    savefile.write('# -> for help with changing fit settings, run ' +
                   '\'./chifit.py --help\'\n')
    names = ['SU3', 'constrained', 'correlated', 'datatype', 'decayname',
             'ensemblesize', 'formfactor', 'hardpiK', 'nexperiments']
    savefile.write('__all__ = [' +
                   ', '.join(["'" + name + "'" for name in names]) + ']\n')
    for name in names:
        value = args.__getattribute__(name)
        if type(value) is str:
            savefile.write("{0} = '{1}'\n".format(name, value))
        else:
            savefile.write('{0} = {1}\n'.format(name, value))
    savefile.close()
    return args


def array2dict(array, keys):
    """
    ----------------------------------------------------------------------------
    Sequentially converts array-like array to dictionary with keys.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    array : numpy.ndarray or list or array-like
        Array-like object to be converted.
    keys : list of strs
        Strings to use as keys in dictionary.
    ----------------------------------------------------------------------------
    Returns
    -------
    dictionary : dict
        Dictionary of keys with values from array.
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function results.
    ----------------------------------------------------------------------------
    """
    dictionary = {}
    index = 0
    for key in keys:
        dictionary[key] = array[index]
        index += 1
    return dictionary


def data(source, xpmtlist, nexperiments_source, nsamples, nsamples_source,
         usecols=(2,)):
    """
    ----------------------------------------------------------------------------
    Reads in data from columns usecols in source according to experiment list
    xpmtlist, sampling desired number of samples nsamples from original data
    whose shape is (nexperiments_source, nsamples_source).
    ----------------------------------------------------------------------------
    Parameters
    ----------
    source : str
        Name of data source.
    xpmtlist : list of ints
        List of experiments to be included.
    nexperiments_source : int
        Number of experiments in data source.
    nsamples : int
        Number of samples to use from data source.
    nsamples_source : int
        Number of samples in data source.
    usecols : tuple of ints (optional; default is (2,))
        Columns of data to use from data source.
    ----------------------------------------------------------------------------
    Returns
    -------
    np.ndarray
        Desired data, with shape ({number of experiments}, {number of samples}).
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    ----------------------------------------------------------------------------
    """
    data = np.loadtxt(source, usecols=usecols)
    data = data.reshape((nexperiments_source, nsamples_source))
    return data[xpmtlist, :nsamples]


def inputs(source, xpmtlist):
    """
    ----------------------------------------------------------------------------
    Reads in inputs from source according to experiment list xpmtlist.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    source : str
        Name of inputs source.
    xpmtlist : list of ints
        List of experiments to be included.
    ----------------------------------------------------------------------------
    Returns
    -------
    inputs : np.ndarray of dicts
        Array of {number of experiments} total dictionaries, where each
        dictionary stores input floats for particular experiment.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    ----------------------------------------------------------------------------
    Notes
    -----
    + Input floats for each experiment are as follows:
        > 'a_fm' : lattice spacing in fm
        > 'a' : lattice spacing in r_1 units
        > 'ml_val' : mass of light valence quark in r_1 units
        > 'ml_sea' : mass of light sea quark in r_1 units
        > 'mh_val' : mass of heavy valence quark in r_1 units
        > 'mh_sea' : mass of heavy sea quark in r_1 units
        > 'E' : energy of pion/Kaon in r_1 units
    ----------------------------------------------------------------------------
    """
    source = np.loadtxt(source)
    inputs = np.array([{'xpmt': xpmt} for xpmt in xpmtlist])
    for i, xpmt in enumerate(xpmtlist):
        inputs[i]['a_fm']   = float(source[xpmt][0])
        inputs[i]['a']      = float(source[xpmt][1])
        inputs[i]['ml_val'] = float(source[xpmt][2])
        inputs[i]['ml_sea'] = float(source[xpmt][3])
        inputs[i]['mh_val'] = float(source[xpmt][4])
        inputs[i]['mh_sea'] = float(source[xpmt][5])
        inputs[i]['E']      = float(source[xpmt][6])
    return inputs


def params():
    """
    ----------------------------------------------------------------------------
    Reads a priori fit parameters from settings/params.py.
    ----------------------------------------------------------------------------
    Returns
    -------
    params_apriori : dict of gvar.gvars or NoneType
        A priori estimates of fit parameters, as Gaussian variables with widths,
        if specified; meant to be used with constrained fit.
    params_initval : dict of floats or NoneType
        Initial values of fit parameters, if specified; meant to be used with
        unconstrained fit.
    ----------------------------------------------------------------------------
    Raises
    ------
    ImportError
        Must know if fit is constrained or not; thus fit parameters may be read
        only after command-line arguments have been parsed. See function
        args_parse.
    ----------------------------------------------------------------------------
    Notes
    -----
    + See settings/params.py for complete descriptions of fit parameters.
    ----------------------------------------------------------------------------
    """
    try:
        from settings.fit import constrained
    except ImportError:
        raise ImportError('fit parameters may be read only after ' +
                          'command-line arguments have been parsed')
    from settings.params import params
    if constrained:
        params_apriori = params
        params_initval = None
    else:
        params_apriori = None
        params_initval = {}
        for param in params.keys():
            params_initval[param] = params[param].mean
    return params_apriori, params_initval


def results(source, delimiter='\t', skipfirst=9, skiplast=None, usecols=(1,)):
    """
    ----------------------------------------------------------------------------
    Reads array of results of fit parameters from pickled source and returns it
    as dictionary with keys from text source.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    source : str
        Name of pickled and text sources, to which '.p' or '.txt' will be
        appended (respectively).
    delimiter : str (optional; default is '\\t')
        String used to separate columns of data in text source.
    skipfirst : int or NoneType (optional; default is 9)
        Number of data values to ignore from beginning of text source. Default
        of 9 will skip values pertaining to fit settings; should not be changed.
    skiplast : int or NoneType (optional; default is None)
        Number of data values to ignore from end of text source. Default of None
        will read every last fit parameter; should not be changed.
    usecols : tuple of ints (optional; default is (1,))
        Columns of data to use from text source.
    ----------------------------------------------------------------------------
    Returns
    -------
    dict of gvar.GVars
        Dictionary of fit parameters as Gaussian variables with widths and their
        attached covariance matrix.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    array2dict : function
    gvar : class, from gvar
    numpy : module, as np
    os : module
    pickle : module
    ----------------------------------------------------------------------------
    Notes
    -----
    + Sorting of fit parameters in pickled source should match that of text
      source. This is accomplished by default if using standard output files
      from previous run of chifit.py.
    + See function fitters.lsq.all.
    ----------------------------------------------------------------------------
    """
    dictkeys = np.loadtxt(source + '.txt', delimiter=delimiter, dtype=str,
                          usecols=usecols)[skipfirst:skiplast]
    dictkeys = sorted([key.strip() for key in dictkeys])
    return array2dict(gvar(pickle.load(open(source + '.p', 'rb'))), dictkeys)

