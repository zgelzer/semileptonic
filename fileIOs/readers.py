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
from sys import stdout
import argparse
import numpy as np
import os
import pickle


def args(maindir):
    """
    ----------------------------------------------------------------------------
    Reads in command-line arguments of main script run from directory maindir.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Passes maindir and arguments to function args_parse. Intended to be called
    before any other function in semileptonic module.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    maindir : str
        Directory of main script.
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
    args = argparse.ArgumentParser(description='Perform chiral fit on ' +
                                                    'form factors data.')
    args.add_argument('decayname',
                      help="name of decay (must be 'B2K' or 'B2pi')")
    args.add_argument('formfactor',
                      help="form factor to be computed (must be 'perp' or " +
                           "'para')")
    args.add_argument('datatype',
                      help="type of data (must be 'bs', 'jk', or 'jk2bs')")
    args.add_argument('-c', '--constrained', dest='constrained',
                      action='store_true',
                      help='uses constrained fit')
    args.add_argument('-C', '--correlated', dest='correlated',
                      action='store_true',
                      help='uses correlated fit')
    args.add_argument('-d', '--datadir', dest='datadir',
                      default='.',
                      help="data directory; default='.'")
    args.add_argument('-e', '--esize', dest='ensemblesize',
                      default=3, type=int,
                      help='ensemble size; default=3')
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
                      help='loads fit results from given LOAD.p and LOAD.txt')
    args.add_argument('-n', '--nsamples', dest='nsamples',
                      default=None,
                      help='number of samples to use from resampled data')
    args.add_argument('-o', '--outputdir', dest='outputdir',
                      default=None,
                      help="data directory; default='./result_{decayname}_" +
                           "{formfactor}'")
    args.add_argument('-p', '--plot', dest='plot',
                      action='store_true',
                      help='displays plot')
    args.add_argument('-s', '--su3', dest='SU3',
                      action='store_true',
                      help='uses SU(3) gauge')
    args.add_argument('-x', '--exclude', dest='exclude',
                      default=None,
                      help='excludes given experiment numbers (must be comma-' +
                           '-separated list)')
    args.add_argument('-X', '--inputs', dest='inputsource',
                      default='X.dat',
                      help="inputs source; default='X.dat'")
    args.add_argument('-Y', '--data', dest='datasource',
                      default='Y.dat',
                      help="data source; default='Y.dat'")
    args = args.parse_args()
    return args_parse(args, maindir)


def args_parse(args, maindir):
    """
    ----------------------------------------------------------------------------
    Parses command-line arguments args of main script run from directory
    maindir.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Applies tests; alters formats if necessary; adds arguments relating to data
    sizes; and saves important arguments to semileptonic/settings/fit.py.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    args : argparse.Namespace
        Strings of command-line arguments along with their values.
    maindir : str
        Directory of main script.
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
        Alteration and addition of arguments, as follows...
        args.exclude or args.include is converted to list of integers.
        args.fitlength is converted to list of floats if specified.
        args.nensembles is added.
        args.nexperiments is added.
        args.nexperiments_source is added.
        args.nsamples is converted to integer if specified.
        args.nsamples_source is added.
        args.outputdir is set to './chifit_{args.decayname}_{args.formfactor}'
                       if not specified.
        args.xpmtlist is added.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    os : module
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'invalid decay name'
        Decay name must be one of: 'B2K' (for B-->K) or 'B2pi' (for B-->pi).
    ValueError : 'invalid form factor'
        Form factor must be one of: 'para' (for parallel) or 'perp' (for
        perpendicular).
    ValueError : 'invalid data type'
        Data type must be one of: 'bs' (for bootstrap), 'jk' (for jackknife), or
        'jk2bs' (for jackknife that will be converted to bootstrap).
    ValueError : 'invalid fit length' [...]
        Fit length must be comma-separated list entered as: {min,max,numpoints}.
    ValueError : 'invalid experiment list' [...]
        Cannot simultaneously specify lists of experiments both to include and
        to exclude. Instead, summarize choices as single list of inclusions or
        of exclusions.
    ValueError : 'invalid number of experiments' [...]
        Number of experiments must be multiple of {ensemble size}.
    ValueError : 'invalid number of samples' [...]
        Number of bootstrap/jackknife samples must be greater than zero and less
        than {number of samples in data source}.
    ----------------------------------------------------------------------------
    """
    if args.decayname not in ['B2K', 'B2pi']:
        raise ValueError('invalid decay name')
    if args.formfactor not in ['para', 'perp']:
        raise ValueError('invalid form factor')
    if args.datatype not in ['bs', 'jk', 'jk2bs']:
        raise ValueError('invalid data type')
    if args.fitlength is not None:
        args.fitlength = [float(n) for n in args.fitlength.split(',')]
        if len(args.fitlength) != 3:
            raise ValueError('invalid length of fit line (must be comma-' +
                             'separated list entered as min,max,numpoints)')
    if args.outputdir is None:
        args.outputdir = './chifit_' + args.decayname + '_' + args.formfactor
    if (args.include is not None) and (args.exclude is not None):
        raise ValueError('invalid experiment list (cannot specify both ' +
                         'include and exclude experiment lists)')
    workdir = os.getcwd()
    os.chdir(args.datadir)
    args.nexperiments_source = len(np.loadtxt(args.inputsource))
    args.nsamples_source = int(len(np.loadtxt(args.datasource)) /
                               args.nexperiments_source)
    os.chdir(workdir)
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
    savefile = open(os.path.join(maindir, 'settings', 'fit.py'), 'w')
    savefile.write('# Fit Settings #\n')
    savefile.write('# -> generated by chifit.py\n')
    savefile.write('# -> not for user input; all manual changes will be ' +
                   'overwritten\n')
    savefile.write('# -> for help with changing fit settings, run \'python ' +
                   'chifit.py --help\'\n')
    savefile.write('\n')
    savefile.write('__all__ = [\'constrained\', \'correlated\', ' +
                   '\'datatype\', \'decayname\', \'ensemblesize\', ' +
                   '\'formfactor\', \'hardpiK\',  \'nexperiments\', \'SU3\']\n')
    savefile.write('\n')
    savefile.write('constrained = {}\n'.format(args.constrained))
    savefile.write('correlated = {}\n'.format(args.correlated))
    savefile.write('datatype = \'{}\'\n'.format(args.datatype))
    savefile.write('decayname = \'{}\'\n'.format(args.decayname))
    savefile.write('ensemblesize = {}\n'.format(args.ensemblesize))
    savefile.write('formfactor = \'{}\'\n'.format(args.formfactor))
    savefile.write('hardpiK = {}\n'.format(args.hardpiK))
    savefile.write('nexperiments = {}\n'.format(args.nexperiments))
    savefile.write('SU3 = {}\n'.format(args.SU3))
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


def results(source, usecols=(1,), delimiter='\t'):
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
    usecols : tuple of ints (optional; default is (1,))
        Columns of data to use from text source.
    delimiter : str (optional; default is '\\t')
        String used to separate columns of data in text source.
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
    stdout : file, from sys
    ----------------------------------------------------------------------------
    Notes
    -----
    + Sorting of fit parameters in pickled source should match that of text
      source. This is accomplished by default if using standard output files
      from previous run of chifit.py.
    + See function fitters.lsq.all.
    ----------------------------------------------------------------------------
    """
    stdout.write('\nReading results from {}\n'.format(
                            os.path.dirname(os.path.realpath(source + '.txt'))))
    dictkeys = np.loadtxt(source + '.txt', usecols=usecols, delimiter=delimiter,
                          dtype=str)
    dictkeys = sorted([key.strip() for key in dictkeys])
    return array2dict(gvar(pickle.load(open(source + '.p', 'rb'))), dictkeys)

