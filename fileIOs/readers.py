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


from datetime import datetime as dt
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
    os : module
    ----------------------------------------------------------------------------
    """
    args = argparse.ArgumentParser(description='Performs chiral fit to ' +
                                               'lattice QCD form factors.')
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
                      default=None,
                      help='ensemble size; default is guessed from INPUTSOURCE')
    args.add_argument('-f', '--formfactor', dest='formfactor',
                      default=None,
                      help="form factor to be computed (must be 'para', " +
                           "'perp', 'scalar', 'tensor', or 'vector')")
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
                      default=None, nargs='*', action='append',
                      help='loads fit parameters from LOAD.p and fit ' +
                           'settings from LOAD.txt (except for ENSEMBLESIZE ' +
                           'and INCLUDE/EXCLUDE)')
    args.add_argument('-n', '--nsamples', dest='nsamples',
                      default=None,
                      help='number of samples to use from resampled data')
    args.add_argument('-o', '--outputdir', dest='outputdir',
                      default=None,
                      help="data directory; default='{decay name}" + os.sep +
                           "{form factor}'")
    args.add_argument('-p', '--plot', dest='plot',
                      action='store_true',
                      help='displays plot')
    args.add_argument('-s', '--save', dest='savename',
                      default=None,
                      help='name to use when saving results; default=' +
                           '{date/time}')
    args.add_argument('-S', '--SU3', dest='SU3',
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
    Applies tests, alters formats if necessary, adds file-related arguments, and
    saves important settings to settings.fit.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments as object attributes.
    ----------------------------------------------------------------------------
    Returns
    -------
    args : argparse.Namespace
        Command-line arguments as object attributes. Arguments are altered or
        added as follows:
        args.fitlength is converted to list of floats if specified;
        args.outputdir is set to './{args.decayname}/{args.formfactor}' if not
        specified;
        args.savename is set to date/time as 'YYYYMMDD-hhmm' if not specified or
        if empty;
        args.workdir is added.
    ----------------------------------------------------------------------------
    Results
    -------
    {outputdir} : directory or directories
        Output directory or directories; default is
        './{decay name}/{form factor}'.
    fit.py : file
        Storage of important settings, located in script directory 'settings'.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    datetime : class, from datetime, as dt
    os : module
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'must load completed runs [...]'
        Scalar or vector form factors require parallel and perpendicular form
        factor results to already have been obtained; these results are then
        loaded via '--load' command-line arguments.
    ValueError : 'must specify form factor'
        Must specify form factor if not loading previous results.
    ValueError : 'invalid form factor'
        Form factor must be one of: 'para' (for parallel), 'perp' (for
        perpendicular), 'scalar', 'tensor', or 'vector'.
    ValueError : 'invalid length of fit line [...]'
        Fit length must be comma-separated list entered as: {min,max,numpoints}.
    ----------------------------------------------------------------------------
    Notes
    -----
    + Divided into two main parsing paths, determined by nature of form factor:
      combined (scalar or vector) or single (parallel, perpendicular, or
      tensor). See functions args_parse_combined, args_parse_single.
    ----------------------------------------------------------------------------
    """
    args.workdir = os.getcwd()
    savenames = ['SU3', 'constrained', 'correlated', 'datatype', 'decayname',
                 'formfactor', 'hardpiK']
    if (args.formfactor == 'scalar') or (args.formfactor == 'vector'):
        if args.load is not None:
            args = args_parse_combined(args, savenames)
        else:
            raise ValueError("must load completed runs via '-L {para} " +
                             "-L {perp}', in any order")
    elif ((args.formfactor in ['para', 'perp', 'tensor']) or
          (args.load is not None)):
        args = args_parse_single(args, savenames)
        savenames = sorted(savenames + ['ensemblesize', 'nexperiments'])
    elif args.formfactor is None:
        raise ValueError('must specify form factor')
    else:
        raise ValueError('invalid form factor')
    if (args.savename is None) or (args.savename == ''):
        args.savename = dt.now().strftime('%Y%m%d-%H%M')
    if args.outputdir is None:
        args.outputdir = args.decayname + os.sep + args.formfactor
    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)
    if args.fitlength is not None:
        args.fitlength = [float(n) for n in args.fitlength.split(',')]
        if len(args.fitlength) != 3:
            raise ValueError('invalid length of fit line (must be comma-' +
                             'separated list entered as min,max,numpoints)')
    savefile = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                   '..', 'settings', 'fit.py'), 'w')
    savefile.write('# Fit Settings #\n')
    savefile.write('# -> generated by semileptonic/fileIOs/readers.py\n')
    savefile.write('# -> not for user input; all manual changes will be ' +
                   'overwritten\n')
    savefile.write('# -> for help with changing fit settings, run ' +
                   '\'./chifit.py --help\'\n')
    savefile.write('__all__ = [' +
                   ', '.join(["'" + savename + "'" for savename in savenames]) +
                   ']\n')
    for savename in savenames:
        value = args.__getattribute__(savename)
        if type(value) is str:
            savefile.write("{0} = '{1}'\n".format(savename, value))
        else:
            savefile.write('{0} = {1}\n'.format(savename, value))
    savefile.close()
    return args


def args_parse_combined(args, savenames):
    """
    ----------------------------------------------------------------------------
    Parses command-line arguments args of main script if using combined (scalar
    or vector) form factor, loading only those fit settings that are contained
    within savenames.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Applies tests, alters formats if necessary, and adds file-related arguments.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments as object attributes.
    savenames : list of strs
        Names of fit settings to be loaded (and eventually saved).
    ----------------------------------------------------------------------------
    Returns
    -------
    args : argparse.Namespace
        Command-line arguments as object attributes. Arguments are altered or
        added as follows:
        args.{savename} (for each savename in savenames parameter, except
        'formfactor') is changed to its value from {args.loadpara}.txt (without
        loss of generality);
        args.load is converted to list of strings;
        args.loadpara is added;
        args.loadperp is added.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    numpy : module, as np
    results : function
    ----------------------------------------------------------------------------
    Raises
    ------
    ValueError : 'invalid form factor in load file [...]'
        May only load results for parallel and perpendicular form factors.
    ValueError : 'must load both parallel and perpendicular form factors'
        Must load results for both parallel and perpendicular form factors.
    ValueError : 'fit settings of both load files must agree [...]'
        Fit settings listed in parameter savenames (except for 'formfactor')
        must agree in both load files (see Notes); conflicts are printed.
    ValueError : 'invalid number of load statements [...]'
        Must specify exactly two load files (one for each parallel and
        perpendicular form factor).
    ----------------------------------------------------------------------------
    Notes
    -----
    + After it has been tested that results for both parallel and perpendicular
      form factors were provided, all fit settings listed in parameter savenames
      (except for 'formfactor') are compared between {args.loadpara}.txt and
      {args.loadperp}.txt to ensure consistency. This requires that parameter
      savenames does not currently contain 'ensemblesize' or 'nexperiments'.
    + See function args_parse.
    ----------------------------------------------------------------------------
    """
    if len(args.load) == 2:
        args.load = [args.load[0][0], args.load[1][0]]
        loaddicts = []
        for loadfile in args.load:
            loaddict = {}
            fitsettings = np.loadtxt(loadfile + '.txt', delimiter='\t',
                                     dtype=str, usecols=(1, 2),
                                     )[:results.func_defaults[1]]
            values = []
            for value in fitsettings[:, 1]:
                try:
                    values.append(eval(value))
                except NameError:
                    values.append(str(value))
            names = np.array([name.strip() for name in fitsettings[:, 0]])
            for savename in savenames:
                loaddict[savename] = values[np.where(names == savename)[0][0]]
                if savename != 'formfactor':
                    args.__setattr__(savename, loaddict[savename])
            if loaddict['formfactor'] not in ['para', 'perp']:
                raise ValueError("invalid form factor in load file (must be " +
                                 "'para' or 'perp')")
            elif loaddict['formfactor'] == 'para':
                args.loadpara = loadfile
            elif loaddict['formfactor'] == 'perp':
                args.loadperp = loadfile
            loaddicts.append(loaddict)
        try:
            assert args.loadpara != args.loadperp
        except AttributeError:
            raise ValueError('must load both parallel and perpendicular form' +
                             'factors')
        conflicts = reduce(lambda S1, S2: S1.difference(S2),
                           (set(loaddict.items()) for loaddict in loaddicts))
        conflicts = sorted([name[0] for name in conflicts])
        conflicts.remove('formfactor')
        if len(conflicts) != 0:
            raise ValueError('fit settings of both load files must agree: ' +
                             ', '.join(conflicts))
    else:
        raise ValueError("invalid number of load statements (must specify " +
                         "exactly two load files via '-L {para} -L {perp}', " +
                         "in any order)")
    return args


def args_parse_single(args, savenames):
    """
    ----------------------------------------------------------------------------
    Parses command-line arguments args of main script if using single (parallel,
    perpendicular, or tensor) form factor, loading only those fit settings that
    are contained within savenames (if such loading is desired).
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Applies tests, alters formats if necessary, and adds data-related arguments.
    ----------------------------------------------------------------------------
    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments as object attributes.
    savenames : list of strs
        Names of fit settings to be loaded (and eventually saved).
    ----------------------------------------------------------------------------
    Returns
    -------
    args : argparse.Namespace
        Command-line arguments as object attributes. Arguments are altered or
        added as follows:
        args.{savename} (for each savename in savenames parameter) is changed to
        its value from {args.load}.txt if loading previous results;
        args.exclude or args.include is converted to list of integers if
        specified;
        args.ensemblesize is converted to integer if specified, otherwise it is
        guessed from args.inputsource (see Notes);
        args.load is converted to string if specified;
        args.nensembles is added;
        args.nexperiments is added;
        args.nexperiments_source is added;
        args.nsamples is converted to integer if specified;
        args.nsamples_source is added;
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
    ValueError : 'invalid form factor'
        Form factor must be one of: 'para' (for parallel), 'perp' (for
        perpendicular), 'scalar', 'tensor', or 'vector'.
    ValueError : 'invalid form factor for two load statements [...]'
        Must specify only one load file, or must specify combined form factor if
        using two load files.
    ValueError : 'invalid number of load statements'
        May not use more than two load files.
    ValueError : 'must specify decay name'
        Must specify decay name if not loading previous results.
    ValueError : 'invalid decay name'
        Decay name must be one of: 'B2K' (for B-->K) or 'B2pi' (for B-->pi).
    ValueError : 'must specify data type'
        Must specify data type if not loading previous results.
    ValueError : 'invalid data type'
        Data type must be 'bs' (for bootstrap).
    ValueError : 'invalid experiment list [...]'
        Cannot simultaneously specify lists of experiments both to include and
        to exclude. Instead, summarize choices as single list of inclusions or
        of exclusions.
    ValueError : 'invalid number of experiments [...]'
        Number of experiments must be multiple of {ensemble size}.
    ValueError : 'invalid number of samples [...]'
        Number of bootstrap/jackknife samples must be greater than two and less
        than {number of samples in data source}.
    ----------------------------------------------------------------------------
    Notes
    -----
    + If loading fit settings from results of previous run (in {args.load}.txt),
      said values of ensemblesize and nexperiments are ignored so that user may
      have full control when importing data during each run. This requires that
      parameter savenames does not currently contain 'ensemblesize' or
      'nexperiments'.
    + If not specified, args.ensemblesize is guessed from args.inputsource by
      comparing ratios of valence light- to heavy-quark masses to that of the
      first experiment.
    + See function args_parse.
    ----------------------------------------------------------------------------
    """
    if args.load is not None:
        if len(args.load) == 1:
            args.load = args.load[0][0]
            fitsettings = np.loadtxt(args.load + '.txt', delimiter='\t',
                                     dtype=str, usecols=(1, 2),
                                     )[:results.func_defaults[1]]
            values = []
            for value in fitsettings[:, 1]:
                try:
                    values.append(eval(value))
                except NameError:
                    values.append(str(value))
            names = np.array([name.strip() for name in fitsettings[:, 0]])
            for savename in savenames:
                args.__setattr__(savename,
                                 values[np.where(names == savename)[0][0]])
        elif len(args.load) == 2:
            if ((args.formfactor is not None) and
                (args.formfactor not in ['para', 'perp', 'tensor'])):
                raise ValueError('invalid form factor')
            else:
                raise ValueError('invalid form factor for two load ' +
                                 'statements (must specify only one load ' +
                                 'file, or must specify combined form factor ' +
                                 'if using two load files)')
        else:
            raise ValueError('invalid number of load statements')
    else:
        if args.decayname is None:
            raise ValueError('must specify decay name')
        elif args.decayname not in ['B2K', 'B2pi']:
            raise ValueError('invalid decay name')
        if args.datatype is None:
            raise ValueError('must specify data type')
        elif args.datatype not in ['bs']:
            raise ValueError('invalid data type')
    if (args.include is not None) and (args.exclude is not None):
        raise ValueError('invalid experiment list (cannot specify both ' +
                         'include and exclude experiment lists)')
    os.chdir(args.datadir)
    Xheader = np.array(open(args.inputsource).readline().split()[1:])
    Xcol_mh = np.where(Xheader == 'mh_val')[0][0]
    Xcol_ml = np.where(Xheader == 'ml_val')[0][0]
    Xsource = np.loadtxt(args.inputsource)
    args.nexperiments_source = len(Xsource)
    args.nsamples_source = int(len(np.loadtxt(args.datasource)) /
                               args.nexperiments_source)
    os.chdir(args.workdir)
    if args.ensemblesize is None:
        Xratios = np.array([Xsource[xpmt][Xcol_ml] / Xsource[xpmt][Xcol_mh]
                            for xpmt in range(args.nexperiments_source)])
        args.ensemblesize = np.sum([np.allclose(Xratios[:1], Xratios[:xpmt])
                                    for xpmt in range(1, len(Xratios))])
    else:
        args.ensemblesize = int(args.ensemblesize)
    if args.include is not None:
        args.include = [int(xpmt) for xpmt in args.include.split(',')]
        args.nexperiments = len(args.include)
        args.xpmtlist = args.include
    elif args.exclude is not None:
        args.exclude = [int(xpmt) for xpmt in args.exclude.split(',')]
        args.nexperiments = args.nexperiments_source - len(args.exclude)
        args.xpmtlist = [xpmt for xpmt in
                         list(range(args.nexperiments_source)) if xpmt not
                         in args.exclude]
    else:
        args.nexperiments = args.nexperiments_source
        args.xpmtlist = list(range(args.nexperiments_source))
    if args.nexperiments % args.ensemblesize == 0:
        args.nensembles = int(args.nexperiments / args.ensemblesize)
    else:
        raise ValueError('invalid number of experiments (must be ' +
                         'multiple of ensemble size)')
    if args.nsamples is not None:
        args.nsamples = int(args.nsamples)
        if (args.nsamples <= 2) or (args.nsamples > args.nsamples_source):
            raise ValueError('invalid number of samples (must be > 2 and ' +
                             '<= number of samples in data source)')
    else:
        args.nsamples = args.nsamples_source
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


def data(source, xpmtlist, nexperiments_source, nsamples, nsamples_source):
    """
    ----------------------------------------------------------------------------
    Reads in data from source according to experiment list xpmtlist, sampling
    desired number of samples nsamples from original data whose shape is
    (nexperiments_source, nsamples_source).
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
    Raises
    ------
    ValueError : 'header in {source} must contain ff'
        Source must contain commented header line containing 'ff' (see Notes);
        this implies that source must also contain float values for 'ff'.
    ----------------------------------------------------------------------------
    Notes
    -----
    + Source should have commented ('#') header line with variable name 'ff'
      (form factor) describing its corresponding column of float values.
    ----------------------------------------------------------------------------
    """
    header = np.array(open(source).readline().split()[1:])
    if ('ff' not in header):
        raise ValueError('header in ' + source + ' must contain ff')
    usecol = (np.where(header == 'ff')[0][0],)
    data = np.loadtxt(source, usecols=usecol)
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
    Raises
    ------
    ValueError : 'header in {source} must contain [...]'
        Source must contain commented header line describing all inputs (see
        Notes); this implies that source must also contain float values for each
        of these variable names.
    ----------------------------------------------------------------------------
    Notes
    -----
    + Input floats for each experiment are as follows:
        > E : energy of pion/Kaon in r_1 units
        > a : lattice spacing in r_1 units
        > a_fm : lattice spacing in fm
        > alpha_V : renormalized QCD coupling in heavy-quark scheme [1]
        > m0 * a : bare mass of clover b-quark on lattice [2]
        > mh_sea : mass of heavy sea quark in r_1 units
        > mh_val : mass of heavy valence quark in r_1 units
        > ml_sea : mass of light sea quark in r_1 units
        > ml_val : mass of light valence quark in r_1 units
    + Source should have commented ('#') header line with variable names
      describing each corresponding column of float values. These variable names
      must exactly match those listed in the previous note (see above).
    ----------------------------------------------------------------------------
    References
    ----------
    [1] G. P. Lepage and P. Mackenzie, "On the Viability of Lattice Perturbation
        Theory", Phys. Rev. D 48, 2250 (1993) [arXiv:hep-lat/9209022].
    [2] J. Bailey, et al. (Fermilab Lattice and MILC Collaborations), "Update of
        |V_{cb}| from the B --> D* l nu form factor at zero recoil with three-
        flavor lattice QCD", Phys. Rev. D 89, 114504 (2014) [arXiv:1403.0635
        [hep-lat]].
    ----------------------------------------------------------------------------
    """
    header = np.array(open(source).readline().split()[1:])
    if (('E' not in header) or ('a' not in header) or ('a_fm' not in header) or
        ('mh_sea' not in header) or ('mh_val' not in header) or
        ('ml_sea' not in header) or ('ml_val' not in header)):
        raise ValueError('header in ' + source + ' must contain ' +
                         'E, a, a_fm, mh_sea, mh_val, ml_sea, ml_val')
    column_E = np.where(header == 'E')[0][0]
    column_a = np.where(header == 'a')[0][0]
    column_afm = np.where(header == 'a_fm')[0][0]
    column_alpha = np.where(header == 'alpha_V')[0][0]
    column_am0 = np.where(header == 'm0*a')[0][0]
    column_mhsea = np.where(header == 'mh_sea')[0][0]
    column_mhval = np.where(header == 'mh_val')[0][0]
    column_mlsea = np.where(header == 'ml_sea')[0][0]
    column_mlval = np.where(header == 'ml_val')[0][0]
    source = np.loadtxt(source)
    inputs = np.array([{'xpmt': xpmt} for xpmt in xpmtlist])
    for i, xpmt in enumerate(xpmtlist):
        inputs[i]['E']      = float(source[xpmt][column_E])
        inputs[i]['a']      = float(source[xpmt][column_a])
        inputs[i]['a_fm']   = float(source[xpmt][column_afm])
        inputs[i]['alpha_V']   = float(source[xpmt][column_alpha])
        inputs[i]['m0 * a']   = float(source[xpmt][column_am0])
        inputs[i]['mh_sea'] = float(source[xpmt][column_mhsea])
        inputs[i]['mh_val'] = float(source[xpmt][column_mhval])
        inputs[i]['ml_sea'] = float(source[xpmt][column_mlsea])
        inputs[i]['ml_val'] = float(source[xpmt][column_mlval])
    return inputs


def params():
    """
    ----------------------------------------------------------------------------
    Reads a priori fit parameters from settings.params.
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
    + See settings.params for complete descriptions of fit parameters.
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

