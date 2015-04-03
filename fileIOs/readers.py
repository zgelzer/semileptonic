"""Defines functions for reading command-line arguments, inputs, data, and
initial fit parameters."""


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


from gvar import gvar
from sys import stdout
import argparse
import numpy as np
import os
import pickle


def args(maindir):
    """Reads in command-line arguments of script run from directory maindir."""
    args = argparse.ArgumentParser(description='Perform chiral fit on ' +
                                                    'form factors data.')
    args.add_argument('decayname',
                      help="name of decay (must be 'B2K' or 'B2pi')")
    args.add_argument('formfactor',
                      help="form factor to be computed (must be 'perp' or " +
                           "'para')")
    args.add_argument('datatype',
                      help="type of data (must be 'bs', 'jk', 'raw', " +
                           "'raw2bs' or 'raw2jk')")
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
    """Parses command-line arguments args of script run from directory maindir:
    applies tests, alters formats if necessary, adds arguments relating to data
    sizes, and saves important arguments to semileptonic/settings/fit.py."""
    if args.decayname not in ['B2K', 'B2pi']:
        raise ValueError('invalid decay name')
    if args.formfactor not in ['para', 'perp']:
        raise ValueError('invalid form factor')
    if args.datatype not in ['bs', 'jk', 'raw2bs', 'raw2jk']:
        raise ValueError('invalid data type')
    if type(args.constrained) is not bool:
        raise ValueError('constrained must be boolean')
    if type(args.correlated) is not bool:
        raise ValueError('correlated must be boolean')
    if args.fitlength is not None:
        args.fitlength = [float(n) for n in args.fitlength.split(',')]
        if len(args.fitlength) != 3:
            raise ValueError('invalid length of fit line (must be comma-' +
                             'separated list entered as min,max,numpoints)')
    else:
        args.fitlength = 'full'
    if args.outputdir is None:
        args.outputdir = './chifit_' + args.decayname + '_' + args.formfactor
    if (args.include is not None) and (args.exclude is not None):
        raise ValueError('cannot specify both include and exclude experiment ' +
                         'lists')
    args.nexperiments_source = len(np.loadtxt(args.inputsource))
    args.nsamples_source = int(len(np.loadtxt(args.datasource)) /
                               args.nexperiments_source)
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
        raise ValueError('number of experiments must be multiple of ensemble' +
                         ' size')
    if args.nsamples is not None:
        args.nsamples = int(args.nsamples)
        if (args.nsamples <= 0) or (args.nsamples > args.nsamples_source):
            raise ValueError('fits require 0 < number of samples <= number ' +
                             'of samples in data source')
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
    """Sequentially converts an array to a dictionary with keys."""
    dictionary = {}
    index = 0
    for key in keys:
        dictionary[key] = array[index]
        index += 1
    return dictionary


def data(source, xpmtlist, nexperiments_source, nsamples, nsamples_source,
         usecols=(2,)):
    """Reads in data from columns usecols in source according to experiment list
    xpmtlist, sampling only the desired number of samples nsamples from the
    original data (whose shape is (nexperiments_source, nsamples_source))."""
    data = np.loadtxt(source, usecols=usecols)
    data = data.reshape((nexperiments_source, nsamples_source))
    return data[xpmtlist, :nsamples]


def inputs(source, xpmtlist):
    """Reads in inputs from source according to experiment list xpmtlist."""
    source = np.loadtxt(source)
    inputs = np.array([{'xpmt': xpmt} for xpmt in xpmtlist])
    for i, xpmt in enumerate(xpmtlist):
        inputs[i]['a_fm'] = float(source[xpmt][0])
        inputs[i]['ml_val'] = float(source[xpmt][1])
        inputs[i]['mh_val'] = float(source[xpmt][2])
        inputs[i]['E'] = float(source[xpmt][3])
        inputs[i]['a'] = float(source[xpmt][4])
        inputs[i]['ml_sea'] = float(source[xpmt][5])
        inputs[i]['mh_sea'] = float(source[xpmt][6])
        inputs[i]['ff_avg'] = float(source[xpmt][7])
        inputs[i]['alpha'] = float(source[xpmt][8])
    return inputs


def params():
    """Reads a priori parameters for use with fitter; must be called after
    readers.args()."""
    try:
        from settings.fit import constrained
    except ImportError:
        raise('fit parameters may be read only after command-line arguments ' +
              'have been read')
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
    """Loads array of parameters from pickled source and returns them as a
    dictionary with keys from text source. Sorting of keys should match that of
    the original source parameters."""
    stdout.write('\nReading results from {}\n'.format(
                            os.path.dirname(os.path.realpath(source + '.txt'))))
    dictkeys = np.loadtxt(source + '.txt', usecols=usecols, delimiter=delimiter,
                          dtype=str)
    dictkeys = sorted([key.strip() for key in dictkeys])
    return array2dict(gvar(pickle.load(open(source + '.p', 'rb'))), dictkeys)

