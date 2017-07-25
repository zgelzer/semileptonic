#!/usr/bin/env python
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
Performs chiral fit to lattice QCD form factors.
    --------        --------        --------        --------        --------    
Meant to be executed from terminal. Run './chifit.py --help' for more info.
--------------------------------------------------------------------------------
Definitions
-----------
main : function
    Runs or loads chiral fits and saves results.
run_combined : function
    Loads chiral fits for combined form factor and saves results.
run_single : function
    Runs or loads chiral fits for single form factor and saves results.
--------------------------------------------------------------------------------
Notes
-----
+ In all file and function docstrings, 'float' describes a 64-bit float whose
  type is either float (built-in) or numpy.float64. No distinction is made
  between these types of floats, since they are identical on all IEEE-754
  compliant systems.
--------------------------------------------------------------------------------
"""


from fileIOs import readers as read
from matplotlib import pyplot as plt
from sys import stdout
import os


def main():
    """
    ----------------------------------------------------------------------------
    Runs or loads chiral fits and saves results.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Default save directory: './{decay name}/{form factor}'
    Default save name: {date/time as 'YYYYMMDD-hhmm'}
    ----------------------------------------------------------------------------
    Results
    -------
    files
        Various files, as determined by nature of form factor (combined or
        single), with root name {save name} in directory {save directory}.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    readers : module, from fileIOs, as read
    run_combined : function
    run_single : function
    ----------------------------------------------------------------------------
    Notes
    -----
    + Meant to be run automatically when script is called from terminal.
    + See functions run_combined, run_single.
    ----------------------------------------------------------------------------
    """

    #~ Read and parse command-line arguments (must be done before importing
    #  other semileptonic submodules). ~#
    args = read.args()

    #~ Use classic style for plots. ~#
    plt.style.use('classic')

    #~ Run combined or single fits, depending on form factor. ~#
    if (args.formfactor == 'scalar') or (args.formfactor == 'vector'):
        run_combined(args)
    else:
        run_single(args)


def run_combined(args):
    """
    ----------------------------------------------------------------------------
    Loads chiral fits for combined form factor (scalar or vector) and saves
    results.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Default save directory: './{decay name}/{form factor}'
    Default save name: {date/time as 'YYYYMMDD-hhmm'}
    ----------------------------------------------------------------------------
    Results
    -------
    {savename}.pdf : file
        Aggregate of all plots.
    {savename}.txt : file
        Fit settings.
    {savename}_fit.dat : file
        Values from plot of continuum fit with errors.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    os : module
    pyplot : module, from matplotlib, as plt
    readers : module, from fileIOs, as read
    stdout : file, from sys
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function main.
    ----------------------------------------------------------------------------
    """

    #~ Import other semileptonic submodules. ~#
    from fileIOs import writers as write

    #~ Load fit parameter results for both parallel and perpendicular form
    #  factors. ~#
    params_para = read.results(args.loadpara)
    params_perp = read.results(args.loadperp)

    #~ Move to output directory; write fit settings to stdout and
    #  '{savename}.txt'. ~#
    os.chdir(args.outputdir)
    write.results(savename=args.savename)

    #~ If plot argument is supplied, plot continuum extrapolation. ~#
    if args.plot:

        #~ Initialize figure for receiving all future plots. ~#
        figure = plt.figure()

        #~ Plot continuum fit for combined form factor, with default color and
        #  transparency in error fills; write values to '{savename}.dat'. ~#
        write.plot_fitcombo(params_para, params_perp, linelength=args.fitlength,
                            savename=args.savename)

        #~ Add axis labels; set bounds of x-axis to that of custom fit line (if
        #  specified). ~#
        write.plot_labels(legendloc=None, xlims=args.fitlength)

        #~ Remove extraneous white space from borders; adjust fit labels so that
        #  they fit within borders; save plots to '{savename}.pdf' ~#
        figure.tight_layout()
        plt.savefig(args.savename + '.pdf')

        #~ Display plots (may be done only after having saved plots). ~#
        plt.show()

    #~ If plot argument is not supplied, warn user. ~#
    else:
        stdout.write("\nIf you would like to plot, use '-p' or '--plot'.\n\n")


def run_single(args):
    """
    ----------------------------------------------------------------------------
    Runs or loads chiral fits for single form factor (parallel, perpendicular,
    or tensor) and saves results.
        ----    ----    ----    ----    ----    ----    ----    ----    ----    
    Default save directory: './{decay name}/{form factor}'
    Default save name: {date/time as 'YYYYMMDD-hhmm'}
    ----------------------------------------------------------------------------
    Results
    -------
    {savename}.dat : file
        Values from plot of input data with error bars.
    {savename}.p : file
        Fit parameter results in pickled binary format.
    {savename}.pdf : file
        Aggregate of all plots.
    {savename}.txt : file
        Fit settings, and results for fit qualities and parameters.
    {savename}_fit.dat : file
        Values from plot of continuum fit with errors.
    {savename}_fits.dat : file
        Values from plot of ensemble fit averages.
    ----------------------------------------------------------------------------
    Requirements
    ------------
    os : module
    pyplot : module, from matplotlib, as plt
    readers : module, from fileIOs, as read
    stdout : file, from sys
    ----------------------------------------------------------------------------
    Notes
    -----
    + See function main.
    ----------------------------------------------------------------------------
    """

    #~ Import other semileptonic submodules. ~#
    from fileIOs import writers as write
    from fitters import lsq as fitlsq

    #~ Move to data directory; read inputs and data; move back to working
    #  directory. ~#
    os.chdir(args.datadir)
    inputs = read.inputs(args.inputsource, args.xpmtlist)
    data = read.data(args.datasource, args.xpmtlist, args.nexperiments_source,
                     args.nsamples, args.nsamples_source)
    os.chdir(args.workdir)

    #~ Read a priori fit parameters (or their initial values) from
    #  settings/params.py; see settings/params.py for editing instructions. ~#
    params_apriori, params_initval = read.params()

    #~ Obtain chi^2, degrees of freedom, and p-value from one-time least squares
    #  fit (discarding the resulting fit parameters). ~#
    _, chi2, dof, p = fitlsq.one(inputs, data, inits=params_initval,
                                 priors=params_apriori)

    #~ If fit parameter results are supplied, load them. ~#
    if args.load is not None:

        #~ Load and write fit parameter results; move to output directory. ~#
        stdout.write('\nReading results from {}.{{p, txt}}\n'.format(
                                                    os.path.abspath(args.load)))
        stdout.write(open(args.load + '.txt').read())
        params_result = read.results(args.load)
        os.chdir(args.outputdir)

    #~ If fit parameter results are not supplied, run least squares fits. ~#
    else:

        #~ Run least squares fit for all bootstrap/jackknife samples. ~#
        stdout.write('\nRunning chiral fit...\n')
        params_result, cvals, cov = fitlsq.all(inputs, data,
                                               inits=params_initval,
                                               priors=params_apriori)

        #~ Move to output directory; write resulting fit parameters
        #  params_result to '{savename}.p'; write all results to stdout and
        #  '{savename}.txt'. ~#
        os.chdir(args.outputdir)
        write.params(cvals, cov, savename=args.savename)
        write.results(chi2=chi2, dof=dof, p=p, params=params_result,
                      savename=args.savename)

    #~ If plot argument is supplied, plot pertinent data and fits. ~#
    if args.plot:

        #~ Initialize figure for receiving all future plots. ~#
        figure = plt.figure()

        #~ Plot input data with error bars; write values to '{savename}.dat'. ~#
        write.plot_data(inputs, data, savename=args.savename)

        #~ Plot averages of fit results for each ensemble; write values to
        #  '{savename}_fits.dat'. ~#
        write.plot_fitavgs(inputs, params_result, linelength=args.fitlength,
                           savename=args.savename)

        #~ Plot continuum extrapolation with dark gray line color and 50%
        #  transparency in error fills; write values to '{savename}_fit.dat'. ~#
        write.plot_fit(inputs, params_result, 'continuum', alphafill=0.5,
                       color='0.25', label='continuum',
                       linelength=args.fitlength, savename=args.savename)

        #~ Add axis labels; create title with chi^2, degrees of freedom, and
        #  p-value from one-time fit; place small legend in upper right corner;
        #  set bounds of x-axis to that of custom fit line (if specified). ~#
        write.plot_labels(chi2=chi2, dof=dof, legendloc='upper right',
                          legendsize='8', p=p, xlims=args.fitlength)

        #~ Remove extraneous white space from borders; adjust fit labels so that
        #  they fit within borders; save plots to '{savename}.pdf' ~#
        figure.tight_layout()
        plt.savefig(args.savename + '.pdf')

        #~ Display plots (may be done only after having saved plots). ~#
        plt.show()


#~ If script is being called from terminal, run main function. ~#
if __name__ == '__main__':
    main()

