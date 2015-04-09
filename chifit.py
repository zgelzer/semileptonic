"""Performs chiral fit to lattice QCD form factors."""


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


from datetime import datetime as dt
from fileIOs import readers as read
from matplotlib import pyplot as plt
from sys import stdout
import os


def main():
    """Runs chiral fit and saves results
    (default save directory: ./results_{decayname}_{formfactor})."""
    
    #~ datetime may be set to False to disable date/time suffixes in results. ~#
    datetime = dt.now().strftime('%Y%m%d-%H%M')
    workdir = os.getcwd()
    scriptdir = os.path.dirname(os.path.realpath(__file__))

    #~ Read and parse command-line arguments. ~#
    args = read.args(scriptdir)

    #~ Import other semileptonic submodules (may be done only after having read
    #  arguments). ~#
    from calculators import chilogs, stats
    from fileIOs import writers as write
    from fitters import lsq as fitlsq

    #~ Move to data directory; read inputs and data; move back to working
    #  directory. ~#
    os.chdir(args.datadir)
    inputs = read.inputs(args.inputsource, args.xpmtlist)
    data = read.data(args.datasource, args.xpmtlist, args.nexperiments_source,
                     args.nsamples, args.nsamples_source)
    os.chdir(workdir)

    #~ Read a priori fit parameters (or their initial values) from
    #  settings/params.py; see settings/params.py for editing instructions. ~#
    params_apriori, params_initval = read.params()

    #~ Obtain chi^2, degrees of freedom, and p-value Q from one-time least
    #  squares fit. ~#
    _, chi2, dof, Q = fitlsq.one(inputs, data, priors=params_apriori,
                                 p0s=params_initval)

    #~ Create output directory if it does not exist. ~#
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)

    #~ If fit parameter results are supplied, load them. ~#
    if args.load is not None:

        #~ Load fit parameter results. ~#
        params_result = read.results(args.load)

        #~ Move to output directory; write results to stdout only. ~#
        os.chdir(args.outputdir)
        write.results(params_result, chi2, dof, Q, fileouts=False)

    #~ If fit parameter results are not supplied, run least squares fits. ~#
    else:

        #~ Run least squares fit for all bootstrap/jackknife samples. ~#
        stdout.write('\nRunning chiral fit...\n')
        params_result, cvals, cov = fitlsq.all(inputs, data,
                                               priors=params_apriori,
                                               p0s=params_initval)

        #~ Move to output directory; write resulting fit parameters
        #  params_result to 'result.p'; write all results to stdout and
        #  'result.txt'. ~#
        os.chdir(args.outputdir)
        write.params(cvals, cov, datetime=datetime)
        write.results(params_result, chi2, dof, Q, datetime=datetime)

    #~ If plot argument is supplied, plot pertinent data and fits. ~#
    if args.plot:

        #~ Initialize figure for receiving all future plots. ~#
        fig = plt.figure()

        #~ Plot input data with error bars; write values to 'result.dat'. ~#
        write.plot_data(inputs, data, datetime=datetime)

        #~ Plot averages of fit results for each ensemble; write values to
        #  'result_fits.dat'. ~#
        write.plot_fitavgs(inputs, params_result, args.fitlength,
                           datetime=datetime)

        #~ Plot continuum extrapolation with error bars, gray line color, and no
        #  transparency in error fills; write values to 'result_fit.dat'. ~#
        write.plot_fit('continuum', inputs, params_result, args.fitlength,
                       color='0.5', label='continuum', datetime=datetime)

        #~ Add axis labels and place small legend in upper right corner. ~#
        write.plot_labels(legendloc='upper right', legendsize='8')

        #~ Remove extraneous white space from borders; adjust fit labels so that
        #  they fit within borders. ~#
        fig.tight_layout()

        #~ Save plots to 'result.pdf'. ~#
        write.plot_save(datetime=datetime)

        #~ Display plots (may be done only after having saved plots). ~#
        plt.show()


if __name__ == '__main__':
    main()

