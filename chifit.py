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
    datetime = dt.now().strftime('%Y%m%d-%H%M')
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    args = read.args(scriptdir)
    from calculators import chilogs, stats
    from fileIOs import writers as write
    from fitters import lsq as fitlsq
    workdir = os.getcwd()
    os.chdir(args.datadir)
    inputs = read.inputs(args.inputsource, args.xpmtlist)
    data = read.data(args.datasource, args.xpmtlist, args.nexperiments_source,
                     args.nsamples, args.nsamples_source)
    params_apriori, params_initval = read.params()
    os.chdir(workdir)
    _, chi2, dof, Q = fitlsq.one(inputs, data, priors=params_apriori,
                                 p0s=params_initval)
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)
    if args.load is not None:
        params_result = read.results(args.load)
        os.chdir(args.outputdir)
        write.results(params_result, chi2, dof, Q, fileouts=False)
    else:
        stdout.write('\nRunning chiral fit...\n')
        params_result, cvals, cov = fitlsq.all(inputs, data,
                                               priors=params_apriori,
                                               p0s=params_initval)
        os.chdir(args.outputdir)
        write.params(cvals, cov, datetime=datetime)
        write.results(params_result, chi2, dof, Q, datetime=datetime)
    if args.plot:
        fig = plt.figure()
        write.plot_data(inputs, data)
        write.plot_fitavgs(inputs, params_result, args.fitlength)
        write.plot_fit('continuum', inputs, params_result, args.fitlength,
                       color='0.5', label='continuum')
        write.plot_labels(legendloc='upper right', legendsize='8')
        fig.tight_layout()
        write.plot_save(datetime=datetime)
        plt.show()


if __name__ == '__main__':
    main()

