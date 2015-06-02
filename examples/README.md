<!---
  Created by Zechariah Gelzer (University of Iowa) on 2015-03-30.
  Copyright (C) 2015 Zechariah Gelzer.
 
  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or any later version (see
  <http://www.gnu.org/licenses/>).
 
  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
-->

Built-In Documentation
----------------------

    $ ./chifit.py --help

Important Notices
-----------------

+ [chifit.py](../chifit.py) can currently handle only bootstrapped data.
  + Future releases may interact with jackknifed data.
+ [lsq.py](../fitters/lsq.py) currently performs each least squares fit in
  series.
  + Older systems may take ~75 seconds to perform ~1,000 fits.
  + Future releases may incorporate parallel least squares fitting.

Usage Suggestions
-----------------

The following steps detail potentially interesting executions of `chifit.py` and
are presented in a sensible order.

1. **Transfer** datasets to current directory. This bypasses the need to specify
a data directory via the `--dir` flag (which is reasonable when executing runs
for a single form factor).
    ```
    $ cp examples/B2K/perp/X.dat examples/B2K/perp/Y.dat .
    ```

2. **Run chiral fit** for *f_perpendicular* of *B* to *K* in the hard-kaon limit
with all constrained NNLO fit parameters (and with correlations), **save**
all outputs with the prefix "chiral_NNLO", and **plot** the results. This also
demonstrates that one may interchangeably use long or short versions of
command-line flags according to one's own preferences.
    ```
    $ ./chifit.py -d B2K -f perp --constrained --correlated --hard \
    --save chiral_NNLO --plot
    ```

3. **Exclude** high-energy simulations from each ensemble, effectively setting
setting **ensemble size** to two, and run chiral fit.
    ```
    $ ./chifit.py -d B2K -f perp -c -C -H -x 2,5,8,11,14,17,20,23,26,29 -e 2 \
    -s chiral_NNLO_sanshighE
    ```

4. **Load results** from previous run, plot (now with all simulations included),
and use custom **fit length** (`min,max,numpoints`). The loading process imports
fit parameters from `{filename}.p` and (most) fit settings from
`{filename}.txt`; it does not, however, import the values in `{filename}.txt`
that pertain to ensemble size and number of simulations. This enables datasets
to be fully controlled during each run. One may see that the fit parameter
results from low-energy simulations accurately predict the high-energy
simulations.
    ```
    $ ./chifit.py --load B2K/perp/chiral_NNLO_sanshighE \
    -s chiral_NNLO_sanshighE -p --length 0.8,2.0,100
    ```

5. Remove NNLO **fit parameters** (see [params.py](../settings/params.py) for
more information), run chiral fit, and plot in the same range as that of Step 4.
One may see that the chiral fit parameters have roughly stabilized as early as
NLO (at least for these particular datasets).
    ```
    $ sed -i "/width_NNLO/s/^/#/g" settings/params.py
    $ ./chifit.py -d B2K -f perp -c -C -H -s chiral_NLO -p -l 0.8,2.0,100
    ```

File Tree
---------

+ [B2K](B2K): Contains examples for *B* to *K* lattice QCD form factors.
  + [para](B2K/para): Contains example inputs for and outputs from lattice QCD
    simulations of *f_parallel*. Contains results from `chifit.py` executed in
    the hard-kaon limit with all constrained NNLO fit parameters:
    ```
    $ ./chifit.py -d B2K -f para -D examples/B2K/para/ -e 4 -c -C -H \
    -s chiral_NNLO -p -l 0.8,2.0,100
    ```
    + [X.dat](B2K/para/X.dat): Inputs for each simulation, organized by
      ensemble.
    + [Y.dat](B2K/para/Y.dat): Bootstrap samples of outputs from each
      simulation, with organization matching that of `X.dat`.
    + [chiral_NNLO.dat](B2K/para/chiral_NNLO.dat): `Y.dat` vs. `X.dat`, with
      error bars.
    + [chiral_NNLO.p](B2K/para/chiral_NNLO.p): Results for fit parameters,
      stored in pickled binary format.
    + [chiral_NNLO.pdf](B2K/para/chiral_NNLO.pdf): Plot of data and fit results.
    + [chiral_NNLO.txt](B2K/para/chiral_NNLO.txt): Fit settings and results from
      standard output.
    + [chiral_NNLO_fit.dat](B2K/para/chiral_NNLO_fit.dat): Results for continuum
      fit, with error bars.
    + [chiral_NNLO_fits.dat](B2K/para/chiral_NNLO_fits.dat): Results for
      ensemble fits, without error bars.
  + [perp](B2K/perp): Contains example inputs for and outputs from lattice QCD
    simulations of *f_perpendicular*. Contains results from `chifit.py` executed
    in the hard-kaon limit with all constrained NNLO fit parameters:
    ```
    $ ./chifit.py -d B2K -f perp -D examples/B2K/perp/ -e 3 -c -C -H \
    -s chiral_NNLO -p -l 0.8,2.0,100
    ```
    + [X.dat](B2K/perp/X.dat): Inputs for each simulation, organized by
      ensemble.
    + [Y.dat](B2K/perp/Y.dat): Bootstrap samples of outputs from each
      simulation, with organization matching that of `X.dat`.
    + [chiral_NNLO.dat](B2K/perp/chiral_NNLO.dat): `Y.dat` vs. `X.dat`, with
      error bars.
    + [chiral_NNLO.p](B2K/perp/chiral_NNLO.p): Results for fit parameters,
      stored in pickled binary format.
    + [chiral_NNLO.pdf](B2K/perp/chiral_NNLO.pdf): Plot of data and fit results.
    + [chiral_NNLO.txt](B2K/perp/chiral_NNLO.txt): Fit settings and results from
      `stdout`.
    + [chiral_NNLO_fit.dat](B2K/perp/chiral_NNLO_fit.dat): Results for continuum
      fit, with error bars.
    + [chiral_NNLO_fits.dat](B2K/perp/chiral_NNLO_fits.dat): Results for
      ensemble fits, without error bars.
+ [README.md](README.md): Details some suggested runs of `chifit.py` and
  explains the included examples.
