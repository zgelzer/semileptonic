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

    $ python chifit.py --help

Important Notices
-----------------

+ [chifit.py](../chifit.py) can currently handle only bootstrapped data.
  + Modules in [calculators/stats](../calculators/stats) are prepped for
    conversions from raw to bootstrapped data, but this is not currently
    implemented in `chifit.main()`.
  + Future releases will interact with jackknifed data.
+ [lsq.py](../fitters/lsq.py) currently performs each least squares fit in
  series.
  + Older systems may take ~60 seconds to perform ~1,000 fits.
  + Future releases will incorporate parallel least squares fitting.

Usage Suggestions
-----------------

+ **transfer** pertinent files to current directory:
    ```
    $ cp examples/X.dat examples/Y.dat .
    ```

+ **run fit** for *f_perp* of *B* to *K*:
    ```
    $ python chifit.py B2K perp bs --constrained --correlated
    ```

+ **load results** from previous run (at datetime=`dt`) and **plot**:
    ```
    $ python chifit.py B2K perp bs -c -C --load chifit_B2K_perp/result_{dt} --plot
    ```

+ **load results** from previous run (at datetime=`dt`) and plot with custom
**fit length** (`min,max,numpoints`):
    ```
    $ python chifit.py B2K perp bs -c -C -L chifit_B2K_perp/result_{dt} -p --length 0.5,4.0,100
    ```

+ **include** only first 100 bootstrap **samples** (as quick test), run fit,
and plot:
    ```
    $ python chifit.py B2K perp bs -c -C -p -n 100
    ```

+ **include** only first and last ensembles of **experiments**, run fit, and
plot:
    ```
    $ python chifit.py B2K perp bs -c -C -p -i 0,1,2,27,28,29
    ```

+ **exclude** high-energy **experiments** from each ensemble (effectively
setting **ensemble size** to two), run fit, and plot:
    ```
    $ python chifit.py B2K perp bs -c -C -p -x 2,5,8,11,14,17,20,23,26,29 -e 2
    ```

List of Files
-------------

All results are from single run of [chifit.py](../chifit.py) for *f_perp* of
bootstrapped *B* to *K* data.

+ [README.md](README.md): Details suggested runs of `chifit.py` and explains the
  included examples.
+ [X.dat](X.dat): Inputs for each experiment, organized by ensemble.
+ [Y.dat](Y.dat): Bootstrap samples of *B* to *K* form factors, with
  organization matching that of `X.dat`.
+ [result.dat](result.dat): Plot results of `Y.dat` vs. `X.dat` with error bars.
+ [result.p](result.p): Fit parameter results from previous run stored in
  pickled binary format.
+ [result.pdf](result.pdf): Plot results from previous run.
+ [result.txt](result.txt): Fit results from stdout of previous run.
+ [result_fit.dat](result_fit.dat): Plot results of continuum fit with error
  bars.
+ [result_fits.dat](result_fits.dat): Plot results of ensemble fit averages.

