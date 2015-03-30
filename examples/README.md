Built-In Documentation
----------------------

    $ python chifit.py --help

Important Notices
-----------------

+ `chifit.py` can currently handle only raw or bootstrapped data.
  + `calculators/stats/*.py` are prepped for conversions from raw to
    bootstrapped data, but this is not currently implemented in `main()` of
    `chifit.py`.
  + Future releases will incorporate jackknife resampling of raw data.
+ `lsq.py` currently performs each least squares fit in series.
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

+ **load results** from previous run (at datetime=`dt`) and plot:
    ```
    $ python chifit.py B2K perp bs -c -C --load chifit_B2K_perp/result_{dt} --plot
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

+ **exclude** high-energy **experiments** from each ensemble, effectively
setting ensemble size to 2, run fit, and plot:
    ```
    $ python chifit.py B2K perp bs -c -C -p -x 2,5,8,11,14,17,20,23,26,29 -e 2
    ```

File Tree
---------

+ **examples**: Contains example inputs and outputs for use with `chifit.py`.
  + `result.p`: Fit parameter results from previous run (for *f_perp* of
    bootstrapped *B* to *K* data) stored in pickled binary form.
  + `result.pdf`: Plot from previous run (for *f_perp* of bootstrapped *B* to
    *K* data).
  + `result.txt`: Fit results from stdout of previous run (for *f_perp* of
    bootstrapped *B* to *K* data).
  + `README.md`: Details suggested runs of `chifit.py` and explains the
    included examples.
  + `X.dat`: Inputs for each experiment, organized by ensemble.
  + `Y.dat`: Bootstrapped *B* to *K* form factors, with organization matching
    that of `X.dat`.

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
