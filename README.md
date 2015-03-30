Semileptonic
============

Python implementation of chiral fits to lattice QCD form factors for
semileptonic decays of *B* mesons.

---

+ To install, see **INSTALL.md**.
+ To view the change history, see **CHANGES.md**.
+ To view examples and usage suggestions, see **examples/README.md**.

---

File Tree
---------

+ **semileptonic**
  + `chifit.py`: Main program; performs chiral fit to lattice QCD form factors.
  + `CHANGES.md`: Lists change history in chronologically descending order.
  + `constants.py`: Defines lattice QCD constants. User should edit if
    necessary.
  + `INSTALL.md`: Details how to install or uninstall this project.
  + `params.py`: Defines chiral fit parameters. User should edit if necessary.
  + `README.md`: Provides an overview of this project.
  + `.gitignore`: Lists files for git to ignore during commits.
  + **calculators**: Contains modules for calculating chiral logs and
    statistics.
    + **chilogs**: Contains modules for calculating chiral logs that are
      separated by decay channel, with their basic functions stored in a
      communal module.
      + `B2K.py`: Calculate pole and loop contributions to *B* to *K* form
        factors.
      + `B2pi.py`: Calculate pole and loop contributions to *B* to *pi* form
        factors.
      + `fcns.py`: Defines basic functions for use in calculating terms related
        to chiral logs. Imported by other modules in this directory.
    + **stats**: Contains modules for resampling data and for reporting
      statistical estimates of unsampled/resampled data.
      + `bootstrap.py`: Defines functions for bootstrapping raw data and
        averaging bootstrapped data.
      + `raw.py`: Defines functions for averaging raw data.
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
  + **fileIOs**: Contains modules for reading inputs and writing outputs.
    + `readers.py`: Defines functions for reading command-line arguments,
      inputs, data, and initial fit parameters.
    + `writers.py`: Defines functions for writing and plotting data.
  + **fitters**: Contains modules
    + `chiral.py`: Defines chiral fit function for B decay form factors.
    + `lsq.py`: Performs least squares fit(s) to data using Prof. G. Peter
      LePage's `lsqfit.nonlinear_fit`.

*Note that all* `__init__.py` *files render their given directories importable
by python. This is required so that the various python files may communicate
with one another.*


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
