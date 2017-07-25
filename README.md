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

[Semileptonic](https://github.com/zgelzer/semileptonic)
=======================================================

Utilities for analyzing lattice-QCD correlators in semileptonic decays.

---

+ To install, see [INSTALL.md](INSTALL.md).
+ To view the history of version changes, see [CHANGES.md](CHANGES.md).
+ To view examples and usage suggestions, see
  [examples/README.md](examples/README.md).

---

File Tree
---------

+ [.gitattributes](.gitattributes): Lists attributes for given paths.
+ [.gitignore](.gitignore): Lists files for git to ignore during commits.
+ [CHANGES.md](CHANGES.md): Lists history of version changes in chronologically
  descending order.
+ [INSTALL.md](INSTALL.md): Details how to install or uninstall this project.
+ [README.md](README.md): Provides an overview of this project.
+ [chifit.py](chifit.py): Main program. Performs chiral fit to lattice QCD form
  factors.
+ [calculators](calculators): Contains modules for calculating chiral logs and
  statistics.
  + [fcns.py](calculators/fcns.py): Defines functions for use in relating
    particle energies.
  + [chilogs](calculators/chilogs): Contains modules for calculating chiral logs
    that are separated by decay channel, with their basic functions stored in a
    communal module.
    + [B2K.py](calculators/chilogs/B2K.py): Calculate pole and loop
      contributions to *B* to *K* form factors.
    + [B2pi.py](calculators/chilogs/B2pi.py): Calculate pole and loop
      contributions to *B* to *pi* form factors.
    + [fcns.py](calculators/chilogs/fcns.py): Defines basic functions for use in
      calculating terms related to chiral logs. Imported by other modules in
      this directory.
  + [stats](calculators/stats): Contains modules for resampling data and for
    reporting statistical estimates of resampled data.
    + [bootstrap.py](calculators/stats/bootstrap.py): Defines functions for
      bootstrapping raw data and reporting averages and errors of bootstrapped
      data.
+ [examples](examples): Contains example inputs for and outputs from lattice QCD
  simulations of *f_parallel* and *f_perpendicular* for *B* to *K* semileptonic
  decays. Contains example results from `chifit.py`.
  + [README.md](examples/README.md): Details some suggested runs of `chifit.py`
    and explains the included examples.
+ [fileIOs](fileIOs): Contains modules for reading inputs and writing outputs.
  + [readers.py](fileIOs/readers.py): Defines functions for reading command-line
    arguments, inputs, data, initial fit parameters, and previous results.
  + [writers.py](fileIOs/writers.py): Defines functions for writing and plotting
    data.
+ [fitters](fitters): Contains modules for fitters and their functions.
  + [chiral.py](fitters/chiral.py): Defines chiral fit functions for *B* decay
    form factors.
  + [lsq.py](fitters/lsq.py): Performs least squares fit(s) to data using Prof.
    G. P. Lepage's [lsqfit.nonlinear_fit](https://github.com/gplepage/lsqfit).
+ [settings](settings) : Contains modules for definitions of settings,
  constants, and parameters.
  + [constants.py](settings/constants.py): Defines lattice QCD constants. User
    should edit if necessary.
  + [params.py](settings/params.py): Defines chiral fit parameters. User should
    edit if necessary.

*Note that all [\_\_init\_\_.py](__init__.py) files render their given
directories importable by python. This is required so that the various python
files may communicate with one another.*

---

> Special thanks is given to Dr. Ran Zhou of Fermilab for providing example
> datasets, assisting with debugging, and offering invaluable advice.
