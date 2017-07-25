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

> + *Possible changes in a future release:*
>   + *incorporate B_s to K decays*
>   + *parallelize least squares fitter*

---

+ Version 0.4.1
  + updated plots for matplotlib v2

---

+ Version 0.4.0
  + added light quark and gluon discretization effects
  + added heavy quark discretization effects
  + added pole option for *f_para*
  + changed poles to true positions

+ Version 0.3.0
  + added *f_0*, *f_+* functionality and examples
  + added *f_T* functionality
  + added *f_para*, *f_perp* examples
  + added constants
  + added energy relation functions
  + changed plot saving and output structure
  + changed loading of previous results
  + changed reading of input sources and detection of ensemble size

+ Version 0.2.2
  + fixed *g_pi* fit usage

+ Version 0.2.1
  + changed default x axis in plots

+ Version 0.2.0
  + added comprehensive info to all function docstrings
  + added ability to save values from plot data/fits to formatted `.dat` files
  + added example `result*.dat` files
  + changed handling of `.dat` headers
  + changed default NNLO fit parameter widths
  + changed default *f_pi* and *Lambda* constants
  + changed default bootstrap errors
  + removed handling of raw data

+ Version 0.1.1
  + added example plot
  + added run suggestion for setting fit line lengths
  + fixed plot save name in `writers.py`

+ Version 0.1.0
  + initial release of beta version
