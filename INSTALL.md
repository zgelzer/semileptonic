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

To Install:
-----------

1. Install **python 2.7**.
  + This guide recommends that you use
    [**Anaconda**](http://http://continuum.io/downloads). Anaconda is a free
    python distribution that has many popular packages for scientific analyses.
    It may be easily updated and expanded since it installs to the user's home
    directory. This also allows for simple management of package versions.
2. Install pertinent python libraries: **gvar**, **lsqfit**, **numpy**,
**matplotlib**.
  + If using Anaconda, numpy and matplotlib are already present.
  + To install any package with Anaconda, use `$ pip install {package_name}`
    from the terminal.
3. Install this **semileptonic** project.
  + *(Option a.)* Download this project as a zip; extract to desired directory.
  + *(Option b.)* Install **git**; clone this project using `$ git clone
    https://github.com/zgelzer/semileptonic` in desired directory.
4. *(Optional)* Examine the provided examples.
  + Consult **examples/README.md**.

---

*(Optional)* To Uninstall:
--------------------------

1. From the parent directory of **semileptonic**, run `$ rm -rf semileptonic`.
2. *(Optional)* From the parent directory of your **Anaconda** distribution,
run `$ rm -rf anaconda`.
