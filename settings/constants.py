"""Defines lattice QCD constants."""


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


from settings.fit import decayname


r1 = 0.31174
hbarc = 197.327
Lambda = 1000 * r1 / hbarc
fpi = 131 * r1 / hbarc
gpi = 0.45
if decayname == 'B2K':
    Delta_B = (5415.4 - 5279.58) * r1 / hbarc
elif decayname == 'B2pi':
    Delta_B = (5325.2 - 5279.58) * r1 / hbarc
tastemults = {'P': 1, 'A': 4, 'T': 6, 'V': 4, 'I': 1}
a_atol = 1e-06
nearzero = 1e-10
r1_continuum = 3.744241
ml_continuum = r1_continuum * 0.0009646
mh_continuum = r1_continuum * 0.02645

