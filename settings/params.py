"""Defines chiral fit parameters. Each parameter is a Gaussian variable gvar
defined as gvar(centralvalue, width). To turn off a parameter, simply comment
out its line. If not using a constrained fit, the widths will be ignored."""


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


from settings.constants import gpi
from gvar import gvar


params = {}
params['C0'] = gvar(0., 2.)
params['Cl'] = gvar(0., 2.)
params['Cl2'] = gvar(0., 2.)
params['CE'] = gvar(0., 2.)
params['CE2'] = gvar(0., 2.)
params['Ca2'] = gvar(0., 2.)
params['Cla2'] = gvar(0., 2.)
params['ClE'] = gvar(0., 2.)
params['ClE2'] = gvar(0., 2.)
params['Ca2E'] = gvar(0., 2.)
params['Ca2E2'] = gvar(0., 2.)
params['CE3'] = gvar(0., 2.)
params['CE4'] = gvar(0., 2.)
params['Ca4'] = gvar(0., 2.)
params['Ch'] = gvar(0., 2.)
params['gpi'] = gvar(gpi, 0.08)

