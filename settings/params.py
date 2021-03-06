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


"""
--------------------------------------------------------------------------------
Defines chiral fit parameters.
    --------        --------        --------        --------        --------    
Each parameter is a Gaussian variable gvar defined with central value and width
as gvar(central_value, width). To turn off a parameter, simply comment out its
line. If not using a constrained fit, all widths will be ignored.
--------------------------------------------------------------------------------
Definitions
-----------
params : dict of gvar.GVars
    Dictionary of chiral fit parameters as Gaussian variables.
width_NLO : float
    Default width of next-to-leading-order (NLO) fit parameters.
width_NNLO : float
    Default width of next-to-next-to-leading-order (NNLO) fit parameters.
--------------------------------------------------------------------------------
Notes
-----
+ Fit parameters are implemented in chiral fit function (see fitters.chiral).
+ Fit parameter names match those from heavy-light meson staggered chiral
  perturbation theory (HMSchiPT) [1].
    > C0 : Parameterizes leading-order (LO) term (1 + {logs}).
    > CE : Parameterizes NLO term chi_E.
    > CE2 : Parameterizes NLO term chi_E^2.
    > CE3 : Parameterizes NNLO term chi_E^3.
    > CE4 : Parameterizes NNLO term chi_E^4.
    > Ca2 : Parameterizes NLO term chi_{a2}.
    > Ca2E : Parameterizes NNLO term (chi_{a2} * chi_E).
    > Ca2E2 : Parameterizes NNLO term (chi_{a2} * chi_E^2).
    > Ca4 : Parameterizes NNLO term chi_{a2}^2.
    > Ch : Parameterizes NLO term chi_{m_h}.
    > Cl : Parameterizes NLO term chi_{m_l}.
    > Cl2 : Parameterizes NNLO term chi_{m_l}^2.
    > ClE : Parameterizes NNLO term (chi_{m_l} * chi_E).
    > ClE2 : Parameterizes NNLO term (chi_{m_l} * chi_E^2).
    > Cla2 : Parameterizes NNLO term (chi_{m_l} * chi_{a2}).
    > D : Parameterizes LO discretization errors for light quarks and gluons.
    > D3 : Parameterizes heavy quark discretization errors from LO term f3.
    > DB : Parameterizes heavy quark discretization errors from LO term fB.
    > DE : Parameterizes heavy quark discretization errors from LO term fE.
    > DX : Parameterizes heavy quark discretization errors from LO term fX.
    > DY : Parameterizes heavy quark discretization errors from LO term fY.
    > gpi : Parameterizes LO term g_pi (see settings.constants for central
      value), whose width is set to 0.08 so as to remain consistent with a
      direct lattice calculation [2], yet liberal enough to accomodate others.
+ All terms are dimensionless; see fitters.chiral.fitfcn for complete
  descriptions.
+ All LO and NLO terms are expected to be of order unity.
+ All NNLO terms are expected to be of order zero.
+ All discretization errors are expected to be minimal, with widths equal to the
  square root of the number of times each term appears.
--------------------------------------------------------------------------------
References
----------
[1] J. Bailey, et al. (Fermilab Lattice and MILC Collaborations), "The B --> pi
    l nu semileptonic form factor from three-flavor lattice QCD: A model-
    independent determination of |V(ub)|", Phys. Rev. D 79, 054507 (2009)
    [arXiv:0811.3640 [hep-lat]].
[2] W. Detmold, C.-J. D. Lin, and S. Meinel, "Calculation of the heavy-hadron
    axial couplings g_1, g_2, and g_3 using lattice QCD", Phys. Rev. D 85,
    114508 (2012) [arXiv:1203.3378 [hep-lat]].
--------------------------------------------------------------------------------
"""


from gvar import gvar
from math import sqrt
from settings.constants import gpi
from settings.fit import formfactor


width_NLO  = 2.
width_NNLO = 1.
params = {}
params['C0']    = gvar(0., 2.)
params['CE']    = gvar(0., width_NLO)
params['CE2']   = gvar(0., width_NLO)
params['CE3']   = gvar(0., width_NNLO)
params['CE4']   = gvar(0., width_NNLO)
params['Ca2']   = gvar(0., width_NLO)
params['Ca2E']  = gvar(0., width_NNLO)
params['Ca2E2'] = gvar(0., width_NNLO)
params['Ca4']   = gvar(0., width_NNLO)
params['Ch']    = gvar(0., width_NLO)
params['Cl']    = gvar(0., width_NLO)
params['Cl2']   = gvar(0., width_NNLO)
params['ClE']   = gvar(0., width_NNLO)
params['ClE2']  = gvar(0., width_NNLO)
params['Cla2']  = gvar(0., width_NNLO)
params['D']     = gvar(0., 1.)
params['D3']    = gvar(0., sqrt(2) * sqrt(2) ** (formfactor != 'para'))
params['DB']    = gvar(0., 1.)
params['DE']    = gvar(0., sqrt(2))
params['DX']    = gvar(0., sqrt(2))
params['DY']    = gvar(0., 1.)
params['gpi']   = gvar(gpi, 0.08)

