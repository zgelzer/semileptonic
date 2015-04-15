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
Defines lattice QCD constants and internal tolerances.
--------------------------------------------------------------------------------
Definitions
-----------
r1 : float
    Physical value of lattice scale r_1, defined [1] by (r_1^2 * F(r_1) = 1.0),
    where F(r) is the static force between quarks; taken to be 0.31174 fm [2].
    Determined by requiring continuum limit of pion decay constant f_pi at
    physical quark masses to match that of PDG [3].
hbarc : float
    Conversion constant (hbar * c) that equates mass or energy with inverse
    length or inverse time; taken from PDG as 197.327 MeV fm [3].
Lambda : float
    Renormalization scale characteristic of heavy-quark expansion; taken to be
    1000 MeV to provide liberal widths for chiral log terms.
fpi : float
    Pion decay constant f_pi that characterizes energy scale of chiral symmetry
    breaking; estimated from PDG as approximately 131 MeV [3].
gpi : float
    Lattice determination [4] of central value of leading order axial coupling
    constant g_pi (i.e., between pions and heavy-light mesons). Hereafter, g_pi
    is referred to simply as the heavy-light coupling constant.
Delta_B : float
    Difference between masses of B*_s and B^0 mesons for B-->K decay, or that
    of B* and B^0 for B-->pi decay.
tastemults : dict of ints
    Multiplicity of tastes of staggered fermions, with tastes defined as 'P' for
    pseudoscalar, 'A' for axial vector, 'T' for tensor, 'V' for vector, and 'I'
    for singlet.
r1_continuum : float
    Continuum value of lattice scale r_1 in lattice units (r_1 / a); taken to be
    3.744241 [2].
ml_continuum : float
    Continuum value of light quark mass m_l in lattice units (a * m_l); taken to
    be 0.0009646 [2].
mh_continuum : float
    Continuum value of heavy quark mass m_h in lattice units (a * m_h); taken to
    be 0.02645 [2].
a_atol : float
    Absolute tolerance used to deterimine equality of two lattice spacings in
    fermi; taken to be 1e-06 fm, since difference between successive lattice
    spacings is typically of order 1e-02 fm. 
nearzero : float
    Number used to represent essentially zero, where exactly zero would be
    problematic in practice; taken to be 1e-10.
--------------------------------------------------------------------------------
Notes
-----
+ Currently, nearzero is used exclusively to calculate hairpin parameters in
  calculators.chilogs.fcns.deltas_calc. Since hairpin parameters reported in
  calculators.chilogs.fcns.deltas were generated with deltas_calc at (nearzero =
  1e-10), user would have to manually run deltas_calc at new nearzero value and
  then update deltas to reflect changes. As such, nearzero is defined here only
  to provide consistency and transparency for this kind of "magic number".
--------------------------------------------------------------------------------
References
----------
[1] C. Bernard, et al., "The static quark potential in three flavor QCD", Phys.
    Rev. D62, 034503 (2000), [arXiv:0002028 [hep-lat]].
[2] C. Bernard, et al. (Fermilab Lattice and MILC Collaborations),
    "UPDATED mass-independent r1 values",
    <http://physics.wustl.edu/~cb/Fermilab-MILC/secure/>, accessed 2015-03-30.
[3] K. Olive, et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014).
[4] W. Detmold, C.-J. D. Lin, S. Meinel, "Calculation of the heavy-hadron axial
    couplings g_1, g_2, and g_3 using lattice QCD", Phys. Rev. D 85, 114508
    (2012) [arXiv:1203.3378 [hep-lat]].
--------------------------------------------------------------------------------
"""


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
r1_continuum = 3.744241
ml_continuum = r1_continuum * 0.0009646
mh_continuum = r1_continuum * 0.02645
a_atol = 1e-06
nearzero = 1e-10

