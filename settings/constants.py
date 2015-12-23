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
hbarc : float
    Conversion constant (hbar * c) that equates mass or energy with inverse
    length or inverse time; taken from PDG as 197.327 MeV fm [1].
r1 : float
    Physical value of lattice scale r_1, defined [2] by (r_1^2 * F(r_1) = 1.0),
    where F(r) is the static force between quarks; taken to be 0.31174 fm [3].
    Determined by requiring continuum limit of pion decay constant f_pi at
    physical quark masses to match that of PDG [1].
Lambda : float
    Renormalization scale characteristic of heavy-quark expansion; taken to be
    700 MeV to provide liberal widths for chiral log terms [4].
fpi : float
    Pion decay constant f_pi that characterizes energy scale of chiral symmetry
    breaking; taken from PDG as 130.41 MeV [1].
mB : float
    Mass of B^0 meson; taken from PDG as 5279.58 MeV [1].
mBstar : float
    Mass of B* meson; taken from PDG as 5325.2 MeV [1].
mBstar_s : float
    Mass of B*_s meson; taken from PDG as 5415.4 MeV [1].
mK : float
    Mass of kaon; taken from PDG as 497.614 MeV [1].
mpi : float
    Mass of pion; taken from PDG as 134.977 MeV [1].
Delta_B : float
    Difference between masses of B*_s and B^0 mesons for B-->K decay, or that
    of B* and B^0 for B-->pi decay.
a_atol : float
    Absolute tolerance used to deterimine equality of two lattice spacings in
    fermi; taken to be 1e-06 fm, since difference between successive lattice
    spacings is typically of order 1e-02 fm.
gpi : float
    Lattice determination [5] of central value of leading order axial coupling
    constant g_pi (i.e., between pions and heavy-light mesons). Hereafter, g_pi
    is referred to simply as the heavy-light coupling constant.
nearzero : float
    Number used to represent essentially zero, where exactly zero would be
    problematic in practice; taken to be 1e-30.
tastemults : dict of ints
    Multiplicity of tastes of staggered fermions, with tastes defined as 'P' for
    pseudoscalar, 'A' for axial vector, 'T' for tensor, 'V' for vector, and 'I'
    for singlet.
r1_a_continuum : float
    Continuum value of ratio (r_1 / a); taken to be 3.744241 [3], by using the
    average of this ratio for fine lattice spacings.
mh_continuum : float
    Continuum value of heavy quark mass m_h in lattice units (a * m_h); taken to
    be 0.02645 [3].
ml_continuum : float
    Continuum value of light quark mass m_l in lattice units (a * m_l); taken to
    be 0.0009646 [3].
--------------------------------------------------------------------------------
Notes
-----
+ nearzero is used in calculators.chilogs.fcns.R31. nearzero is defined here
  only to provide consistency and transparency for this kind of "magic number".
--------------------------------------------------------------------------------
References
----------
[1] K. Olive, et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014).
[2] C. Bernard, et al., "The static quark potential in three flavor QCD", Phys.
    Rev. D62, 034503 (2000), [arXiv:0002028 [hep-lat]].
[3] C. Bernard, et al. (Fermilab Lattice and MILC Collaborations),
    "UPDATED mass-independent r1 values",
    <http://physics.wustl.edu/~cb/Fermilab-MILC/secure/>, accessed 2015-03-30.
[4] A. Kronfeld and J. Simone, "Computation of Lambda-bar and lambda_1 with
    Lattice QCD", Phys. Lett. B 490, 228 (2000) [Erratum-ibid. B 495, 441
    (2000)] [arXiv:0006345 [hep-ph]].
[5] W. Detmold, C.-J. D. Lin, S. Meinel, "Calculation of the heavy-hadron axial
    couplings g_1, g_2, and g_3 using lattice QCD", Phys. Rev. D 85, 114508
    (2012) [arXiv:1203.3378 [hep-lat]].
--------------------------------------------------------------------------------
"""


from settings.fit import decayname


hbarc = 197.327
r1 = 0.31174
Lambda = 700 / hbarc * r1
fpi = 130.41 / hbarc * r1
mB = 5279.58 / hbarc * r1
mBstar = 5325.2 / hbarc * r1
mBstar_s = 5415.4 / hbarc * r1
mK = 497.614 / hbarc * r1
mpi = 134.977 / hbarc * r1
if decayname == 'B2K':
    Delta_B = mBstar_s - mB
elif decayname == 'B2pi':
    Delta_B = mBstar - mB
a_atol = 1e-06
gpi = 0.45
nearzero = 1e-30
tastemults = {'P': 1, 'A': 4, 'T': 6, 'V': 4, 'I': 1}
r1_a_continuum = 3.744241
mh_continuum = r1_a_continuum * 0.02645
ml_continuum = r1_a_continuum * 0.0009646

