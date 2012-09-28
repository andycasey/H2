"""
From Prochaska et al 1999 (http://adsabs.harvard.edu/abs/1999ApJ...511L..71P)

The detection of the fine-structure C ii* l1335 transition allows
further insight into the physical characteristics of this system. The
fine-structure line can be excited by several processes, but for a
highly ionized system at this redshift, electron collisions dominate
(Morris, S. L., Weymann, R. J., Foltz, C. B., Turnshek, D. A.,
Shectman, S., Price, C., & Boronson, T. A. 1986, ApJ, 310,
40). Therefore, the observed N(C ii*)/N(C ii) ratio provides a measure
of the electron density via

N(CII*)/N(CII) = 3.9e-2 * n_e * [1 + (0.22 * n_p/n_e)]

n_e is electron density, n_p is proton density. Estimate n_p/n_e from
the expected ratio of HI to HII. (If all in HI, n_e is very small).

Unfortunately, the N(C ii*)/N(C ii) ratio cannot be measured directly
because the C ii l1334 profile is saturated over the velocity region
where the C ii* absorption is detected. We can accurately estimate the
N(C ii*)/N(C ii) ratio, however, provided the assumption that the C ii
l1334 and Fe ii l1608 profiles track one another in velocity
space. Low ion profiles always track one another in the damped Lya
systems (e.g., Prochaska & Wolfe 1996; Lu et al. 1996), and one notes
that all of the transitions-irrespective of ionization state-trace one
another in this system. Our approach, then, is to (1) measure N(C
ii)/N(Fe ii) in Reg4, where the effects of saturation are minimal, and
(2) correct a N(Fe ii) measurement in Reg2 by the N(C ii)/N(Fe ii)
ratio to estimate N(C ii) in this velocity region.
"""
import numpy as np

# Assume NCII as given by Cloudy model (~10^15.2)

# NCIIstar ~ 14.

# uplim
NCIIstar = 10**14.6

NCII = 10**14.6

# Assume 30% HI, 70% HII (as in best fitting IGM cloudy model)
np_on_ne  = 1 / 0.70

#N(CII*)/N(CII) = 3.9e-2 * n_e * [1 + (0.22 * n_p/n_e)]

ne = NCIIstar / NCII * (3.9e-2 * (1 + (0.22 + np_on_ne)))**-1

print ne, np.log10(ne)
