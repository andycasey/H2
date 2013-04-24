# solar O abundance = 8.73 +/- 0.07

# where the abundance is defined A(O) = log N(O)/N(H) + 12

from barak.abundances import Asolar, calc_abund, cond_temp

import numpy as np

def calc_Zlim(X, Y, NX, NY):
    """
    X = atom X
    Y = atom Y
    NX = lower lim, best, upper lim
    NY = lower lim, best, upper lim
    """
    Zmin = calc_abund(X, Y, NX[0], NY[2])
    Z = calc_abund(X, Y, NX[1], NY[1])
    Zmax = calc_abund(X, Y, NX[2], NY[0])
    return np.array([Zmin, Z, Zmax])

# Asolar has values  log10 n(el)/n(H) + 12
# OI = 15.39, 15.43, 15.6 #15.47
# HI = 19.3, 19.5, 19.7
# SiII = 14.15, 14.79, 15.02
# SIII = 10, 14.8, 14.8   # no constraint on SII, SI is lower (< 14)
# CII = 15.1, 15.1, 18  # CIII, CI are both lower (< 14)
# FeII = 14.02, 14.06, 14.10 
# MgII = 13.98, 14.01, 14.04  # MgI much lower (~12.3)

OI = 15.27, 15.54, 15.78 #15.47
HI = 19.3, 19.5, 19.7
SiII = 14.15, 14.79, 15.02
SIII = 10, 14.8, 14.8   # no constraint on SII, SI is lower (< 14)
CII = 15.1, 15.1, 18  # CIII, CI are both lower (< 14)
FeII = 14.02, 14.06, 14.10 
MgII = 13.98, 14.01, 14.04  # MgI much lower (~12.3)


print 'OI   %5.2f %5.2f %5.2f' % tuple(calc_Zlim('O', 'H', OI, HI)      ) 
print 'OI/HI sigma = %5.2f' % np.hypot(0.25, 0.2)
print '%.2f %.2f %.2f' % (10**(-0.72 - 0.32), 10**-0.72,10**(-0.72 + 0.32))
print 'SiII %5.2f %5.2f %5.2f' % tuple(calc_Zlim('Si','H', SiII, HI)    ) 
print 'SIII %5.2f %5.2f %5.2f' % tuple(calc_Zlim('S', 'H', SIII, HI)    ) 
print '[CII/O]  %5.2f %5.2f %5.2f' % tuple(calc_Zlim('C', 'O', CII,  OI)) 
print '[FeII/O] %5.2f %5.2f %5.2f' % tuple(calc_Zlim('Fe','O', FeII, OI)) 
print 'FeII/OI sigma = %5.2f' % np.hypot(0.04, 0.2)
print 10**(-0.26 - 0.20), 10**-0.26,10**(-0.26 + 0.20) 
print '[MgII/O] %5.2f %5.2f %5.2f' % tuple(calc_Zlim('Mg','O', MgII, OI)) 
print 'MgII/II sigma = %5.2f' % np.hypot(0.03, 0.2)



