# solar O abundance = 8.73 +/- 0.07

# where the abundance is defined A(O) = log N(O)/N(H) + 12

#Measured N(OI) = 15.3 +/- 0.04  
#Measured N(HI) = 19.5 +/ 0.1

from barak.abundances import Asolar, calc_abund, cond_temp

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
OI = 15.39, 15.43, 15.6 #15.47
HI = 19.3, 19.5, 19.7
SiII = 14.15, 14.79, 15.02
SIII = 10, 14.8, 14.8   # no constraint on SII, SI is lower (< 14)
CII = 15.1, 15.1, 18  # CIII, CI are both lower (< 14)
FeII = 14.02, 14.06, 14.10 
MgII = 13.98, 14.01, 14.04  # MgI much lower (~12.3)

print 'OI   %5.2f %5.2f %5.2f' % tuple(calc_Zlim('O', 'H', OI, HI)      ) 
print 'SiII %5.2f %5.2f %5.2f' % tuple(calc_Zlim('Si','H', SiII, HI)    ) 
print 'SIII %5.2f %5.2f %5.2f' % tuple(calc_Zlim('S', 'H', SIII, HI)    ) 
print '[CII/O]  %5.2f %5.2f %5.2f' % tuple(calc_Zlim('C', 'O', CII,  OI)) 
print '[FeII/O] %5.2f %5.2f %5.2f' % tuple(calc_Zlim('Fe','O', FeII, OI)) 
print '[MgII/O] %5.2f %5.2f %5.2f' % tuple(calc_Zlim('Mg','O', MgII, OI)) 



