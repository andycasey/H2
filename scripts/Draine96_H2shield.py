""" Analytic expressions for H2 self-shielding from Draine and Bertoldi 1996
"""
import numpy as np
import pylab as pl

def calc_H2_shield(NH2):
    """ fshield is the fractional drop in pumping/dissociating UV
    photons due to H2 self-shielding.

    NH2 in cm^-2

    Equation 36 from Draine and Bertoldi 1996 ApJ, 468, 269.
    """
    NH2 = np.atleast_1d(NH2)
    logNH2 = np.log10(NH2)
    f = np.ones(NH2.shape, float)
    c0 = logNH2 > 14 
    f[c0] = (10**(logNH2[c0] - 14))**-0.75
    if len(f) == 1:
        return f[0]
    return f

def calc_H2_shield2(NH2, b):
    """ fshield is the fractional drop in pumping/dissociating UV
    photons due to H2 self-shielding.

    NH2 in cm^-2
    b in km/s

    Equation 36 from Draine and Bertoldi 1996 ApJ, 468, 269.
    """
    x = NH2 / 5e14
    temp = (1 + x)**0.5
    f = 0.965 / (1 + x/b)**2 + 0.035 / temp * np.exp(-8.5e-4*temp)
    return f

if 0:
    NH2 = 10**np.arange(13, 22, 0.1)
    pl.semilogy(log10(NH2), calc_H2_shield(NH2))
    pl.semilogy(log10(NH2), calc_H2_shield2(NH2, 1),label='1')
    pl.semilogy(log10(NH2), calc_H2_shield2(NH2, 3),label='3')
    pl.semilogy(log10(NH2), calc_H2_shield2(NH2, 5),label='5')
    pl.semilogy(log10(NH2), calc_H2_shield2(NH2, 10),label='10')
    pl.legend(frameon=0, loc='best')
    pl.show()
