from __future__ import division
import numpy as np
from H2utils import read_H2
from Draine96_H2shield import calc_H2_shield, calc_H2_shield2

c5,c5lo,c5hi,c5a,c6,c6lo,c6hi,c6a = read_H2(sigma=1)

# gas temperature (Jnu and n dependence on this is weak)
T = 100.
# dust grain radius in um (solar is 0.1). n dependence is a^-1
a = 0.1
# dust grain density (solar is 2). n dependence is d^-1
delta = 2.

# HM-only 5.34e-23
# ism-like 1.54e-20
# starburst 6.69e-21

#JnuLW = 5.34e-23

# from Cloudy MW ISM models
#JnuLW = 3.77e-20

# Cloudy IGM Haardt-Madau models
JnuLW = 3.44e-22

# metallicity (use as proxy for gas to dust ratio kappa)
Zhi = 10 ** (-0.86 + 0.24)
Zlo = 10 ** (-0.86 - 0.24)
Z = 10 ** -0.86
kappa = Z
kappalo = Zlo
kappahi = Zhi

# kappa = 10**-1.4
# kappahi = 10**-1.0
# kappalo = 10**-1.9

NH =   10**19.5
NHhi = 10**(19.5 + 0.2)
NHlo = 10**(19.5 - 0.2)

print 'T %.3g K' % T
print 'log10 NH %.3g' % np.log10(NH)
print 'JnuLW %.3g erg/s/cm^2/Hz/sr' % JnuLW
print 'kappa %.3g' % kappa
print ''

C5 = adict()
# comp 6, cm^-2

C5.NH2 = sum(10**c5a)
C5.NH2hi = sum(10**(c5 + c5hi))
C5.NH2lo = sum(10**(c5 - c5lo))

C5.NJ0 = 10**16.65
C5.NJ0hi = 10**(16.65 + 0.34)
C5.NJ0lo = 10**(16.65 - 0.34)
C5.NJ1 = 10**17.20
C5.NJ1hi = 10**(17.20 + 0.36)
C5.NJ1lo = 10**(17.20 - 0.36)

# upper limits
C5.NJ4hi = 10**14.5
C5.NJ5hi = 10**14.5

C6 = adict()
# component 7
C6.NH2 = sum(10**c6a)
C6.NH2hi = sum(10**(c6+c6hi))
C6.NH2lo = sum(10**(c6-c6lo))

C6.NJ0 = 10**15.34
C6.NJ0hi = 10**(15.34 + 0.23)
C6.NJ0lo = 10**(15.34 - 0.23)
C6.NJ1 = 10**16.45
C6.NJ1hi = 10**(16.45 + 0.37)
C6.NJ1lo = 10**(16.45 - 0.37)

# upper limits
C6.NJ4hi = 10**14.5
C6.NJ5hi = 10**14.3

C = adict()
C.NH2 = sum(10**(0.5*(2*c5 + c5hi - c5lo))) +\
        sum(10**(0.5*(2*c6 + c6hi - c6lo)))
C.NH2lo = sum(10**(c5-c5lo)) + sum(10**(c6-c6lo))
C.NH2hi = sum(10**(c5+c5hi)) + sum(10**(c6+c6lo))


# can only put upper limit on beta values. Therefore want NJ0lo,
# NJ1lo, not sure about NH2...

# eqn A7 and A8 from Jorgenson et al. 2010 ApJ 722, 460
def calc_betaJ40(NJ4, NJ0, NH2):
    A42 = 2.75e-9    # s^-1     
    p40 = 0.26
    return  A42 * NJ4/NH2 / (p40 * NJ0/NH2 + 0.021)

def calc_betaJ51(NJ5, NJ1, NH2):
    A53 = 9.9e-9     # s^-1     
    p51 = 0.12
    return A53 * NJ5/NH2 / (p51 * NJ1/NH2 + 0.049)

def calc_tauUV(kappa, NH, a, delta):
    """ UV shielding from dust
    """
    return 0.879 * (a / 0.1)**-1 * (delta/2.)**-1 * kappa * (NH / 10**21)

def calc_Sshield(NH2,  NH, kappa, a, delta):
    """ eqn A13 from Jorgenson et al.

    kappa is the ratio of the dust to gas ratio relative to solar
    (with is ~0.01). Can use the metallicity as a proxy.

    If using more detailed Draine and Bertoldi shielding expression, b
    is assumed to be 5 km/s.
    """
    s1 = calc_H2_shield2(NH2, 5)
    return s1 * np.exp(-calc_tauUV(kappa, NH, a, delta))


def calc_Td(JnuLW, QUV=1, a=0.1, A=3.20e-3):
    """ Find dust temperature using eqn A19 from Jorgenson et al.
    
    JnuLW : float
      
      Radiation field at 1000 Ang in erg/s/cm^2/Hz/sr.
      
    A : float

      Constant (cm) that depends on th optical properties of the
      grains. Default value for carbonaceous grains, for silicate
      grains it should be 1.34e-3.

    a : float

      Grain diameter in micrometers.
    """
    # 3.2e-20 is JnuLW in the solar neighbourhood.
    chi = JnuLW / 3.2e-20
    
    Td = 12 * (chi * QUV)**(1./6) * (A/3.2e-3)**-(1./6) * (a/0.1)**-(1./6)
    return Td

def calc_Sd(T, Td):
    """ Dust sticking coefficient.

    T is gas temperature

    Td is dust temperature
    """
    Sd = (1 + 0.04*(T + Td)**0.5 + 2e-3*T + 8e-6*T**2)**-1 * \
         (1 + np.exp( 7.5e2* (1./75 - 1./Td) ) )**-1
    return Sd
    
def calc_R(Sd, kappa=1., a=0.1, T=100., delta=2.):
    """ Eqn A17 from Jorgenson et al.

    Sd : float

     Dust sticking coefficient

    delta : float

      Grain material density in g cm^-3

    kappa : float

      Ratio of Dust-to-gas mass ratio to solar neighbourhood
      value. (Can use metallicity as a proxy.)

    """
    R = 4.1e-17 * Sd * (a / 0.1)**-1 * kappa * (T/100)**0.5 * (delta / 2.)**-1
    return R

def calc_approx_fH2(NH, n, JnuLW, kappa, a, delta, A, QUV, T):
    """ Neglect dust shielding, assume NH >> NH2, following Hirashita
    and Ferrara 2005.

    n = H density in cm^-3
    JnuLW = radiation strength at LW bands (1000 ang) in erg/s/cm^2/Hz/sr

    kappa = dust to gas ratio / (solar neighbourhodd gas to dust ratio)


    a = dust grain radius in um (solar is 0.1).

    delta = dust grain density (solar is 2).


    T = gas temperature K (fH2 has only ^(1/2) dependence on this)

    A = Constant (in cm) that depends on th optical properties of the
        grains. Default value for carbonaceous grains, for silicate
        grains it should be 1.34e-3. only ^(1/6) dependence.

    QUV = used to calculate dust temperature
        (fHs has only ^(1/6) dependence on this)
    """
    # find R
    Td = calc_Td(JnuLW, QUV=QUV, a=a, A=A)
    Sd = calc_Sd(T, Td)
    R = calc_R(Sd, kappa=kappa, a=a, T=T, delta=delta)

    # without shielding
    fH2 = np.ones_like(NH) * 0.01 * (n / 92.) * \
          (3.98e-20 / JnuLW) * (R / 3e-17)
    c0 = fH2 * NH  > 2*10**14
    # need shielding
    fH2s = ((n / 92.) * (3.98e-20 / JnuLW) * (R / 3e-17))**4 \
           * (NH[c0] / 10**19.63)**3
    fH2[c0] = np.where(fH2[c0] > fH2s, fH2[c0], fH2s)

    return fH2


c = C

if 0:

    # upper limits
    beta0 = calc_betaJ40(c.NJ4hi, c.NJ0lo, c.NH2lo)
    beta1 = calc_betaJ51(c.NJ5hi, c.NJ1lo, c.NH2lo)
     
    # upper limit
    Rdiss0 = 0.11 * beta0
    Rdiss1 = 0.11 * beta1
     
    # JnuLW erg/s/cm^2/ster
    print 'Jnu upper limits at 1000 Ang from J 4-0, J 5-1'
    print ' ', Rdiss0 / (4 * np.pi * 1.1e8 * Slo)
    print ' ', Rdiss1 / (4 * np.pi * 1.1e8 * Slo)
     
    #JnuLW = Rdiss0 / (4 * np.pi * 1.1e8 * Slo)
    print 'using Jnu=%.2g' % JnuLW 
     
    print 'n(HI)'
     

if 1:
    Slo = calc_Sshield(c.NH2lo,  NHhi, kappahi, a, delta)
    Shi = calc_Sshield(c.NH2hi,  NHlo, kappalo, a, delta)

    Td = calc_Td(JnuLW, QUV=1, a=a, A=3.20e-3)
    Sd = calc_Sd(T, Td)
    R = calc_R(Sd, kappa=kappa, a=a, T=T, delta=delta)
    Rlo = calc_R(Sd, kappa=kappalo, a=a, T=T, delta=delta)
    Rhi = calc_R(Sd, kappa=kappahi, a=a, T=T, delta=delta)
     
    print '  R/4e-17 = %.2f'  % (R / 4e-17)
    #R = 4.1e-17 * sqrt(2)
    fH2 = 2*c.NH2 / (NH + 2*c.NH2)
     
    Rdisshi = 4 * pi * 1.1e8 * JnuLW *  Slo
    Rdisslo = 4 * pi * 1.1e8 * JnuLW *  Shi
    #RdissnoS = 4 * pi * 1.1e8 * JnuLW 
    from barak.constants import Gyr
    print 't_Rdisslo in Myr', 1 / Rdisslo / Gyr * 1e3
    print 't_Rdisshi in Myr', 1 / Rdisshi / Gyr * 1e3

    nmin = 0.5 * Rdisslo / Rhi * fH2
    nmax = 0.5 * Rdisshi / Rlo * fH2
    print ' n min in cm^-3', nmin
    print ' n max in cm^-3', nmax

    print 't_formlo in Gyr', 1 / (Rlo * nmin) / Gyr
    print 't_formhi in Gyr', 1 / (Rhi * nmax) / Gyr

    #print ' ', 0.5 * RdissnoS / R * fH2

    # NH2lo seems to give the most conservative upper limit, stick with that
    # print 'C5'
    # print 'beta0 max =', 
    # print 'beta1 max =', 
     
    # print 'C6'
    # print 'beta0 max =', calc_betaJ40(C6.NJ4hi, C6.NJ0lo, C6.NH2lo)
    # print 'beta1 max =', calc_betaJ51(C6.NJ5hi, C6.NJ1lo, C6.NH2lo)

if 0:
    # Some sample N-fH2 curves
    Ntot = 10**np.arange(18,23, 0.05)

    # typical of IGM
    JnuLW = 5.34e-23
    #JnuLW = 3.98e-20            # ISM


    # dust grain radius in um (solar is 0.1). n dependence is a^-1
    a = 0.1
    # dust grain density (solar is 2). n dependence is d^-1
    delta = 2.
    # Constant (in cm) that depends on th optical properties of the
    # grains. Default value for carbonaceous grains, for silicate
    # grains it should be 1.34e-3.
    A=3.20e-3
    # QUV ?
    QUV = 1
    # gas temperature (Jnu and n dependence on this is weak)
    T = 100.
    print kappa, n, JnuLW
    Z = 0.1
    kappa = 10**-1.4
    # particle density in cm^-3
    nvals = [0.5]
    for n in nvals:
        fH2 = calc_approx_fH2(Ntot, n, JnuLW, kappa, a, delta, A, QUV, T)
        NH2 = 0.5 * Ntot * fH2
        pl.figure(1)
        plot(log10(Ntot), np.log10(fH2))
        pl.figure(2)
        plot(log10(NH2), np.log10(fH2))

    show()
