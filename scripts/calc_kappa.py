# kappa is a measure of the dust to gas ratio. See Appendix of Wolfe
# et al. 2003: http://adsabs.harvard.edu/abs/2003ApJ...593..215W.
# Note that we assume the total (gas + dust) ratio Fe/Si in the DLA is
# the same as for the MW, so the first term inside the brackets is 1.
#
# kappa_X = 10^[X/H] * (1 - 10^[Fe/X])
#
# X is a reference element that is not depleted onto dust. (Fe is
# heavily depleted). Hows the dust ratio defined? Is it proportional
# to the number of grains, or the mass? People seem to use any dust
# measure, i.e. A(V), or in this case N(Fe depleted onto grains). 
#
# Calculate depletion using Fe/O (Assume all O is in OI, all Fe in
# FeII).

from barak.absorb import split_trans_name
from barak.utilities import adict
from cloudy.cloudy_utils import read_observed, get_ratio

def calc_abund_range(k1, k2, obs):
    el1, el2 = (split_trans_name(k)[0] for k in (k1, k2))
    best = calc_abund(el1, el2, obs[k1][0], obs[k2][0])
    lo = calc_abund(el1, el2, obs[k1][0]-obs[k1][1], obs[k2][0]+obs[k2][1])
    hi = calc_abund(el1, el2, obs[k1][0]+obs[k1][1], obs[k2][0]-obs[k2][1])
    return np.array([best, lo, hi])

def calc_kappa(idep, inodep, obs):
    eldep, elnodep = (split_trans_name(k)[0] for k in (idep, inodep))
    Ndep = obs[idep][0]
    Nnodep = obs[inodep][0]
    NH = obs['HI'][0]
    Ndeplo = obs[idep][0] - obs[idep][1]
    Nnodeplo = obs[inodep][0] - obs[inodep][1]
    NHlo = obs['HI'][0] - obs['HI'][1]
    Ndephi = obs[idep][0] + obs[idep][2]
    Nnodephi = obs[inodep][0] + obs[inodep][2]
    NHhi = obs['HI'][0] + obs['HI'][2]

    t1 = 10**calc_abund(elnodep, 'H', Nnodep, NH)
    t2 = (1 - 10**calc_abund(eldep, elnodep, Ndep, Nnodep))
    kappa = t1 * t2

    # upper limit
    t1 = 10**calc_abund(elnodep, 'H', Nnodephi, NHlo)
    t2 = (1 - 10**calc_abund(eldep, elnodep, Ndeplo, Nnodephi))

    kappahi = t1 * t2

    # lower limit
    t1 = 10**calc_abund(elnodep, 'H', Nnodeplo, NHhi)
    t2 = (1 - 10**calc_abund(eldep, elnodep, Ndephi, Nnodeplo))

    kappalo = t1 * t2
    return np.array([kappa, kappalo, kappahi])

obs = adict(read_observed('../cloudy/observed_logN/total'))

from barak.abundances import calc_abund

if 1:

    print np.log10(calc_kappa('FeII', 'OI',  obs))
    print np.log10(calc_kappa('FeII', 'SiII',obs))
    print np.log10(calc_kappa('MgII', 'OI',  obs))

     
    print '[Fe/O]', calc_abund_range('FeII', 'OI', obs)
    print '[Mg/O]', calc_abund_range('MgII', 'OI',   obs)
    print '[Fe/Si]', calc_abund_range('FeII', 'SiII',obs)
    print '[Mg/Si]', calc_abund_range('MgII', 'SiII',obs)
    print '[Si/Fe]', calc_abund_range('SiII', 'FeII',obs)
    print '[O/Fe]', calc_abund_range('OI', 'FeII',   obs)
    print '[Si/H]', calc_abund_range('SiII', 'HI',   obs)
    print '[O/H]', calc_abund_range('OI', 'HI',      obs)
    print '[Fe/H]', calc_abund_range('FeII', 'HI',   obs)

