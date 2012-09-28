from barak.interp import trilinear_interp
from barak.utilities import adict, indexnear
from barak.absorb import split_trans_name
from barak.plot import arrplot, puttext, draw_arrows
from barak.io import loadobj, parse_config
from cloudy.cloudy_utils import read_observed, get_ratio

from glob import glob
import os
import numpy as np

prefix = '../cloudy/'
simname = 'igm'
#simname = 'agn'
#simname = 'igm_dust'
#simname = 'igm_dust_const_press'  # also cosmic rays
#simname = 'ism_with_grains'
#simname = 'starburst'

gridname = os.path.join(prefix, simname, 'grid.cfg')

roman_map = {'I':0, 'II':1, 'III':2, 'IV':3, 'V':4, 'VI':5,
             'VII':6, 'VIII':7, 'IX':8, 'X':9, '2':2}

#COLORS = dict(Si='y', C='k', Al='c',O='r',N='g', Fe='orange', Ne='pink',
#               Mg='b', Ca='0.7', Zn='purple', Cr='m')
COLORS = dict(Si='orange', C='k', Al='c',O='g', N='c', Fe='r', Ne='pink',
              Mg='b', Ca='m', Zn='purple', Cr='m', H='y')

LS = ['solid','dashed','dotted', '-.']

def plot_trans_vs_NHI(iZ, inH, grid, trans, ax):

    atoms,nums = zip(*[split_trans_name(t) for t in trans])

    count = dict((k, 0) for k in COLORS)
    for atom, num in zip(atoms, nums):
        col = COLORS[atom]
        N = grid.N[atom][:, inH, iZ, roman_map[num]]
        pl.plot(grid.NHI, N, lw=2, ls=LS[count[atom] % len(LS)],
                color=col,label=atom+num)
        count[atom] += 1
 
    pl.legend(frameon=0, ncol=3)
    pl.ylabel('$\log_{10}\ N (\mathrm{cm}^{-2}$')
    pl.xlabel('$\log_{10}\ N_\mathrm{HI} (\mathrm{cm}^{-2})$')
    puttext(0.05, 0.95, '$n_{H}=%.3g\  Z=%.3g$' % (
        grid.nH[inH], grid.Z[iZ]), ax)

def plot_trans_vs_U(iZ, iNHI, grid, trans, ax):

    atoms,nums = zip(*[split_trans_name(t) for t in trans])

    count = dict((k, 0) for k in COLORS)
    for atom, num in zip(atoms, nums):
        col = COLORS[atom]
        N = grid.N[atom][iNHI, :, iZ, roman_map[num]]
        pl.plot(grid.nH, N, lw=3, ls=LS[count[atom] % len(LS)],
                color=col,label=atom+num)
        count[atom] += 1
 
    ax.legend(frameon=0, ncol=2)
    ax.set_ylabel('$\log_{10}\ N\ (\mathrm{cm}^{-2})$')
    ax.set_xlabel('$\log_{10}\ n_\mathrm{H}\ (\mathrm{cm}^{-3})$')

    puttext(0.05, 0.05, '$\log_{10}\ N_\mathrm{HI}=%.3g,\ Z=%.3g$' % (
        grid.NHI[iNHI],grid.Z[iZ]),ax)
    ax.set_xlim(grid.nH[0]+0.01, grid.nH[-1]-0.01)

    ax1 = ax.twiny()
    x0,x1 = ax.get_xlim()
    const = (grid.U + grid.nH)[0]
    assert np.allclose(const, grid.U + grid.nH)
    ax1.set_xlim(const - x0, const - x1)
    ax1.set_xlabel('$\log_{10}\ U$')


if 1:
    ##############################################
    # Read the model
    ##############################################
    
    cfg = parse_config(gridname)

    M = loadobj(os.path.join(prefix, simname, cfg.prefix + '_grid.sav.gz'))
    M = adict(M)

    # A finer grid of parameter values for interpolation below
    NHI = np.linspace(M.NHI[0], M.NHI[-1], 100)
    nH = np.linspace(M.nH[0], M.nH[-1], 101)
    Z = np.linspace(M.Z[0], M.Z[-1], 102)

    dNHI = NHI[1] - NHI[0]
    dnH = nH[1] - nH[0]
    dZ = Z[1] - Z[0]


if 1:
    ##############################################
    # Read the observed column densities
    ##############################################

    obs = read_observed('../data/observed_logN/total')
    ratios = ('MgI/MgII FeI/FeII CaI/CaII NI/NII SiII/SiIII').split()

    obs_ratios = []
    for ratio in ratios:
        numer, denom = ratio.split('/')
        low, best, high = get_ratio(obs[numer], obs[denom])
        obs_ratios.append( (low, best, high) )

    r = obs_ratios[0]
    MgI_MgII = r[1]
    MgI_MgII_sig = max(r[2] - r[1], r[1] - r[0])
    FeI_FeII_max = obs_ratios[1][2]
    CaI_CaII_max = obs_ratios[2][2]
    NI_NII_max = obs_ratios[3][2]
    SiII_SiIII_max = obs_ratios[4][2]


if 0:
    #############################################
    # plot the incident radiation field
    #############################################
    pl.figure()
    pl.loglog(M.cont.ryd, M.cont.fnu)
    pl.xlabel(r'$Energy (Rydbergs)$')
    pl.ylabel(r'$F_{\nu} (ergs/s/cm^2/Hz)$')

    pl.draw()

if 1:
    #############################################
    # grid interpolation
    #############################################

    minterp = {}

    print 'interpolating...'
    # MgI / MgII
    ratio = M.N['Mg'][:,:,:,0] - M.N['Mg'][:,:,:,1]
    minterp['MgI/MgII'] = trilinear_interp(
        NHI, nH, Z, M.NHI, M.nH, M.Z, ratio)
    # FeI / FeII
    ratio = M.N['Fe'][:,:,:,0] - M.N['Fe'][:,:,:,1]
    minterp['FeI/FeII'] = trilinear_interp(
        NHI, nH, Z, M.NHI, M.nH, M.Z, ratio)
    # CaI / CaII
    ratio = M.N['Ca'][:,:,:,0] - M.N['Ca'][:,:,:,1]
    minterp['CaI/CaII'] = trilinear_interp(
        NHI, nH, Z, M.NHI, M.nH, M.Z, ratio)
    # NI / NII
    ratio = M.N['N'][:,:,:,0] - M.N['N'][:,:,:,1]
    minterp['NI/NII'] = trilinear_interp(
        NHI, nH, Z, M.NHI, M.nH, M.Z, ratio)
    # SiII / SiIII
    ratio = M.N['Si'][:,:,:,1] - M.N['N'][:,:,:,2]
    minterp['SiII/SiIII'] = trilinear_interp(
        NHI, nH, Z, M.NHI, M.nH, M.Z, ratio)
    # HI
    #minterp['HI'] = trilinear_interp(
    #    NHI, nH, Z, M.NHI, M.nH, M.Z, M.N['H'][:,:,:,0])

    print 'done'

    # from barak.plot import arrplot
    # arrplot(nH, Z, MgII[0,:,:])
    # arrplot(M.nH, M.Z, M.N['Mg'][0,:,:,1])
    # pl.show()

if 1:
    #############################################
    # find the log likelihood.
    #############################################

    # gaussian probability around MgI/MgII value

    ln_prob = np.zeros_like(minterp['MgI/MgII'])
    ln_prob += -0.5 * ((minterp['MgI/MgII'] - MgI_MgII) / MgI_MgII_sig)**2

    invalid = np.zeros(minterp['MgI/MgII'].shape, bool)
    invalid |= minterp['NI/NII'] > NI_NII_max
    invalid |= minterp['SiII/SiIII'] > SiII_SiIII_max

    # these don't provide any constraint, so can skip them
    invalid |= minterp['CaI/CaII'] > CaI_CaII_max
    invalid |= minterp['FeI/FeII'] > FeI_FeII_max

    # skip this, it's not very constraining
    #invalid |= minterp['HI'] > 19.7

    ln_prob[invalid] = -np.inf

    # create an array of Z values to apply Z prior
    Zvals = np.ones_like(minterp['MgI/MgII']) * Z
    ln_prob +=  -0.5 * ((Zvals - (-0.86)) / 0.24)**2

    prob = np.exp(ln_prob)
    # normalise
    prob /= (prob.sum() * dNHI * dZ * dnH)

if 0:
    #############################################
    # plot marginalised posteriors
    #############################################

    pl.figure()
    ax = pl.subplot(311)
    ax.plot(NHI, prob.sum(1).sum(1) * dZ * dnH)
    ax.set_xlabel('NHI')
    ax = pl.subplot(312)
    ax.plot(nH, prob.sum(0).sum(1) * dNHI * dZ)
    ax.set_xlabel('nH')
    ax = pl.subplot(313)
    ax.plot(Z, prob.sum(0).sum(0) * dNHI * dnH)
    ax.set_xlabel('Z')
    pl.show()

if 0:
    #############################################
    # marginalised posteriors for paper
    #############################################

    pl.rc('font', size=13)
    fig = pl.figure(figsize=(4.3, 2.8))
    fig.subplots_adjust(bottom=0.16, top=0.85, left=0.08, right=0.98,
                        wspace=1e-5)
    ax0 = pl.subplot(121)
    ax0.plot(NHI, prob.sum(1).sum(1) * dZ * dnH, lw=1.5)
    ax0.set_xlabel('$\log_{10}\ N_\mathrm{HI}$')
    ax1 = pl.subplot(122)
    ax1.plot(nH, prob.sum(0).sum(1) * dNHI * dZ, lw=1.5)
    ax1.set_xlabel('$\log_{10}\ n_\mathrm{H}$')
    ax1a = ax1.twiny()
    ax1a.set_xlabel('$\log_{10}\ U$')
    x0,x1 = ax1.get_xlim()
    # get the conversion constant from nH to U
    const = (M.U + M.nH)[0]
    assert np.allclose(const, M.U + M.nH)
    ax1a.set_xlim(const - x0, const - x1)
    ax1.set_yticklabels('')
    ax0.set_ylim(*ax1.get_ylim())
    ax0.set_xlim(18.01, 20.1)
    ax1.set_xlim(-3.9, 1.2)
    pl.savefig('NHI_U.pdf')
    pl.show()

    # Conclusions: Ratios of transitions in the same species (Si, C)
    # are very sensitive to the ionization parameter (thus density),
    # but mostly insensitive to the HI column density over the range
    # logN 14 -> 18 and metallicity over the range logZ -2 -> 0.


if 1:
    ##################################################################
    # plot the predicted column densities as a function of NHI, and
    # compare to observed values for the maximum likelihood model
    ################################################################## 

    # for solar abundances, plotting vs logU/lognH
    
    # OVI, CaI, AlI, FeI are very low (logN = 9-10)
    # SiIV, CIV consistent with upper limits.

    imax = where(prob == prob.max())
    nHbest = nH[imax[1][0]]    
    iZ = indexnear(M.Z, Z[imax[2][0]])
    inH = indexnear(M.nH, nH[imax[1][0]])
    #iNHI = indexnear(M.NHI, 19.5)
    iNHI = indexnear(M.NHI, NHI[imax[0][0]])

    pl.rc('font', size=14)
    pl.rc('legend', fontsize=10.5)

if 0:
    trans1 = ('FeII MgI MgII OI OVI NI NII').split()

    fig = pl.figure(figsize=(4.3, 6))
    fig.subplots_adjust(right=0.97, top=0.92,bottom=0.09)
    ax = pl.subplot(111)

    plot_trans_vs_U(iZ, iNHI, M, trans1, ax)
    ax.set_ylim(11.1, 16.3)

    # observed values
    val,sig = obs['MgII'][:2]
    ax.plot(nHbest, val, 'o', mec=COLORS['Mg'], mew=2, mfc='w', ms=6)
    val,sig = obs['MgI'][:2]
    ax.plot(nHbest, val, 'o', color=COLORS['Mg'], ms=9, mew=2, mec='w')
    val,sig = obs['NI'][:2]
    ax.plot(nHbest, val, 'o', color=COLORS['N'], ms=9, mew=2, mec='w')
    draw_arrows(nHbest, val, ax=ax, c=COLORS['N'],ms=1.5,capsize=5,
                lw=2, direction='down', zorder=6)
    val,sig = obs['NII'][:2]
    ax.plot(nHbest, val, 'o',mec=COLORS['N'], mew=2, mfc='w', ms=6)
    val,siglo,sighi = obs['OI']
    ax.plot(nHbest, val, 'o', color=COLORS['O'], ms=9, mew=2, mec='w')
    val,sig = obs['FeII'][:2]
    ax.plot(nHbest-0.2, val, 'o', color=COLORS['Fe'], ms=9, mew=2, mec='w')
    #pl.savefig('cloudy1.pdf')
    #pl.savefig('cloudy1.png')

if 0:
    trans2 = ('HIII CaII SiII SiIII SiIV CI CII CIII').split()
    fig = pl.figure(figsize=(4.3, 6))
    fig.subplots_adjust(right=0.97, top=0.92,bottom=0.09)
    ax = pl.subplot(111)

    plot_trans_vs_U(iZ, iNHI, M, trans2, ax)
    ax.set_ylim(10.1, 17.3)

    # Now plot observed column densities
    val,sig = obs['CaII'][:2]
    ax.plot(nHbest, val, 'o', c=COLORS['Ca'], ms=9, mew=2,mec='w')
    val,sig = obs['CIII'][:2]
    draw_arrows(nHbest+0.2, val, ax=ax, c=COLORS['C'],ms=1.5,capsize=5,
                lw=2, zorder=5)#, linestyle=[(0, (3, 2))])
    col = ax.scatter(nHbest+0.2, val, marker='o', edgecolors=COLORS['C'],
                     lw=2, s=40, c='w', zorder=10)
    col.set_linestyles([(0, (2, 1))])
    val,siglo,sighi = obs['CI']
    ax.plot(nHbest, val, 'o', color=COLORS['C'],ms=9, mew=2,mec='w')
    yerr = np.atleast_2d([siglo, sighi]).T
    _,cap,err = ax.errorbar(nHbest, val, yerr=yerr, fmt=None,
                            lw=2, ecolor=COLORS['C'], mew=2, capsize=4)
    err[0].set_zorder(10)
    val,sig = obs['CII'][:2]
    ax.plot(nHbest, val, 'o',mec=COLORS['C'], mew=2, mfc='w', ms=6,zorder=6)
    draw_arrows(nHbest, val, ax=ax, c=COLORS['C'],ms=1.5,capsize=5,
                lw=2)
    val,sig = obs['SiIII'][:2]
    ax.plot(nHbest-0.2, val, 'o',mec=COLORS['Si'], mew=2, mfc='w', ms=6)
    draw_arrows(nHbest-0.2, val, ax=ax, c=COLORS['Si'],ms=1.5,capsize=5,
                lw=2)
    val,siglo,sighi = obs['SiII']
    ax.plot(nHbest-0.2, val, 'o', color=COLORS['Si'],ms=9, mew=2,mec='w')
    yerr = np.atleast_2d([siglo, sighi]).T
    _,cap,err = ax.errorbar(nHbest-0.2, val, yerr=yerr, fmt=None,
                            lw=2, ecolor=COLORS['Si'], mew=2, capsize=4)
    err[0].set_zorder(10)
    #pl.savefig('cloudy2.pdf')
    #pl.savefig('cloudy2.png')

    pl.show()

show()


if 0:
    #######################################################################
    # Plots marginalised along only one parameter, can be good for
    # diagnostics
    #######################################################################
    pl.figure()
    # collapse the likelihood along each dimension
    prob_NHI_nH = prob.sum(axis=2)
    prob_NHI_Z = prob.sum(axis=1)
    prob_Z_nH = prob.sum(axis=0).T
    pl.figure()

    ax = pl.subplot(221)
    arrplot(NHI, nH, prob_NHI_nH, ax=ax, cmap=pl.cm.hot)#, vmin=0, vmax=1)
    ax.set_ylabel('nH')
    ax = pl.subplot(222)
    arrplot(Z, nH, prob_Z_nH, ax=ax, cmap=pl.cm.hot)#, vmin=0, vmax=1)
    ax.set_xlabel('Z')
    ax = pl.subplot(223)
    arrplot(NHI, Z, prob_NHI_Z, ax=ax, cmap=pl.cm.hot)#, vmin=0, vmax=1)
    ax.set_xlabel('NHI')
    ax.set_ylabel('Z')

    pl.show()


##################################
# Debugging
#################################
if 0:
    pl.figure()
    ax = pl.gca()
    for i in range(len(NHI)):
        ax.cla()
        c = arrplot(Z, nH, prob[i,:,:].T, ax=ax, cmap=pl.cm.hot, vmin=0, vmax=1)
        ax.set_xlabel('Z')
        ax.set_ylabel('nH')
        ax.set_title('NHI %.2f' % NHI[i])
        pl.show()
        pl.waitforbuttonpress()

if 0:
    pl.figure()
    ax = pl.gca()
    for i in range(len(nH)):
        ax.cla()
        c = arrplot(NHI, Z, prob[:,i,:], ax=ax, cmap=pl.cm.hot, vmin=0, vmax=1)
        ax.set_xlabel('NHI')
        ax.set_ylabel('Z')
        ax.set_title('nH %.2f' % nH[i])
        pl.show()
        pl.waitforbuttonpress()

if 0:
    pl.figure()
    ax = pl.gca()
    for i in range(len(Z)):
        ax.cla()
        c = arrplot(NHI, nH, prob[:,:,i], ax=ax, cmap=pl.cm.hot, vmin=0, vmax=1)
        ax.set_xlabel('NHI')
        ax.set_ylabel('nH')
        ax.set_title('Z %.2f' % Z[i])
        pl.show()
        pl.waitforbuttonpress()
