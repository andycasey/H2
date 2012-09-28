from scipy.constants import Boltzmann as kboltz
from scipy.constants import eV
from H2utils import read_H2, calc_OPR
from barak.plot import shade_to_line, puttext, draw_arrows, shade_to_line_vert
import asciitable
import atpy
import os
from astropy.cosmology import WMAP7
from model_Krumholz import calc_fH2

def detlim(NHI, NH2=13.5):
    Ntot = 10**NHI + 2*10**NH2
    return np.log10(Ntot), np.log10(10**NH2 / Ntot)

J = 0,1,2,3,4,5
gJ = 1, 9, 5, 21, 9, 33
EJ = 0, 0.01469, 0.04394, 0.08747, 0.14491, 0.21575

# Note for LMC and SMC only, negative column densities are upper
# limits.

MWplane = atpy.Table('../data/Savage77_tab1.tbl')

LMC = Table('../data/Welty12_LMC.tbl')
LMC = LMC.where((LMC.NHI > 0) & (LMC.NH2 > 0) &
                (LMC.lNHI != '<') & (LMC.lNH2 != '<'))
SMC = Table('../data/Welty12_SMC.tbl')
SMC = SMC.where((SMC.NHI > 0) & (SMC.NH2 > 0) & (SMC.lNHI != '<'))

IVC = Table('../data/Richter03_t2.tbl')
IVC = IVC.where((IVC.lNH2 == '') & ~np.isnan(IVC.NH2))


# this is just used for fH2 measurements
DLAs = atpy.Table('../data/f_H2_DLA.tbl').data.view(np.recarray)

c5,c5lo,c5hi,c5a,c6,c6lo,c6hi,c6a = read_H2(sigma=1)
c5_NH2hi = log10(sum(10**(c5+c5hi)))
c5_NH2lo = log10(sum(10**(c5-c5lo)))
c5_NH2 = log10(sum(10**c5a))
c6_NH2hi = log10(sum(10**(c6+c6hi)))
c6_NH2lo = log10(sum(10**(c6-c6lo)))
c6_NH2 = log10(sum(10**c6a))

c_NH2 = log10(sum(10**(0.5*(2*c5 + c5hi - c5lo))) +
              sum(10**(0.5*(2*c6 + c6hi - c6lo))) )
c_NH2lo = log10(sum(10**(c5-c5lo)) + sum(10**(c6-c6lo)))
c_NH2hi = log10(sum(10**(c5+c5hi)) + sum(10**(c6+c6lo)))

NHI = 19.5
NHIlo = NHI-0.2
NHIhi = NHI+0.2

c_NHI = NHI
c_NHIlo = c_NHI - 0.2
c_NHIhi = c_NHI + 0.2

fval = 2*10**c_NH2 / (10**c_NHI + 2*10**c_NH2)

# calculate error assuming 2*NH2 is negligible compared to NHI, then
#
# logf = np.log10(2) + c_NH2 - c_NHI
#
# sigma_logf = np.sqrt(sigmac_NH2**2 + sigmac_NHI**2)
#
# sigma_logf = np.sqrt(0.3**2 + 0.2**2) = 0.36


if 1:
    pl.rc('text', usetex=True)

    ###################################################################
    # f = 2*NH2 / (NHI + 2*NH2)
    ###################################################################
    F = adict()
    xaxis = 'Ntot'   # 'NHI', 'Ntot'

    pl.rc('font', size=13)
    pl.rc('legend', fontsize=8)
    #fig = pl.figure(figsize=(4.3, 4))
    #fig.subplots_adjust(top=0.99, right=0.99, left=0.13, bottom=0.115)
    #ax = pl.gca()

    fig = pl.figure(figsize=(7.6, 4))
    fig.subplots_adjust(top=0.99, right=0.99, left=0.07, bottom=0.115,
                        wspace=0.0001)
    ax = fig.add_subplot(121)

    m = MWplane[(MWplane.lNH2 == '-') & (MWplane.lNHI == '-')]
    d = DLAs[DLAs.lNH2 == '']

    # conservative lower limit for our sub-DLA
    cN2lolim = 16.46

    p02_NHI = 19.4
    p02_NH2 = np.log10(2.39e16 + 3.46e16)

    # sembach 2001 leading arm Z~0.3
    s01_NHI = 19.9
    s01_NH2 = 16.8

    # Richter 2001 central area of MS  Z~0.3
    #
    #5.4e-4 = 2*NH2 / (NHI + 2*NH2)
    #NHI = (2*NH2) / 5.4e-4  - 2*NH2
    r01_NHI = 20.0
    r01_NH2 = 16.4

    # Richter 1999 HVC in front of LMC? Z~0.5
    r99_NH2 = 15.45        # pm 0.1
    r99_NHI = 19.08


    if xaxis == 'NHI':
        F.MWplane = m['NHI'], 2*10**m['NH2'] / (10**m['NHI'] + 2*10**m['NH2'])
        F.LMC = LMC.NHI, 2*10**LMC.NH2 / (10**LMC.NHI + 2*10**LMC.NH2)
        F.SMC = SMC.NHI, 2*10**SMC.NH2 / (10**SMC.NHI + 2*10**SMC.NH2)
        F.C5 = c5_NHI, 2*10**c5_NH2 / (10**c5_NHI + 2*10**c5_NH2)
        F.C6 = c6_NHI, 2*10**c6_NH2 / (10**c6_NHI + 2*10**c6_NH2)
        F.Clo = NHI, 2*10**cN2lolim / (10**NHI + 2*10**cN2lolim)
        F.IVC = IVC.NHI, 2*10**IVC.NH2 / (10**IVC.NHI + 2*10**IVC.NH2)
        F.DLA = d.NHI, 2*10**d.NH2 / (10**d.NHI + 2*10**d.NH2)
        F.p02 = p02_NHI, 2*10**p02_NH2 / (10**p02_NHI + 2*10**p02_NH2)
        F.s01 = s01_NHI, 2*10**s01_NH2 / (10**s01_NHI + 2*10**s01_NH2)
        F.r01 = r01_NHI, 2*10**r01_NH2 / (10**r01_NHI + 2*10**r01_NH2)
    elif xaxis == 'Ntot':
        F.MWplane = np.log10(10**m['NHI'] + 2*10**m['NH2']), 2*10**m['NH2'] / (10**m['NHI'] + 2*10**m['NH2']),MWplane.EBmV
        F.LMC =     np.log10(10**LMC.NHI + 2*10**LMC.NH2),2*10**LMC.NH2 / (10**LMC.NHI + 2*10**LMC.NH2)      ,LMC['E(BV)LMC'] 
        F.SMC =     np.log10(10**SMC.NHI + 2*10**SMC.NH2),2*10**SMC.NH2 / (10**SMC.NHI + 2*10**SMC.NH2)      ,SMC['E(BV)SMC'] 
        #F.C5 =      np.log10(10**c5_NHI + 2*10**c5_NH2)    ,2*10**c5_NH2 / (10**c5_NHI + 2*10**c5_NH2)           
        #F.C6 =      np.log10(10**c6_NHI + 2*10**c6_NH2)    ,2*10**c6_NH2 / (10**c6_NHI + 2*10**c6_NH2)
        F.C =      np.log10(10**c_NHI + 2*10**c_NH2)    ,2*10**c_NH2 / (10**c_NHI + 2*10**c_NH2)
        #F.Clo = np.log10(10**NHI + 2*10**cN2lolim), 2*10**cN2lolim / (10**NHI + 2*10**cN2lolim)            
        F.DLA = np.log10(10**d.NHI + 2*10**d.NH2), 2*10**d.NH2 / (10**d.NHI + 2*10**d.NH2)                   
        F.IVC = np.log10(10**IVC.NHI + 2*10**IVC.NH2), 2*10**IVC.NH2 / (10**IVC.NHI + 2*10**IVC.NH2)         
        F.p02 = np.log10(10**p02_NHI + 2*10**p02_NH2), 2*10**p02_NH2 / (10**p02_NHI + 2*10**p02_NH2)         
        F.s01 = np.log10(10**s01_NHI + 2*10**s01_NH2), 2*10**s01_NH2 / (10**s01_NHI + 2*10**s01_NH2)         
        F.r01 = np.log10(10**r01_NHI + 2*10**r01_NH2), 2*10**r01_NH2 / (10**r01_NHI + 2*10**r01_NH2)
        F.r99 = np.log10(10**r99_NHI + 2*10**r99_NH2), 2*10**r99_NH2 / (10**r99_NHI + 2*10**r99_NH2)         

    if xaxis == 'Ntot':
        ax.scatter(F.MWplane[0],np.log10(F.MWplane[1]), marker='v', c='c', s=35, linewidths=0.5, label='$\mathrm{In MW plane}$')
        ax.scatter(F.IVC[0],    np.log10(F.IVC[1]    ), marker='^', c='b', s=35, linewidths=0.5, label='$\mathrm{IVC}$')
        ax.scatter(F.s01[0],    np.log10(F.s01[1]    ), marker='*', c='m', s=85, linewidths=0.5, label='$\mathrm{Magellanic stream}$')
        ax.scatter(F.r01[0],    np.log10(F.r01[1]    ), marker='*', c='m', s=75, linewidths=0.5) 
        ax.scatter(F.r99[0],    np.log10(F.r99[1]    ), marker='p', c='m', s=55, linewidths=0.5, label='$\mathrm{HVC}$')
        ax.scatter(F.SMC[0],    np.log10(F.SMC[1]    ), marker='D', c='r', s=20, linewidths=0.5, label='$\mathrm{SMC}$')
        ax.scatter(F.LMC[0],    np.log10(F.LMC[1]    ), marker='s', c='g', s=20, linewidths=0.5, label='$\mathrm{LMC}$')
        ax.scatter(F.DLA[0],    np.log10(F.DLA[1]    ), marker='o', c='none',s=25, linewidths=0.5, label=r'$z > 1.5\ \mathrm{DLA}$')
        #draw_arrows(F.p02[0],    np.log10(F.p02[1]), ax, ms=1.5, capsize=4)
        ax.scatter(F.p02[0],    np.log10(F.p02[1]    ), marker='o', c='w',s=25, linewidths=0.5)
        #ax.scatter(F.C5[0],     np.log10(F.C5[1]     ), marker='o', edgecolors='w', c='k',s=75, zorder=10,linewidths=0.5)
        #ax.scatter(F.C6[0],     np.log10(F.C6[1]     ), marker='o', edgecolors='w', c='k',s=75, zorder=10,linewidths=0.5, label='This paper')
        ax.scatter(F.C[0],     np.log10(F.C[1]     ), marker='o', c='k',s=70, zorder=10,linewidths=0.5, label='$\mathrm{This paper}')

        # error bars for measurements from this paper
        # NHI error (diagonal)
        x = np.log10(10**c_NHIlo + 2*10**c_NH2), np.log10(10**c_NHIhi + 2*10**c_NH2)
        y = (np.log10(2*10**c_NH2 / (10**c_NHIlo + 2*10**c_NH2)),
                 np.log10(2*10**c_NH2 / (10**c_NHIhi + 2*10**c_NH2)))
        ax.plot(x, y, '-k',zorder=1)
        # NH2 error (vertical)
        x = 2 * [np.log10(10**c_NHI + 2*10**c_NH2)]
        y = [np.log10(F.C[1]) - 0.36, np.log10(F.C[1]) + 0.36]
        ax.plot(x, y, '-k',zorder=1)


        # x = np.log10(10**c5_NHIlo + 2*10**c5_NH2), np.log10(10**c5_NHIhi + 2*10**c5_NH2)
        # y = (np.log10(2*10**c5_NH2 / (10**NHIlo + 2*10**c5_NH2)),
        #          np.log10(2*10**c5_NH2 / (10**NHIhi + 2*10**c5_NH2)))
        # ax.plot(x, y, '-k',zorder=1)
        # x = np.log10(10**c6_NHIlo + 2*10**c6_NH2), np.log10(10**c6_NHIhi + 2*10**c6_NH2)
        # y = (np.log10(2*10**c6_NH2 / (10**c6_NHIlo + 2*10**c6_NH2)),
        #          np.log10(2*10**c6_NH2 / (10**c6_NHIhi + 2*10**c6_NH2)))
        # ax.plot(x, y, '-k',zorder=1)
        # x = 2 * [np.log10(10**c5_NHI + 2*10**c5_NH2)]
        # y = [np.log10(2*10**c5_NH2lo / (10**c5_NHI + 2*10**c5_NH2lo)),
        #      np.log10(2*10**c5_NH2hi / (10**c5_NHI + 2*10**c5_NH2hi))]
        # ax.plot(x, y, '-k',zorder=1)
        # x = 2 * [np.log10(10**c6_NHI + 2*10**c6_NH2)]
        # y = [np.log10(2*10**c6_NH2lo / (10**c6_NHI + 2*10**c6_NH2lo)),
        #      np.log10(2*10**c6_NH2hi / (10**c6_NHI + 2*10**c6_NH2hi))]
        # ax.plot(x, y, '-k',zorder=1)

    else:
        ax.plot(F.MWplane[0], np.log10(F.MWplane[1]), 'cv', ms=5,  label='$\mathrm{In MW plane}$')
        ax.plot(F.IVC[0],     np.log10(F.IVC[1]    ), 'b^', ms=5, label= '$\mathrm{IVC}$')
        ax.plot(F.s01[0],     np.log10(F.s01[1]    ), 'mp', ms=7)
        ax.plot(F.r01[0],     np.log10(F.r01[1]    ), 'mp', ms=7)
        ax.plot(F.r99[0],     np.log10(F.r99[1]    ), 'mp', ms=7,   label='$\mathrm{HVC}$')
        ax.plot(F.SMC[0],     np.log10(F.SMC[1]    ), 'rD', ms=4.5, label='$\mathrm{SMC}$')
        ax.plot(F.LMC[0],     np.log10(F.LMC[1]    ), 'gs', ms=5,   label='$\mathrm{LMC}$')
        ax.plot(F.DLA[0],     np.log10(F.DLA[1]    ), 'ok', ms=6, mfc='none', label='\mathrm{z>1.5 DLA}$')
        draw_arrows(F.p02[0],     np.log10(F.p02[1]), ax, ms=1.5, capsize=4)
        ax.plot(F.p02[0],     np.log10(F.p02[1]    ), 'ok', ms=6, mfc='w')
        ax.plot(F.C5[0],      np.log10(F.C5[1]     ), '*k', ms=12)
        ax.plot(F.C6[0],      np.log10(F.C6[1]     ), '*k', ms=12,label='$\mathrm{z=0.56 sub-DLA}$')
        #col = ax.scatter(F.Clo[0],   np.log10(F.Clo[1]     ), marker='o', c='w',s=40, 
        #           label='z=0.56 sub-DLA (lower limit)')
        #col.set_linestyles([(0, (2, 1))])

    if 0:#xaxis != 'NH2':
        bins = np.arange(-9.1, 2.01, 0.6)
        y,_ = np.histogram(np.log10(np.concatenate(
            [np.atleast_1d(F[k][1]) for k in F])), bins)
        b = np.repeat(bins, 2)
        Y = np.concatenate([[0], np.repeat(y,2), [0]])
        Y = 0.2 * Y / Y.max()
        import matplotlib.transforms as mtransforms
        trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
        ax.plot(Y, b, color='0.5', transform=trans)

    if xaxis == 'NHI':
        ax.set_xlabel('$N(\mathrm{HI})')
    elif xaxis == 'Ntot':      
        ax.set_xlabel('$\log_{10}\ [N_\mathrm{HI}\ +\ 2N_{\mathrm{H}_2}]\ (\mathrm{cm}^{-2})$')

    if xaxis == 'Ntot':
        Nvals, fvals = detlim(np.linspace(17,23), NH2=15.)
        im = shade_to_line(Nvals, fvals, blend=2., y0=-10, color='0.7')
 
        puttext(0.05, -6.8, '$\mathrm{Detection\ limit}$', pl.gca(), rotation=-19, fontsize=14,
                ycoord='data') #'$N$(H$_2)\ < \ 14.5$'

        #pl.plot(Nvals, fvals, '-', label='Nh2=12.5')
        x0,y0 = NHIlo, 0.7
        
        # typical error bars for other points
        # NHI
        pl.plot([x0-0.21, x0+0.21], [y0+0.21, y0-0.21], 'k', lw=0.5)
        pl.plot([x0, x0], [y0-0.3, y0+0.3], 'k', lw=0.5)

    #pl.ylabel('$N$(H$_2$)/($N$(HI) $+\ 2N$(H$_2$)')
    pl.ylabel('$\log_{10}\ f_{\mathrm{H}_2}$')
    if xaxis == 'Ntot':      
        pl.xlim(18.7, 22.15)
    ax.set_ylim(-8.2, 1.5)
    #leg = ax.legend(loc='upper left', ncol=2, frameon=0, scatterpoints=1)
    #if xaxis == 'Ntot':
    #    for h in leg.legendHandles:
    #        s = h._sizes[0]
    #        h._sizes = (s,s,s)

    ls = ['solid']*3
    colors= ['0.3']*3
    if not os.path.lexists('model_fH2a.sav'):
        Ntot = 10**np.arange(16, 23, 0.0005)
        tosave = []
        for i,zfrac in enumerate([0.2, 1]):
            Z = np.log10(zfrac)
            f_H2 = calc_fH2(Z, Ntot, sigma_d_m21=1, R_m16p5=1, phiCNM=3, Z_SN=0.0204)
            c0 = f_H2 > 0
            NH2 = 0.5 * Ntot[c0] * f_H2[c0]
            NHI = Ntot[c0] - 2*NH2
            x = [np.log10(Ntot[c0])[0]] + list(np.log10(Ntot[c0]))
            y = [-4] + list(np.log10(f_H2[c0]))
            ax.plot(x, y, color=colors[i], ls=ls[i], lw=0.5)
            tosave.append((x, y))
        saveobj('model_fH2a.sav', tosave)
    else:
        model = loadobj('model_fH2a.sav')
        for i,(x,y) in enumerate(model):
            x = np.array(x)
            y = np.array(y)
            c0 = np.array(y) > -4
            ax.plot(x[c0], y[c0], color=colors[i], ls=ls[i], lw=0.5)
            #ax.plot(x, y, color='0.3', ls=ls[i], lw=0.5)

    fH2_models = loadobj('models.sav')
    for nkJ, NH2, Ntot, fH2 in fH2_models:
        x = np.log10(NH2)
        c0 = x < 17.
        ax.plot(np.log10(Ntot)[c0], np.log10(fH2)[c0], '0.3', lw=0.5)

    ax.text(21.5, 0.5, '$Z_{\odot}$',ha='right',fontsize=10)
    ax.plot([21.48, 21.6], [0.36, -0.1], color='0.3', lw=0.5)
    ax.text(22.0, 0.5, '$0.2 Z_{\odot}$', ha='right', fontsize=10)
    ax.plot([21.78, 21.9], [0.36, -0.3], color='0.3', lw=0.5)
        
    #leg.get_frame().set_lw(0.5)
    
    #if xaxis == 'Ntot':
    #    pl.savefig('pics/f_vs_ntot.png', dpi=300)
    #    pl.savefig('pics/f_vs_ntot.pdf')
    #else:
    #    pl.savefig('pics/f_vs_HI.png', dpi=300)
    #    pl.savefig('pics/f_vs_HI.pdf')

    
    pl.show()


if 1:
    F.MWplane = m['NH2'], 2*10**m['NH2'] / (10**m['NHI'] + 2*10**m['NH2'])
    F.LMC = LMC.NH2, 2*10**LMC.NH2 / (10**LMC.NHI + 2*10**LMC.NH2)
    F.SMC = SMC.NH2, 2*10**SMC.NH2 / (10**SMC.NHI + 2*10**SMC.NH2)
    # F.C5 = c5_NH2, 2*10**c5_NH2 / (10**c5_NHI + 2*10**c5_NH2)
    # F.C5hi = c5_NH2, 2*10**c5_NH2 / (10**c5_NHIlo + 2*10**c5_NH2)
    # F.C5lo = c5_NH2, 2*10**c5_NH2 / (10**c5_NHIhi + 2*10**c5_NH2)
    # F.C6 = c6_NH2, 2*10**c6_NH2 / (10**c6_NHI + 2*10**c6_NH2)
    # F.C6hi = c6_NH2, 2*10**c6_NH2 / (10**c6_NHIlo + 2*10**c6_NH2)
    # F.C6lo = c6_NH2, 2*10**c6_NH2 / (10**c6_NHIhi + 2*10**c6_NH2)
    F.C = c_NH2, 2*10**c_NH2 / (10**c_NHI + 2*10**c_NH2)
    F.Chi = c_NH2, 2*10**c_NH2 / (10**c_NHIlo + 2*10**c_NH2)
    F.Clo = c_NH2, 2*10**c_NH2 / (10**c_NHIhi + 2*10**c_NH2)
    F.IVC = IVC.NH2, 2*10**IVC.NH2 / (10**IVC.NHI + 2*10**IVC.NH2)
    F.DLA = d.NH2, 2*10**d.NH2 / (10**d.NHI + 2*10**d.NH2)
    F.p02 = p02_NH2, 2*10**p02_NH2 / (10**p02_NHI + 2*10**p02_NH2)
    F.s01 = s01_NH2, 2*10**s01_NH2 / (10**s01_NHI + 2*10**s01_NH2)
    F.r01 = r01_NH2, 2*10**r01_NH2 / (10**r01_NHI + 2*10**r01_NH2)
    F.r99 = r99_NH2, 2*10**r99_NH2 / (10**r99_NHI + 2*10**r99_NH2)

    # fit line to MW plane data
    x0,y0 = F.MWplane
    c0 = ~np.isnan(y0)
    x,y = x0[c0], y0[c0]
    i = x.argsort()
    coeff = np.polyfit(x[i], y[i], 1)
    
    #ax = pl.gca()
    ax = fig.add_subplot(122)

    Nvals, fvals = detlim(np.linspace(17,23), NH2=15.)
    Nvals = np.ones(50)*15.
    fvals = np.linspace(-9, 1.5)
    im = shade_to_line_vert(fvals, Nvals, blend=2., x0=-2, color='0.7')

    ax.plot(F.MWplane[0], np.log10(F.MWplane[1]), 'cv', ms=5,  label='$\mathrm{In\ MW\ plane}$')
    ax.plot(F.IVC[0],     np.log10(F.IVC[1]    ), 'b^', ms=5,  label='$\mathrm{IVC}$')
    ax.plot(F.s01[0],     np.log10(F.s01[1]    ), 'm*', ms=8,  label='$\mathrm{Magellanic\ \ Stream}$')
    ax.plot(F.r01[0],     np.log10(F.r01[1]    ), 'm*', ms=8)
    ax.plot(F.r99[0],     np.log10(F.r99[1]    ), 'mp', ms=7,  label='$\mathrm{HVC}$')
    ax.plot(F.SMC[0],     np.log10(F.SMC[1]    ), 'rD', ms=4.5,label='$\mathrm{SMC}$')
    ax.plot(F.LMC[0],     np.log10(F.LMC[1]    ), 'gs', ms=5,  label='$\mathrm{LMC}$')
    ax.plot(F.DLA[0],     np.log10(F.DLA[1]    ), 'ok', ms=5, mfc='none', label='$z>1.5\ \ \mathrm{DLA/sub\,DLA}$')
    #draw_arrows(F.p02[0],     np.log10(F.p02[1]), ax, ms=1.5, capsize=4)
    ax.plot(F.p02[0],     np.log10(F.p02[1]    ), 'ok', ms=5, mfc='w')
    # ax.scatter(F.C5[0],     np.log10(F.C5[1]     ), marker='o', c='k',s=75, edgecolors='w', zorder=10,linewidths=0.5)
    # ax.scatter(F.C6[0],     np.log10(F.C6[1]     ), marker='o', c='k',s=75, edgecolors='w', zorder=10,linewidths=0.5, label='This paper')        
    ax.scatter(F.C[0],     np.log10(F.C[1]     ), marker='o', c='k',s=75, zorder=10,linewidths=0.5, label='$\mathrm{This\ paper}$')        
    #col = ax.scatter(F.Clo[0],   np.log10(F.Clo[1]     ), marker='o', c='w',s=40, 
    #           label='z=0.56 sub-DLA (lower limit)')
    #col.set_linestyles([(0, (2, 1))])

    # just shows the errors from NHI, because Nh2 errors are along
    # direction parallel to relation.

    ax.plot((c_NH2lo, c_NH2hi),
            (np.log10(2*10**c_NH2lo / (10**c_NHI + 2*10**c_NH2lo)),
             np.log10(2*10**c_NH2hi / (10**c_NHI + 2*10**c_NH2hi))), '-k')
    ax.plot((c_NH2, c_NH2),
            [np.log10(F.C[1]) - 0.36, np.log10(F.C[1]) + 0.36], '-k')

    # ax.plot((c5_NH2lo, c5_NH2hi),
    #         (np.log10(2*10**c5_NH2lo / (10**c5_NHI + 2*10**c5_NH2lo)),
    #          np.log10(2*10**c5_NH2hi / (10**c5_NHI + 2*10**c5_NH2hi))), '-k')
    # ax.plot((c6_NH2lo, c6_NH2hi),
    #         (np.log10(2*10**c6_NH2lo / (10**c6_NHI + 2*10**c6_NH2lo)),
    #          np.log10(2*10**c6_NH2hi / (10**c6_NHI + 2*10**c6_NH2hi))), '-k')
    # ax.plot((c5_NH2, c5_NH2),
    #         (np.log10(F.C5hi[1]), np.log10(F.C5lo[1])), '-k')
    # ax.plot((c6_NH2, c6_NH2),
    #         (np.log10(F.C6hi[1]), np.log10(F.C6lo[1])), '-k')

    ls = ['solid']*3

    if not os.path.lexists('model_fH2b.sav'):
        Ntot = 10**np.concatenate([np.arange(20.9, 22.8, 0.0001)])
        tosave = []
        for i,zfrac in enumerate([0.2, 1]):
            Z = np.log10(zfrac)
            f_H2 = calc_fH2(Z, Ntot, sigma_d_m21=1, R_m16p5=1,
                            phiCNM=3, Z_SN=0.0204)
            c0 = f_H2 > 0
            NH2 = 0.5 * Ntot[c0] * f_H2[c0]
            NHI = Ntot[c0] - 2*NH2
            x = np.log10(NH2)
            y = np.log10(f_H2[c0])
            X = np.arange(17, 21.5, 0.05)
            Y = np.interp(X, x, y)
            ax.plot(X, Y, color=colors[i], ls=ls[i], lw=0.5)
            tosave.append((X, Y))
        saveobj('model_fH2b.sav', tosave)
    else:
        model = loadobj('model_fH2b.sav')
        for i,(x,y) in enumerate(model):
            c0 = x > 17.8
            ax.plot(x[c0],y[c0], color=colors[i], ls=ls[i], lw=0.5)

    for nkJ, NH2, Ntot, fH2 in fH2_models:
        x = np.log10(NH2)
        c0 = x < 17.
        ax.plot(np.log10(NH2)[c0], np.log10(fH2)[c0], '0.3', lw=0.5)


    ax.set_xlabel('$\log_{10}\, N_{\mathrm{H}_2}\ (\mathrm{cm}^{-2})$')

    #pl.ylabel('$\log_{10}\ f\ (\mathrm{H}_2$)')
    leg = ax.legend(loc='upper left', frameon=0, scatterpoints=1)
    #leg.get_frame().set_lw(0.5)

    #leg = ax.legend(loc='upper left', ncol=2, frameon=0, scatterpoints=1)
    #if xaxis == 'Ntot':
    #    for h in leg.legendHandles:
    #        s = h._sizes[0]
    #        h._sizes = (s,s,s)

    ax.set_xlim(12.21, 21.49)
    #ax.set_ylim(-7.2, 0.499)
    ax.set_ylim(-8.2, 1.5)
    ax.set_yticklabels([])
        
    leg.get_frame().set_lw(0.5)

    x0,y0 = 20, 0.6

    ax.plot([x0-0.21, x0+0.21], [y0-0.21, y0+0.21], 'k', lw=0.5)
    # NHI
    ax.plot([x0, x0], [y0-0.3, y0+0.3], 'k', lw=0.5)

    #ax.text(13.4, 0.5, '$Z_{\odot}$',ha='right',fontsize=10)
    #ax.text(14.3, 0.5, '$0.2 Z_{\odot}$', fontsize=10)

    pl.savefig('pics/f_vs_ntot.png', dpi=300)
    pl.savefig('pics/f_vs_ntot.pdf')

    pl.show()



if 1:
    ###################################################################
    # f(z)
    ###################################################################
    F = adict()

    m = MWplane[(MWplane.lNH2 == '-') & (MWplane.lNHI == '-')]
    d = DLAs[DLAs.lNH2 == ''].view(np.recarray)

    F.MWplane = [0]*len(m), 10**m['NH2'] / (10**m['NHI'] + 2*10**m['NH2'])
    F.LMC = [0]*len(LMC), 10**LMC.NH2 / (10**LMC.NHI + 2*10**LMC.NH2)
    F.SMC = [0]*len(SMC), 10**SMC.NH2 / (10**SMC.NHI + 2*10**SMC.NH2)
    # F.C5 = 0.567, 10**c5_NH2 / (10**NHI + 2*10**c5_NH2)
    # F.C6 = 0.567, 10**c6_NH2 / (10**NHI + 2*10**c6_NH2)
    F.C = 0.567, 2*10**c_NH2 / (10**c_NHI + 2*10**c_NH2)
    F.DLA = d.zabs, 10**d.NH2 / (10**d.NHI + 2*10**d.NH2)

    pl.rc('font', size=14)
    #pl.rc('legend', fontsize=12)
    fig = pl.figure(figsize=(3.5, 3.7))
    fig.subplots_adjust(top=0.88, right=0.97, left=0.165, bottom=0.12)
    ax = pl.gca()
    fvals = np.array(list(np.log10(F.MWplane[1])) + list(np.log10(F.SMC[1])) + \
                              list(np.log10(F.LMC[1])))
    f0,f1,f2 = np.percentile(fvals[fvals > -20], [10,50, 90])

    ax.plot(0, f1, 'r^', ms=8, mew=0, label='$\mathrm{Local\ group}$',)
    ax.errorbar(0, f1, yerr=np.transpose([(f1-f0, f2-f1)]),
                ecolor='r', capsize=3, mew=1, fmt=None) 
    
    ax.plot(WMAP7.lookback_time(F.DLA[0]),np.log10(F.DLA[1]), 'ok', ms=5, mfc='none', mew=0.5, label='$z>1.5\ \ \mathrm{DLA/sub\,DLA}$')
    x = WMAP7.lookback_time(F.C[0])
    y = np.log10(F.C[1])
    ax.scatter(x,y, marker='o', c='k',s=70, zorder=10,linewidths=0.5,label='$\mathrm{This\ paper}$')
    plot([x,x], [y-0.36, y+0.36], 'k')


    ax.set_xlabel('$\mathrm{Lookback\  time\ (Gyr)}$')
    ax.set_ylabel('$\log_{10}\ f_{\mathrm{H}_2}$')
    ax.set_xlim(-0.9, 13.5)
    ax.set_ylim(-8.99, 0.49)

    leg = ax.legend(loc='lower left', frameon=0, scatterpoints=1)
    leg.get_frame().set_lw(0.5)

    ax1 = pl.twiny(ax)
    zvals = np.linspace(0,10,1000)
    times = WMAP7.lookback_time(zvals)
    zminorticks = list(np.arange(0,1.05,0.1)) + list(np.arange(2., 10, 1.))
    minorGyrticks = np.interp(zminorticks, zvals, times)
    zmajorticks = [0, 0.5, 1, 5]
    majorGyrticks = np.interp(zmajorticks, zvals, times)
    ax1.xaxis.set_ticks(minorGyrticks, minor=1)
    zlabels = [('$%.1g$' % t if '%.1g' % t in '0 0.5 1 2 5'.split()
                else '' ) for t in zminorticks ]
    ax1.xaxis.set_ticklabels(zlabels, minor=1)
    ax1.xaxis.set_ticks(majorGyrticks)
    ax1.xaxis.set_ticklabels('')
    ax1.set_xlim(*ax.get_xlim())

    ax1.set_xlabel('$\mathrm{Redshift}$')

    pl.savefig('pics/fz.png', dpi=300)
    pl.savefig('pics/fz.pdf')
    pl.show()
