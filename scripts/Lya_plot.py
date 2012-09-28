import barak.spec
from barak.absorb import find_tau, readatom
from barak.convolve import convolve_psf

lines =  [('HI', 0.55725616, 20, 19.5)]
lineslo =  [('HI', 0.55725616, 20, 19.3)]
lineshi =  [('HI', 0.55725616, 20, 19.7)]
atom = readatom()

pl.rc('font', size=12.5)

if 1:
    sp = barak.spec.read('fos_colo2.txt')
    sp.fl *= 1e16
    sp.co *= 1e16
    sp.er *= 1e16
    #sp = barak.spec.read('../vpfit_COS/qsoc_G190H.txt')
    wa = sp.wa + 0.2
    dw = np.median(np.diff(sp.wa))

    # resolution seems to be 0.8 expected value? Or sign of lower N
    # comoponents in wings?
    npix = 1.39 / dw * 0.8

    tau, ticks = find_tau(wa, lines, atom)
    taulo, ticks = find_tau(wa, lineslo, atom)
    tauhi, ticks = find_tau(wa, lineshi, atom)

    fig = pl.figure(figsize=(4.3, 3))
    fig.subplots_adjust(bottom=0.16, left=0.16, right=0.98, top=0.96)


    ax = fig.add_subplot(111)

    ax.plot(wa, sp.fl,color='0.3', lw=0.5,drawstyle='steps-mid')
    ax.plot(wa, sp.co, '0.7', ls='dashed')
    ax.plot(wa, sp.co*convolve_psf(np.exp(-taulo), npix), 'r', lw=0.3)
    ax.plot(wa, sp.co*convolve_psf(np.exp(-tau), npix), 'r')
    ax.plot(wa, sp.co*convolve_psf(np.exp(-tauhi), npix), 'r', lw=0.3)
    ax.axhline(0, color='k', lw=0.3)
    ax.minorticks_on()
    ax.set_xlim(1881, 1904.9)
    #ax.set_xlim(1861, 1924.9)
    ymax = 3.49
    resid = (sp.fl - sp.co*convolve_psf(np.exp(-tau), npix)) / sp.er
    ax.axhline((-0.12-0.03)*ymax, color='0.7', lw=0.5)
    ax.axhline((-0.12+0.03)*ymax, color='0.7', lw=0.5)
    ax.scatter(wa, -0.12*ymax + 0.03*ymax*resid, c='g', marker='.', s=20,
               faceted=False)
    ax.set_ylim(-0.22*ymax, ymax)
    ax.set_xlabel('$\mathrm{Wavelength}\ (\AA)$')
    pl.figtext(0.012, 0.55,
               r'$\mathrm{Flux} \times\, 10^{16}\ \ (\mathrm{erg}\,\mathrm{s}^{-1}\mathrm{cm}^{-2}\AA^{-1})$',
               rotation=90, va='center')
    puttext(0.95, 0.25, r'$\mathrm{HI\ Ly}\alpha$', ax, fontsize=16, ha='right')
    pl.savefig('pics/NHI.png', dpi=300)
    pl.savefig('pics/NHI.pdf')
    pl.show()

