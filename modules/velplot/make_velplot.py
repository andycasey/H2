import velplot
reload(velplot)
from velplot import make_models, make_plot
from barak.io import parse_config
from barak.constants import c_kms

cfg = """\
f26name = HIRES.f26
spfilename = q0107c_HIRES.txt
trfilename = transitions/general
vmax = 299
wadiv = None
Rfwhm = 6.6
osc = False
residuals = True
redshift = 0.5571530371
"""

cfgname = 'velplot.cfg'

with open(cfgname, 'w') as fh:
    fh.write(cfg)

opt = parse_config(cfgname)

pl.rc('xtick',labelsize=11)
pl.rc('ytick',labelsize=11)

if 1:
    sp, transitions, model, models, ticks = make_models()
    fig, axes = make_plot(sp, transitions, model, models, [])

    mg2796 = 2796.3542699
    waobs = (opt.redshift + 1) * mg2796
    tickpos = (ticks[ticks.wa0 == mg2796].wa / waobs  - 1) * c_kms

    for i,ax in enumerate(axes):
        if i > 3:
            ax.set_yticklabels([])
        ax.minorticks_on()
        if i in (3, 6, 9):
            if i in (6,):
                ax.set_xlabel('Velocity offset (km s$^{-1}$)', fontsize=12)
        else:
            ax.set_xticklabels([])
        for j,v in enumerate(tickpos):
            if j in (3, 4):
                ax.plot([v, v], [1.1, 1.3], color='k', lw=1.5)
            else:
                ax.plot([v, v], [1.1, 1.3], color='0.5')

    pl.savefig('HIRES.pdf')
