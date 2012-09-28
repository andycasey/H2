import velplot.velplot
reload(velplot.velplot)
from velplot.velplot import make_models, make_plot
from barak.io import parse_config
from barak.constants import c_kms

cfg = """\
f26name = HIRES.f26
spfilename = q0107c_HIRES.txt
trfilename = transitions/general
vmax = 299
vmin = None
wadiv = None
Rfwhm = 6.6
osc = False
residuals = True
redshift = 0.5572884888
"""

cfgname = 'velplot.cfg'

with open(cfgname, 'w') as fh:
    fh.write(cfg)

opt = parse_config(cfgname)

pl.rc('xtick',labelsize=11)
pl.rc('ytick',labelsize=11)
pl.rc('xtick.major', size=3)      # major tick size in points
pl.rc('xtick.minor', size=1.5)      # minor tick size in points
pl.rc('ytick.major', size=3)      # major tick size in points
pl.rc('ytick.minor', size=1.5)      # minor tick size in points

if 1:
    cols = dict(model='r', resid='g', data='0.5')
    sp, transitions, model, models, ticks = make_models(opt)
    fig, axes = make_plot(sp, transitions, model, models, [], opt, width=8,
                          cols=cols)
    fig.subplots_adjust(bottom=0.12,left=0.075)
    pl.figtext(0.014, 0.5, 'Transmission', rotation=90, va='center',fontsize=12)

    mg2796 = 2796.3542699
    waobs = (opt.redshift + 1) * mg2796
    tickpos = (ticks[ticks.wa0 == mg2796].wa / waobs  - 1) * c_kms
    print repr(tickpos)

    for i,ax in enumerate(axes):
        if i > 2:
            ax.set_yticklabels([])
        ax.minorticks_on()
        if i in (2, 5, 8):
            if i in (5,):
                ax.set_xlabel('Velocity offset (km s$^{-1}$)', fontsize=12)
        else:
            ax.set_xticklabels([])
        for j,v in enumerate(tickpos):
            if j in (4, 5):
                ax.plot([v, v], [1.1, 1.3], color='k', lw=1.5)
            else:
                ax.plot([v, v], [1.1, 1.3], color='b')

    pl.savefig('HIRES.pdf')
    pl.savefig('HIRES.png', dpi=400)
    show()
