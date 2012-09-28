import velplot.velplot
reload(velplot.velplot)
from velplot.velplot import make_plot, make_models
from barak.io import parse_config
from barak.constants import c_kms
from barak.convolve import convolve_psf

cfg = """\
f26name = all.f26
spfilename = Cjoin.txt
trfilename = transitions/paper
vmax = 299
vmin = None
wadiv = None
Rfwhm = convolve_with_COS_FOS
osc = False
residuals = 0
redshift = 0.5572884888
"""

unrelated =[
    #(1512.490, 1512.698  ),
    (1513.431,1515.614),
  ( 1608.062  , 1608.772  ), 
  (  1605.440 ,   1606.218), 
  ( 1688.709  , 1689.313  ), 
  ( 1690.611  , 1691.020  ), 
  ( 1472.989  , 1474.508  ), 
  (  1470.121 ,   1471.004), 
  (  1587.195 ,   1588.116), 
  ( 1590.794  , 1591.765  ), 
  (  1618.654 ,1620.563   ), 
  (   1615.927, 1617.836  ),
   ( 1614.409, 1615.820), 
   ( 1519.939, 1520.820), 
   ( 1539.532, 1541.674), 
   ( 1542.015, 1544.012), 
    ] 


cfgname = 'velplot.cfg'

with open(cfgname, 'w') as fh:
    fh.write(cfg)


opt = parse_config(cfgname)

pl.rc('xtick',labelsize=8.5)
pl.rc('ytick',labelsize=8.5)
pl.rc('xtick.major', size=3)      # major tick size in points
pl.rc('xtick.minor', size=1.5)      # minor tick size in points
pl.rc('xtick.major', pad=3 )     # distance to major tick label in points
pl.rc('xtick.minor', pad=3 )     # distance to the minor tick label in points
pl.rc('ytick.major', size=3)      # major tick size in points
pl.rc('ytick.minor', size=1.5)      # minor tick size in points
pl.rc('ytick.major', pad=3 )     # distance to major tick label in points
pl.rc('ytick.minor', pad=3 )     # distance to the minor tick label in points

if 1:
    sp, transitions, model, models, ticks = make_models(opt)
    #sp.fl = convolve_psf(sp.fl, 1.5)
    cols = dict(model='r', resid='g', data='0.5')

    fig, axes = make_plot(sp, transitions, model, models, [], opt, width=4.2,
                          unrelated=unrelated, ncols=2, nrows=5,cols=cols,
                          aspect=0.6)

    fig.subplots_adjust(bottom=0.065,left=0.11)
    pl.figtext(0.018, 0.5, 'Transmission', rotation=90, va='center',fontsize=11)
    tickpos =  array([-110.26510868,  -93.01127355,  -72.24951539,
                      -49.79515804,
                      -26.07570683,    0.        ,   44.06937431,
                      68.34348321,
                      100.02608005,  124.18839882])

    for i,ax in enumerate(axes):
        ax.set_yticks([0., 0.5, 1])
        if i > 4:
            ax.set_yticklabels([])
        else:
            ax.set_yticklabels(['0', '0.5', '1'])

        if i not in (4, 9):
            #if i in (4,):
            #    ax.set_xlabel('Velocity offset (km s$^{-1}$)', fontsize=11)
            ax.set_xticklabels([])

        if i == 8:
            ax.text(120, 0.1,'(FOS)', fontsize=13)
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        ax.set_ylim(-0.2, 1.5)
        for j,v in enumerate(tickpos):
            if j in (4, 5):
                ax.plot([v, v], [1.1, 1.3], color='k', lw=1.5)
            else:
                ax.plot([v, v], [1.1, 1.3], color='b')

        ax.minorticks_on()

    fig.text(0.55, 0.01, 'Velocity offset (km s$^{-1}$)', fontsize=11,ha='center')
    pl.savefig('COS.pdf')
    pl.show()
