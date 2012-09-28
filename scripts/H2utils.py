import numpy as np
from barak.io import readtxt
from StringIO import StringIO

"""
z1   0.5571716       3.3328989e-06  2.9235276e-06  6.1819434e-06 5.8875732e-06
z2   0.55731413      6.9815339e-06  5.6234853e-07  1.0778247e-05 3.907031e-06
"""

vals = readtxt(StringIO("""\
N1   15.980229       0.051283581    0.44173713     0.25458389    0.78321895
b1   7.1385398       1.0740454      0.099691753    1.6881329     0.62992256
N2   15.401429       0.15744999     0.61836035     0.28957997    1.4460613
b2   4.4966749       0.88646246     0.50331688     2.491544      1.4914935
N3   16.86705        0.093238909    0.46774437     0.31209545    0.72591071
N4   16.467646       0.45228955     0.35720928     0.66221697    0.96223281
N5   16.184167       0.18668953     0.19701228     0.24673393    0.54682153
N6   15.640367       0.2410867      0.25817481     0.31882251    1.6768936
N7   15.686983       0.030926233    0.19910952     0.11180025    0.39515113
N8   15.398693       0.10942594     0.24660601     0.18772311    1.7104885
"""),
               names='name,logN,siglo,sighi,sig2lo,sig2hi')

# old
# vals = readtxt(StringIO("""\
# N1   16.92    0.6101    0.07884   0.9654    0.3432
# b1   4.314    0.1769    0.7664    0.4878    1.713
# N2   15.26    0.1464    0.3157    0.3125    1.7
# b2   7.421    1.328     0.004784  4.99      0.5069
# N3   17.52    0.6729    0.042     1.078     0.1888
# N4   16.22    0.1401    0.6071    0.2793    1.8
# N5   17.03    0.8219    -0.1849   1.161     0.06227
# N6   15.46    0.01716   0.2225    0.07043   2.17
# N7   16.16    0.3642    0.0795    0.6157    0.4
# N8   15.32    0.05069   0.1174    0.119     2.154"""),
#                names='name,logN,siglo,sighi,sig2lo,sig2hi')


def read_H2(sigma=1):
    """ Read the info about the H2 limits for the z=0.56 system

    for 2 sigma limits, set keyword sigma=2"""
    

    # components 6 and 7

    Jnames = 'N0 N1 N2 N3'.split()
    c6 = np.array([vals[i][1] for i in [0,4,6,8]])
    c6lo = np.array([vals[i][2] for i in [0,4,6,8]])
    c6hi = np.array([vals[i][3] for i in [0,4,6,8]])
    c6a = 0.5 * ((c6 + c6hi) + (c6 - c6lo))
    if sigma == 2:
        c6lo = np.array([vals[i][4] for i in [0,4,6,8]])
        c6hi = np.array([vals[i][5] for i in [0,4,6,8]])
    c7 = np.array([vals[i][1] for i in [2,5,7,9]])
    c7lo = np.array([vals[i][2] for i in [2,5,7,9]])
    c7hi = np.array([vals[i][3] for i in [2,5,7,9]])
    c7a = 0.5 * (c7 + c7hi + (c7 - c7lo))
    if sigma == 2:
        c7lo = np.array([vals[i][4] for i in [2,5,7,9]])
        c7hi = np.array([vals[i][5] for i in [2,5,7,9]])

    return c6,c6lo,c6hi,c6a, c7, c7lo,c7hi,c7a


def H2_table_paper():
    c6,c6lo,c6hi,c6a,c7,c7lo,c7hi,c7a = read_H2(sigma=1)
    #J = 0
    Nfmt = '%5.2f'
    bfmt = '%4.1f'

    fmt1 = [Nfmt]*2 + [bfmt]*2 
    fmt1 = tuple(['%i'] + 2*fmt1)
    fmt2 = tuple(['%i'] + [Nfmt]*4)

    hd1 = (r'& \multicolumn{2}{c}{Component 5} & '
             r'\multicolumn{2}{c}{Component 6} \\')
    hd2 = (r'$J$ & $\log\,N$ & $b$ & $\log\,N$ & $b$ \\')
    f1 = r'%s & $%s\pm%s$ & $%s\pm%s$ & $%s\pm%s$ & $%s\pm%s$ \\' % fmt1
    f2 = r'%s & $%s\pm%s$ &           & $%s\pm%s$ &           \\' % fmt2
    print r"""\begin{table}
\addtolength{\tabcolsep}{-2.5pt}
\begin{center}
\begin{tabular}{ccccc}"""
    print r'\hline'
    print hd1
    print hd2
    print r'\hline'
    b6 = 0.5 * (vals.logN[1] - vals.siglo[1] + vals.sighi[1] + vals.logN[1])
    b6sig = 0.5 * (vals.siglo[1] + vals.sighi[1])
    b7 = 0.5 * (vals.logN[3] - vals.siglo[3] + vals.sighi[3] + vals.logN[3])
    b7sig = 0.5 * (vals.siglo[3] + vals.sighi[3])
    print f1 % (0, c6a[0], 0.5*(c6hi[0] + c6lo[0]), b6, b6sig,
                c7a[0], 0.5*(c7hi[0] + c7lo[0]), b7, b7sig) 
    for i in range(1,4):
        print f2 % (i, c6a[i], 0.5*(c6hi[i]+c6lo[i]),
                    c7a[i], 0.5*(c7hi[i]+c7lo[i]))
    print r"""\hline
\end{tabular}
\caption{\label{t_logN_Htwo} Column densities for the two components
  showing H$_2$.}
\end{center}
\end{table}
"""


def calc_OPR(N, Nlo=None, Nhi=None):
    """ Find the otho-para ratio (OPR).

    Parameters
    ----------
    N : array_like
      Array of log10(H2) column densities, sorted from J=0 to higher J.

    Nhi, Nlo : array_like, optional
      Maximum and minimum log10 column densities.

    Returns
    -------
    OPR : float

    or

    OPR, OPRlo, OPRhi : float
    """

    N = np.array(N)
    if Nhi is not None and Nlo is not None:
        Nhi = np.array(Nhi)
        Nlo = np.array(Nlo)
        assert len(N) == len(Nhi) == len(Nlo)

    # we must have the first two levels!
    assert not np.isnan(N[:2]).any()
    assert (N[:2] > 0).all()

    # remove nan values by setting the column densities to a very low
    # value.
    N[np.isnan(N)] = 0 
    # handle upper limits - set the value equal to the limit
    N = np.where(N > 0, N, -N)

    Neven = N[::2]
    Nodd = N[1::2]

    if Nhi is not None and Nlo is not None:
        Nevenhi = Nhi[::2]
        Nevenlo = Nlo[::2]
        Noddhi = Nhi[1::2]
        Noddlo = Nlo[1::2]
        OPR = sum(10**Nodd) / sum(10**Neven)
        OPRhi = sum(10**Noddhi) / sum(10**Nevenlo)
        OPRlo = sum(10**Noddlo) / sum(10**Nevenhi)
        return OPR, OPRlo, OPRhi

    return sum(10**Nodd) / sum(10**Neven)
    
