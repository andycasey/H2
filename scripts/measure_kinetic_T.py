from scipy.optimize import curve_fit
from H2utils import read_H2

# CI and H_2 are normally observed together. This is because neutral
# carbon and molecular hydrogen are photoionized and photodissociated
# respectively, by photons of the same energies and therefore, they
# are usually found together.

# below value is for vibrational ground state

Bv = 85.36    # K

# g for level J is defined
#
# g(J) = (2*J + 1) * (2*I + 1)
#
# where I is 0 for even J and 1 for odd J.

g = np.array([1, 9, 5, 21, 9, 33])

# According to the Boltzmann distribution (see Equation (8) in
# Levshakov et al. 2002),

#N[J] / N[0] = g[J] / g[0] * exp(-Bv * J * (J+1) / Tex)

# The excitation diagram is typically plotted as log(N_J /g_J) versus
# the relative energy between the level J and J = 0. The excitation
# temperature, defined in Equation (A1), is inversely proportional to
# the negative slope of the line connecting the excitation diagram
# points, i.e.,

# Tex_0J = -Bv * J*(J+1) / np.log ((g[0] * N[J]) / (g[J] * N[0]))

c5,c5lo,c5hi,c5a,c6,c6lo,c6hi,c6a = read_H2(sigma=1)

def straight_line(x, a, b):
    return a*x + b

if 0:
    #########################################
    # Measure OPR assuming LTE
    #########################################

    # OPR = 3 * sum(odd) / sum(even)
    
    Nhi, Nlo, N = c5+c5hi, c5-c5lo, c5a

    Nevenlo = Nlo[::2]
    Noddlo = Nlo[1::2]
    Neven = N[::2]
    Nodd = N[1::2]
    Nevenhi = Nhi[::2]
    Noddhi = Nhi[1::2]

    OPRhi = sum(10**Noddhi) / sum(10**Nevenlo)
    OPRlo =  sum(10**Noddlo) / sum(10**Nevenhi)
    OPR =  sum(10**Nodd) / sum(10**Neven)
    print 'OPR Comp. 6: %.2f %.2f %.2f' % (OPRlo, OPR, OPRhi)

    Nhi, Nlo, N = c6+c6hi, c6-c6lo, c6a

    Nevenlo = Nlo[::2]
    Noddlo = Nlo[1::2]
    Neven = N[::2]
    Nodd = N[1::2]
    Nevenhi = Nhi[::2]
    Noddhi = Nhi[1::2]

    OPRhi = sum(10**Noddhi) / sum(10**Nevenlo)
    OPRlo = sum(10**Noddlo) / sum(10**Nevenhi)
    OPR = sum(10**Nodd) / sum(10**Neven)
    print 'OPR Comp. 7: %.2f %.2f %.2f' % (OPRlo, OPR, OPRhi)

if 0:
    ###################################################################
    # measure NH2 tot vs log10(N[2] / N[0]) and NH2 tot vs log10(N[3]
    # / N[1]) to compare with fig 5 and 6 from Srianand et al 2005.
    ###################################################################
    Nhi, Nlo, N = c5+c5hi, c5-c5lo, c5a
    print 'Comp. 6'
    print '  total log10 N(H_2) %.2f' % log10(sum(10**N))
    print '  log 10 N[2]/N[0]: %.2f %.2f %.2f' % (
        Nlo[2] - Nhi[0], N[2] - N[0], Nhi[2] - Nlo[0])
    print '  log 10 N[3]/N[1]: %.2f %.2f %.2f' % (
        Nlo[3] - Nhi[1], N[3] - N[1], Nhi[3] - Nlo[1])
    Nhi, Nlo, N = c6+c6hi, c6-c6lo, c6a
    print 'Comp. 7'
    print '  total log10 N(H_2) %.2f' % log10(sum(10**N))
    print '  log 10 N[2]/N[0]: %.2f %.2f %.2f' % (
        Nlo[2] - Nhi[0], N[2] - N[0], Nhi[2] - Nlo[0])
    print '  log 10 N[3]/N[1]: %.2f %.2f %.2f' % (
        Nlo[3] - Nhi[1], N[3] - N[1], Nhi[3] - Nlo[1])

if 1:
    #############################################
    # Make the excitation plot for the paper
    #############################################

    pl.rc('text', usetex=True)
    pl.rc('font', size=14)

    # N[J] / N[0] = g[J] / g[0] * exp(-Bv * J * (J+1) / Tex)
    
    # The excitation diagram is typically plotted as log(N_J /g_J)
    # versus the relative energy between the level J and J = 0.

    # so
    #
    # x = Bv*J*(J+1)
    # y = log10(N[J] / g[J])

    # log10(NJ / gJ) = log10( N0 / g0 * exp(-x / Tex) )
    # log10(NJ/gJ) = log10(N0/g0) + log10(exp(-x / Tex))
    # log10(NJ/gJ) = log10(N0/g0) - x / (ln(10) * Tex)

    # So slope is given by -(ln(10) * Tex)^-1
 
    # level 1 - 0

    J = np.arange(6)
    xvals = Bv * J[:4] * (J[:4] + 1)
    xvals2 = np.linspace(-100, 4000)
    
    fig = pl.figure(figsize=(4.3, 4.5))
    fig.subplots_adjust(bottom=0.1, left=0.13, top=0.91, right=0.99)
    ax = pl.gca()

    # use half way between instead of maximum likelihood values.
    Nhi, Nlo, N = c5+c5hi, c5-c5lo, c5a
    yhi = Nhi - log10(g[:4])
    ylo = Nlo - log10(g[:4])
    y =   N - log10(g[:4])
    yerr = yhi - y
    popt, pcov = curve_fit(straight_line, xvals, y, sigma=yerr)
    perr = np.sqrt(np.diag(pcov))
    T = -1 / (np.log(10) * popt[0])
    Thi = -1 / (np.log(10) * popt[0] + perr[0])
    Tlo = -1 / (np.log(10) * popt[0] - perr[0])
    #label = 'T=%.0f +/- %.0f K' % (T, dT) 
    label = r'$T=240K$'

    ax.errorbar(xvals-10, y, yerr=(y-ylo, yhi-y), fmt=None, ecolor='k', lw=0.5, mew=0.5,capsize=0)
    ax.plot(xvals-10, y, 'ob', label='$\mathrm{Component\ 5}$')
    ax.plot(xvals2, straight_line(xvals2, *popt), '-', color='0.5', lw=0.5,
            zorder=1, label=label)
 
    Nhi, Nlo, N = c6+c6hi, c6-c6lo, c6a
    yhi = Nhi - log10(g[:4])
    ylo = Nlo - log10(g[:4])
    y =   N - log10(g[:4])
    yerr = yhi - y
    popt, pcov = curve_fit(straight_line, xvals, y, sigma=yerr)
    perr = np.sqrt(np.diag(pcov))
    T = -1 / (np.log(10) * popt[0])
    Thi = -1 / (np.log(10) * popt[0] + perr[0])
    Tlo = -1 / (np.log(10) * popt[0] - perr[0])
    label = r'$T=290K$'
    #label = 'T=%.0f +/- %.0f K' % (T, dT) 

    #ax.fill_between(xvals, ylo, yhi, lw=0, color='0.9')
    ax.errorbar(xvals+10, y, yerr=(y-ylo, yhi-y), fmt=None, capsize=0, lw=0.5, ecolor='k', mew=0.5)
    ax.plot(xvals+10, y, 'r^', label='$\mathrm{Component\ 6}$')
    ax.plot(xvals2, straight_line(xvals2, *popt), '--', color='0.5', lw=0.5,
            label=label)
    
    ax1 = pl.twiny(pl.gca())
    ax1.set_xticks(Bv * J * (J+1))
    ax1.set_xticklabels(['$0$','$1$','$2$','$3$','$4$','$5$'])
    ax1.set_xlabel('$J$')
    ax.set_ylabel('$\log_{10}\ (N_J\ [\mathrm{cm}^{-2}] /g_J)$')
    ax.set_xlabel('$B_v\ J(J+1)\ \ (\mathrm{K})$')
    ax.set_xlim(-40, 2690)
    ax1.set_xlim(-40, 2690)
    ax.set_ylim(11.5, 17.3)

    ax.legend(frameon=0)

    # upper limits for J=4 and J=5
    # both components, J=4
    # uplim = 14.5 - log10(g[4])
    # J=5
    # uplim6 = 14.5 - log10(g[5])
    # uplim7 = 14.3 - log10(g[5])
    Tvals = Bv * J * (J+1)
    arrow_len = 10.
    capsize = 2
    capsize = min(capsize, arrow_len) 
    arrowdown_verts = [[0.,0.], [0, -arrow_len],
                       [0.5*capsize, -(arrow_len-capsize)], [0, -arrow_len],
                       [-0.5*capsize, -(arrow_len-capsize)], [0, -arrow_len]]
    off = 30
    pl.scatter(Tvals[-2]-off, 14.5 - log10(g[4]), s=1000, marker=None, 
               verts=arrowdown_verts, linewidths=0.5)
    pl.plot(Tvals[-2]-off, 14.5 - log10(g[4]), 'ob')
            
    pl.scatter(Tvals[-2]+off, 14.5 - log10(g[4]), s=1000, marker=None, 
               verts=arrowdown_verts, linewidths=0.5)
    pl.plot(Tvals[-2]+off, 14.5 - log10(g[4]), '^r')

    pl.scatter(Tvals[-1]-off, 14.5 - log10(g[5]), s=1000, marker=None, 
               verts=arrowdown_verts, linewidths=0.5)
    pl.plot(Tvals[-1]-off, 14.5 - log10(g[5]), 'ob')
    pl.scatter(Tvals[-1]+off, 14.3 - log10(g[5]), s=1000, marker=None, 
               verts=arrowdown_verts, linewidths=0.5)
    pl.plot(Tvals[-1]+off, 14.3 - log10(g[5]), '^r')

    pl.show()
    fig.canvas.print_figure('Tex.pdf')
    fig.canvas.print_figure('Tex.png', dpi=300)

def Tex_0J(J, logNJ, logN0):
    """ Find the excitation temperature given log10 of two column
    densities, for level J and and the J=0 level.

    returns temperature in Kelvin
    """
    # note we are taking log to base e here.
    ratio = log(10**logN0 * g[J] / (g[0] * 10**logNJ))
    if ratio < 0:
        return np.inf
    return Bv*J*(J+1) / ratio

if 1:
    #####################################
    # measure individual Tex_0J values
    #####################################

    # log(N[J] / g[J]) = log(N[0]/g[0]) - Bv*J*(J+1)/Tex

    # Tex = Bv*J*(J+1) / (log(N[0] * g[J] / (g[0] * N[J])
    print  'comp. 5'
    Nhi, Nlo, N = c5+c5hi, c5-c5lo, c5a
    for i in range(1,4):
        Thi = Tex_0J(i, Nhi[i], Nlo[0])
        T = Tex_0J(i, N[i], N[0])
        Tlo = Tex_0J(i, Nlo[i], Nhi[0])
        print '  T (K): J%i-0  %5.0f %5.0f %5.0f' % (i, Tlo, T, Thi)

    print  'comp. 6'
    Nhi, Nlo, N = c6+c6hi, c6-c6lo, c6a
    for i in range(1,4):
        Thi = Tex_0J(i, Nhi[i], Nlo[0])
        T = Tex_0J(i, N[i], N[0])
        Tlo = Tex_0J(i, Nlo[i], Nhi[0])
        print '  T (K): J%i-0  %5.0f %5.0f %5.0f' % (i, Tlo, T, Thi)

