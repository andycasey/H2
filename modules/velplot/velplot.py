from barak.pyvpfit import readf26
from barak.absorb import findtrans, readatom, find_tau
from barak.constants import c_kms
from barak.sed import make_constant_dv_wa_scale
from barak.convolve import convolve_constant_dv
from COS import convolve_with_COS_FOS
import barak.spec
from barak.utilities import between, adict
from barak.plot import puttext

import pylab as pl
import numpy as np
import matplotlib.transforms as mtransforms

atomdat = readatom(molecules=1)

def get_fig_axes(nrows, ncols, npar, width=11.7, height=None, aspect=0.5):
    """ Generate a figure with a number of subplots.

    Parameters
    ----------
    nrows : int
      Number of rows
    ncols : int
      Number of columns
    npar : int
      Number of plots in total
    width : float
      Width of the figure in inches
    aspect : float
      width / height of each sub-plot.

    Returns
    -------
    fig : matplotlib figure object
    axes : list of matplotlib axes objects
      ordered from top left going down coluns.


    Notes
    """
    if height is None:
        height = width*aspect*nrows/ncols
    fig = pl.figure(figsize=(width, height))    

    axes = [fig.add_subplot(nrows, ncols, i+1) for i in range(npar)]

    # reorder axes so they go from top left down columns instead of across
    axes1 = []
    ileft = []
    ibottom = []
    for j in range(ncols):
        for i in range(nrows):
            ind = j + i*ncols
            if ind > npar - 1:
                continue
            axes1.append(axes[ind])
    # find the indices of the left and bottom plots (used to set axes
    # labels)
    ileft = range(nrows)
    ibottom = [i*nrows - 1 for i in range(1, ncols+1)]
    for i in range(ncols*nrows - npar):
        ibottom[-(i+1)] -= ncols*nrows - npar - i

    velplot = adict(axes=axes1, nrows=nrows, ncols=ncols,
                    ileft=ileft, ibottom=ibottom)
    return fig, velplot

def get_nrows_ncols(npar):
    """ Get the number of rows and columns
    """ 
    nrows = max(int(np.sqrt(npar)), 1)
    ncols = nrows
    while npar > (nrows * ncols):
        nrows += 1

    return nrows, ncols

def plot_tick_vel(ax, vpos, offset, t,  col='b', alpha=1, ls='-',
                  ticklabels=False, lw=1, height=0.2):
    """ Plot a single velocity tick."""

    T = ax.vlines([vpos], 1.1 + offset, 1.1 + height + offset,
                  color=col, lw=lw,
                  alpha=alpha, linestyle=ls)

    if ticklabels:
        label = '%s %.0f %.2f' % (t.name, t.wa0, t.z)
        label = label.replace('NeVII', 'NeVIII')

        T.append(ax.text(vpos, 1.05 + offset, label, rotation=60,
                    fontsize=8, va='bottom',alpha=0.7, color=col))
    return T


def plot_velocity_regions(wmin, wmax, w0, w1, obswa, ax, offset):
    """ wmin, wmax is minimum and maximum wavelengths of the plot.

    w0 and w1 are the min and max wavelengths of the fitting
    regions. obswa is the wavelength of the transition for this plot.
    """
    cond = ((w1 >= wmax) & (w0 < wmax)) | \
           ((w1 <= wmax) & (w0 >= wmin)) | \
           ((w1 > wmin) & (w0 <= wmin))
    regions = []

    if not cond.any():
        return regions

    vel0 = (w0[cond] / obswa - 1) * c_kms
    vel1 = (w1[cond] / obswa - 1) * c_kms
    for v0,v1 in zip(vel0, vel1):
        yoff = 1.1 + offset
        R, = ax.plot([v0, v1], [yoff,yoff],'r',lw=3, alpha=0.7)
        regions.append(R)

    return regions

def process_Rfwhm(Rfwhm, wa, model, models):
    """ Convolve the input models using the Rfwhm option

    Return the new models.

    wa:  wavelength array, shape (N,)
    model: model flux array, shape (N,)
    models: list of model flux arrays each with shape (N,)

    Rfwm is one of:

      'convolve_with_COS_FOS'
      a float

    Returns
    -------
    new_model, new_models
    """

    model_out = None
    models_out = []

    if Rfwhm is None:
        return model, models

    elif Rfwhm == 'convolve_with_COS_FOS':
        print 'convolving with COS/FOS instrument profile'
        model_out = convolve_with_COS_FOS(model, wa, wa[1] - wa[0])
        models_out = [convolve_with_COS_FOS(m, wa, wa[1] - wa[0]) for m
                      in models]

    elif isinstance(Rfwhm, float):
        print 'Convolving with fwhm %.2f km/s' % Rfwhm
        # use a pixel velocity width 4 times smaller than the FWHM
        ndiv = 4.
        wa_dv = make_constant_dv_wa_scale(wa[0], wa[-1], Rfwhm / ndiv)
        model_out = convolve_constant_dv(wa, model, wa_dv, ndiv)
        # do the same to every model if there's more than one
        for m in models:
            models_out.append(convolve_constant_dv(wa, m, wa_dv, ndiv))
    else:
        raise ValueError('Unknown value for Rfwhm option')

    return model_out, models_out

def read_transitions(fh, atomdat):
    """ Read transitions from a file

    Parameters
    ----------
    fh : str or file object
    atomdat : atom.dat object
      Read with `barak.absorb.readatom()`
    
    Returns
    -------
    linelist : list of atom.dat entries
      Atom.dat entry for each transition

    Notes
    -----
    Example file format:

      HI 1215
      CIV 1549
      CIV 1550
      MgII 2796
      MgII 2803

    """
    if isinstance(fh, basestring):
        print 'Reading transitions from', fh
        fh = open(fh, 'rt')

    linelist = []
    for tr in fh:
        tr = tr.strip()
        if tr and not tr.startswith('#'):
            name, t = findtrans(tr, atomdat=atomdat)
            linelist.append((name, t))

    fh.close()
    return linelist

def get_colours(transitions):
    """ Generate a single color for each ion in a list of
    transitions.

    Returns a dictionary mapping the ion name to a colour.
    """
    colours = 'b', 'r', 'g', 'orangered', 'c', 'purple'
    #colours = ['0.5']

    ions = [tr[0].split()[0] for tr in transitions]
    # want an ordered set
    ionset = []
    for ion in ions:
        if ion not in ionset:
            ionset.append(ion)

    colour = dict(zip(ionset, colours * (len(ions) // len(colours) + 1)))
    return colour

def make_models(opt, path=None):
    if opt.f26name is not None:
        n = (path + opt.f26name if path is not None else opt.f26name)
        f26 = readf26(n)

    n = (path + opt.trfilename if path is not None else opt.trfilename)
    transitions = read_transitions(n, atomdat)
    ntrans = len(transitions)
    if isinstance(opt.spfilename, basestring):
        n = (path + opt.spfilename if path is not None else opt.spfilename)
        sp = barak.spec.read(n)
    else:
        sp = opt.spfilename
    wa = sp.wa
    edges = barak.spec.find_bin_edges(wa)     # pixel edges
    dwa = edges[1:] - edges[:-1]             # width per pixel
    ticks = []
    model = None
    models = []

    if opt.f26name is not None:
        dw = np.median(dwa)
        lines = []
        for l in f26.lines:
            lines.append((l['name'].replace(' ',''),l['z'],l['b'],l['logN']))

        if opt.wadiv is not None:
            dw1 = dw / opt.wadiv
            wa1 = np.arange(wa[0], wa[-1]+0.5*dw1, dw1)
            tau, ticks = find_tau(wa1, lines, atomdat)
        else:
            tau, ticks = find_tau(wa, lines, atomdat)

        model = np.exp(-tau)
        models = []
         
        #models = [np.exp(-t) for t in taus]
         
        if opt.wadiv is not None:
            model, models = process_Rfwhm(opt.Rfwhm, wa1, model, models)
        else:
            model, models = process_Rfwhm(opt.Rfwhm, wa, model, models)
         
        if opt.wadiv is not None:
            model = np.interp(wa, wa1, model)
            models = [np.interp(wa, wa1, m) for m in models]

    return sp, transitions, model, models, ticks

def make_plot(sp, transitions, model, models, ticks, opt, width=None,
              height=None, aspect=0.5, textloc='lower left',
              unrelated=[], cols=None, flw=0.5, mlw=0.5, nrows=None, ncols=None,
              unrelated_alpha=0.5, trans_fontsize='medium', zerovel_lw=0.5,
              tlw=1):
    """
    Parameters
    ----------
    cols : dict
      An optional dictionary of colours. The keys 'data', 'resid',
      'ticks' and 'model' can be defined.

    """
    if cols is None:
        cols = {}
    if 'data' not in cols:
        cols['data'] = None
    if 'resid' not in cols:
        cols['resid'] = 'g'
    if 'model' not in cols:
        cols['model'] = 'r'
    if 'ticks' not in cols:
        cols['ticks'] = 'b'

    wa = sp.wa
    wextra = np.diff(wa)[0] * 2
    nfl = sp.fl / sp.co
    ner = sp.er / sp.co
    if opt.f26name is not None:
        resid = (nfl - model) / ner

    # actual plotting
    ntrans = len(transitions)
    if None in (nrows, ncols):
        nrows, ncols = get_nrows_ncols(ntrans)
    assert ntrans <= nrows*ncols
    fig, vplot = get_fig_axes(nrows, ncols, ntrans, width=width, height=height,
                              aspect=aspect)

    colour = get_colours(transitions)

    zp1 = opt.redshift + 1
    vmax = opt.vmax
    if opt.vmin is None:
        vmin = -abs(vmax)
    else:
        vmin = opt.vmin

    betamin = vmin / c_kms
    betamax = vmax / c_kms

    fig.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.99,
                        wspace=0.0001, hspace=0.00001)

    # plot top down, so we need reversed()
    for i,trans in enumerate(transitions):
        ax = vplot.axes[i]
        ion = trans[0].split()[0]
        watrans = trans[1][0]
        obswa = watrans * zp1
        wmin = obswa * (1 - 3*abs(betamin))
        wmax = obswa * (1 + 3*abs(betamax))
        cond = between(wa, wmin, wmax)
        wa1 = wa[cond]

        fl = nfl[cond]
        vel = (wa1 / obswa - 1) * c_kms
        ax.axhline(0, color='gray', lw=0.5)
        c = cols['data']
        if c is None:
            c = colour[ion]
        ax.plot(vel, fl, color=c, lw=flw, ls='steps-mid')
        for w0,w1 in unrelated:
            c0 = between(wa1, w0, w1)
            c1 = between(wa1, w0-wextra, w1+wextra)
            if c0.any():
                ax.plot(vel[c0], fl[c0], 'w', lw=flw+2, ls='-',
                        drawstyle='steps-mid')
                #ax.plot(vel[c0], fl[c0], '--', lw=flw, color=c,
                #        drawstyle='steps-mid')
                ax.plot(vel[c1], fl[c1], '-', lw=flw, color=c,
                        alpha=unrelated_alpha, drawstyle='steps-mid')

        ax.axhline(1, color='gray', lw=0.5, ls='dashed', zorder=3)

        for m in models:
            ax.plot(vel, m[cond], cols['model'], lw=0.2)
        if opt.f26name is not None:
            if opt.modelfill:
                ax.fill_between(vel, model[cond], 1, lw=0,
                                color=cols['model'], alpha=0.4, zorder=0)
            else:
                ax.plot(vel, model[cond], cols['model'], lw=mlw)

            for w0,w1 in unrelated:
                c0 = between(wa1, w0, w1)
                if c0.any():
                    ax.plot(vel[c0], model[cond][c0], 'w', lw=flw+2, ls=':')

        if opt.residuals and opt.f26name is not None:
            mult = 0.05
            offset = -0.2
            ax.scatter(vel, offset + mult*resid[cond], marker='.', s=5,
                       faceted=False, c=cols['resid'])
            ax.axhline(offset - mult, color='0.5', lw=0.3)
            ax.axhline(offset + mult, color='0.5', lw=0.3)

        if len(ticks) > 0:
            tickwmin = obswa * (1 + betamin)
            tickwmax = obswa * (1 + betamax)
            wticks = ticks.wa
            cond = between(wticks, tickwmin, tickwmax)
            if cond.any():
                # should really check ion name here too...
                c0 = np.abs(trans[1][0] - ticks.wa0[cond]) < 1e-2
                vel = (wticks[cond][c0] / obswa - 1) * c_kms
                #import pdb; pdb.set_trace()
                for j,t in enumerate(ticks[cond][c0]):
                    T = plot_tick_vel(ax, vel[j], 0, t, lw=tlw, col=cols['ticks'])
                vel = (wticks[cond][~c0] / obswa - 1) * c_kms
                for j,t in enumerate(ticks[cond][~c0]):
                    T = plot_tick_vel(ax, vel[j], 0, t, lw=tlw, col=cols['ticks'],
                                      ls=':')

        #bbox = dict(facecolor='w', edgecolor='None')
        transf = mtransforms.blended_transform_factory(
            ax.transAxes, ax.transData)
        name = trans[0]
        if opt.osc:
            name = name + ' %.3g' % trans[1]['osc']

        if textloc == 'lower left':
            ax.text(0.03, 0.02, name, fontsize=trans_fontsize, transform=transf) # bbox=bbox,
        elif textloc =='upper right':
            puttext(0.97, 0.97, name, ax, fontsize=trans_fontsize, va='top',
                    ha='right') # bbox=bbox,
        #ax.text(0.98, 0.02, name, ha='right', fontsize=13, transform=transf) # bbox=bbox,
    
    for i,ax in enumerate(vplot.axes):
        ax.axvline(0, color='0.7', lw=zerovel_lw, zorder=0)
        ax.set_xlim(vmin, vmax)
        ax.set_ylim(-0.49, 1.49)

    return fig, vplot
