"""
This model must define the following objects:

- a dictionary P with keys. The value of every key is a tuple with the
same length (the number of model parameters)

    name  : parameter names
    min   : minimum allowed parameter values
    max   : maximum allowed parameter values
    guess : parameter values to use to generate initial walker positions

- array of values x

- array of data values ydata

- array of data one sigma errors ysigma

- a ymodel(x, par) function that generates the model of the data given
  an array of parameter values

- a ln_likelihood(x, ydata, ysigma) function

"""
from __future__ import division
from barak.absorb import calc_iontau, readatom, findtrans
from barak.pyvpfit import  readf26
from barak.utilities import adict
from COS import convolve_with_COS_FOS
from barak.constants import c_kms
import barak.spec
from barak.plot import puttext
import numpy as np
from barak.io import writetxt

import pylab as pl

pl.rc('xtick', labelsize=11)
pl.rc('xtick.major', size=3, pad=5)
pl.rc('xtick.minor', size=1.5)
pl.rc('ytick', labelsize=11)
pl.rc('ytick.major', size=3, pad=4)
pl.rc('ytick.minor', size=1.5)

# divide wavelength pixels into this many sub-pixels
nwadiv = 10
# cache for spectra
SPECTRA = {}
atom = readatom(molecules=True)
options = {}
options['wa_vary'] = True

###################################################################
# Read initial model parameters along with fitting regions from vpfit
# f26-style file
###################################################################


def readf26file(filename):
    """ Read a f26 style file, and return with regions sorted by
    minimum wavelength.
    """
    vpin = readf26(filename)
    isort = np.argsort([r.wmin for r in vpin.regions])
    vpin.regions = vpin.regions[isort]
    return vpin

def find_tau(wa, lines):
    """ generate a tau array given a list of transitions and a
    wavelength array.

    The transition list is in a vp.lines format
    """ 
    tau = np.zeros_like(wa)
    for l in lines:
        trans = atom[l.name.strip()]
        tau += calc_iontau(wa, trans, l.z + 1, l.logN, l.b, maxdv=None) 

    return tau


def calc_region_info(vp, lines_fixed, npad=10):
    """read in the spectra for each region.

    return a list with one entry for each region giving the spectrum,
    and start and end indices of the fitting region in that spectrum.
    """
    regions = []
    
    for reg in vp.regions:
        try:
            sp = SPECTRA[reg['filename']]
        except KeyError:
            sp = barak.spec.read(reg['filename'])
            SPECTRA[reg['filename']] = sp

        #print reg['filename'], reg['wmin'], reg['wmax']
        # cut out the fitting region
        i0,i1 = sp.wa.searchsorted([reg['wmin'], reg['wmax']])
        i1 += 1
        # generate tau for the non-varying transitions
        j0 = max(0, i0 - npad)
        j1 = min(len(sp.wa) - 1, i1 + npad)
        # split into number of divisions
        wa = np.linspace(sp.wa[j0], sp.wa[j1], (j1-j0) * nwadiv)
        # get non-varying lines
        tau = find_tau(wa, lines_fixed)
        temp = reg.resolution.split()
        trans = ' '.join(temp[1:3])
        regions.append( (wa, tau, sp, (i0,i1), trans) )
    # for reg in vp.regions:
    #     try:
    #         sp = SPECTRA[reg['filename']]
    #     except KeyError:
    #         sp = barak.spec.read(reg['filename'])
    #         SPECTRA[reg['filename']] = sp

    #     #print reg['filename'], reg['wmin'], reg['wmax']
    #     # cut out the fitting region
    #     i0,i1 = sp.wa.searchsorted([reg['wmin'], reg['wmax']])
    #     i1 += 1
    #     # generate tau for the non-varying transitions
    #     j0 = max(0, i0 - npad)
    #     j1 = min(len(sp.wa) - 1, i1 + npad)
    #     # split into number of divisions
    #     wa = np.linspace(sp.wa[j0], sp.wa[j1], (j1-j0) * nwadiv)
    #     # get non-varying lines
    #     tau = find_tau(wa, lines_fixed)
    #     regions.append( (wa, tau, sp, (i0,i1)) )

    return regions

def get_par_from_lines(lines, blim=(2, 15), Nlim=(14, 19.), dz=5e-5):
    """ Get all the parameters we want to vary from a list of lines.
    """
    P = adict(names=[], min=[], max=[], tied={})

    i = 1
    for l in lines:
        if l.zpar == '':
            P.names.append('z%i' % i)
            P.min.append(l.z - dz)
            P.max.append(l.z + dz)
        elif l.zpar.islower():
            if l.zpar not in P.tied:
                # entry gives the index of the parameter corresponding
                # to this tied character. Note all tied lines must
                # have the same lower case character (this is
                # different to VPFIT, where one character is lower
                # case, and all the others are the same letter but
                # upper case)
                
                P.tied[l.zpar] = len(P.names)
                P.names.append('z%i' % i)
                P.min.append(l.z - dz)
                P.max.append(l.z + dz)

        if l.logNpar == '':
            P.names.append('N%i' % i)
            P.min.append(Nlim[0])
            P.max.append(Nlim[1])
        if l.bpar == '':
            P.names.append('b%i' % i)
            P.min.append(blim[0])
            P.max.append(blim[1])
        elif l.bpar.islower():
            if l.bpar not in P.tied:

                # entry gives the index of the parameter corresponding
                # to this tied character. Note all tied lines must
                # have the same lower case character (this is
                # different to VPFIT, where one character is lower
                # case, and all the others are the same letter but
                # upper case)
                
                P.tied[l.bpar] = len(P.names)
                P.names.append('b%i' % i)
                P.min.append(blim[0])
                P.max.append(blim[1])


        if '' in (l.zpar, l.logNpar, l.bpar) or l.bpar.islower() or \
               l.zpar.islower():
            i += 1

    return P

def copy_par_to_lines(par, lines):
    """ given a list of parameters, populate a list of lines with the
    correct value. lines can then be passed to ymodel.
    """
    lines = lines.copy()
    i = 0
    for l in lines:
        if l.zpar == '':
            l.z = par[i]
            i += 1
        elif l.zpar.islower():
            # multiple b parameters share a single parameter
            l.z = par[P.tied[l.zpar]]
            
            # if we are at the tied paramter index, move to next
            # parameter, otherwise we don't want to skip any
            # parameters.

            if i == P.tied[l.zpar]:
                i += 1

        if l.logNpar == '':
            l.logN = par[i]
            i += 1
        if l.bpar == '':
            l.b = par[i]
            i += 1
        elif l.bpar.islower():
            # multiple b parameters share a single parameter
            l.b = par[P.tied[l.bpar]]
            
            # if we are at the tied paramter index, move to next
            # parameter, otherwise we don't want to skip any
            # parameters.

            if i == P.tied[l.bpar]:
                i += 1

    assert i == len(par)

    return lines

 
def make_data_arrays(regions):
    """ Make 1-d wavelength, normalised flux and error arrays from all
    the regions to pass to the ln_likelihood function
    """ 
    x = []
    ydata = []
    ysigma = []
    for wa,tau,sp,(i,j),tname in regions:
        x.extend(sp.wa[i:j])
        ydata.extend(sp.fl[i:j] / sp.co[i:j])
        ysigma.extend(sp.er[i:j] / sp.co[i:j])

    return map(np.array, (x, ydata, ysigma))
 

def ymodel(wa, pars):
    """Generate the model normalised flux array for every region.

    Note this needs the global variables `lines_vary` and `regions` to
    be defined.
    """
    linepars = pars
    dz = 0
    if options['wa_vary']:
        linepars, regpars = pars[:-len(regions)], pars[-len(regions):]
    lines = copy_par_to_lines(linepars, lines_vary)
    fl = []
    for ireg, (wa, tau0, sp, (i, j),tname) in enumerate(regions):
        tau = tau0.copy()
        for l in lines:
            if options['wa_vary']:
                dz = regpars[ireg] / wa[len(wa)//2] * (l.z + 1)

            trans = atom[l.name.strip()]
            tau += calc_iontau(wa, trans, l.z+1 + dz, l.logN, l.b)

        dwpix = wa[1] - wa[0]
        flmodel0 = convolve_with_COS_FOS(np.exp(-tau), wa, dwpix)
        flmodel = np.interp(sp.wa[i:j], wa, flmodel0)
        fl.extend(flmodel)
        
    return fl

     
def ln_likelihood(pars, x, y, ysigma):
    # if we are outside the allowable parameter ranges, return 0
    # likelihood.
    for i,p in enumerate(pars):
        if not (P.min[i] < p < P.max[i]):
            return -np.inf

    ymod = ymodel(x, pars)
    #pl.plot(x, ymod)
    #pl.draw()
    #pl.waitforbuttonpress()
    resid = (y - ymod) / ysigma
    return -0.5 * np.dot(resid, resid)


def get_initial_positions(nwalkers):
    # Get initial parameter positions (guesses!) for each walker
    
    p0 = np.random.uniform(size=(nwalkers, len(P.names)))
    for i in range(len(P.names)):
        p0[:, i] = P.min[i] + p0[:, i] * (P.max[i] - P.min[i])

    return p0

def plot_model(pars):
    """
    `regions`, `x` must be defined in the model.py module namespace
    """
    from run_emcee.plot_mcmc import get_fig_axes, get_nrows_ncols

    linepars = pars
    dz = 0
    if options['wa_vary']:
        linepars, regpars = pars[:-len(regions)], pars[-len(regions):]
    lines = copy_par_to_lines(linepars, lines_vary)
    
    nrows, ncols = get_nrows_ncols(len(regions))
    fig, axes = get_fig_axes(nrows, ncols, len(regions), width=8.2)
    fig.subplots_adjust(wspace=1e-6, hspace=1e-6, right=1, top=1,
                        left=0.07, bottom=0.06)
    z0 = lines.z[1]
    zp1 = lines.z.mean()
    colours = dict(J0='purple',J1='g',J2='r',J3='b')

    for ind, (wa, tau0, sp, (i, j), tname) in enumerate(regions):

        tau = tau0.copy()
        for l in lines:
            if options['wa_vary']:
                #import pdb; pdb.set_trace()
                dz = regpars[ind] / wa[len(wa)//2] * (l.z + 1)
                #print l.z, dz, l.logN
            
            trans = atom[l.name.strip()]
            tau += calc_iontau(wa, trans, l.z+1 + dz, l.logN, l.b)

        dwpix = wa[1] - wa[0]
        flmodel0 = convolve_with_COS_FOS(np.exp(-tau), wa, dwpix)
        ymod = np.interp(sp.wa[i:j], wa, flmodel0)


        ax = axes[ind]
        i1 = max(i-50, 0)
        j1 = min(j+50, len(sp.wa)-1)
        nfl = sp.fl[i1:j1] / sp.co[i1:j1]
        ner = sp.er[i1:j1] / sp.co[i1:j1]

        tstr, t = findtrans(tname, atom)
        obswa = t.wa * (1 + z0)
        dv = c_kms * (sp.wa[i:j] / obswa - 1)
        dv1 = c_kms * (sp.wa[i1:j1] / obswa - 1)

        tstr = tstr[2:]
        
        ax.axhline(0, color='0.5',lw=0.5)
        ax.axhline(1, color='k', lw=0.5, ls='--')
        ax.fill_between([-150, 150], -0.35, -0.15, color='0.9')
        ax.axhline(-0.25, color='k', lw=0.25)
        ax.plot(dv1, nfl, color='0.7', lw=0.5, drawstyle='steps-mid')
        ax.plot(dv, sp.fl[i:j]/sp.co[i:j], 'k', lw=0.5, drawstyle='steps-mid')
        ax.plot(dv, ymod, color='0.3')
        #ax.axvline(0, color='0.5', lw=0.5)

        ax.plot([-23]*2, [1.1, 1.4], '0.5')
        ax.plot([0]*2, [1.1, 1.4], '0.5')
        
        puttext(0.95, 0.95, tstr, ax, fontsize=9, va='top', ha='right',
                color=colours[tstr[:2]])

        resid = (sp.fl[i:j]/sp.co[i:j] - ymod) / sp.er[i:j] * sp.co[i:j]

        ax.plot(dv, -0.25 + resid*0.1, '.', color='0.3', ms=2)
        ax.set_ylim(-0.5, 1.9)
        ax.set_xlim(-120, 120)
        ax.set_yticks([0, 0.5, 1, 1.5])
        ax.set_yticklabels(['0.0', '0.5', '1.0', '1.5'])
        ax.set_xticks([-100, -50, 0, 50, 100])
        ax.set_xticklabels(['', '-50', '0', '50', ''])
        if ind+1 < (ncols*(nrows-1)):
            ax.set_xticklabels([])
        if ind % ncols:
            ax.set_yticklabels([])

    fig.text(0.5, 0.00, 'Velocity offset (km/s)', fontsize=12, ha='center',va='bottom')
    fig.text(0.015, 0.5, 'Transmission', fontsize=12, rotation=90,va='center')
    ndf = len(ydata) - len(pars)
    print (((ydata - ymodel(x, pars))/ ysigma) **2).sum() / ndf
    pl.savefig('fig/model.pdf')
    pl.savefig('fig/model.png',dpi=300)

    return fig

def print_par(par):
    """ Print the maximum likelihood parameters and their
    uncertainties.
    """
    rec = []
    for i in range(len(P.names)):
        p = P.ml[i]
        m1 = P.p1sig[i]
        m2 = P.p2sig[i]
        j1 = P.p1sig_joint[i]
        j2 = P.p2sig_joint[i]
        rec.append( (P.names[i], p,  p - j1[0], j1[1] - p,
                     p - j2[0], j2[1] - p, p - m1[0], m1[1] - p,
                     p - m2[0], m2[1] - p) )

    names = 'name,ml,j1l,j1u,j2l,j2u,m1l,m1u,m2l,m2u'
    rec = np.rec.fromrecords(rec, names=names)

    hd = """\
# name : parameter name
# ml   : maximum likelihood value
# j1l  : 1 sigma lower error (joint with all other parameters) 
# j1u  : 1 sigma upper error (joint)
# j2l  : 2 sigma lower error (joint) 
# j2u  : 2 sigma upper error (joint) 
# m1l  : 1 sigma lower error (marginalised over all other parameters)
# m1u  : 1 sigma upper error (marginalised)
# m2l  : 2 sigma lower error (marginalised) 
# m2u  : 2 sigma upper error (marginalised) 
"""
    #pdb.set_trace()
    writetxt('parameters.txt', rec, header=hd, fmt_float='.8g', overwrite=1)


if 1:
    print 'Reading H2J0123.f26'
    vp = readf26file('H2J0123.f26')
    vary = np.array(['' in (l.logNpar, l.bpar, l.zpar) or l.bpar.islower() \
                     or l.zpar.islower() for l in vp.lines])
    lines_vary = vp.lines[vary]
    lines_fixed = vp.lines[~vary]
    regions = calc_region_info(vp, lines_fixed)

    # parameter info.
    P = get_par_from_lines(lines_vary, blim=(0.5,15), Nlim=(14.0, 19.0))
    if options['wa_vary']:
        print 'Adding possible wavelength shift to each fitting region'
        P.names.extend(['dw%i' % i for i in range(len(regions))])
        P.min.extend([-0.02] * len(regions))
        P.max.extend([0.02] * len(regions))

    print len(P.names), 'parameters'
    assert len(P.names) == len(P.min) == len(P.max)
    x, ydata, ysigma = make_data_arrays(regions)

if 0:
    plot(x, ydata, '.', x, ysigma, '.')
    y = ymodel(x, P.guess)
    plot(x, y)
    show()
