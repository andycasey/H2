from __future__ import division
import numpy as np
""" Determine f_H2 from metallicity from Krumholz' model (as described
by Kuhlen, Krumholz et al. 2012)

Assumes:

  - formation-dissociation balance

  - based on a radiative transfer calculation of an idealized
    spherical giant atomic-molecular complex, subject to a uniform and
    isotropic Lyman-Werner (LW) radiation field.

  - two-phase equilibrium between a cold neutral medium (CNM) and a
    warm neutral medium (WNM) (Wolfire et al. 2003) as in ISM.  The
    assumption of pressure balance between these two ISM components
    forces the minimum CNM density to be linearly proportional to the
    intensity of the LW radiation field, with only a weak dependence
    on metallicity.

""" 
# Need: column density, metallicity

def calc_fH2(Z, Ntot, sigma_d_m21=1, R_m16p5=1, phiCNM=3, Z_SN=0.0204):
    """ Find the molecular mass fraction 

    find the mass fraction f_H2 = m(H2) / m(Htotal) for a metallciity
    using the relation from Krumholz. These equations are taken from
    Kuhlen, Krumholz et al. (2012).

    sigma_d_m21 : float (1)
      Dust cross section per H nucleus to 1000A radiation normalised
      to 10^-21 cm^-2 as in Krumholz and Gnedin (2011).

    R_m16p5 : float (1)
      Rate coefficient for H2 formation on dust grains normalised to
      Milky Way value of 10^-16.5 cm^3 s^-1 (Wolfire et al. 2008) as
      in Krumholz and Gnedin (2011).

    phiCNM : float (3)
      Multiplier to the minimum density nmin as in Krumholz and Gnedin
      (2011).

    Z_SN : float (0.0204)
      Metallicity in solar neighbourhood (Rodriguez and Delgao-Inglada
      2011).
    """

    # dust cross section per H nucleus (cm**2), from Kuhlen et al 2012.
    sigma_d = sigma_d_m21 * 1e-21 * 10**(Z - Z_SN)
     
    # dust optical depth of the cloud
    tau_c = sigma_d * Ntot
     
    chi = 2.3 * (sigma_d_m21 / R_m16p5) * \
          (1 + 3.1 * (10**(Z - Z_SN))**0.365) / phiCNM

    s = np.log(1 + 0.6*chi + 0.01*chi**2) / (0.6 * tau_c)
    # note that for f to be positive, s must be <= 2
    f_H2 = 1 - 3/4. * s / (1 + 0.25*s)

    return f_H2


if 0:
    # pressure balance assumption lets us eliminate dependence on nH0 and
    # Gdash0, but we can check the values we expect to find.
     
    # ambient UV radiation field intensity, normalized to the Draine
    # (1978) value for the Milky Way.
    #Gdash0 = 1
     
    # volume density of H nuclei (cm^-3)
    #nmin = 31 / (1 + 3.1 * (10**(Z - Z_SN))**0.365) * Gdash0
     
    #n = phiCNM * nmin
     
    Ntot = 10**np.arange(20, 23, 0.001)
     
    #pl.figure(figsize=(10,5))
    for zfrac in (0.12, 0.2, 0.5, 1):
        Z = np.log10(zfrac)
        f_H2 = calc_fH2(Z, Ntot, sigma_dm21=1, R_m16p5=1, phiCNM=3, Z_SN=0.0204)
        c0 = f_H2 > 0
        NH2 = 0.5 * Ntot[c0] * f_H2[c0]
        NHI = Ntot[c0] - 2*NH2
        plot([np.log10(Ntot[c0])[0]] + list(np.log10(Ntot[c0])),
             [12] + list(np.log10(2*NH2)), 'r')
     
    pl.xlabel('log10 NHtot')
    pl.ylabel('log10 NH2')
    pl.xlim(18.8, 22.2)
    pl.ylim(12, 22)
    pl.minorticks_on()
    pl.show()
