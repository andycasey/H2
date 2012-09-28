from Draine96_H2shield import calc_H2_shield, calc_H2_shield2
"""
following Jorgenson Appendix

eqn A11
"""

def calc_Td(JnuLW, QUV, a=0.1, A=3.20e-3):
    """ Find dust temperature using eqn A19 from Jorgenson et al.
    
    JnuLW : float
      Radiation field at 1000 Ang in erg/s/cm^2/Hz/sr.
      
    A : float
      Constant (cm) that depends on th optical properties of the
      grains. Default value for carbonaceous grains, for silicate
      grains it should be 1.34e-3.

    a : float
      Grain diameter in micrometers.
    """
    # 3.2e-20 is JnuLW in the solar neighbourhood.
    chi = JnuLW / 3.2e-20
    
    Td = 12 * (chi * QUV)**(1./6) * (A/3.2e-3)**-(1./6) * (a/0.1)**-(1./6)
    return Td

def calc_Sd(T, Td):
    """ Dust sticking coefficient.
    
    T is gas temperature

    Td is dust temperature
    """
    Sd = (1 + 0.04*(T + Td)**0.5 + 2e-3*T + 8e-6*T**2)**-1 * \
         (1 + np.exp( 7.5e2* (1./75 - 1./Td) ) )**-1
    return Sd
    
def calc_R2(Sd, kappa, a, T, delta):
    """ Eqn A17 from Jorgenson et al.

    Sd : float

     Dust sticking coefficient

    delta : float

      Grain material density in g cm^-3

    kappa : float

      Ratio of Dust-to-gas mass ratio to solar neighbourhood
      value. (Can use metallicity as a proxy.)

    """
    R1 = 4.1e-17 * Sd * (a/0.1)**-1 * (delta/2.)**-1
    return R1



def approx_fH2(NH2, n, JnuLW, kappa, a, T, delta, b):
    """ Find approximate fH2 using taylor expansion of exppnential
    dust shielding term.
    """
    Td = calc_Td(JnuLW, 1, a, 3.20e-3)
    Sd = calc_Sd(T, Td)
    R1 = (T/100.)**0.5 * calc_R2(Sd, kappa, a, T, delta)

    fH2 = kappa * (0.0136 * (JnuLW / 3.2e-20)**-1 * (R1/3e-17) * \
                    (n/100.) * \
                    (calc_H2_shield(NH2)/0.01)**-1 + \
                    #(calc_H2_shield2(NH2,b)/0.01)**-1 + \
                    1.758 * NH2 / 10**21)
        
    return fH2

# gas temperature (Jnu and n dependence on this is weak)
T = 100.
# dust grain radius in um (solar is 0.1). n dependence is a^-1
a = 0.1
# dust grain density (solar is 2). n dependence is d^-1
delta = 2.

# Solar neighbourhood, Habing et al. 1968. ergs/s/cm^2/Hz/sr
JnuLW = 3.2e-20   # ^-1
# dust to gas ratio / solar neighbourhood dust to gas ratio
kappa = 1         # ^1
b = 5.
n = 100.           # ^1

# cm^-2
NH2 = 10**np.arange(11, 22, 0.02)

models = []

for nkJ in [0.01, 0.1, 1, 10, 100, 1000]:
    n1 = n * nkJ
    fH2 = np.array(approx_fH2(NH2, n1, JnuLW, kappa, a, T, delta, b), float,)
    fH2[fH2 >1] = 1
    Ntot = 2*NH2 / fH2
    models.append([nkJ, NH2, Ntot, fH2])
    pl.figure(1)
    #pl.clf()
    plot(np.log10(Ntot), np.log10(fH2),
         label=r'$n\kappa/J^{LW}_{\nu} %g$' % nkappa_on_J)

    pl.figure(2)
    #pl.clf()
    plot(np.log10(NH2), np.log10(fH2))
    pl.xlabel('NH2')

saveobj('models.sav', models)

pl.figure(1)
pl.xlabel('Ntot')
pl.legend(loc=0)
show()

# note R1 = R * kappa

# calc n

#n = 2*pi * 1.1e8 * JnuLW * (T/100.)**-0.5 *  R1**-1 * \
#    calc_H2_shield(NH2, b) * (kappa**-1 * fH2 - 0.879 * 2*NH2/10**21)

