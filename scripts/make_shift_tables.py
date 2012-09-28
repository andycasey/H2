from barak.interp import interp_spline, InterpCubicSpline

if 1:
    wa = np.arange(1425, 1810, 25)

    # spline fitted to the H2 shifts (see
    # /home/nhmc/data/COS/Crighton11585/scripts/combine.py and
    # /home/nhmc/data/COS/Crighton11585/scripts/H2_shift.py)
    x = [1467.39713564, 1541.25417348, 1637.0329482, 1700.77445989]
    y = [-0.07306245, -0.02829985, -0.00470822,  0.00099504]
    spl = InterpCubicSpline(x,y)

    plot(wa, spl(wa), '.-', label='H2 shift')

    # these are polynomial coefficients. FUV A is the higher wavelength range
    G160M_FUVA_1627_to_1589 = np.array([ 1.23279433e-03, -2.11545155])
    G160M_FUVB_1589_to_1627 = np.array([-9.51673753e-04,  1.46794556])

    wa1 = np.linspace(1425, 1600)
    wa2 = np.linspace(1600, 1800)
    plot(wa1, np.polyval(G160M_FUVB_1589_to_1627, wa1), label='FUVB shift')
    plot(wa2, np.polyval(G160M_FUVA_1627_to_1589, wa2), label='FUVA shift')
    pl.legend()
    show()

    for i in range(len(wa)):
        print '%.0f &  %.4f \\\\' % (wa[i], spl(wa[i]))
    
