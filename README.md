# rcSLW

This code sets up a class for getting local gray gas coefficients κ's, and weights a's for the Rank Correlated SLW model.
The following initializes the code:

    import rcslw
    nGG = 3      # number of GG not including the clear j=0 gstrong6494
    P = 1.5        # atm
    rad = rcslw.rcslw(P, nGG)

Then, at any give gas state, (here written for co2), the κ's and a's are computed as:

    Tg = 1200   # K
    Yco2=0.5    # mol frac
    Nconc = 8   # mol/m3
    sp = 'co2'
    k, a = rad.get_k_a(sp, Tg, Yco2, Nconc)

This would then be used in conjunction with an ODE solver for the dI_j/ds equations, which we need to write separately.
In such a solver, the composition and temperature may vary spatially, so we would set Tg, Yco2, Nconc as above and redo the call to get_k_a
locally as needed.

The code can be obtained from github:

git clone https://github.com/BYUignite/rcSLW.git
