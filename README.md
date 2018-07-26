# Rank Correlated SLW model

This code sets up a class for getting local gray gas coefficients κ's, and weights a's for the Rank Correlated SLW model.
* See the paper by [Solovjov et al., J. Quantitative Spectroscopy and Radiative Transfer 197:26-44(2017)](https://www.sciencedirect.com/science/article/pii/S0022407316306434?via%3Dihub).
* The following initializes the code:

    import rcslw
    nGG = 3                   # number of GG not including the clear j=0 gstrong6494
    P = 1.5                   # atm
    rad = rcslw.rcslw(P, nGG)

* Then, at any give gas state, the κ's and a's are computed as:

    Tg = 1200    # K
    Yco2 =0.5    # mol frac
    Yco  =0.1    # mol frac
    Yh2o =0.1    # mol frac
    Nconc = 8    # mol/m3

    k, a = rad.get_k_a(Tg, Nconc, Yco2, Yco, Yh2o)

* This would then be used in conjunction with an ODE solver for the dI_j/ds equations, which we need to write separately.
* In such a solver, the composition and temperature may vary spatially, so we would set Tg, Nconc, Yco2, Yco, Yh2o as above
and redo the call to get_k_a locally as needed.
