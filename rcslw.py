

import numpy as np
from   scipy.interpolate         import RegularGridInterpolator
from   scipy.optimize            import fsolve
from   numpy.polynomial.legendre import leggauss

#------------------------------------------------------------------------------

class rcslw():

    #--------------------------------------------------------------------------

    def __init__(self, P, nGG):

        s = self

        s.P    = P              # pressure (atm)
        s.Tref = 1000.0         # reference temperature (Tb in Falbdf) (K)
        s.nGG  = nGG            # number of grey gases not including the clear gas 
        s.Fmin = 0.02           
        s.Fmax = 0.98

        s.nGGa = s.nGG+1        # number of grey gases not including the clear gas 

        s.P_table    = np.array([0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 30, 50])
        s.C_table    = 0.0001*(1000/0.0001)**(np.arange(0,71)/70.0)
        s.Tg_table   = np.linspace(300.0, 3000.0, 28)
        s.Tb_table   = np.linspace(300.0, 3000.0, 28)
        s.Yh2o_table = np.array([0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0])

        s.nP     = len(s.P_table)     #10
        s.nC     = len(s.C_table)     #71
        s.nTg    = len(s.Tg_table)    #28
        s.nTb    = len(s.Tb_table)    #28
        s.ny_h2o = len(s.Yh2o_table)  #9

        s.set_Falbdf_co2_co_h2o_at_P()
        s.set_interpolating_functions()
        s.set_Fpts()

    #--------------------------------------------------------------------------

    def get_k_a(self, sp, Tg, Y, Nconc):
        '''
        return the local gray gas coefficients (k) and the local weights (a).
        Function only works for co2 or co or h2o, not mixtures.
        sp:    input; string; species name: one of 'co2', 'co', 'h2o'
        Tg:    input; float; gas temperature
        Y:     input; float; mole fraction of co2 or co or h2o
        Nconc: input; float; molar concentration: mol/m3
        '''

        s = self

        Yh2o = Y if sp=='h2o' else None

        C  = np.empty(s.nGG)       # C = C(F, \phi_{loc}, Tref)
        Ct = np.empty(s.nGGa)     # \tilde{C} = C(\tilde{F}, \phi_{loc}, Tref)

        for j in range(s.nGG):
            C[j] = s.get_FI_albdf(sp, s.F_pts[j], Tg, s.Tref, Yh2o)

        for j in range(s.nGGa):
            Ct[j] = s.get_FI_albdf(sp, s.Ft_pts[j], Tg, s.Tref, Yh2o)

        k = np.empty(s.nGGa)
        k[0] = 0.0
        k[1:] = Nconc * Y_coco2h2o * C

        FCt = np.empty(s.nGGa)
        for j in range(s.nGGa):
            FCt[j] = s.get_F_albdf(sp, Ct[j], Tg, Tg)

        a = np.empty(s.nGGa)
        a[0]  = FCt[0]
        a[1:] = FCt[1:] - FCt[0:-1]

        return k, a

    #--------------------------------------------------------------------------

    def get_FI_albdf(self, sp, F, Tg, Tb, Yh2o=None ):
        '''
        Inverse F_albdf: pass in F and get out C.
        sp:   input; string; one of 'co2', 'co', 'h2o'
        C:    input; float; cross section
        Tg:   input; float; gas temperature
        Tb:   input; float; black temperature
        Yh2o: input; float; h2o mole fraction (only used for the h2o table)
        returns C.
        '''

        s = self

        if Tg < s.Tg_table[0] : Tg = s.Tg_table[0]
        if Tg > s.Tg_table[-1]: Tg = s.Tg_table[-1]
        if Tb < s.Tb_table[0] : Tb = s.Tb_table[0]
        if Tb > s.Tb_table[-1]: Tb = s.Tb_table[-1]
        if Yh2o:
            if Yh2o < s.Yh2o_table[0] : Tb = s.Yh2o_table[0]
            if Yh2o > s.Yh2o_table[-1]: Tb = s.Yh2o_table[-1]

        def Func(C):
            if C  < s.C_table[0]  : C  = s.C_table[0]
            if C  > s.C_table[-1] : C  = s.C_table[-1]
            pt = np.array([Yh2o, Tg, Tb, C]) if Yh2o else np.array([Tg, Tb, C])
            return s.interp_F_albdf[sp](pt) - F

        return fsolve(Func, s.C_table[30])[0]


    #--------------------------------------------------------------------------

    def get_F_albdf(self, sp, C, Tg, Tb, Yh2o=None ):
        '''
        sp:   input; string; one of 'co2', 'co', 'h2o'
        C:    input; float; cross section
        Tg:   input; float; gas temperature
        Tb:   input; float; black temperature
        Yh2o: input; float; h2o mole fraction (only used for the h2o table)
        returns the albdf function F
        '''

        s = self

        if C  < s.C_table[0]  : C  = s.C_table[0]
        if C  > s.C_table[-1] : C  = s.C_table[-1]
        if Tg < s.Tg_table[0] : Tg = s.Tg_table[0]
        if Tg > s.Tg_table[-1]: Tg = s.Tg_table[-1]
        if Tb < s.Tb_table[0] : Tb = s.Tb_table[0]
        if Tb > s.Tb_table[-1]: Tb = s.Tb_table[-1]
        if Yh2o:
            if Yh2o < s.Yh2o_table[0] : Tb = s.Yh2o_table[0]
            if Yh2o > s.Yh2o_table[-1]: Tb = s.Yh2o_table[-1]

        pt = np.array([Yh2o, Tg, Tb, C]) if Yh2o else np.array([Tg, Tb, C])

        return s.interp_F_albdf[sp](pt)

    #--------------------------------------------------------------------------

    def set_Fpts(self):
        '''
        Set grid of F points based on Gauss-Legendre Quadrature.
        '''

        s = self

        x,w = leggauss(s.nGG)

        s.Ft_pts = np.empty(s.nGGa)          # \tilde{F} grid
        s.Ft_pts[0] = s.Fmin
        s.Ft_pts[1:] = s.Fmin + (s.Fmax-s.Fmin)*np.cumsum(w)

        s.F_pts = s.Fmin + x*(s.Fmax-s.Fmin)  # F grid  (values between \tild{F} pnts

    #--------------------------------------------------------------------------

    def set_interpolating_functions(self):
        ''' 
        The table is read using multi-linear interpolation. 
        The scipy RegularGridInterpolator returns a function that is then
            called with the desired grid point and returns the ALBDF at 
            that point.
        This function just sets those interpolation functions as a dictionary
            over radiating species.
        '''

        s = self

        s.eval_F_albdf = {}

        sp = 'co'
        s.interp_F_albdf[sp] = RegularGridInterpolator((s.Tg_table, s.Tb_table, s.C_table), s.Falbdf[sp])
        sp = 'co2'
        s.interp_F_albdf[sp] = RegularGridInterpolator((s.Tg_table, s.Tb_table, s.C_table), s.Falbdf[sp])
        sp = 'h2o'
        s.interp_F_albdf[sp] = RegularGridInterpolator((s.Yh2o_table, s.Tg_table, s.Tb_table, s.C_table), s.Falbdf[sp])

    #--------------------------------------------------------------------------

    def set_Falbdf_co2_co_h2o_at_P(self):
        ''' 
        Read the albdf table for co2, co, and h2o.
        Separate files are given for various pressures.
        Interpolate the files to the desired pressure: s.P.
        Reshape from the given 1D to the 3D (co/co2) or 4D (h2o) tables.
        '''

        s = self

        if not (s.Pfiles[0] <= s.P <= s.Pfiles[-1]):
            raise ValueError('Pressure = ', s.P, 'atm is out of range')

        if s.P==s.P_table[0]:
            i1 = 0
            i2 = 1
            P1 = s.P_table[i1]
            P2 = s.P_table[i2]
        elif s.P==s.P_table[-1]:
            i1 = s.nP-2 
            i2 = s.nP-1
            P1 = s.P_table[i1]
            P2 = s.P_table[i2]
        else:
            P1 = s.P_table[P_table <= s.P][-1]
            P2 = s.P_table[P_table >= s.P][0]
            i1 = np.where(s.P_table==P1)[0][0]
            i2 = np.where(s.P_table==P2)[0][0]

        f = (s.P-P1)/(P2-P1)          # interpolation factor

        s.Falbdf = {}

        #------------- CO2

        file1 = 'co2_p' + str(P1).replace('.', '_') + '.txt'
        file1 = 'co2_p' + str(P2).replace('.', '_') + '.txt'

        F1 = np.loadtxt('ALBDF_Tables/'+file1)
        F2 = np.loadtxt('ALBDF_Tables/'+file2)

        s.Falbdf['co2'] = F1*(1-f) + F2*(f)
        s.Falbdf['co2'] = np.reshape(np.Falbdf['co2'], (s.nTg, s.nTb, s.nC))

        #------------- CO 

        file1 = 'co_p' + str(P1).replace('.', '_') + '.txt'
        file1 = 'co_p' + str(P2).replace('.', '_') + '.txt'

        F1 = np.loadtxt('ALBDF_Tables/'+file1)
        F2 = np.loadtxt('ALBDF_Tables/'+file2)

        s.Falbdf['co'] = F1*(1-f) + F2*(f)
        s.Falbdf['co'] = np.reshape(np.Falbdf['co'], (s.nTg, s.nTb, s.nC))

        #------------- H2O

        file1 = 'h2o_p' + str(P1).replace('.', '_') + '.txt'
        file1 = 'h2o_p' + str(P2).replace('.', '_') + '.txt'

        F1 = np.loadtxt('ALBDF_Tables/'+file1)
        F2 = np.loadtxt('ALBDF_Tables/'+file2)

        s.Falbdf['h2o'] = F1*(1-f) + F2*(f)
        s.Falbdf['h2o'] = np.reshape(np.Falbdf['h2o'], (S.ny_h2o, s.nTg, s.nTb, s.nC))

    #--------------------------------------------------------------------------

