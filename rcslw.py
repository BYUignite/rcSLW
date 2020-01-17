

import numpy as np
from   scipy.interpolate         import RegularGridInterpolator
from   scipy.optimize            import fsolve
from   numpy.polynomial.legendre import leggauss
import sys

#------------------------------------------------------------------------------

class rcslw():
    '''
    The Rank Correlated SLW model.
    V.P. Solovjov, F. Andre, D. Lemonnier, B.W. Webb
    http://www.sciencedirect.com/science/article/pii/S0022407316306434
    @author David O. Lignell
    Note: This was verified by comparing to Solovjov's code. 
        The code here is independent of his code, except the inversion function 
        follows the same approach.
    '''

    #--------------------------------------------------------------------------

    def __init__(self, nGG, P, Tg, Yco2, Yco, Yh2o, fvsoot):

        s = self

        s.P    = P              # pressure (atm)
        s.Tref = 1000.0         # reference temperature (Tb in Falbdf) (K)
        s.nGG  = nGG            # number of grey gases not including the clear gas
        s.Cmin = 0.0001
        s.Cmax = 100.0

        s.nGGa = s.nGG+1        # number of grey gases including the clear ga

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

        s.Fmin = s.get_F_albdf(s.Cmin, Tg, Tg, Yco2, Yco, Yh2o, fvsoot)
        s.Fmax = s.get_F_albdf(s.Cmax, Tg, Tg, Yco2, Yco, Yh2o, fvsoot)

        s.set_Fpts()
        #sys.exit() # doldb

    #--------------------------------------------------------------------------

    def get_k_a(self, Tg, Yco2, Yco, Yh2o, fvsoot):
        '''
        THIS IS THE CLASS INTERFACE FUNCTION
        return the local gray gas coefficients (k) and the local weights (a).
        Tg:     input; float; gas temperature
        Yco2:   input; float; mole fraction co2
        Yco:    input; float; mole fraction co
        Yh2o:   input; float; mole fraction h2o
        fvsoot: input; float; soot volume fraction = rho*Ysoot/rhoSoot
        returns k, a arrays
        '''

        s = self

        C  = np.empty(s.nGG)      # C = C(F, \phi_{loc}, Tref)
        Ct = np.empty(s.nGGa)     # \tilde{C} = C(\tilde{F}, \phi_{loc}, Tref)

        for j in range(s.nGG):
            C[j] =  s.get_FI_albdf(s.F_pts[j],  Tg, s.Tref, Yco2, Yco, Yh2o, fvsoot)

        for j in range(s.nGGa):
            Ct[j] = s.get_FI_albdf(s.Ft_pts[j], Tg, s.Tref, Yco2, Yco, Yh2o, fvsoot)

        Nconc = s.P*101325/8.31446/Tg    # mol/m3

        k = np.empty(s.nGGa)
        k[0] = 0.0
        k[1:] = Nconc * C

        FCt = np.empty(s.nGGa)
        for j in range(s.nGGa):
            FCt[j] = s.get_F_albdf(Ct[j], Tg, Tg, Yco2, Yco, Yh2o, fvsoot)

        a = np.empty(s.nGGa)
        a[0]  = FCt[0]
        a[1:] = FCt[1:] - FCt[0:-1]

        return k, a

    #--------------------------------------------------------------------------

    def get_F_albdf(self, C, Tg, Tb, Yco2, Yco, Yh2o, fvsoot):
        '''
        C:    input; float; cross section
        Tg:   input; float; gas temperature
        Tb:   input; float; black temperature
        Yco2: input; float; mole fraction co2
        Yco:  input; float; mole fraction co
        Yh2o: input; float; mole fraction h2o
        returns the albdf function F
        '''

        s = self

        if Yco2 <= 1E-20: Yco2 = 1E-20
        if Yco  <= 1E-20: Yco  = 1E-20
        if Yh2o <= 1E-20: Yh2o = 1E-20

        if Tg   < s.Tg_table[0]   : Tg   = s.Tg_table[0]
        if Tg   > s.Tg_table[-1]  : Tg   = s.Tg_table[-1]
        if Tb   < s.Tb_table[0]   : Tb   = s.Tb_table[0]
        if Tb   > s.Tb_table[-1]  : Tb   = s.Tb_table[-1]
        if Yh2o < s.Yh2o_table[0] : Yh2o = s.Yh2o_table[0]
        if Yh2o > s.Yh2o_table[-1]: Yh2o = s.Yh2o_table[-1]

        CYco2 = C/Yco2
        CYco  = C/Yco
        CYh2o = C/Yh2o

        if CYco2 < s.C_table[0]  : CYco2 = s.C_table[0]
        if CYco2 > s.C_table[-1] : CYco2 = s.C_table[-1]
        if CYco  < s.C_table[0]  : CYco  = s.C_table[0]
        if CYco  > s.C_table[-1] : CYco  = s.C_table[-1]
        if CYh2o < s.C_table[0]  : CYh2o = s.C_table[0]
        if CYh2o > s.C_table[-1] : CYh2o = s.C_table[-1]

        F_co2  = s.interp_F_albdf['co2'](np.array([      Tg, Tb, CYco2 ]))[0]
        F_co   = s.interp_F_albdf['co']( np.array([      Tg, Tb, CYco  ]))[0]
        F_h2o  = s.interp_F_albdf['h2o'](np.array([Yh2o, Tg, Tb, CYh2o ]))[0]

        return F_co2 * F_co * F_h2o * s.F_albdf_soot(C, Tg, Tb, fvsoot)

    #--------------------------------------------------------------------------

    def get_FI_albdf(self, F, Tg, Tb, Yco2, Yco, Yh2o, fvsoot):
        '''
        Inverse F_albdf: pass in F and get out C.
        C:      input; float; cross section
        Tg:     input; float; gas temperature
        Tb:     input; float; black temperature
        Yco2:   input; float; mole fraction co2
        Yco:    input; float; mole fraction co
        Yh2o:   input; float; mole fraction h2o
        fvsoot: input; float; soot volume fraction = rho*Ysoot/rhoSoot
        returns C.
        '''

        s = self

        if F <= s.get_F_albdf(s.Cmin, Tg, Tb, Yco2, Yco, Yh2o, fvsoot):
            return s.Cmin
        if F >= s.get_F_albdf(s.Cmax, Tg, Tb, Yco2, Yco, Yh2o, fvsoot):
            return s.Cmax

        #----------- find table location

        iLo = 0
        iHi = s.nC - 1

        maxit = 100
        for it in range(maxit):
            if iHi-iLo == 1:
                break
            iMi = int((iLo+iHi)/2)
            if F > s.get_F_albdf(s.C_table[iMi], Tg, Tb, Yco2, Yco, Yh2o, fvsoot):
                iLo = iMi
            else:
                iHi = iMi
        if(it==maxit):
            print("WARNING, NO CONVERGENCE IN", maxit, "iterations")
            return s.C_table[iMi]
        #print("# iterations", it)

        #----------- now interpolate to the solution

        Flo = s.get_F_albdf(s.C_table[iLo], Tg, Tb, Yco2, Yco, Yh2o, fvsoot)
        Fhi = s.get_F_albdf(s.C_table[iHi], Tg, Tb, Yco2, Yco, Yh2o, fvsoot)
        cc = s.C_table[iLo] + (F-Flo)*(s.C_table[iHi]-s.C_table[iLo])/(Fhi-Flo)

        #print("relative error:", (F - s.get_F_albdf(cc, Tg, Tb, Yco2, Yco, Yh2o, fvsoot))/F)
        #return cc #doldb

        #----------- but F(cc) is not equal to F! So give it another interpolation:
        # rel error ~ 0.1% --> 1E-6
        
        FF = s.get_F_albdf(cc, Tg, Tb, Yco2, Yco, Yh2o, fvsoot)
        ccc = s.C_table[iLo] + (F-Flo)*(cc-s.C_table[iLo])/(FF-Flo)

        #----------- and one more for fun...
        # rel error ~ 1E-6 --> 1E-10

        FFF = s.get_F_albdf(ccc, Tg, Tb, Yco2, Yco, Yh2o, fvsoot)
        cccc = cc + (F-FF)*(ccc-cc)/(FFF-FF)

        #print("relative error:", (F - s.get_F_albdf(cccc, Tg, Tb, Yco2, Yco, Yh2o, fvsoot))/F)

        return cccc

    #--------------------------------------------------------------------------

    def get_FI_albdf_old(self, F, Tg, Tb, Yco2, Yco, Yh2o, fvsoot):
        '''
        Inverse F_albdf: pass in F and get out C.
        C:    input; float; cross section
        Tg:   input; float; gas temperature
        Tb:   input; float; black temperature
        Yco2: input; float; mole fraction co2
        Yco:  input; float; mole fraction co
        Yh2o: input; float; mole fraction h2o
        fvsoot: input; float; soot volume fraction = rho*Ysoot/rhoSoot
        returns C.
        '''

        s = self

        if F < s.get_F_albdf(s.Cmin, Tg, Tb, Yco2, Yco, Yh2o, fvsoot):
            return s.Cmin
        if F > s.get_F_albdf(s.Cmax, Tg, Tb, Yco2, Yco, Yh2o, fvsoot):
            return s.Cmax

        def Func(C):
            return s.get_F_albdf(C, Tg, Tb, Yco2, Yco, Yh2o, fvsoot) - F
        cc = s.bisection(Func,s.Cmin,s.Cmax,1E-8,100)
        return cc


    #--------------------------------------------------------------------------

    def set_Fpts(self):
        '''
        Set grid of F points based on Gauss-Legendre Quadrature.
        '''

        s = self

        x,w = leggauss(int(s.nGG * 2))
        x   = x[s.nGG:]
        w   = w[s.nGG:]

        s.Ft_pts = np.empty(s.nGGa)          # \tilde{F} grid
        s.Ft_pts[0] = s.Fmin
        s.Ft_pts[1:] = s.Fmin + (s.Fmax-s.Fmin)*np.cumsum(w)


        s.F_pts = s.Fmin + x*(s.Fmax-s.Fmin)  # F grid (vals bet. \tild{F} pnts

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

        s.interp_F_albdf = {}

        sp = 'co'
        s.interp_F_albdf[sp] = RegularGridInterpolator(
                               (s.Tg_table, s.Tb_table, s.C_table), s.Falbdf[sp])
        sp = 'co2'
        s.interp_F_albdf[sp] = RegularGridInterpolator(
                               (s.Tg_table, s.Tb_table, s.C_table), s.Falbdf[sp])
        sp = 'h2o'
        s.interp_F_albdf[sp] = RegularGridInterpolator(
                               (s.Yh2o_table, s.Tg_table, s.Tb_table, s.C_table), s.Falbdf[sp])

    #--------------------------------------------------------------------------

    def set_Falbdf_co2_co_h2o_at_P(self):
        '''
        Read the albdf table for co2, co, and h2o.
        Separate files are given for various pressures.
        Interpolate the files to the desired pressure: s.P.
        Reshape from the given 1D to the 3D (co/co2) or 4D (h2o) tables.
        '''

        s = self

        if not (s.P_table[0] <= s.P <= s.P_table[-1]):
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
            P1 = s.P_table[s.P_table <  s.P][-1]
            P2 = s.P_table[s.P_table >= s.P][0]
            i1 = np.where(s.P_table==P1)[0][0]
            i2 = np.where(s.P_table==P2)[0][0]

        f = (s.P-P1)/(P2-P1)          # interpolation factor

        s.Falbdf = {}

        #------------- CO2

        file1 = 'co2_p' + str(P1).replace('.', '_') + '.txt'
        file2 = 'co2_p' + str(P2).replace('.', '_') + '.txt'

        F1 = np.loadtxt('ALBDF_Tables/'+file1)
        F2 = np.loadtxt('ALBDF_Tables/'+file2)

        s.Falbdf['co2'] = F1*(1-f) + F2*(f)
        s.Falbdf['co2'] = np.reshape(s.Falbdf['co2'], (s.nTg, s.nTb, s.nC))


        #------------- CO

        file1 = 'co_p' + str(P1).replace('.', '_') + '.txt'
        file2 = 'co_p' + str(P2).replace('.', '_') + '.txt'

        F1 = np.loadtxt('ALBDF_Tables/'+file1)
        F2 = np.loadtxt('ALBDF_Tables/'+file2)

        s.Falbdf['co'] = F1*(1-f) + F2*(f)
        s.Falbdf['co'] = np.reshape(s.Falbdf['co'], (s.nTg, s.nTb, s.nC))

        #------------- H2O

        file1 = 'h2o_p' + str(P1).replace('.', '_') + '.txt'
        file2 = 'h2o_p' + str(P2).replace('.', '_') + '.txt'

        F1 = np.loadtxt('ALBDF_Tables/'+file1)
        F2 = np.loadtxt('ALBDF_Tables/'+file2)

        s.Falbdf['h2o'] = F1*(1-f) + F2*(f)
        s.Falbdf['h2o'] = np.reshape(s.Falbdf['h2o'],(s.ny_h2o,s.nTg,s.nTb,s.nC))

    #--------------------------------------------------------------------------

    def F_albdf_soot(self, c, Tg, Tb, fvsoot):
        '''
        C:    input; float; cross section
        Tg:   input; float; gas temperature
        Tb:   input; float; black temperature
        returns the albdf function F for soot
        Note: This comes from Solovjov 2001 and Chang 1984
        http://heattransfer.asmedigitalcollection.asme.org/article.aspx?articleid=1444909
        http://www.sciencedirect.com/science/article/pii/0735193384900514
        '''

        if fvsoot < 1E-12:
            return 1.0

        s = self

        csoot = 7.0                 # general soot parameter
        hCokb = 0.01438777354       # m*K = h*Co/kb = Planck*lightSpeed/Boltzmann (for soot)

        Nconc = s.P*101325/8.31446/Tg # mol/m3

        x = hCokb*c*Nconc / (csoot*fvsoot*Tb)

        n = np.array([1,2,3])       # sum from n=1 to oo, but n=1 to 3 is enough.
        return 1.0 - 15.0/np.pi**4*np.sum(np.exp(-n*x)/n**4*(n*x*(n*x*(n*x+3)+6)+6))

    #--------------------------------------------------------------------------
    
    def bisection(self, f,a,b,tol,maxit):

        s = self
    
        if f(a)*f(b) > 0 :
            print("Error: a, b, don't bracket the root")
            return np.nan,-1
    
        for nit in range(1,maxit+1) :
        
            c = 0.5*(a+b)
            if    f(a)*f(c) < 0.0 : 
                b=c
            else : 
                a=c
        
            if np.abs((b-a)/(0.5*(b+a))) <= tol : 
                print('nit=', nit)
                return 0.5*(a+b)
            if nit==maxit : 
                print('nit=', nit)
                print("No convergence in ", maxit, " iterations")
                return 0.5*(a+b)

#------------------------------------------------------------------------------

P      = 1.0       # atm
Tg     = 1000      # gas temperature
Yco2   = 0.8       # mole fraction co2
Yco    = 0.1       # mole fraction co
Yh2o   = 0.10      # mole fraction h2o
nGG    = 3         # number of gray gases (not including the clear gas)
fvsoot = 0.0      # soot volume fraction (=rho*Ysoot/rhoSoot, where Ysoot = mass frac)

slw   = rcslw(nGG, P, Tg, Yco2, Yco, Yh2o, fvsoot)

k, a = slw.get_k_a(Tg, Yco2, Yco, Yh2o, fvsoot)

print('\n\nk =', k)
print('a =', a)





