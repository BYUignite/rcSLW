/**
 * @file rcslw.h
 * Source file for class rcslw
 */

#include "multi-D.h"
#include "rcslw.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>    // exit

////////////////////////////////////////////////////////////////////////////////
/** Constructor function
 */

rcslw::rcslw(const int    p_nGG, 
             const double p_P,
             const double Tg, 
             const double Yco2, 
             const double Yco, 
             const double Yh2o,
             const double fvsoot){

    P     = p_P;
    nGG   = p_nGG;
    nGGa  = nGG + 1;
    Tref  = 1000.0; 
    Cmin  = 0.0001;
    Cmax  = 100.0;

    P_table  = vector<double>{0.1, 0.25, 0.5, 1,2,4,8,15,30,50};
    C_table  = vector<double>(71);
    for(int  i=0; i<71; i++)
        C_table[i] = 0.0001*pow(1000/0.0001, i/70.0);

    Tg_table = vector<double>(28);
    for(int i=0; i<28; i++)
        Tg_table[i] = 300 + i*(3000.0-300.0)/27;

    Tb_table = vector<double>(28);
    for(int i=0; i<28; i++)
        Tb_table[i] = 300 + i*(3000.0-300.0)/27;

    Yh2o_table = vector<double>{0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0};

    nP     =  P_table.size();        // 10 
    nC     =  C_table.size();        // 71
    nTg    =  Tg_table.size();       // 28
    nTb    =  Tb_table.size();       // 28
    ny_h2o =  Yh2o_table.size();     //  9

    set_Falbdf_co2_co_h2o_at_P();

    Fmin = get_F_albdf(Cmin, Tg, Tg, Yco2, Yco, Yh2o, fvsoot);
    Fmax = get_F_albdf(Cmax, Tg, Tg, Yco2, Yco, Yh2o, fvsoot);

    set_Fpts();

}

////////////////////////////////////////////////////////////////////////////////
/** THIS IS THE CLASS INTERFACE FUNCTION
 *  Given the gas state, set the k and a vectors.
 *  These can then be accessed by the user.
 *  return the local gray gas coefficients (k) and the local weights (a).
 *  Tg:      input; double; gas temperature
 *  Yco2:    input; double; mole fraction co2
 *  Yco:     input; double; mole fraction co
 *  Yh2o:    input; double; mole fraction h2o
 *  fvsoot:  input; double; soot volume fraction = rho*Ysoot/rhosoot
 *  returns k, a vectors
 */

void rcslw::get_k_a(double Tg, double Yco2, double Yco, double Yh2o, double fvsoot, 
                    vector<double> &k, vector<double> &a){

    vector<double> C(nGG);
    vector<double> Ct(nGGa);
   
    vector<double> FCt(nGGa);


    for(int j=0; j < nGG; j++)
        C[j] = get_FI_albdf(F_pts[j], Tg, Tref, Yco2, Yco, Yh2o, fvsoot);
    
    
    for(int j=0; j < nGGa; j++)
        Ct[j] = get_FI_albdf(Ft_pts[j], Tg, Tref, Yco2, Yco, Yh2o, fvsoot);
    
    double Nconc = P*101325/8.31446/Tg;    // mol/m3
    
    k.resize(nGGa);
    k[0] = 0.0;
    for(int j=1; j < nGGa; j++)
        k[j] = Nconc * C[j-1];

    for(int j=0; j < nGGa; j++)
        FCt[j] = get_F_albdf(Ct[j], Tg, Tg, Yco2, Yco, Yh2o, fvsoot);

    a.resize(nGGa);
    a[0] = FCt[0];
    for(int j=1; j < nGGa; j++)
        a[j] = FCt[j] - FCt[j-1];

}

/////////////////////////////////////////////////////////////
/** returns the albdf function F
 *  C:      input; double; cross section
 *  Tg:     input; double; gas temperature  
 *  Tb:     input; double; black temperature
 *  Yco2:   input; double; mole fraction co2
 *  Yco:    input; double; mole fraction co
 *  Yh2o:   input; float; mole fraction h2o
 */

double rcslw::get_F_albdf(double C, double Tg, double Tb, double Yco2, double Yco, double Yh2o,
                        double fvsoot){


    if (Yco2 <= 1E-20)               Yco2 = 1E-20;
    if (Yco <= 1e-20)                Yco = 1e-20;
    if (Yh2o <= 1e-20)               Yh2o = 1e-20;

    if (Tg < Tg_table[0])            Tg = Tg_table[0];
    if (Tg > Tg_table[nTg-1])        Tg = Tg_table[nTg-1];

    if (Tb < Tb_table[0])            Tb = Tb_table[0];
    if (Tb > Tb_table[nTb-1])        Tb = Tb_table[nTb-1];

    if (Yh2o < Yh2o_table[0])        Yh2o = Yh2o_table[0];
    if (Yh2o > Yh2o_table[ny_h2o-1]) Yh2o = Yh2o_table[ny_h2o-1];

    double CYco2 = C/Yco2;
    double CYco  = C/Yco;
    double CYh2o = C/Yh2o;

    if (CYco2  < C_table[0])         CYco2 = C_table[0];
    if (CYco2  > C_table[nC-1])      CYco2 = C_table[nC-1];
    if (CYco   < C_table[0])         CYco  = C_table[0];
    if (CYco   > C_table[nC-1])      CYco  = C_table[nC-1];
    if (CYh2o  < C_table[0])         CYh2o = C_table[0];
    if (CYh2o  > C_table[nC-1])      CYh2o = C_table[nC-1];

    double F_co2, F_co, F_h2o;

    F_co2 = LI_3D(            Tg_table, Tb_table, C_table, Falbdf_co2,        Tg, Tb, CYco2);
    F_co  = LI_3D(            Tg_table, Tb_table, C_table, Falbdf_co,         Tg, Tb, CYco);
    F_h2o = LI_4D(Yh2o_table, Tg_table, Tb_table, C_table, Falbdf_h2o,  Yh2o, Tg, Tb, CYh2o);

    return F_co2 * F_co * F_h2o * F_albdf_soot(C, Tg, Tb, fvsoot);
}


///////////////////////////////////////////////////////
/** Inverse F_albdf: pass in F and get out C.
 *  F:       input; double; cross section
 *  Tg:      input; double; gas temperature
 *  Tb:      input; double; black temperature
 *  Yco2:    input; double; mole fraction co2
 *  Yco:     input; double; mole fraction co
 *  Yh2o:    input; double; mole fraction h2o
 *  fvsoot:  input; double; soot volume fraction = rho*Ysoot/rhoSoot
 *  returns C.
 */

double rcslw::get_FI_albdf(double F, double Tg, double Tb, double Yco2, double Yco, double Yh2o, 
                         double fvsoot){

   if (F <= get_F_albdf(Cmin, Tg, Tb, Yco2, Yco, Yh2o, fvsoot))
       return(Cmin);
   if (F >= get_F_albdf(Cmax, Tg, Tb, Yco2, Yco, Yh2o, fvsoot))
        return(Cmax);

   // Find table location

   int iLo = 0;
   int iHi = nC - 1;      

   int maxit = 100;
   int it;
   for (it=0; it < maxit; it++){
       if (iHi-iLo == 1)
           break;
       int iMi = (iLo+iHi)/2;
       if (F > get_F_albdf(C_table[iMi], Tg, Tb, Yco2, Yco, Yh2o, fvsoot))
           iLo = iMi;
       else
           iHi = iMi;
   }

   int iMi = (iLo+iHi)/2;
   if(it == maxit){
       cout << "WARNING, NO CONVERGENCE IN" << maxit << "iterations" << "\n";
       return C_table[iMi] ;
   }

   // Now interpolate to the solution

   double Flo = get_F_albdf(C_table[iLo], Tg, Tb, Yco2, Yco, Yh2o, fvsoot);
   double Fhi = get_F_albdf(C_table[iHi], Tg, Tb, Yco2, Yco, Yh2o, fvsoot);
   double cc = C_table[iLo] + (F-Flo) * (C_table[iHi]-C_table[iLo])/(Fhi-Flo);

   // cout << "relative error:" << (F - get_F_albdf(cc, Tg, Tb, Yco2, Yco, Yh2o, fvsoot))/F << "\n";
   // return cc; //doldb

   //----------- but F(cc) is not equal to F! So give it another interpolation:
   // rel error ~ 0.1% --> 1E-6

   double FF = get_F_albdf(cc, Tg, Tb, Yco2, Yco, Yh2o, fvsoot);
   double ccc = C_table[iLo] + (F-Flo)*(cc-C_table[iLo])/(FF-Flo);
   
   //----------- and one more for fun...
   // rel error ~ 1E-6 --> 1E-10

   double FFF = get_F_albdf(ccc, Tg, Tb, Yco2, Yco, Yh2o, fvsoot);
   double cccc = cc + (F-FF)*(ccc-cc)/(FFF-FF);

   //cout << "relative error:" << (F - get_F_albdf(ccc, Tg, Tb, Yco2, Yco, Yh2o, fvsoot))/F << "\n";

   return cccc;

}


////////////////////////////////////////////////////////////////
/** Set grid of F points based on Gauss-Legendre Quadrature.
 */

void rcslw::set_Fpts(){

    vector<vector<double> > W(21);
    
    W[0]  = vector<double> {0.0,0.0};
    W[1]  = vector<double> {1.0000000000000000E+00, 1.0000000000000000E+00};
    W[2]  = vector<double> {3.4785484513745368E-01, 6.5214515486254621E-01, 6.5214515486254621E-01, 3.4785484513745368E-01};
    W[3]  = vector<double> {1.7132449237916975E-01, 3.6076157304813894E-01, 4.6791393457269137E-01, 4.6791393457269137E-01, 3.6076157304813894E-01, 1.7132449237916975E-01};
    W[4]  = vector<double> {1.0122853629037669E-01, 2.2238103445337434E-01, 3.1370664587788705E-01, 3.6268378337836177E-01, 3.6268378337836177E-01, 3.1370664587788705E-01, 2.2238103445337434E-01, 1.0122853629037669E-01};
    W[5]  = vector<double> {6.6671344308688069E-02, 1.4945134915058036E-01, 2.1908636251598201E-01, 2.6926671930999652E-01, 2.9552422471475298E-01, 2.9552422471475298E-01, 2.6926671930999652E-01, 2.1908636251598201E-01, 1.4945134915058036E-01, 6.6671344308688069E-02}; 
    W[6]  = vector<double> {4.7175336386512022E-02, 1.0693932599531888E-01, 1.6007832854334611E-01, 2.0316742672306565E-01, 2.3349253653835464E-01, 2.4914704581340269E-01, 2.4914704581340269E-01, 2.3349253653835464E-01, 2.0316742672306565E-01, 1.6007832854334611E-01, 1.0693932599531888E-01, 4.7175336386512022E-02};
    W[7]  = vector<double> {3.5119460331752374E-02, 8.0158087159760305E-02, 1.2151857068790296E-01, 1.5720316715819341E-01, 1.8553839747793763E-01, 2.0519846372129555E-01, 2.1526385346315766E-01, 2.1526385346315766E-01, 2.0519846372129555E-01, 1.8553839747793763E-01, 1.5720316715819341E-01, 1.2151857068790296E-01, 8.0158087159760305E-02, 3.5119460331752374E-02};
    W[8]  = vector<double> {2.7152459411754037E-02, 6.2253523938647706E-02, 9.5158511682492591E-02, 1.2462897125553403E-01, 1.4959598881657676E-01, 1.6915651939500262E-01, 1.8260341504492361E-01, 1.8945061045506859E-01, 1.8945061045506859E-01, 1.8260341504492361E-01, 1.6915651939500262E-01, 1.4959598881657676E-01, 1.2462897125553403E-01, 9.5158511682492591E-02, 6.2253523938647706E-02, 2.7152459411754037E-02};
    W[9]  = vector<double> {2.1616013526484130E-02, 4.9714548894969221E-02, 7.6425730254889246E-02, 1.0094204410628699E-01, 1.2255520671147836E-01, 1.4064291467065063E-01, 1.5468467512626521E-01, 1.6427648374583273E-01, 1.6914238296314363E-01, 1.6914238296314363E-01, 1.6427648374583273E-01, 1.5468467512626521E-01, 1.4064291467065063E-01, 1.2255520671147836E-01, 1.0094204410628699E-01, 7.6425730254889246E-02, 4.9714548894969221E-02, 2.1616013526484130E-02};
    W[10] = vector<double> {1.7614007139153273E-02, 4.0601429800386217E-02, 6.2672048334109443E-02, 8.3276741576704671E-02, 1.0193011981724026E-01, 1.1819453196151825E-01, 1.3168863844917653E-01, 1.4209610931838187E-01, 1.4917298647260366E-01, 1.5275338713072578E-01, 1.5275338713072578E-01, 1.4917298647260366E-01, 1.4209610931838187E-01, 1.3168863844917653E-01, 1.1819453196151825E-01, 1.0193011981724026E-01, 8.3276741576704671E-02, 6.2672048334109443E-02, 4.0601429800386217E-02, 1.7614007139153273E-02}; 
    W[11] = vector<double> {1.4627995298274705E-02, 3.3774901584815178E-02, 5.2293335152682870E-02, 6.9796468424520197E-02, 8.5941606217067396E-02, 1.0041414444288072E-01, 1.1293229608053883E-01, 1.2325237681051199E-01, 1.3117350478706188E-01, 1.3654149834601478E-01, 1.3925187285563156E-01, 1.3925187285563156E-01, 1.3654149834601478E-01, 1.3117350478706188E-01, 1.2325237681051199E-01, 1.1293229608053883E-01, 1.0041414444288072E-01, 8.5941606217067396E-02, 6.9796468424520197E-02, 5.2293335152682870E-02, 3.3774901584815178E-02, 1.4627995298274705E-02}; 
    W[12] = vector<double> {1.2341229799987091E-02, 2.8531388628933743E-02, 4.4277438817419551E-02, 5.9298584915436742E-02, 7.3346481411080411E-02, 8.6190161531953288E-02, 9.7618652104114065E-02, 1.0744427011596561E-01, 1.1550566805372561E-01, 1.2167047292780342E-01, 1.2583745634682830E-01, 1.2793819534675221E-01, 1.2793819534675221E-01, 1.2583745634682830E-01, 1.2167047292780342E-01, 1.1550566805372561E-01, 1.0744427011596561E-01, 9.7618652104114065E-02, 8.6190161531953288E-02, 7.3346481411080411E-02, 5.9298584915436742E-02, 4.4277438817419551E-02, 2.8531388628933743E-02, 1.2341229799987091E-02}; 
    W[13] = vector<double> {1.0551372617343395E-02, 2.4417851092631938E-02, 3.7962383294363120E-02, 5.0975825297148079E-02, 6.3274046329574674E-02, 7.4684149765659763E-02, 8.5045894313485068E-02, 9.4213800355914160E-02, 1.0205916109442532E-01, 1.0847184052857647E-01, 1.1336181654631956E-01, 1.1666044348529646E-01, 1.1832141527926213E-01, 1.1832141527926213E-01, 1.1666044348529646E-01, 1.1336181654631956E-01, 1.0847184052857647E-01, 1.0205916109442532E-01, 9.4213800355914160E-02, 8.5045894313485068E-02, 7.4684149765659763E-02, 6.3274046329574674E-02, 5.0975825297148079E-02, 3.7962383294363120E-02, 2.4417851092631938E-02, 1.0551372617343395E-02}; 
    W[14] = vector<double> {9.1242825930943974E-03, 2.1132112592771271E-02, 3.2901427782304517E-02, 4.4272934759003985E-02, 5.5107345675716936E-02, 6.5272923966999755E-02, 7.4646214234568811E-02, 8.3113417228900935E-02, 9.0571744393032852E-02, 9.6930657997929923E-02, 1.0211296757806078E-01, 1.0605576592284637E-01, 1.0871119225829413E-01, 1.1004701301647524E-01, 1.1004701301647524E-01, 1.0871119225829413E-01, 1.0605576592284637E-01, 1.0211296757806078E-01, 9.6930657997929923E-02, 9.0571744393032852E-02, 8.3113417228900935E-02, 7.4646214234568811E-02, 6.5272923966999755E-02, 5.5107345675716936E-02, 4.4272934759003985E-02, 3.2901427782304517E-02, 2.1132112592771271E-02, 9.1242825930943974E-03};
    W[15] = vector<double> {7.9681924961695228E-03, 1.8466468311091087E-02, 2.8784707883322873E-02, 3.8799192569626793E-02, 4.8402672830594434E-02, 5.7493156217619093E-02, 6.5974229882180324E-02, 7.3755974737704802E-02, 8.0755895229419811E-02, 8.6899787201082698E-02, 9.2122522237785789E-02, 9.6368737174643990E-02, 9.9593420586794934E-02, 1.0176238974840521E-01, 1.0285265289355848E-01, 1.0285265289355848E-01, 1.0176238974840521E-01, 9.9593420586794934E-02, 9.6368737174643990E-02, 9.2122522237785789E-02, 8.6899787201082698E-02, 8.0755895229419811E-02, 7.3755974737704802E-02, 6.5974229882180324E-02, 5.7493156217619093E-02, 4.8402672830594434E-02, 3.8799192569626793E-02, 2.8784707883322873E-02, 1.8466468311091087E-02, 7.9681924961695228E-03}; 
    W[16] = vector<double> {7.0186100094692984E-03, 1.6274394730905965E-02, 2.5392065309262427E-02, 3.4273862913021626E-02, 4.2835898022226426E-02, 5.0998059262376244E-02, 5.8684093478535704E-02, 6.5822222776361752E-02, 7.2345794108848449E-02, 7.8193895787070311E-02, 8.3311924226946846E-02, 8.7652093004403908E-02, 9.1173878695763863E-02, 9.3844399080804566E-02, 9.5638720079274833E-02, 9.6540088514727812E-02, 9.6540088514727812E-02, 9.5638720079274833E-02, 9.3844399080804566E-02, 9.1173878695763863E-02, 8.7652093004403908E-02, 8.3311924226946846E-02, 7.8193895787070311E-02, 7.2345794108848449E-02, 6.5822222776361752E-02, 5.8684093478535704E-02, 5.0998059262376244E-02, 4.2835898022226426E-02, 3.4273862913021626E-02, 2.5392065309262427E-02, 1.6274394730905965E-02, 7.0186100094692984E-03};
    W[17] = vector<double> {6.2291405559100326E-03, 1.4450162748594548E-02, 2.2563721985495038E-02, 3.0491380638445777E-02, 3.8166593796386906E-02, 4.5525611523353778E-02, 5.2507414572678060E-02, 5.9054135827524654E-02, 6.5111521554076457E-02, 7.0629375814255727E-02, 7.5561974660031797E-02, 7.9868444339771819E-02, 8.3513099699845522E-02, 8.6465739747035669E-02, 8.8701897835693738E-02, 9.0203044370640612E-02, 9.0956740330259786E-02, 9.0956740330259786E-02, 9.0203044370640612E-02, 8.8701897835693738E-02, 8.6465739747035669E-02, 8.3513099699845522E-02, 7.9868444339771819E-02, 7.5561974660031797E-02, 7.0629375814255727E-02, 6.5111521554076457E-02, 5.9054135827524654E-02, 5.2507414572678060E-02, 4.5525611523353778E-02, 3.8166593796386906E-02, 3.0491380638445777E-02, 2.2563721985495038E-02, 1.4450162748594548E-02, 6.2291405559100326E-03}; 
    W[18] = vector<double> {5.5657196642477837E-03, 1.2915947284064104E-02, 2.0181515297735174E-02, 2.7298621498568355E-02, 3.4213810770307482E-02, 4.0875750923645232E-02, 4.7235083490266047E-02, 5.3244713977759678E-02, 5.8860144245324549E-02, 6.4039797355015429E-02, 6.8745323835736297E-02, 7.2941885005653004E-02, 7.6598410645870627E-02, 7.9687828912071559E-02, 8.2187266704339651E-02, 8.4078218979661792E-02, 8.5346685739338499E-02, 8.5983275670394627E-02, 8.5983275670394627E-02, 8.5346685739338499E-02, 8.4078218979661792E-02, 8.2187266704339651E-02, 7.9687828912071559E-02, 7.6598410645870627E-02, 7.2941885005653004E-02, 6.8745323835736297E-02, 6.4039797355015429E-02, 5.8860144245324549E-02, 5.3244713977759678E-02, 4.7235083490266047E-02, 4.0875750923645232E-02, 3.4213810770307482E-02, 2.7298621498568355E-02, 2.0181515297735174E-02, 1.2915947284064104E-02, 5.5657196642477837E-03}; 
    W[19] = vector<double> {5.0028807496389329E-03, 1.1613444716468455E-02, 1.8156577709613410E-02, 2.4579739738232007E-02, 3.0839500545175719E-02, 3.6894081594024859E-02, 4.2703158504674703E-02, 4.8228061860758530E-02, 5.3432019910332196E-02, 5.8280399146997154E-02, 6.2740933392133172E-02, 6.6783937979140381E-02, 7.0382507066898942E-02, 7.3512692584743411E-02, 7.6153663548446479E-02, 7.8287844658210981E-02, 7.9901033243527819E-02, 8.0982493770597061E-02, 8.1525029280385797E-02, 8.1525029280385797E-02, 8.0982493770597061E-02, 7.9901033243527819E-02, 7.8287844658210981E-02, 7.6153663548446479E-02, 7.3512692584743411E-02, 7.0382507066898942E-02, 6.6783937979140381E-02, 6.2740933392133172E-02, 5.8280399146997154E-02, 5.3432019910332196E-02, 4.8228061860758530E-02, 4.2703158504674703E-02, 3.6894081594024859E-02, 3.0839500545175719E-02, 2.4579739738232007E-02, 1.8156577709613410E-02, 1.1613444716468455E-02, 5.0028807496389329E-03};
    W[20] = vector<double> {4.5212770985300181E-03, 1.0498284531151609E-02, 1.6421058381907345E-02, 2.2245849194166653E-02, 2.7937006980023528E-02, 3.3460195282547678E-02, 3.8782167974472377E-02, 4.3870908185673324E-02, 4.8695807635072405E-02, 5.3227846983937115E-02, 5.7439769099391892E-02, 6.1306242492929319E-02, 6.4804013456601486E-02, 6.7912045815234398E-02, 7.0611647391287169E-02, 7.2886582395804478E-02, 7.4723169057968677E-02, 7.6110361900626741E-02, 7.7039818164248389E-02, 7.7505947978425332E-02, 7.7505947978425332E-02, 7.7039818164248389E-02, 7.6110361900626741E-02, 7.4723169057968677E-02, 7.2886582395804478E-02, 7.0611647391287169E-02, 6.7912045815234398E-02, 6.4804013456601486E-02, 6.1306242492929319E-02, 5.7439769099391892E-02, 5.3227846983937115E-02, 4.8695807635072405E-02, 4.3870908185673324E-02, 3.8782167974472377E-02, 3.3460195282547678E-02, 2.7937006980023528E-02, 2.2245849194166653E-02, 1.6421058381907345E-02, 1.0498284531151609E-02, 4.5212770985300181E-03};
    
    
    vector<vector<double> > X(21);
    
    X[0]  = vector<double> {0.0, 0.0};
    X[1]  = vector<double> {-5.7735026918962573E-01,  5.7735026918962573E-01};
    X[2]  = vector<double> {-8.6113631159405257E-01, -3.3998104358485626E-01,  3.3998104358485626E-01,  8.6113631159405257E-01};
    X[3]  = vector<double> {-9.3246951420315205E-01, -6.6120938646626448E-01, -2.3861918608319693E-01,  2.3861918608319693E-01,  6.6120938646626448E-01,  9.3246951420315205E-01};
    X[4]  = vector<double> {-9.6028985649753618E-01, -7.9666647741362673E-01, -5.2553240991632899E-01, -1.8343464249564978E-01,  1.8343464249564978E-01,  5.2553240991632899E-01,  7.9666647741362673E-01,  9.6028985649753618E-01};
    X[5]  = vector<double> {-9.7390652851717174E-01, -8.6506336668898454E-01, -6.7940956829902444E-01, -4.3339539412924721E-01, -1.4887433898163122E-01,  1.4887433898163122E-01,  4.3339539412924721E-01,  6.7940956829902444E-01,  8.6506336668898454E-01,  9.7390652851717174E-01};
    X[6]  = vector<double> {-9.8156063424671924E-01, -9.0411725637047480E-01, -7.6990267419430469E-01, -5.8731795428661748E-01, -3.6783149899818018E-01, -1.2523340851146891E-01,  1.2523340851146891E-01,  3.6783149899818018E-01,  5.8731795428661748E-01,  7.6990267419430469E-01,  9.0411725637047480E-01,  9.8156063424671924E-01};
    X[7]  = vector<double> {-9.8628380869681231E-01, -9.2843488366357352E-01, -8.2720131506976502E-01, -6.8729290481168548E-01, -5.1524863635815410E-01, -3.1911236892788974E-01, -1.0805494870734367E-01,  1.0805494870734367E-01,  3.1911236892788974E-01,  5.1524863635815410E-01,  6.8729290481168548E-01,  8.2720131506976502E-01,  9.2843488366357352E-01,  9.8628380869681231E-01};
    X[8]  = vector<double> {-9.8940093499164994E-01, -9.4457502307323260E-01, -8.6563120238783176E-01, -7.5540440835500300E-01, -6.1787624440264377E-01, -4.5801677765722737E-01, -2.8160355077925892E-01, -9.5012509837637454E-02,  9.5012509837637454E-02,  2.8160355077925892E-01,  4.5801677765722737E-01,  6.1787624440264377E-01,  7.5540440835500300E-01,  8.6563120238783176E-01,  9.4457502307323260E-01,  9.8940093499164994E-01};
    X[9]  = vector<double> {-9.9156516842093090E-01, -9.5582394957139782E-01, -8.9260246649755570E-01, -8.0370495897252314E-01, -6.9168704306035322E-01, -5.5977083107394754E-01, -4.1175116146284263E-01, -2.5188622569150548E-01, -8.4775013041735292E-02,  8.4775013041735292E-02,  2.5188622569150548E-01,  4.1175116146284263E-01,  5.5977083107394754E-01,  6.9168704306035322E-01,  8.0370495897252314E-01,  8.9260246649755570E-01,  9.5582394957139782E-01,  9.9156516842093090E-01};
    X[10] = vector<double> {-9.9312859918509488E-01, -9.6397192727791381E-01, -9.1223442825132584E-01, -8.3911697182221878E-01, -7.4633190646015080E-01, -6.3605368072651502E-01, -5.1086700195082713E-01, -3.7370608871541955E-01, -2.2778585114164510E-01, -7.6526521133497338E-02,  7.6526521133497338E-02,  2.2778585114164510E-01,  3.7370608871541955E-01,  5.1086700195082713E-01,  6.3605368072651502E-01,  7.4633190646015080E-01,  8.3911697182221878E-01,  9.1223442825132584E-01,  9.6397192727791381E-01,  9.9312859918509488E-01};
    X[11] = vector<double> {-9.9429458548239924E-01, -9.7006049783542869E-01, -9.2695677218717398E-01, -8.6581257772030018E-01, -7.8781680597920811E-01, -6.9448726318668275E-01, -5.8764040350691160E-01, -4.6935583798675706E-01, -3.4193582089208424E-01, -2.0786042668822130E-01, -6.9739273319722211E-02,  6.9739273319722211E-02,  2.0786042668822130E-01,  3.4193582089208424E-01,  4.6935583798675706E-01,  5.8764040350691160E-01,  6.9448726318668275E-01,  7.8781680597920811E-01,  8.6581257772030018E-01,  9.2695677218717398E-01, 9.7006049783542869E-01, 9.9429458548239924E-01};
    X[12] = vector<double> {-9.9518721999702131E-01, -9.7472855597130947E-01, -9.3827455200273280E-01, -8.8641552700440096E-01, -8.2000198597390295E-01, -7.4012419157855436E-01, -6.4809365193697555E-01, -5.4542147138883956E-01, -4.3379350762604513E-01, -3.1504267969616340E-01, -1.9111886747361631E-01, -6.4056892862605630E-02,  6.4056892862605630E-02,  1.9111886747361631E-01,  3.1504267969616340E-01,  4.3379350762604513E-01,  5.4542147138883956E-01,  6.4809365193697555E-01,  7.4012419157855436E-01,  8.2000198597390295E-01, 8.8641552700440096E-01, 9.3827455200273280E-01, 9.7472855597130947E-01, 9.9518721999702131E-01};
    X[13] = vector<double> {-9.9588570114561692E-01, -9.7838544595647092E-01, -9.4715906666171423E-01, -9.0263786198430707E-01, -8.4544594278849805E-01, -7.7638594882067880E-01, -6.9642726041995728E-01, -6.0669229301761807E-01, -5.0844071482450570E-01, -4.0305175512348629E-01, -2.9200483948595690E-01, -1.7685882035689018E-01, -5.9230093429313208E-02,  5.9230093429313208E-02,  1.7685882035689018E-01,  2.9200483948595690E-01,  4.0305175512348629E-01,  5.0844071482450570E-01,  6.0669229301761807E-01,  6.9642726041995728E-01, 7.7638594882067880E-01, 8.4544594278849805E-01, 9.0263786198430707E-01, 9.4715906666171423E-01, 9.7838544595647092E-01, 9.9588570114561692E-01};
    X[14] = vector<double> {-9.9644249757395442E-01, -9.8130316537087281E-01, -9.5425928062893817E-01, -9.1563302639213207E-01, -8.6589252257439497E-01, -8.0564137091717913E-01, -7.3561087801363179E-01, -6.5665109403886501E-01, -5.6972047181140173E-01, -4.7587422495511827E-01, -3.7625151608907870E-01, -2.7206162763517810E-01, -1.6456928213338079E-01, -5.5079289884034266E-02,  5.5079289884034266E-02,  1.6456928213338079E-01,  2.7206162763517810E-01,  3.7625151608907870E-01,  4.7587422495511827E-01,  5.6972047181140173E-01, 6.5665109403886501E-01, 7.3561087801363179E-01, 8.0564137091717913E-01, 8.6589252257439497E-01, 9.1563302639213207E-01, 9.5425928062893817E-01, 9.8130316537087281E-01, 9.9644249757395442E-01};
    X[15] = vector<double> {-9.9689348407464951E-01, -9.8366812327974729E-01, -9.6002186496830755E-01, -9.2620004742927431E-01, -8.8256053579205263E-01, -8.2956576238276836E-01, -7.6777743210482619E-01, -6.9785049479331585E-01, -6.2052618298924289E-01, -5.3662414814201986E-01, -4.4703376953808915E-01, -3.5270472553087812E-01, -2.5463692616788985E-01, -1.5386991360858354E-01, -5.1471842555317698E-02,  5.1471842555317698E-02,  1.5386991360858354E-01,  2.5463692616788985E-01,  3.5270472553087812E-01,  4.4703376953808915E-01, 5.3662414814201986E-01, 6.2052618298924289E-01, 6.9785049479331585E-01, 7.6777743210482619E-01, 8.2956576238276836E-01, 8.8256053579205263E-01, 9.2620004742927431E-01, 9.6002186496830755E-01, 9.8366812327974729E-01, 9.9689348407464951E-01};
    X[16] = vector<double> {-9.9726386184948157E-01, -9.8561151154526838E-01, -9.6476225558750639E-01, -9.3490607593773967E-01, -8.9632115576605220E-01, -8.4936761373256997E-01, -7.9448379596794239E-01, -7.3218211874028971E-01, -6.6304426693021523E-01, -5.8771575724076230E-01, -5.0689990893222936E-01, -4.2135127613063533E-01, -3.3186860228212767E-01, -2.3928736225213706E-01, -1.4447196158279649E-01, -4.8307665687738310E-02,  4.8307665687738310E-02,  1.4447196158279649E-01,  2.3928736225213706E-01,  3.3186860228212767E-01, 4.2135127613063533E-01, 5.0689990893222936E-01, 5.8771575724076230E-01, 6.6304426693021523E-01, 7.3218211874028971E-01, 7.9448379596794239E-01, 8.4936761373256997E-01, 8.9632115576605220E-01, 9.3490607593773967E-01, 9.6476225558750639E-01, 9.8561151154526838E-01, 9.9726386184948157E-01};
    X[17] = vector<double> {-9.9757175379084195E-01, -9.8722781640630952E-01, -9.6870826253334430E-01, -9.4216239740510710E-01, -9.0780967771832455E-01, -8.6593463833456441E-01, -8.1688422790093362E-01, -7.6106487662987299E-01, -6.9893911321626290E-01, -6.3102172708052851E-01, -5.5787550066974667E-01, -4.8010654519032703E-01, -3.9835927775864594E-01, -3.1331108133946323E-01, -2.2566669161644948E-01, -1.3615235725918298E-01, -4.5509821953102533E-02,  4.5509821953102533E-02,  1.3615235725918298E-01,  2.2566669161644948E-01, 3.1331108133946323E-01, 3.9835927775864594E-01, 4.8010654519032703E-01, 5.5787550066974667E-01, 6.3102172708052851E-01, 6.9893911321626290E-01, 7.6106487662987299E-01, 8.1688422790093362E-01, 8.6593463833456441E-01, 9.0780967771832455E-01, 9.4216239740510710E-01, 9.6870826253334430E-01, 9.8722781640630952E-01, 9.9757175379084195E-01};
    X[18] = vector<double> {-9.9783046248408580E-01, -9.8858647890221230E-01, -9.7202769104969799E-01, -9.4827298439950758E-01, -9.1749777451565906E-01, -8.7992980089039707E-01, -8.3584716699247530E-01, -7.8557623013220657E-01, -7.2948917159355664E-01, -6.6800123658552102E-01, -6.0156765813598057E-01, -5.3068028592624517E-01, -4.5586394443342027E-01, -3.7767254711968923E-01, -2.9668499534402826E-01, -2.1350089231686559E-01, -1.2873610380938480E-01, -4.3018198473708608E-02,  4.3018198473708608E-02,  1.2873610380938480E-01, 2.1350089231686559E-01, 2.9668499534402826E-01, 3.7767254711968923E-01, 4.5586394443342027E-01, 5.3068028592624517E-01, 6.0156765813598057E-01, 6.6800123658552102E-01, 7.2948917159355664E-01, 7.8557623013220657E-01, 8.3584716699247530E-01, 8.7992980089039707E-01, 9.1749777451565906E-01, 9.4827298439950758E-01, 9.7202769104969799E-01, 9.8858647890221230E-01, 9.9783046248408580E-01};
    X[19] = vector<double> {-9.9804993053568758E-01, -9.8973945426638554E-01, -9.7484632859015352E-01, -9.5346633093352962E-01, -9.2574133204858433E-01, -8.9185573900463222E-01, -8.5203502193236214E-01, -8.0654416760531689E-01, -7.5568590375397071E-01, -6.9979868037918436E-01, -6.3925441582968168E-01, -5.7445602104780713E-01, -5.0583471792793111E-01, -4.3384716943237650E-01, -3.5897244047943500E-01, -2.8170880979016527E-01, -2.0257045389211670E-01, -1.2208402533786741E-01, -4.0785147904578239E-02,  4.0785147904578239E-02, 1.2208402533786741E-01, 2.0257045389211670E-01, 2.8170880979016527E-01, 3.5897244047943500E-01, 4.3384716943237650E-01, 5.0583471792793111E-01, 5.7445602104780713E-01, 6.3925441582968168E-01, 6.9979868037918436E-01, 7.5568590375397071E-01, 8.0654416760531689E-01, 8.5203502193236214E-01, 8.9185573900463222E-01, 9.2574133204858433E-01, 9.5346633093352962E-01, 9.7484632859015352E-01, 9.8973945426638554E-01, 9.9804993053568758E-01};
    X[20] = vector<double> {-9.9823770971055925E-01, -9.9072623869945708E-01, -9.7725994998377430E-01, -9.5791681921379168E-01, -9.3281280827867652E-01, -9.0209880696887434E-01, -8.6595950321225956E-01, -8.2461223083331170E-01, -7.7830565142651942E-01, -7.2731825518992710E-01, -6.7195668461417957E-01, -6.1255388966798030E-01, -5.4946712509512818E-01, -4.8307580168617870E-01, -4.1377920437160498E-01, -3.4199409082575849E-01, -2.6815218500725369E-01, -1.9269758070137111E-01, -1.1608407067525521E-01, -3.8772417506050816E-02, 3.8772417506050816E-02, 1.1608407067525521E-01, 1.9269758070137111E-01, 2.6815218500725369E-01, 3.4199409082575849E-01, 4.1377920437160498E-01, 4.8307580168617870E-01, 5.4946712509512818E-01, 6.1255388966798030E-01, 6.7195668461417957E-01, 7.2731825518992710E-01, 7.7830565142651942E-01, 8.2461223083331170E-01, 8.6595950321225956E-01, 9.0209880696887434E-01, 9.3281280827867652E-01, 9.5791681921379168E-01, 9.7725994998377430E-01, 9.9072623869945708E-01, 9.9823770971055925E-01};


    vector<double> cumsum_w(nGG, W[nGG][nGG]);
    for(int i=1; i<nGG; i++)
        cumsum_w[i] = cumsum_w[i-1] + W[nGG][i+nGG];

    Ft_pts.resize(nGGa);                            // \tilde{F} grid
    Ft_pts[0] = Fmin;
    for(int i=1; i<nGGa; i++)
        Ft_pts[i] = Fmin + (Fmax-Fmin)*cumsum_w[i-1];   

    F_pts.resize(nGG);
    for(int i=nGG; i<X[nGG].size(); i++)
        F_pts[i-nGG] = Fmin + X[nGG][i]*(Fmax - Fmin);    // F grid (vals bet. \tild{F} pnts)

}

////////////////////////////////////////////////////////////////
/** Read the albdf table for co2, co, and h2o.
 *  Separate files are given for various pressures.
 *  Interpolate the files to the desired pressure: P.
 *  Reshape from the given 1D to the 3D (co/co2) or 4D (h2o) tables.
 */

void rcslw::set_Falbdf_co2_co_h2o_at_P(){

    double P1;
    double P2;
    double f;

    if(P < P_table[0] || P > P_table.back() ) {
        cout << "Pressure = " << P << " atm is out of range\n";
        exit(0);
    }
    else if(P==P_table[0]){
        P1 = P_table[0];
        P2 = P_table[1];
    }
    else if(P==P_table.back()){
        P1 = P_table[nP-2];
        P2 = P_table[nP-1];
    }
    else{
        int i;
        for(i=0; P_table[i] < P; i++)
            P1 = P_table[i];
        P2 = P_table[i];
    }

    f  = (P-P1)/(P2-P1);               // Interpolation factor

    string Pres_1 = to_string(P1);
    string Pres_2 = to_string(P2);
    Pres_1.erase(Pres_1.find_last_not_of('0')+1, string::npos);
    Pres_2.erase(Pres_2.find_last_not_of('0')+2, string::npos);
    if(Pres_1.back() == '.') Pres_1.push_back('0');
    if(Pres_2.back() == '.') Pres_2.push_back('0');
    replace(Pres_1.begin(), Pres_1.end(), '.', '_');
    replace(Pres_2.begin(), Pres_2.end(), '.', '_');

    //--------------------- CO2
    
    string co2_file1 = "ALBDF_Tables/co2_p" + Pres_1 + ".txt";
    string co2_file2 = "ALBDF_Tables/co2_p" + Pres_2 + ".txt";

    vector<vector<vector<double> > > co2_F1(nTg, vector<vector<double> >(nTb, vector<double>(nC))); 
    vector<vector<vector<double> > > co2_F2(nTg, vector<vector<double> >(nTb, vector<double>(nC))); 
    Falbdf_co2.resize(nTg, vector<vector<double> >(nTb, vector<double>(nC)));

    get_FI_albdf_tables(co2_file1, nTg, nTb, nC, co2_F1);
    get_FI_albdf_tables(co2_file2, nTg, nTb, nC, co2_F2);

    for (int i=0; i< nTg; i++){
        for(int j=0; j<nTb; j++){
            for(int k=0; k<nC; k++){
                Falbdf_co2[i][j][k] = co2_F1[i][j][k]*(1-f) + co2_F2[i][j][k]*f;
            }
        }
    }


    //---------------------CO

    string co_file1 = "ALBDF_Tables/co_p" + Pres_1 + ".txt";
    string co_file2 = "ALBDF_Tables/co_p" + Pres_2 + ".txt";

    vector<vector<vector<double> > > co_F1(nTg, vector<vector<double> >(nTb, vector<double>(nC))); 
    vector<vector<vector<double> > > co_F2(nTg, vector<vector<double> >(nTb, vector<double>(nC))); 
    Falbdf_co.resize(nTg, vector<vector<double> >(nTb, vector<double>(nC)));

    get_FI_albdf_tables(co_file1, nTg, nTb, nC, co_F1);
    get_FI_albdf_tables(co_file2, nTg, nTb, nC, co_F2);

    for (int i=0; i< nTg; i++){
        for(int j=0; j<nTb; j++){
            for(int k=0; k<nC; k++){
                Falbdf_co[i][j][k] = co_F1[i][j][k]*(1-f) + co_F2[i][j][k]*f;
            }
        }
    } 

    //---------------------H2O
   
    string h2o_file1 = "ALBDF_Tables/h2o_p" + Pres_1 + ".txt";
    string h2o_file2 = "ALBDF_Tables/h2o_p" + Pres_2 + ".txt";

    vector<vector<vector<vector<double> > > > h2o_F1(ny_h2o, vector<vector<vector<double> > >(nTg, vector<vector<double> >(nTb, vector<double>(nC)))); 
    vector<vector<vector<vector<double> > > > h2o_F2(ny_h2o, vector<vector<vector<double> > >(nTg, vector<vector<double> >(nTb, vector<double>(nC)))); 
    Falbdf_h2o.resize(ny_h2o,vector<vector<vector<double> > >(nTg, vector<vector<double> >(nTb, vector<double>(nC))));

    get_FI_albdf_tables(h2o_file1, ny_h2o, nTg, nTb, nC, h2o_F1);
    get_FI_albdf_tables(h2o_file2, ny_h2o, nTg, nTb, nC, h2o_F2);

    for(int i=0; i< ny_h2o; i++){
        for (int j=0; j< nTg; j++){
            for(int k=0; k<nTb; k++){
                for(int n=0; n<nC; n++){
                    Falbdf_h2o[i][j][k][n] = h2o_F1[i][j][k][n]*(1-f) + h2o_F2[i][j][k][n]*f;
                }
            }
        }
    }
  
}
    
/////////////////////////////////////////////////////////////////////////
/** Import ALBDF files to matrices
 */

void rcslw::get_FI_albdf_tables(string Ptable_file_name, int nx, int ny, int nz, 
        vector<vector<vector<double> > > &myarray){

    ifstream ifile;
    ifile.open(Ptable_file_name);
    if(!ifile){
        cout << endl << "error opening file: " << Ptable_file_name << endl;
        exit(0);
    }
    for (int i=0; i< nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                ifile >> myarray[i][j][k];
            }
        }
    }
    ifile.close();
    
}

/////////////////////////////////////////////////////////////////////////
/** Import ALBDF files to matrices
 */

void rcslw::get_FI_albdf_tables(string Ptable_file_name, int nx, int ny, int nz, int nw,
        vector<vector<vector<vector<double> > > > &myarray){
 
    ifstream ifile;
    ifile.open(Ptable_file_name);
    if(!ifile){
        cout << endl << "error opening file: " << Ptable_file_name << endl;
        exit(0);
    }
    for (int i=0; i< nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                for(int n=0; n<nw; n++){
                    ifile >> myarray[i][j][k][n];
                }
            }
        }
    }
    ifile.close();
    
}
////////////////////////////////////////////////////////////////
/** C;   input; double; cross section
 *  Tg;  input; double; gas temperature
 *  Tb;  input; double; black temperature
 *  returns the albdf function F for soot
 *  Note: This comes from Solovjov 2001 and Chang 1984
 *  http://heattransfer.asmedigitalcollection.asme.org/article.aspx?articleid=1444909
 *  http://www.sciencedirect.com/science/article/pii/0735193384900514
 */

double rcslw::F_albdf_soot(double c, double Tg, double Tb, double fvsoot){

    if (fvsoot < 1E-12)
        return 1.0;
    
    double csoot = 7.0;                // general soot parameter
    double hCokb = 0.01438777354;      // m*K = h*Co/kb = Planc*lightSpeed/Boltzmann (for soot)
    
    double Nconc = P*101325/8.31446/Tg;      // mol/m3

    double x = hCokb*c*Nconc / (csoot*fvsoot*Tb);

    vector<double> n = {1, 2, 3};      // sum from n=1 to oo, but n=1 to 3 is enough

    double theSum = 0.0;
    double nx;
    for(int i=0; i<n.size(); i++){
        nx = n[i]*x;
        theSum += exp(-nx)/pow(n[i],4)*(nx*(nx*(nx+3)+6)+6);
    }
    return 1.0 - 15.0/pow(M_PI,4)*theSum;
}

////////////////////////////////////////////////////////

int main(){

    double P      = 1.0;                     // atm
    double Tg     = 1000.0;                  // gas temperature
    double Yco2   = 0.9;                     // mole fraction co2
    double Yco    = 0.1;                     // mole fraction co
    double Yh2o   = 0.00;                    // mole fraction h2o
    int    nGG    = 3;                       // number of gray gases (not including the clear gas)
    double fvsoot = 0.0;                     // soot volume fraction (=rho*Ysoot/rhoSoot, where Ysoot = mass frac)

    rcslw slw(nGG, P, Tg, Yco2, Yco, Yh2o, fvsoot);

    vector<double> k;
    vector<double> a; 

    slw.get_k_a(Tg, Yco2, Yco, Yh2o, fvsoot, k, a);

    cout << "k = ";
    for(int j=0; j<k.size(); j++)
        cout << k[j] << " ";
    
    cout << "\n\n" << "a = ";
    for(int j=0; j<a.size(); j++)
        cout << a[j] << " "; 

    return (0);
}
