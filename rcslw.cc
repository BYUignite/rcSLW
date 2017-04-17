/**
 * @file rcslw.cc
 * Source file for class rcslw
 */

#include "rcslw.h"
#include <cmath>       // pow


using namespace std;

////////////////////////////////////////////////////////////////////////////////
/** Constructor
 */

rcslw::rcslw(const int    p_nGG,
             const double p_Tg,
             const double p_Nconc,
             const double p_Yco2,
             const double p_Yco,
             const double p_Yh2o){

    nGG   = p_nGG;
    Tg    = p_Tg;
    Nconc = p_Nconc;
    Yco2  = p_Yco2;
    Yco   = p_Yco;
    Yh2o  = p_Yh2o;

    P_table    = vector<double>(10);
    C_table    = vector<double>(71);

    for(int i=0; i<71; i++) 
        C_table[i] = 0.0001*pow(1000/0.0001, i/70.0);

    P_table[0] = 0.1;
    P_table[1] = 0.25;
    P_table[2] = 0.5;
    P_table[3] = 1;
    P_table[4] = 2;
    P_table[5] = 4;
    P_table[6] = 8;
    P_table[7] = 15;
    P_table[8] = 30;
    P_table[9] = 50;


}

////////////////////////////////////////////////////////////////////////////////
/** User interface function. 
 *  Given the gas state, set the k and a vectors.
 *  These can then be accessed by the user.
 */

void rcslw::set_k_a(){

    vector<double> C(nGG);
    vector<double> Ct(nGGa);

    for(int j=0; j<nGGA; j++)
        C[j] = get_FI_albdf(F_pts[j], Tg, Tref, Yco2, Yco, Yh2o)

}


