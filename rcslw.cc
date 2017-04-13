/**
 * @file rcslw.cc
 * Source file for class rcslw
 */

#include "rcslw.h"


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
}

////////////////////////////////////////////////////////////////////////////////
/** User interface function. 
 *  Given the gas state, set the k and a vectors.
 *  These can then be accessed by the user.
 */

void rcslw::set_k_a(){


}


