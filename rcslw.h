/**
 * @file rcslw.h
 * Header file for class rcslw
 */

#pragma once

#include "multi-D.h"
//#include "leggauss_data.cc"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector> 

using namespace std;
////////////////////////////////////////////////////////////////////////////////

/** Class implementing rcslw object
 *
 *  @author ignite.byu.edu
 *  todo: add some description and references, etc.
 */

class rcslw {


    //////////////////// DATA MEMBERS //////////////////////

    public:

        double P;         ///<   pressure (atm)
        double Tref;      ///<   reference temperature (Tb in Falbdf) (K)
        int    nGG;       ///<   number of gray gases, not including the clear gas
        double Cmin;
        double Cmax;
        int    nGGa;      ///<   number of grey gases including the clear gas


        std::vector<double> P_table;
        std::vector<double> C_table;
        std::vector<double> Tg_table;
        std::vector<double> Tb_table;
        std::vector<double> Yh2o_table;

        std::vector<double> F_pts;
        std::vector<double> Ft_pts;

        std::vector<vector<vector<double> > >          Falbdf_co2;
        std::vector<vector<vector<double> > >          Falbdf_co;
        std::vector<vector<vector<vector<double> > > > Falbdf_h2o;

        int nP;   
        int nC;
        int nTg;
        int nTb;
        int ny_h2o;

        double Fmin;
        double Fmax;

        double co;
        double co2;
        double h2o;
         
        // not needed: double Tg;
        // not needed: double Tb;
        // not needed: double Nconc;
        // not needed: double Yco2; 
        // not needed: double Yco; 
        // not needed: double Yh2o;

        // not needed: std::vector<double> k; 
        // not needed: std::vector<double> a; 

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        //void rcslw(const double p_P, const int p+)
        void get_k_a(double Tg, double Yco2, double Yco, double Yh2o, double fvsoot,
                     vector<double> &k, vector<double> &a);
        // not needed: vector<double> get_k();
        // not needed: vector<double> get_a();
        double get_F_albdf(double C, double Tg, double Tb, double Yco2, double Yco, double Yh2o,
                         double fvsoot);
        double get_FI_albdf(double F, double Tg, double Tb, double Yco2, double Yco, double Yh2o,
                         double fvsoot);
        void get_FI_albdf_tables(string P_table, int nx, int ny, int nz,
                vector<vector<vector<double> > > &myarray);
        void get_FI_albdf_tables(string P_table, int nx, int ny, int nz, int nw,
                vector<vector<vector<vector<double> > > > &myarray);
        void set_Fpts(void);
        void set_Falbdf_co2_co_h2o_at_P(void);
        double F_albdf_soot(double c, double Tg, double Tb, double fvsoot);
        

    private:

                     
                     
                     


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////
    
    public: 
    
        rcslw(const double p_P,
              const int    p_nGG,
              const double Tg, 
              const double Yco2, 
              const double Yco, 
              const double Yh2o,
              const double fvsoot);
        //~rcslw(){}


};

