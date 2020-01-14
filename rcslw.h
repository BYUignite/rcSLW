/**
 * @file rcslw.h
 * Header file for class rcslw
 */

#pragma once

//#include "leggauss_data.cc"

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

        int    nGG;       ///<   number of gray gases, not including the clear gas
        int    nGGa;      ///<   number of grey gases including the clear gas
        double P;         ///<   pressure (atm)
        double Tref;      ///<   reference temperature (Tb in Falbdf) (K)
        double Cmin;
        double Cmax;


        vector<double> P_table;
        vector<double> C_table;
        vector<double> Tg_table;
        vector<double> Tb_table;
        vector<double> Yh2o_table;

        vector<double> F_pts;
        vector<double> Ft_pts;

        vector<vector<vector<double> > >          Falbdf_co2;
        vector<vector<vector<double> > >          Falbdf_co;
        vector<vector<vector<vector<double> > > > Falbdf_h2o;

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
         
    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        void get_k_a(double Tg, double Yco2, double Yco, double Yh2o, double fvsoot,
                     vector<double> &k, vector<double> &a);

    private:

        double get_F_albdf(double C, double Tg, double Tb, double Yco2, double Yco, double Yh2o,
                         double fvsoot);
        double get_FI_albdf(double F, double Tg, double Tb, double Yco2, double Yco, double Yh2o,
                         double fvsoot);
        void get_FI_albdf_tables(std::string Ptable_file_name, int nx, int ny, int nz,
                vector<vector<vector<double> > > &myarray);
        void get_FI_albdf_tables(string Ptable_file_name, int nx, int ny, int nz, int nw,
                vector<vector<vector<vector<double> > > > &myarray);
        void set_Fpts(void);
        void set_Falbdf_co2_co_h2o_at_P();
        double F_albdf_soot(double c, double Tg, double Tb, double fvsoot);
        
    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////
    
    public: 
    
        rcslw(const int    p_nGG,
              const double p_P,
              const double Tg, 
              const double Yco2, 
              const double Yco, 
              const double Yh2o,
              const double fvsoot);


};

