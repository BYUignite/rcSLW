/**
 * @file rcslw.h
 * Header file for class rcslw
 */

#ifndef RCSLW_H
#define RCSLW_H

#include <vector>


////////////////////////////////////////////////////////////////////////////////

/** Class implementing rcslw object
 *
 *  @author ignite.byu.edu
 *  todo: add some description and references, etc.
 */

class rcslw {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        int    nGG;               ///<   number of gray gases, not including the clear gas
        double P;                 ///<   pressure (atm)
        double Tref = 1000.0;     ///<   reference temperature (Tb in Falbdf) (K)
        double Cmin = 0.0001;
        double Cmax = 100.0;
        
        
        

        int    nGGa = nGG+1;       ///<   number of grey gases including the clear gas
        
        std::vector<double> P_table;
        std::vector<double> C_table;
        std::vector<double> Tg_table;
        std::vector<double> Tb_table;
        std::vector<double> Yh2o_table;

         
        double Tg; 
        double Nconc;
        double Yco2; 
        double Yco; 
        double Yh2o;

   public:

        std::vector<double> k; 
        std::vector<double> a; 

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        void set_k_a();
        /**
         * THIS IS THE CLASS INTERFACE FUNCTION
         * return the local gray gas coefficients (k) and the local weights (a).
         * Tg:    input; float; gas temperature
         * Nconc: input; float; molar concentration: mol/m3
         * Yco2:  input; float; mole fraction co2
         * Yco:   input; float; mole fraction co
         * Yh2o:  input; float; mole fraction h2o
         */

        std::vector<int> C;
        std::vector<int> Ct;

        

        for()



    private:

                     
                     
                     


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////
    
    public: 
    
        rcslw(const int    p_nGG,
              const double p_Tg, 
              const double p_Nconc, 
              const double p_Yco2, 
              const double p_Yco, 
              const double p_Yh2o);
        rcslw(){}
        ~rcslw(){}


};

