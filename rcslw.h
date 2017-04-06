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

        int    nGG;      ///<   number of gray gases, not including the clear gas

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

    private:

                     
                     
                     


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////
    
    public: 
    
        rcslw(const double p_Tg, 
              const double p_Nconc, 
              const double p_Yco2, 
              const double p_Yco, 
              const double p_Yh2o);
        rcslw(){}
        ~rcslw(){}


};

