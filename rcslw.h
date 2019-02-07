#ifndef RCSLW_H
#define RCSLW_H
#include <vector>


class rcslw {
	/*
	The Rank Correlated SLW model.
	V.P.Solovjov, F.Andre, D.Lemonnier, B.W.Webb
	http ://www.sciencedirect.com/science/article/pii/S0022407316306434
	@author David O.Lignell
	Note : This was verified by comparing to Solovjov's code.
	The code here is the C++ modified code of David O. Lignell's python code.
	follows the same approach.
	*/
public:
	void setData(double Pressure = 0, int nGG = 0, double Tg = 0, double Yco2 = 0, double Yco = 0, double Yh2o = 0, double fvsoot = 0);
	void set_Falbdf_co2_co_h2o_at_P();
	void set_interpolating_functions();
	void get_k_a(double Tg, double Yco2, double Yco, double Yh2o, double fvsoot, vector<double>& k, vector<double>& a);
	double get_F_albdf(double C, double Tg, double Tb, double Yco2, double Yco, double Yh2o, double fvsoot);
	void get_FI_albdf();
	void get_FI_albdf_old();
	void set_Fpts();
	double F_albdf_soot(double c, double Tg, double Tb, double fvsoot);
	void bisection();

private:
	double P;
	int nGG;
	double Tg;
	double Yco2;
	double Yco;
	double Yh2o;
	double fvsoot;
	int nGGa;
	double Cmin;
	double Cmax;
	double Tref;
	vector<double> Tg_table;
	vector<double> Tb_table;
	vector<double> C_table;		// Table of C values
	vector<double> Yh2o_table;
	vector<double> P_table;		// atm
	vector<string> Pstr;
	int numYh2o;
	double numTg;
	double numTb;
	int numC;
	int numP;
	double Fmin;
	double Fmax;
	vector<double> Ft_pts;
	vector<double> F_pts;
};

#endif
