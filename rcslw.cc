#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;

#include "rcslw.h"

void rcslw::set_Falbdf_co2_co_h2o_at_P() {
	/*
	Read the albdf table for co2, co, and h2o.
	Separate files are given for various pressures.
	Interpolate the files to the desired pressure: s.P.
	Reshape from the given 1D to the 3D (co/co2) or 4D (h2o) tables.
	*/

	int i1, i2;
	double P1, P2;
	double f = 0.0;

	if (P <= P_table[0] || P >= P_table[P_table.size() - 1]) {
		cout << "Pressure = " << P << " atm is out of range" << endl;
	}
	else if (P == P_table[0]) {
		i1 = 0;
		i2 = 1;
		P1 = P_table[i1];
		P2 = P_table[i2];
	}
	else if (P == P_table[P_table.size() - 1]) {
		i1 = P_table.size() - 2;
		i2 = P_table.size() - 1;
		P1 = P_table[i1];
		P2 = P_table[i2];
	}
	else {
		for (int i = 0; i <P_table.size(); i++) {
			if (P > P_table[i] && P < P_table[i + 1]) {
				i1 = i;
				i2 = i + 1;
				P1 = P_table[i1];
				P2 = P_table[i2];
			}
		}
	}
	f = (P - P1) / (P2 - P1);			// interpolation factor

										// Need to read in ALBDF files for CO2, CO, H2O

	ifstream file_h2o_plo("h2o_" + Pstr[i1] + ".dat");
	ifstream file_h2o_phi("h2o_" + Pstr[i2] + ".dat");
	ifstream file_co2_plo("co2_" + Pstr[i1] + ".dat");
	ifstream file_co2_phi("co2_" + Pstr[i2] + ".dat");
	ifstream file_co_plo("co_" + Pstr[i1] + ".dat");
	ifstream file_co_phi("co_" + Pstr[i2] + ".dat");

	const int ntot = numYh2o*numC*numTg*numTb;
	vector<double> F_h2o_lo(ntot);
	vector<double> F_h2o_hi(ntot);
	vector<double> F_h2o(ntot);

	vector<double> F_co2_lo(ntot);
	vector<double> F_co2_hi(ntot);
	vector<double> F_co2(ntot);

	vector<double> F_co_lo(ntot);
	vector<double> F_co_hi(ntot);
	vector<double> F_co(ntot);

	for (int i = 0; i<ntot; i++) {
		file_h2o_plo >> F_h2o_lo[i];
		file_h2o_phi >> F_h2o_hi[i];
		F_h2o[i] = F_h20_lo[i] + f * (F_h2o_hi - F_h2o_lo);
	}

	for (int i = 0; i<ntot; i++) {
		file_co2_plo >> F_co2_lo[i];
		file_co2_phi >> F_co2_hi[i];
		F_co2[i] = F_co2_lo[i] + f * (F_co2_hi - F_co2_lo);
	}

	for (int i = 0; i<ntot; i++) {
		file_co_plo >> F_co_lo[i];
		file_co_phi >> F_co_hi[i];
		F_co[i] = F_co_lo[i] + f * (F_co_hi - F_co_lo);
	}


	vector<vector<vector<vector<double> > > > F_H2O(numYh2o, vector<vector<vector<double> > >(numTg, vector<vector<double> >(numTb, vector<double>(numC, 0.0))));
	vector<vector<vector<double> > >          F_CO2(numTg, vector<vector<double> >(numTb, vector<double>(numC, 0.0)));
	vector<vector<vector<double> > >          F_CO(numTg, vector<vector<double> >(numTb, vector<double>(numC, 0.0)));

	for (int iY = 0; iY < numYh2o; iY++) {
		for (int iTg = 0; iTg < numTg; iTg++) {
			for (int iTb = 0; iTb < numTb; iTb++) {
				for (int iC = 0; iC < numC; iC++) {
					F_h2o >> F_H2O[iY][iTg][iTb][iC];
				}
			}
		}
	}

	for (int iTg = 0; iTg < numTg; iTg++) {
		for (int iTb = 0; iTb < numTb; iTb++) {
			for (int iC = 0; iC < numC; iC++) {
				F_co2 >> F_CO2[iTg][iTb][iC];
			}
		}
	}

	for (int iTg = 0; iTg < numTg; iTg++) {
		for (int iTb = 0; iTb < numTb; iTb++) {
			for (int iC = 0; iC < numC; iC++) {
				F_co >> F_CO[iTg][iTb][iC];
			}
		}
	}







	file_h2o_plo.close();
	file_h2o_phi.close();
	file_co2_plo.close();
	file_co2_phi.close();
	file_co_plo.close();
	file_co_phi.close();



	cout << "NEEDS WORK" << endl;
	return;
}

void rcslw::set_interpolating_functions() {
	/*
	The table is read using multi-linear interpolation.
	The scipy RegularGridInterpolator returns a function that is then
	called with the desired grid point and returns the ALBDF at
	that point.
	This function just sets those interpolation functions as a dictionary
	over radiating species.
	*/
	cout << "NEEDS WORK" << endl;
	return;
}


void rcslw::setData(double Pressure, int nGG, double Tg, double Yco2, double Yco, double Yh2o, double fvsoot) {
	P = Pressure;
	nGG = nGG;						 // number of grey gases not including the clear gas
	Tg = Tg;
	Yco2 = Yco2;
	Yco = Yco;
	Yh2o = Yh2o;
	fvsoot = fvsoot;
	Tref = 1000.0;					// reference temperature (Tb in Falbdf) (K)
	Cmin = 0.0001;
	Cmax = 100.0;
	numYh2o = 9;					// size of array for Yh2o
	numTg = 28;						// size of array for Temperature Tg
	numTb = 28;						// size of array for Temperature Tb
	numC = 71;						// size of array for C
	numP = 10;						// size of array for Pressure
	double Fmin, Fmax;

	nGGa = nGG + 1;					// number of grey gases including the clear gas


									// Table of the pressures
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

	// Table that identifies the pressure
	Pstr[0] = "p0.1";
	Pstr[1] = "p0.25";
	Pstr[2] = "p0.5";
	Pstr[3] = "p1";
	Pstr[4] = "p2";
	Pstr[5] = "p4";
	Pstr[6] = "p8";
	Pstr[7] = "p15";
	Pstr[8] = "p30";
	Pstr[9] = "p50";


	for (double i = 0; i < numC; ++i) {	// Fill table of C values
		C_table[i] = 0.0001 * pow((1000.0 / 0.0001), (i / 70.0));

	}

	for (int i = 0; i < numTg; ++i) {
		if (i == 0) {
			Tg_table[i] = 300.0;
		}
		else {
			Tg_table[i] = Tg_table[i - 1] + (3000.0 - 300.0) / (numTg - 1);
		}
	}

	for (int i = 0; i < numTb; ++i) {
		if (i == 0) {
			Tb_table[i] = 300.0;
		}
		else {
			Tb_table[i] = Tb_table[i - 1] + (3000.0 - 300.0) / (numTb - 1);
		}
	}

	Yh2o_table[0] = 0.0;
	Yh2o_table[1] = 0.05;
	Yh2o_table[2] = 0.1;
	Yh2o_table[3] = 0.2;
	Yh2o_table[4] = 0.3;
	Yh2o_table[5] = 0.4;
	Yh2o_table[6] = 0.6;
	Yh2o_table[7] = 0.8;
	Yh2o_table[8] = 1.0;

	set_Falbdf_co2_co_h2o_at_P();
	set_interpolating_functions();

	Fmin = get_F_albdf(Cmin, Tg, Tg, Yco2, Yco, Yh2o, fvsoot);
	Fmax = get_F_albdf(Cmax, Tg, Tg, Yco2, Yco, Yh2o, fvsoot);

	cout << "Fmin, Fmax = " << Fmin << " " << Fmax << endl;
	set_Fpts();
	return;
}





void rcslw::get_k_a(double Tg, double Yco2, double Yco, double Yh2o, double fvsoot, vector<double>& k, vector<double>& a) {
	/*
	THIS IS THE CLASS INTERFACE FUNCTION
	return the local gray gas coefficients(k) and the local weights(a).
	Tg:     input; float; gas temperature
	Yco2 : input; float; mole fraction co2
	Yco : input; float; mole fraction co
	Yh2o : input; float; mole fraction h2o
	fvsoot : input; float; soot volume fraction = rho*Ysoot / rhoSoot
	returns k, a arrays
	*/
	vector<double> C(nGG);			// C = C(F, \phi_{loc}, Tref)
	vector<double> Ct(nGGa);		// \tilde{C} = C(\tilde{F}, \phi_{loc}, Tref)
	for (int i = 0; i < nGGa; ++i) {
		Ct[i] = get_FI_albdf(Ft_pts[i], Tg, Tref, Yco2, Yco, Yh2o, fvsoot);
		if (i == 0) {
			cout << "Ct =    [  " << C[i] << "   ";
		}
		else if (i > 0 && i != nGG - 1) {
			cout << Ct[i] << "   ";
		}
		else {
			cout << Ct[i] << "]" << endl;
		}
	}

	for (int i = 0; i < nGG; ++i) {
		C[i] = get_FI_albdf(F_pts[i], Tg, Tref, Yco2, Yco, Yh2o, fvsoot);
		if (i == 0) {
			cout << "C =    [  " << C[i] << "   ";
		}
		else if (i > 0 && i != nGG - 1) {
			cout << C[i] << "   ";
		}
		else {
			cout << C[i] << "]" << endl;
			cout << endl;
		}
	}

	double Nconc = P * 101325 / 8.31445 / Tg;

	k[0] = 0.0;

	for (int i = 0; i < nGG; ++i) {
		k[i + 1] = Nconc * C[i];
	}

	vector<double> FCt(nGGa);       // doldb
	for (int i = 0; i < nGGa; ++i) {
		FCt[i] = get_F_albdf(Ct[i], Tg, Tg, Yco2, Yh2o, fvsoot);

		if (i == 0) {
			cout << "FCt =    [  " << C[i] << "   ";
		}
		else if (i > 0 && i != nGG - 1) {
			cout << FCt[i] << "   ";
		}
		else {
			cout << FCt[i] << "]" << endl;
			cout << endl;
		}
	}

	a[0] = FCt[0];
	for (int i = 1; i < nGGa; ++i) {
		a[i] = FCt[i] - FCt[i - 1];
	}
	cout << "NEEDS WORK" << endl;
	return;
}

double rcslw::get_F_albdf(double C, double Tg, double Tb, double Yco2, double Yco, double Yh2o, double fvsoot) {
	/*
	C:    input; float; cross section
	Tg:   input; float; gas temperature
	Tb:   input; float; black temperature
	Yco2: input; float; mole fraction co2
	Yco:  input; float; mole fraction co
	Yh2o: input; float; mole fraction h2o
	returns the albdf function F
	*/
	double CYco2, CYco, CYh2o;
	double F_co2, F_co, F_h2o;

	if (Yco2 <= pow(10, -20)) {
		Yco2 = pow(10, -20);
	}
	if (Yco <= pow(10, -20)) {
		Yco = pow(10, -20);
	}
	if (Yh2o <= pow(10, -20)) {
		Yh2o = pow(10, -20);
	}
	if (Tg < Tg_table[0]) {
		Tg = Tg_table[0];
	}
	if (Tg > Tg_table[-1]) {
		Tg = Tg_table[-1];
	}
	if (Tb < Tb_table[0]) {
		Tb = Tb_table[0];
	}
	if (Tb > Tb_table[-1]) {
		Tb = Tb_table[-1];
	}
	if (Yh2o < Yh2o_table[0]) {
		Yh2o = Yh2o_table[0];
	}
	if (Yh2o > Yh2o_table[-1]) {
		Yh2o = Yh2o_table[-1];
	}
	CYco2 = C / Yco2;
	CYco = C / Yco;
	CYh2o = C / Yh2o;

	if (CYco2 < C_table[0]) {
		CYco2 = C_table[0];
	}
	if (CYco2 > C_table[-1]) {
		CYco2 = C_table[-1];
	}
	if (CYco < C_table[0]) {
		CYco = C_table[0];
	}
	if (CYco > C_table[-1]) {
		CYco = C_table[-1];
	}
	if (CYh2o < C_table[0]) {
		CYh2o = C_table[0];
	}
	if (CYh2o > C_table[-1]) {
		CYh2o = C_table[-1];
	}


	//	F_co2 = interp_F_albdf;
	//	F_co = interp_F_albdf;
	//	F_h2o = interp_F_albdf;

	cout << "NEEDS WORK" << endl;
	return F_co2 * F_co * F_h2o * F_albdf_soot(C, Tg, Tb, fvsoot);
}

void rcslw::get_FI_albdf() {
	/*
	Inverse F_albdf: pass in F and get out C.
	C:      input; float; cross section
	Tg:     input; float; gas temperature
	Tb:     input; float; black temperature
	Yco2:   input; float; mole fraction co2
	Yco:    input; float; mole fraction co
	Yh2o:   input; float; mole fraction h2o
	fvsoot: input; float; soot volume fraction = rho*Ysoot/rhoSoot
	returns C.
	*/
	cout << "NEEDS WORK" << endl;
	return;
}

void rcslw::get_FI_albdf_old(F, ) {
	/*  LOW PRIORITY
	Inverse F_albdf: pass in F and get out C.
	C:    input; float; cross section
	Tg:   input; float; gas temperature
	Tb:   input; float; black temperature
	Yco2: input; float; mole fraction co2
	Yco:  input; float; mole fraction co
	Yh2o: input; float; mole fraction h2o
	fvsoot: input; float; soot volume fraction = rho*Ysoot/rhoSoot
	returns C.
	*/
	cout << "NEEDS WORK" << endl;
	return;
}

void rcslw::set_Fpts() {
	/*
	Set grid of F points based on Gauss - Legendre Quadrature.
	*/

	// Need to figure out x and w for Gauss-Legendre Quadrature

	Ft_pts[0] = Fmin;
	double w_cum_sum = w[0];
	for (int i = 1; i < nGGa; ++i) {
		w_cum_sum = w_cum_sum + w[i];
		Ft_pts[i] = Fmin + (Fmax - Fmin) * w_cum_sum;
		if (i == 0) {
			cout << "Ft_pts =    [  " << Ft_pts[i] << "   ";
		}
		else if (i > 0 && i != nGG - 1) {
			cout << Ft_pts[i] << "   ";
		}
		else {
			cout << Ft_pts[i] << "]" << endl;
		}
	}
	for (int i = 0; i < nGG; ++i) {
		F_pts[i] = Fmin + x[i] * (Fmax - Fmin);
		if (i == 0) {
			cout << "F_pts =    [  " << F_pts[i] << "   ";
		}
		else if (i > 0 && i != nGG - 1) {
			cout << F_pts[i] << "   ";
		}
		else {
			cout << F_pts[i] << "]" << endl;
			cout << endl;
		}
	}

	cout << "NEEDS WORK" << endl;
	return;
}

double rcslw::F_albdf_soot(double c, double Tg, double Tb, double fvsoot) {
	/*
	C:    input; float; cross section
	Tg : input; float; gas temperature
	Tb : input; float; black temperature
	returns the albdf function F for soot
	Note : This comes from Solovjov 2001 and Chang 1984
	http ://heattransfer.asmedigitalcollection.asme.org/article.aspx?articleid=1444909
	http ://www.sciencedirect.com/science/article/pii/0735193384900514
	*/
	double FunctionF;
	double csoot = 7.0;			  // general soot parameter
	double hCokb = 0.01438777354; // m*K = h*Co/kb = Planck*lightSpeed/Boltzmann (for soot)
	double Nconc, x, n;
	const double PI = 3.141592653589793238463;
	vector<double> n(3);          // sum from n = 1 to oo, but n=1 to 3 is enough
	double sumValue;

	Nconc = P * 101325 / 8.31446 / Tg;   // mol/m3

	x = hCokb * c * Nconc / (csoot * fvsoot * Tb);

	if (fvsoot < pow(1, -12)) {
		return 1.0;
	}
	else {
		for (int i = 0; i < 3; ++i) {
			n = i + 1;
			FunctionF = (exp(-n*x) / pow(n, 4) * (n*x*(n*x*(n*x + 3) + 6) + 6));
			sumValue = sumValue + FunctionF;
		}
	}


	cout << "NEEDS TO BE CHECKED" << endl;
	return 1.0 - 15.0 / pow(PI, 4) * sumValue;
}
void rcslw::bisection() {
	// LOW PRIORITY
	cout << "NEEDS WORK" << endl;
	return;
}

rcslw.cpp
Displaying rcslw.cpp.