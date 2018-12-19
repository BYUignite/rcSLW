/* Mikel Zaitzeff
Multi Dimensional Interpolator
*/

#include<iostream>
#include <vector>
using namespace std;

double LI_1D(vector <double> &x, vector<double> &f, const double xi) {
	double dx = x.at(1) - x.at(0);
	int ilo = xi / dx;
	double interp_value;

	if (ilo < 0) {
		ilo = 0.0;
	}
	else if (ilo >= x.size() - 1) {
		ilo = x.size() - 2;
	}
	int ihi = ilo + 1;

	interp_value = f[ilo] + (xi - x.at(ilo)) * (f[ihi] - f[ilo]) / (x.at(ihi) - x.at(ilo));
	return interp_value;
}

double LI_2D(vector <double> &x, vector <double> &y, vector <vector <double> > &f, const double xi, const double yi) {
	vector<double> newX(2);
	vector <double> newF(2);
	double dx = x.at(1) - x.at(0);
	int ilo = xi / dx;
	if (ilo < 0) {
		ilo = 0.0;
	}
	else if (ilo >= x.size() - 1) {
		ilo = x.size() - 2;
	}
	int ihi = ilo + 1;
	newX[0] = x.at(ilo);
	newX[1] = x.at(ihi);
	// --------- interpolate grid to yi
	newF[0] = LI_1D(y, f[ilo], yi);
	newF[1] = LI_1D(y, f[ihi], yi);

	// --------- interpolate final x direction
	return LI_1D(newX, newF, xi);
}

double LI_3D(vector <double> &x, vector <double> &y, vector <double> &z, vector<vector<vector<double> > > &f, const double xi, const double yi, const double zi) {
	vector<double> newX(2);
	vector <double> newF(2);
	double dx = x.at(1) - x.at(0);
	int ilo = xi / dx;
	if (ilo < 0) {
		ilo = 0.0;
	}
	else if (ilo >= x.size() - 1) {
		ilo = x.size() - 2;
	}
	int ihi = ilo + 1;
	newX[0] = x.at(ilo);
	newX[1] = x.at(ihi);
	// --------- interpolate grid to yi
	newF[0] = LI_2D(y, z, f[ilo], yi, zi);
	newF[1] = LI_2D(y, z, f[ihi], yi, zi);

	// --------- interpolate final x direction
	return LI_1D(newX, newF, xi);

}

double LI_4D(vector <double> &x, vector <double> &y, vector <double> &z, vector <double> &w, vector<vector<vector<vector<double> > > > &f, const double xi, const double yi, const double zi, const double wi) {

	vector<double> newX(2);
	vector <double> newF(2);
	double dx = x.at(1) - x.at(0);
	int ilo = xi / dx;
	if (ilo < 0) {
		ilo = 0.0;
	}
	else if (ilo >= x.size() - 1) {
		ilo = x.size() - 2;
	}
	int ihi = ilo + 1;
	newX[0] = x.at(ilo);
	newX[1] = x.at(ihi);
	// --------- interpolate grid to yi
	newF[0] = LI_3D(y, z, w, f[ilo], yi, zi, wi);
	newF[1] = LI_3D(y, z, w, f[ihi], yi, zi, wi);

	// --------- interpolate final x direction
	return LI_1D(newX, newF, xi);
}

double LI_5D(vector <double> &x, vector <double> &y, vector <double> &z, vector <double> &w, vector <double> &a, vector<vector<vector<vector<vector<double> > > > > &f, const double xi, const double yi, const double zi, const double wi, const double ai) {

	vector<double> newX(2);
	vector <double> newF(2);
	double dx = x.at(1) - x.at(0);
	int ilo = xi / dx;
	if (ilo < 0) {
		ilo = 0.0;
	}
	else if (ilo >= x.size() - 1) {
		ilo = x.size() - 2;
	}
	int ihi = ilo + 1;
	newX[0] = x.at(ilo);
	newX[1] = x.at(ihi);
	// --------- interpolate grid to yi
	newF[0] = LI_4D(y, z, w, a, f[ilo], yi, zi, wi, ai);
	newF[1] = LI_4D(y, z, w, a, f[ihi], yi, zi, wi, ai);

	// --------- interpolate final x direction
	return LI_1D(newX, newF, xi);
}


int main() {
	int Nx = 8.0;
	int Ny = 6.0;
	int Nz = 3.0;
	int Nw = 4.0;
	int Na = 10.0;
	vector <double> x(Nx);
	vector <double> y(Ny);
	vector <double> z(Nz);
	vector <double> w(Nw);
	vector <double> a(Na);
	vector <vector<vector<vector<vector<double>>>>> f(Nx, vector<vector<vector<vector<double> > > >(Ny, vector<vector<vector<double> > >(Nz, vector<vector<double> >(Nw,vector<double> (Na,0.0)))));
	vector <double> value(Ny);
	double test;

	double dx = 10.0 / (Nx - 1.0);
	double dy = 5.0 / (Ny - 1.0);
	double dz = 7.0 / (Nz - 1.0);
	double dw = 8.0 / (Nw - 1.0);
	double da = 2.0 / (Na - 1.0);
	for (int i = 0; i < Nx; i++) {
		x.at(i) = i * dx;
	}
	for (int i = 0; i < Ny; i++) {
		y.at(i) = i * dy;
	}
	for (int i = 0; i < Nz; i++) {
		z.at(i) = i * dz;
	}
	for (int i = 0; i < Nw; i++) {
		w.at(i) = i * dw;
	}
	for (int i = 0; i < Na; i++) {
		a.at(i) = i * da;
	}

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int l = 0; l < Nz; l++) {
				for (int k = 0; k < Nw; k++) {
					for (int p = 0; p < Na; p++) {
						f[i][j][l][k][p] = sin(x.at(i))*cos(y.at(j))*sin(z.at(l))*cos(w.at(k))*sin(a.at(p));
					}
				}
			}
		}
	}
	test = LI_5D(x, y, z, w, a, f, 2, 4, 5, 7, 1);
	cout << test << endl;
	system("pause");
	return 0;
}