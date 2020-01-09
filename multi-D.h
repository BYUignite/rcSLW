// COMPILE AS g++ -std=c++11 multi-D-interpolation.cc
/**
 * A multi-dimensional interpolator
 * Interpolates up to 5-D
 */ 

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

double LI_1D(vector<double> x, vector<double> f, double xi){
    double dx = x[1] - x[0];
    int ilo = int(xi/dx);
    if (ilo < 0){
        ilo = 0;
    }
    else if (ilo >= x.size() - 1){
        ilo = x.size() - 2;
    }
    int ihi = ilo + 1;

    return (f[ilo] + (xi - x[ilo]) * (f[ihi]-f[ilo])/(x[ihi]-x[ilo]));

}


double LI_2D(vector<double> x, vector<double> y, vector<vector<double> > f, double xP, double yP){
    double dx = x[1] - x[0];

    int ilo = int(xP/dx);
    if (ilo < 0){
        ilo = 0;
    }

    else if (ilo >= x.size() - 1){
        ilo = x.size() - 2;
    }

    int ihi = ilo + 1;

    ////////////// interpolate grid to yP ////////////////////////

    double flo = LI_1D(y, f[ilo], yP);
    double fhi = LI_1D(y, f[ihi], yP);

    ////////////// interpolate final x direction /////////////////

    return LI_1D( vector<double> {x[ilo], x[ihi]}, vector<double> {flo, fhi}, xP);
}


double LI_3D(vector<double> x, vector<double> y, vector<double> z, vector<vector<vector<double> > > f, double xP,
        double yP, double zP){
    double dx = x[1] - x[0];

    int ilo = int(xP/dx);
    if (ilo < 0){
        ilo = 0;
    }
    else if (ilo >= x.size() - 1){
        ilo = x.size() - 2;
    }
    int ihi = ilo + 1;


    ///////////////////// interpolate grid to yP, zP /////////////////

    double flo = LI_2D(y, z, f[ilo], yP, zP);
    double fhi = LI_2D(y, z, f[ihi], yP, zP);

    //////////////////// interpolate final x direction //////////////

    return LI_1D(vector<double> {x[ilo], x[ihi]}, vector<double> {flo, fhi}, xP);

}


double LI_4D(vector<double> x, vector<double> y, vector<double> z, vector<double> w, vector<vector<vector<vector<double> > > > f,
        double xP, double yP, double zP, double wP){
    double dx = x[1] - x[0];

    int ilo = int(xP/dx);
    if (ilo <= 0){
        ilo = 0;
    }
    else if (ilo >= x.size() - 1){
        ilo = x.size() - 2;
    }
    int ihi = ilo +1;

    ////////////////// interpolate grid to yP, zP ////////////////////

    double flo = LI_3D(y, z, w, f[ilo], yP, zP, wP);
    double fhi = LI_3D(y, z, w, f[ihi], yP, zP, wP);
    
    //////////////// interpolate final x direction //////////////////

    return LI_1D( vector<double> {x[ilo], x[ihi]}, vector<double> {flo, fhi}, xP);

    }


double LI_5D(vector<double> x, vector<double> y, vector<double> z, vector<double> w, vector<double> a,
        vector<vector<vector<vector<vector<double> > > > > f, double xP, double yP, double zP, double wP, double aP){
    double dx = x[1] - x[0];

    int ilo = int(xP/dx);
    if (ilo <= 0){
        ilo = 0;
    }
    else if (ilo >= x.size() - 1){
        ilo = x.size() - 2;
    }
    int ihi = ilo +1;

    ////////////////// interpolate grid to yP, zP ////////////////////

    double flo = LI_4D(y, z, w, a,  f[ilo], yP, zP, wP, aP);
    double fhi = LI_4D(y, z, w, a,  f[ihi], yP, zP, wP, aP);
    
    //////////////// interpolate final x direction //////////////////

    return LI_1D( vector<double> {x[ilo], x[ihi]}, vector<double> {flo, fhi}, xP);

    }


/*int main(){
    int Nx = 8;
    int Ny = 6;

    vector<double> x(Nx);
    vector<double> y(Ny);
    vector<vector<double> > f(Nx, vector<double>(Ny));



    for (int i=0; i<Nx; i++){
        x[i] = (i)*10.0/(Nx-1);
        cout << x[i] << " ";
    }
    cout << "\n\n";


    for(int i=0; i<Ny; i++){
        y[i] = (i)*5.0/(Ny-1);
        cout << y[i] << " ";
    }
    cout << "\n\n";

    for (int i=0; i<Nx; i++){
        for (int j=0; j<Ny; j++){
            f[i][j] = sin(x[i])*cos(y[j]);
            cout << setw(12) << f[i][j] << " ";
        }
        cout << "\n";
    }
    
    double interp = LI_2D(x,y,f,0.5,0.2);
    cout << interp << "\n";

    return(0);

}*/


