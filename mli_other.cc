#include <vector>
#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int get_location_dim(const vector<double> &x, const double xP) {
    int ihi;
    if(xP <= x[0])
        ihi = 1;
    else if(xP >= x.back())
        ihi = x.size()-1;
    else {
        vector<const double>::iterator itHi = lower_bound(x.begin(), x.end(), xP); // lower_bound gives values >= xP
        ihi = itHi - x.begin();
    }
    return ihi - 1;
}

////////////////////////////////////////////////////////////////////////////////

double interp3D(const vector<double> &x,
                const vector<double> &y,
                const vector<double> &z,
                const vector<vector<vector<double>>> &f,
                const double xP,
                const double yP,
                const double zP){

    constexpr int ndim = 3;        // number of dimensions
    constexpr int nptbx = 8;       // number of points in the cubical bounding box

    int i = get_location_dim(x,xP);
    int j = get_location_dim(y,yP);
    int k = get_location_dim(z,zP);

    int ip = i+1;
    int jp = j+1;
    int kp = k+1;

    double xWta = (xP - x[i])/(x[i+1]-x[i]);
    double yWta = (yP - y[j])/(y[j+1]-y[j]);
    double zWta = (zP - z[k])/(z[k+1]-z[k]);

    double xWtb = 1-xWt;
    double yWtb = 1-yWt;
    double zWtb = 1-zWt;

    return f[i ,j ,k ]*xWtb*yWtb*zWtb +
           f[ip,j ,k ]*xWta*yWtb*zWtb +
           f[i ,jp,k ]*xWtb*yWta*zWtb +
           f[ip,jp,k ]*xWta*yWta*zWtb +
           f[i ,j ,kp]*xWtb*yWtb*zWta +
           f[ip,j ,kp]*xWta*yWtb*zWta +    
           f[i ,jp,kp]*xWtb*yWta*zWta +
           f[ip,jp,kp]*xWta*yWta*zWta;


}


////////////////////////////////////////////////////////////////////////////////

double interp3D(const vector<double> &x,
                const vector<double> &y,
                const vector<double> &z,
                const vector<vector<vector<double>>> &f,
                const double xP,
                const double yP,
                const double zP){

    constexpr int ndim = 3;        // number of dimensions
    constexpr int nptbx = 8;       // number of points in the cubical bounding box

    vector<double> indLocs(ndim);
    indLocs[0] = get_location_dim(x, xP);
    indLocs[1] = get_location_dim(y, yP);
    indLocs[2] = get_location_dim(z, zP);

    vector<double> dimWts(ndim);
    dimWts[0] = (xP - x[indLocs[0]])/(x[indLocs[0]+1]-x[indLocs[0]]);
    dimWts[1] = (yP - y[indLocs[1]])/(y[indLocs[1]+1]-y[indLocs[1]]);
    dimWts[2] = (zP - z[indLocs[2]])/(z[indLocs[2]+1]-z[indLocs[2]]);

    vector<double> oneMdimWts(ndim);
    oneMdimWts[0] = 1.0 - dimWts[0];
    oneMdimWts[1] = 1.0 - dimWts[1];
    oneMdimWts[2] = 1.0 - dimWts[2];

    int indx;
    double fP = 0.0;
    for(int ipt=0; ipt<nptbx; ++ipt){
        ptWt = 1.0;
        for(int jdim=0; jdim<ndim; ++jdim)
            ptWt *= (ipt & (1<<jdim)) ? dimWts[jdim] : oneMdimWts[jdim];
        indx = get_fpt_indx(ipt, indLocs); 
        fP += ptWt*f[]
    }
    


}

