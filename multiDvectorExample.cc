#include <iostream>
#include <vector>

using namespace std;

int main() {

    int nx = 5;
    int ny = 4;
    int nz = 3;

    vector<vector<vector<double> > > w(nx, vector<vector<double> >(ny, vector<double>(nz, 0.0)));

    vector<vector<vector<double> > > v(nx);
    for(int i=0; i<nx; i++) {
        v[i].resize(ny);
        for(int j=0; j<ny; j++) {
            v[j].resize(nz);
            for(int k=0; k<nz; k++) {
                v[k] = 0.0;
            }
        }
    }
        
    return 0;

}

