#include <iostream>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

void read_table_data(const string fname, vector<double> &data);

int main() {

    vector<double> data(100);
      
    //////////////////////////////////////////

    ifstream ifile("data.dat"); 

    for (int i=0; i<100; i++)
        ifile >> data[i];

    for (int i=0; i<100; i++)
        cout << endl << data[i];

    ifile.close();


    //////////////////////////////////////////

    fstream file;

    file.open("data.dat", fstream::in);

    int i=0;
    while(!file.eof())
        file >> data[i++];

    for (int i=0; i<100; i++)
        cout << endl << "\t" << data[i];

    file.close();
    
    //////////////////////////////////////////
    
    vector<double> mydata;
    read_table_data("data.dat", mydata);

    for (int i=0; i<mydata.size(); i++)
        cout << endl << "\t\t" << mydata[i];
    
    return 0;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

void read_table_data(const string fname, vector<double> &data){

    fstream file;
    file.open(fname.c_str(), fstream::in);

    double dmb;
    while(!file.eof()) {
        file >> dmb;
        data.push_back(dmb);
    }

    file.close();

}
