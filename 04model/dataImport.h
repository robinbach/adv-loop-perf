// eecs545 final project
// data improt header
// dataImport.h

// import data set from file
// the data set is float numbers

#ifndef dataImport_H
#define dataImport_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include <cassert>

using namespace std;

class RawDataSet
{

public:
    // the array stores the dictionary strings
	vector< vector<double> > rawDataTable;
    
    // the data file has to be within the same folder of binary executable file.
    // will generate an datamatrix inside
    // Optional second parameter: choose how many instances to read.
    void readDataFromFile(string filename, const int instances = 1e6);

    // show the current dataset size in console
    void printDataMatrixSize() const;

    // show the current data metrix in console
    // WARNING: can be extra large
    // try to use:
    // ./program > outputfile.txt
    // to redirect the output to a file
    void printDataMatrixContent() const;

};


#endif