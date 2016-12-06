
// wordMap.cpp
// inplementation of map methods

#include "dataImport.h"
#include <ctype.h>
#include <sstream>
using namespace std;

// will generate an datamatrix inside
void RawDataSet::readDataFromFile(const string filename, const int instances)
{
    ifstream ifs (filename);
    

    if (ifs.is_open()) {
        
        rawDataTable.clear();
        
        string oneLine;
        getline(ifs, oneLine);
        while (ifs.good())
        {
            for(int i=0; i < oneLine.size(); ++i)
            {
                if(oneLine[i] == ',') oneLine[i] = ' ';
            }
            
            stringstream oneLineStream(oneLine);
            vector<double> oneDataInstance;


            while(oneLineStream.good())
            {
                double oneNumber;
                oneLineStream >> oneNumber;
                oneDataInstance.push_back(oneNumber);
            }

            oneDataInstance.pop_back();


            rawDataTable.push_back(oneDataInstance);
            if(rawDataTable.size() >= instances)
                break;

            getline(ifs, oneLine);
        }

        if(rawDataTable.size() != instances)
            rawDataTable.pop_back();
    }
    else {
        // show message:
        std::cout << "Error opening file" << endl;
    }
}
    
// show the current dataset size in console
void RawDataSet::printDataMatrixSize() const
{
    if(rawDataTable.empty())
    {
        cout << "Current data matrix is empty." << endl;
        return;
    }

    cout << "Current data matrix size: (n, d) = {" 
        << rawDataTable.size() << ", " 
        << rawDataTable[0].size() << "}." << endl;
}

// show the current data metrix in console
// WARNING: can be extra large
// try to use:
// ./program > outputfile.txt
// to redirect the output to a file
void RawDataSet::printDataMatrixContent() const
{
    printDataMatrixSize();

    for (int i = 0; i < rawDataTable.size(); i++)
    {
        for (int j = 0; j < rawDataTable[i].size(); j++)
        {
            cout << rawDataTable[i][j] << " ";
        }
        cout << endl;
    }
}



