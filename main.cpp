
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <array>
#include <stdio.h>
#include <algorithm>

/**
 * read data from file
 */

int rowA = 0;
int colA = 0;
using namespace std;
int main ()
{
    string lineA, lineB;
    int x, xx, yy;
    int arrayA[500000][2] = {0};
    int arrayB[500000] = {0};

    //int **tmp;
    //tmp = new int*[200000];

    // Open the input files
    fstream  file_ext, file_num;
    file_ext.open("extents.txt", ios::in);

    file_num.open("numbers.txt", ios_base::in);

    cout << "after opening filename" << endl;

    // Check for error
    if (file_ext.fail() && file_num.fail()){
        cerr << "Error opening file" << endl;
        exit(1);
    }

    if (file_ext.good()){
        while(getline(file_ext, lineA)){
            istringstream streamA(lineA);
            colA=0;
            while(streamA >> x) {
                arrayA[rowA][colA] = x;
                colA++;
            }
            rowA++;
        }
    }

    //Read a file until you have reached the end
    int count = 0;
    while (file_num >> xx) {
        arrayB[count] = xx;
        count ++;
    }

    if (false) {
        for (int i=0; i<rowA; i++){
            for (int j=0; j<colA; j++){
                cout << left << setw(3) << arrayA[i][j] << " " ;
            }
            cout << endl;
        }
        for (int j=0; j<count; j++){
            cout << left << setw(2) << arrayB[j] << " " ;
        }
        cout << endl;
    }

    //Now we call the sort function

    //sort(arrayB, arrayB + count);
    //cout << "Sorted Array looks like this." << endl;
    //for (size_t i = 0; i != count; ++i)
    //    cout << arrayB[i] << " ";

    cout << "rowA:" << rowA << endl;
    cout << "colA:" << colA << endl;
    cout << "Number:" <<  count << endl;
    ofstream outfile;// declaration of file pointer named outfile
    outfile.open("output.txt", ios::out); // opens file named "filename" for output
    for ( int j = 0; j < 10; j++ ) //loop numbers
    {
        int output=0;
        for ( int i = 0; i < 10; i++ ) // loop intervals
        {
            if (arrayB[j] <= arrayA[i][1] && arrayB[j]>=arrayA[i][0] )
            {
                output++;
            }
        }
        arrayA[j][0]=output;
        if (j >= count-10) {
            cout << "This number " << arrayB[j] << " appears " << arrayA[j][0] << " times in the list of ranges" << endl;
        }
        outfile << setw(2) << arrayA[j][0] << " ";

    }

    outfile.close();// closes

    cout << " Done " << endl;
    
    // close opning file
    
    file_ext.close();
    file_num.close();
    
    return 0;
}
