#include <iostream>
#include "myexception.h"
#include "matrix.h"

using namespace std;

int main(int argc, char *argv[])
{

    try
    {
        double c;
        Matrix<double> mat(50, 50);

        for(int i=0; i<mat.rows(); ++i)
        {
            for(int j=0; j<mat.cols(); ++j)
            {
                cin >> c;
                mat(i, j) = c;
            }
        }
        c = det<double>(mat);
        cout << c << endl;
    }
    catch(const MyException& e)
    {
        cout << e.what() << endl;
    }

    return 0;
}
