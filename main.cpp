#include <iostream>
#include "myexception.h"
#include "matrix.h"

using namespace std;

int main(int argc, char *argv[])
{

    try
    {
        double c;
        Matrix<double> mat(4, 4), L(4, 4), U(4, 4);

        for(int i=0; i<mat.rows(); ++i)
        {
            for(int j=0; j<mat.cols(); ++j)
            {
                cin >> c;
                mat(i, j) = c;
            }
        }
       if(LU<double>(mat, L, U))
       {
           cout << endl;
           cout << L.toString() << endl;
           cout << U.toString() << endl;
       }
       else
       {
           cout << "Matrix is singular!" << endl;
       }
    }
    catch(const MyException& e)
    {
        cout << e.what() << endl;
    }

    return 0;
}
