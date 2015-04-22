#ifndef SOLVER
#define SOLVER

#include "vector.h"
#include "matrix.h"


template<typename T>
Matrix<T> solver(const Matrix<T>& mat, const Matrix<T>& vec)
{
    if(mat.rows() != vec.rows())
        throw MyException("Error solve linear system: Dimension not match!");

    int n = vec.rows();
    Matrix<T> L(n, n), U(n, n), P(n, n);
    Matrix<T> res(n, 1), tmp(n, 1);

    LUP<T>(mat, L, U, P);
    tmp = P * vec;

    //LUx=b => Ld=b
    res(0, 0) = tmp(0, 0) / L(0, 0);
    for(int i=1; i<n; ++i)
    {
        res(i, 0) = tmp(i, 0);
        for(int j=0; j<i; ++j)
        {
            res(i, 0) -= L(i, j) * res(j, 0);
        }
        res(i, 0) /= L(i, i);
    }

    res(n - 1) /= U(n - 1, n - 1);
    for(int i=n-2; i>=0; --i)
    {
        for(int j=n-1; j>i; --j)
        {
            res(i, 0) -= U(i, j) * res(j, 0);
        }
        res(i, 0) /= U(i, i);
    }

    return res;
}

template<typename T>
Matrix<T> solverLS(const Matrix<T>& mat, const Matrix<T>& vec)
{
    Matrix<T> AT = Transpos<T>(mat);
    Matrix<T> res = Inv<T>(AT * mat) * (AT * vec);

    return res;
}

#endif // SOLVER

