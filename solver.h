#ifndef SOLVER
#define SOLVER

#include "vector.h"
#include "matrix.h"

template<typename T>
Vector<T> solver(const Matrix<T>& mat, const Vector<T>& vec)
{
    if(mat.rows() != vec.getDim())
        throw MyException("Error solve linear system: Dimension not match!");

    int n = mat.rows();
    Matrix<T> L(n, n), U(n, n), P(n, n);
    Vector<T> res(n), tmp(n);

    LUP<T>(mat, L, U, P);
    tmp = P * vec;

    //LUx=b => Ld=b
    res(0) = tmp(0) / L(0, 0);
    for(int i=1; i<n; ++i)
    {
        res(i) = tmp(i);
        for(int j=0; j<i; ++j)
        {
            res(i) -= L(i, j) * res(j);
        }
        res(i) /= L(i, i);
    }

    res(n - 1) /= U(n - 1, n - 1);
    for(int i=n-2; i>=0; --i)
    {
        for(int j=n-1; j>i; --j)
        {
            res(i) -= U(i, j) * res(j);
        }
        res(i) /= U(i, i);
    }

    return res;
}

#endif // SOLVER

