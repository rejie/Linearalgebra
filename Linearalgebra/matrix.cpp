#include "matrix.h"

template<typename T>
Matrix<T>::Matrix(int row, int col, const T& init_var)
{

    if(row > 0 && col > 0)
    {
        data.resize(row);
        mdata.resize(row * col, init_var);
        for(int i=0; i<row; ++i)
        {
            data[i] = &mdata[i * col];
        }

        _row = row;
        _col = col;
    }
    else
        throw MyException("Matrix Initialize Error!");



}

template<typename T>
Matrix<T>::Matrix(int row, int col)
{

    if(row > 0 && col > 0)
    {
        data.resize(row);
        mdata.resize(row * col, 0);
        for(int i=0; i<row; ++i)
        {
            data[i] = &mdata[i * col];
        }

        _row = row;
        _col = col;
    }
    else
        throw MyException("Matrix Initialize Error!");
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& mat)
{
    /*
    for(int i=0; i<_row; ++i)
    {
        for(int j=0; j<_col; ++j)
        {
            data[i][j] = mat(i, j);
        }
    }
    _row = mat.rows();
    _col = mat.cols();
    */
    *this = mat;
}

template<typename T>
Matrix<T>::~Matrix(){}

//get data=================================================================================
template<typename T>
void Matrix<T>::resize(int new_row, int new_col)
{
    if(new_row > 0 && new_col > 0)
    {
        std::vector<T> tmp(mdata.begin(), mdata.end());

        data.resize(new_row);
        mdata.resize(new_row * new_col);
        for(int i=0; i<new_row; ++i)
        {
            data[i] = &mdata[i * new_col];
        }

        for(int i=0; i<_row && tmp.size()>0; ++i)
        {
            for(int j=0; j<_col; ++j)
            {
                mdata[i*new_col + j] = tmp[i*_col +j];
            }
        }

        _row = new_row;
        _col = new_col;
    }
    else
        throw MyException("Matrix size error!");
}

template<typename T>
void Matrix<T>::swapRows(int ri, int rj)
{
    if(!checkr(ri) || !checkr(rj))
        throw MyException("Row index out of range!");

    if(ri != rj)
    {
        T* tmp = data[ri];
        data[ri] = data[rj];
        data[rj] = tmp;
    }
}

template<typename T>
void Matrix<T>::swapCols(int ci, int cj)
{
    if(!checkc(ci) || !checkc(cj))
        throw MyException("Column index out of range!");

    if(ci != cj)
    {
        T tmp;

        for(int i=0; i<_row; ++i)
        {
            tmp = data[i][ci];
            data[i][ci] = data[i][cj];
            data[i][cj] = tmp;
        }
    }
}

template<typename T>
std::string Matrix<T>::toString()
{
    std::stringstream ss;

    for(int i=0; i<_row; ++i)
    {
        for(int j=0; j<_col; ++j)
        {
            ss << std::fixed << data[i][j] << " ";
        }
        ss << "\n";
    }

    return ss.str();
}

//overload operator========================================================================
template<typename T>
T& Matrix<T>::operator()(int row, int col)
{
    if(check(row, col))
        return data[row][col];
    else
        throw MyException("Index out of range!");
}

template<typename T>
const T& Matrix<T>::operator()(int row, int col) const
{
    if(check(row, col))
        return data[row][col];
    else
        throw MyException("Index out of range!");
}


template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& mat)
{
    if(&mat == this)
        return *this;

     int new_row = mat.rows(),
                 new_col = mat.cols();

    resize(new_row, new_col);
    for(int i=0; i<new_row; ++i)
    {
        for(int j=0; j<new_col; ++j)
        {
            data[i][j] = mat(i, j);
        }
    }

    _row = new_row;
    _col = new_col;

    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& mat)
{
    Matrix<T> res(_row, _col);

    if( _row == mat.rows() && _col == mat.cols() )
    {
        for(int i=0; i<_row; ++i)
        {
            for(int j=0; j<_col; ++j)
            {
               res(i, j) = data[i][j] + mat(i, j);
            }
        }
    }
    else
        throw MyException("Matrix addition error: dimension not match!");

    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const T& var)
{
    Matrix<T> res(_row, _col);

    for(int i=0; i<_row; ++i)
    {
        for(int j=0; j<_col; ++j)
        {
           res(i, j) = data[i][j] + var;
        }
    }

    return res;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& mat)
{
    return *this = *this + mat;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const T& var)
{
    return *this = *this + var;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-()
{
    for(int i=0; i<_row; ++i)
    {
        for(int j=0; j<_col; ++j)
        {
            data[i][j] = -data[i][j];
        }
    }

    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& mat)
{
    Matrix<T> tmp_mat = mat;
    return *this + (-tmp_mat);
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const T& var)
{
    T tmp_var = var;
    return *this + (-tmp_var);
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& mat)
{
    Matrix<T> tmp_mat = mat;
    return *this = *this + (-tmp_mat);
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const T& var)
{
    T tmp_var = var;
    return *this = *this + (-tmp_var);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& mat)
{
    Matrix<T> res(_row, mat.cols());

    if(_col == mat.rows())
    {
        for(int i=0; i<_row; ++i)
        {
            for(int j=0; j<mat.cols(); ++j)
            {
                for(int k=0; k<_col; ++k)
                {
                    res(i, j) += data[i][k] * mat(k, j);
                }
            }
        }
    }
    else
        throw MyException("Matrix multiplication error: dimension not match!");

    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& var)
{
    Matrix<T> res(_row, _col);

    for(int i=0; i<_row; ++i)
    {
        for(int j=0; j<_col; ++j)
        {
            res(i, j) = data[i][j] * var;
        }
    }

    return res;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& mat)
{
    return *this = *this * mat;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T& var)
{
    return *this = *this * var;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(const T& var)
{
    Matrix res(_row, _col);

    if(var != 0)
    {
        for(int i=0; i<_row; ++i)
        {
            for(int j=0; j<_col; ++j)
            {
                res(i, j) = data[i][j] / var;
            }
        }
    }
    else
        throw MyException("Matrix division error: divide by zero!");

    return res;
}

template<typename T>
Matrix<T>& Matrix<T>::operator/=(const T& var)
{
    return *this = *this / var;
}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T>& mat)
{
    if(_row != mat.rows() || _col != mat.cols())
        return false;

    for(int i=0; i<_row; ++i)
    {
        for(int j=0; j<_col; ++j)
        {
            if(data[i][j] != mat(i, j))
                return false;
        }
    }

    return true;
}

//matrix operation=========================================================================
template<typename T>
Matrix<T> Identity(int dim)
{
    if(dim > 0)
    {
        Matrix<T> res(dim, dim, 0);

        for(int i=0; i<dim; ++i)
        {
            res(i, i) = (T)1;
        }

        return res;
    }
    else
        throw MyException("Matrix size error!");
}

template<typename T>
Matrix<T> Transpos(Matrix<T>& mat)
{
   Matrix<T> res(mat.cols(), mat.rows());

   for(int i=0; i<mat.cols(); ++i)
   {
       for(int j=0; j<mat.rows(); ++j)
       {
           res(i, j) = mat(j, i);
       }
   }

   return res;
}

template<typename T>
T Det(const Matrix<T>& mat)
{
    Matrix<T> tmp(mat);

    if(tmp.rows() != tmp.cols())
        throw MyException("Det(mat): Matrix is not square!");

    if(tmp.rows() == 1)
        return tmp(0, 0);

    if(tmp.rows() == 2)
        return tmp(0, 0) * tmp(1, 1) - tmp(1, 0) * tmp(0, 1);

    T n, d = 1, size = tmp.rows();
    int maxr, first;
    for(int i=0; i<tmp.rows()-1; ++i)
    {
        maxr = i;
        first = i;
        while(first < size)
        {

            for(int j=i+1; j<size; ++j)
            {
                if(std::abs(tmp(maxr, first)) < std::abs(tmp(j, first)))
                    maxr = j;
            }

            if(std::abs(tmp(maxr, first)) < ERR)
                first++;
            else
                break;
        }

        if(first >= tmp.cols())
            return 0;

        if(maxr != i)
            tmp.swapRows(maxr, i);

        for(int j=i+1; j<tmp.rows(); ++j)
        {
            n = tmp(j, i) / tmp(i, i);

            for(int k=i+1; k<tmp.cols(); ++k)
            {
                tmp(j, k) -= tmp(i, k) * n;
            }
            tmp(j, i) = 0;
        }
        d *= tmp(i, i);
    }
    d *= tmp(tmp.rows()-1, tmp.cols()-1);

    return d;
}

template<typename T>
bool LU(const Matrix<T>& mat, Matrix<T>& L, Matrix<T>& U)
{
    if(mat.rows() != mat.cols())
        throw MyException("LU(mat): Matrix is not square!");


    T n;
    int maxr, first, size = mat.rows();
    Matrix<T> l(size, size), u(mat);

    for(int i=0; i<size-1; ++i)
    {
        first = i;
        maxr = i;
        while(first < size)
        {

            for(int j=i+1; j<size; ++j)
            {
                if(std::abs(u(maxr, first)) < std::abs(u(j, first)))
                    maxr = j;
            }

            if(std::abs(u(maxr, first)) < ERR)
                first++;
            else
                break;
        }

        if(first >= size)
            break;

        if(maxr != i)
        {
            u.swapRows(maxr, i);
            l.swapRows(maxr, i);
        }

        for(int j=i+1; j<size; ++j)
        {
            l(j, first) = u(j, first) / u(i, first);
            l(j, first) = (l(j, first) == 0.0) ? 0.0 : l(j, first); // avoid -0
        }

        for(int j=i+1; j<size; ++j)
        {
            n = u(j, first) / u(i, first);
            for(int k=i+1; k<size; ++k)
            {
                u(j, k) -= u(i, k) * n;
            }
            u(j, first) = 0;
        }
    }

    for(int i=0; i<size; ++i)
        l(i, i) = 1;

    L = l;
    U = u;

    return true;
}

template<typename T>
Matrix<T> Minor(const Matrix<T>& mat, int r, int c)
{
    int n = mat.rows();

    if(mat.rows() != mat.cols())
        throw MyException("Minor(mat): Matrix is not suqare!");

    if(n == 1)
        throw MyException("Minor(mat): Matrix size <= 1*1");

    Matrix<T> m(n-1, n-1);
    for(int i=0, ii=0; i<n; ++i)
    {
        if(i == r)continue;

        for(int j=0, jj=0; j<n; ++j)
        {
            if(j == c)continue;

            m(ii, jj) = mat(i, j);
            jj++;
        }
        ii++;
    }

    return m;
}

template<typename T>
Matrix<T> Adj(const Matrix<T>& mat)
{
    int n = mat.rows();

    if(mat.rows() != mat.cols())
        throw MyException("Adj(mat): Matrix is not square!");

    if(n == 1)
        throw MyException("Adj(mat): Matrix size <= 1*1");

    Matrix<T> res(n, n);
    for(int i=0, k=1; i<n; ++i)
    {
        for(int j=0; j<n; ++j, k*=-1)
        {
           res(i, j) = k * Det<T>(Minor<T>(mat, i, j));
        }
    }

    return Transpos<T>(res);
}

template<typename T>
int Rank(const Matrix<T>& mat)
{
    T n;
    int maxr, first, size = mat.rows(), d = 0;
    Matrix<T> tmp(mat);

    for(int i=0; i<size-1; ++i)
    {
        first = i;
        maxr = i;
        while(first < size)
        {

            for(int j=i+1; j<size; ++j)
            {
                if(std::abs(tmp(maxr, first)) < std::abs(tmp(j, first)))
                    maxr = j;
            }

            if(std::abs(tmp(maxr, first)) < 1e-7)
                first++;
            else
                break;
        }

        if(first >= size)
            break;

        if(maxr != i)
            tmp.swapRows(maxr, i);

        for(int j=i+1; j<size; ++j)
        {
            n = tmp(j, first) / tmp(i, first);
            for(int k=i+1; k<size; ++k)
            {
                tmp(j, k) -= tmp(i, k) * n;
            }
            tmp(j, first) = 0;
        }
        d++;
    }

    return d;
}

template<typename T>
Matrix<T> Inv(const Matrix<T>& mat)
{
    if(mat.rows() != mat.cols())
        throw MyException("Inv(mat): Matrix is not square!");

    T n;
    int maxr, size = mat.rows();
    Matrix<T> tmp(mat), invm(Identity<T>(size));

    for(int i=0; i<size-1; ++i)
    {
        maxr = i;

        for(int j=i+1; j<size; ++j)
        {
            if(std::abs(tmp(maxr, i)) < std::abs(tmp(j, i)))
                maxr = j;
        }

        if(std::abs(tmp(maxr, i)) < ERR)
           throw MyException("Inv(mat): Matrix is singular!");

        if(maxr != i)
        {
            tmp.swapRows(maxr, i);
            invm.swapRows(maxr, i);
        }

        for(int j=i+1; j<size; ++j)
        {
            n = tmp(j, i) / tmp(i, i);
            for(int k=0; k<size; ++k)
            {
                tmp(j, k) -= tmp(i, k) * n;
                invm(j, k) -= invm(i, k) * n;
            }
        }
    }

    for(int i=0; i<size; ++i)
    {
        n = tmp(i, i);
        for(int j=0; j<size; ++j)
        {
            invm(i, j) /= n;
            tmp(i, j) /= n;
        }
    }

    for(int i=size-1; i>0; --i)
    {
        for(int j=i-1; j>=0; --j)
        {
           for(int k=size-1; k>=0; --k)
           {
               invm(j, k) -= invm(i, k) * tmp(j, i);
           }
        }
    }

    return invm;
}
