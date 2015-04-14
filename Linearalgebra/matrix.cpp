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
            ss << data[i][j] << " ";
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
Matrix<T> identity(int dim)
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
Matrix<T> transpos(Matrix<T>& mat)
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
