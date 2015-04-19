#include "vector.h"
#include "math.h"


template<typename T>
Vector<T>::Vector(int d)
{

    if(d > 0){
        dimensions = d;
        data.resize(dimensions);
    }


    else
        throw MyException("Matrix Initialize Error!");

}
//get data===================================================================================
template<typename T>
std::string Vector<T>::toString()
{
    std::stringstream ss;

    for(int i=0; i<dimensions; ++i)
    {
        ss << data[i] << " ";
    }

    return ss.str();
}
//overload operation=========================================================================
template<typename T>
T& Vector<T>::operator()(int site)
{
    if(check(site))
        return data[site];
    else
        throw MyException("Index out of range!");
}

template<typename T>
const T& Vector<T>::operator()(int site) const
{
    if(check(site))
        return data[site];
    else
        throw MyException("Index out of range!");
}

template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& vec)
{
    if(&vec == this)
        return *this;

    int new_dimensions = vec.dimensions;

    data.resize(new_dimensions);
    for(int i=0; i<new_dimensions; ++i)
    {
        data[i] = vec(i);
    }

    dimensions = new_dimensions;

    return *this;
}

template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& vec)
{
    Vector<T> res(dimensions);

    if( dimensions == vec.dimensions )
    {
        for(int i=0; i<dimensions; ++i)
        {
            res(i) = data[i] + vec(i) ;
        }
    }
    else
        throw MyException("Vector addition error: dimension not match!");

    return res;
}

template<typename T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& vec)
{
    return *this = *this + vec;
}

template<typename T>
Vector<T>& Vector<T>::operator-()
{
    for(int i=0; i<dimensions; ++i)
    {
        data[i] = -data[i];
    }

    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator-(const Vector<T>& vec)
{
    Vector<T> tmp = vec;

    return *this + (-tmp);
}

template<typename T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& vec)
{
    return *this = *this - vec;
}

template<typename T>
T Vector<T>::operator*(const Vector<T>& vec)
{
    T dot;

    if( dimensions == vec.getDim() )
    {
        for(int i=0; i<dimensions; ++i)
        {
            dot += data[i]*vec(i) ;
        }
    }
    else
        throw MyException("Vector multiplication error: dimension not match!");

    return dot;
}

template<typename T>
Vector<T> Vector<T>::operator*(const T& var)
{
    Vector<T> res(dimensions);

    for(int i=0; i<dimensions; ++i)
    {
        res(i) = data[i] * var;
    }

    return res;
}

template<typename T>
Vector<T&> Vector<T>::operator*=(const T& var)
{
    return *this = this * var;
}

template<typename T>
Vector<T> Vector<T>::operator/(const T& var)
{
    Vector<T> res(dimensions);

    if( var!=0 )
    {
        for(int i=0; i<dimensions; ++i)
        {
            res(i) = data[i] / var ;
        }
    }
    else
        throw MyException("Vector division error: divide by zero!");

    return res;
}

template<typename T>
Vector<T>& Vector<T>::operator/=(const T& var)
{
    return *this = this / var;
}

template<typename T>
bool Vector<T>::operator==(const Vector<T>& vec)
{
    if(dimensions != vec.dimensions)
        return false;

    for(int i=0; i<dimensions; ++i)
    {
        if(data[i]!=vec(i))
            return false;
    }

    return true;
}

template<typename T>
bool Vector<T>::operator!=(const Vector<T>& vec)
{
    return !(*this==vec);
}

//vector operation=========================================================================
template<typename T>
T Vector<T>::norm(const Vector<T>& vec)
{
    T sum;

    sum = vec * vec;

    return sqrt(sum);
}

template<typename T>
Vector<T> Vector<T>::normal(const Vector<T>& vec)
{
    Vector<T> res(vec.dimensions);

    res = vec / norm(vec);

    return res;
}

template<typename T>
Vector<T> Vector<T>::cross(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    Vector<T> res(vec1.dimensions);

    res(0) = (vec1(1)*vec2(2)) - (vec2(2)*vec1(1));
    res(1) = -((vec2(0)*vec1(2)) - (vec2(2)*vec1(0)));
    res(3) = (vec2(0)*vec1(1)) - (vec2(1)*vec1(0));

    return res;
}

template<typename T>
T Vector<T>::com(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    T com;

    com = (vec1 * vec2) / norm(vec2);

    return com;
}

template<typename T>
Vector<T> Vector<T>::proj(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    Vector<T> res(vec1.dimensions);

    res = com(vec1,vec2) * normal(vec2);

    return res;
}

template<typename T>
T Vector<T>::area(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    T area;

    area = norm(cross(vec1 , vec2))/2;

    return area;
}

template<typename T>
bool Vector<T>::isParallel(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    Vector<T> res(vec1.dimensions);

    for(int i=0 ; i<vec1.size() ;++i)
    {
        if(vec1[i]!=0)
            res(i) = vec2(i) / vec1(i);
        else
            res(i) = vec1(i) / vec2(i);
    }
    for(int i=0 ; i<vec1.size()-1 ;++i)
    {
        if(abs(res(i)-res(i+1))>1e-7)
            return false;
    }

    return true;
}

template<typename T>
bool Vector<T>::isOrthogonal(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    T dot;

    dot = vec1 * vec2;
    if(dot==0)
        return true;
    else
        return false;
}

template<typename T>
double Vector<T>::angle(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    double angle;

    angle = acos((vec1 * vec2)/(norm(vec1)*norm(vec2))) * 180.0 / M_PI;

    return angle;
}

template<typename T>
Vector<T>  Vector<T>::pn(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    Vector<T> res(vec1.dimensions);

    res = cross(vec1 , vec2);

    return res;
}

template<typename T>
std::vector< Vector<T> > Vector<T>::ob(const std::vector< Vector<T> > v)
{
    std::vector< Vector<T> > res;
    res = v;

    for(int i=1 ; i<v.size() ; ++i){
        for(int j=0 ; j<i ; ++j){
            res[i] -= pj(res[i] , v[j]);
        }
        res[i] = normal(res[i]);
    }

    return res;
}





