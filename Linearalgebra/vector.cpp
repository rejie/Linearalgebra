#include "vector.h"

template<typename T>
Vector<T>::Vector()
{

   dimensions = 1;
   data.resize(dimensions);

}

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

template<typename T>
Vector<T>::Vector(int d, int init_var)
{

    if(d > 0){
        dimensions = d;
        data.resize(dimensions, init_var);
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

    int new_dimensions = vec.getDim();

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

    if( dimensions == vec.getDim() )
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
Vector<T> Vector<T>::operator-(const Vector<T>& vec)
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
    T dot = 0;

    if( dimensions == vec.getDim() )
    {
        for(int i=0; i<dimensions; ++i)
        {
            dot += data[i]*vec(i) ;
        }
        return dot;
    }
    else
        throw MyException("Vector multiplication error: dimension not match!");

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
Vector<T>& Vector<T>::operator*=(const T& var)
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
    return *this = *this / var;
}

template<typename T>
bool Vector<T>::operator==(const Vector<T>& vec)
{
    if(dimensions != vec.getDim())
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
T norm(const Vector<T>& vec)
{
    T sum = 0;
    Vector<T> vec2 = vec;
    sum = vec2 * vec2;

    return sqrt(sum);
}

template<typename T>
Vector<T> normal(const Vector<T>& vec)
{
    Vector<T> res(vec.getDim());
    Vector<T> vec2 = vec;

    res = vec2 / norm(vec2);

    return res;
}

template<typename T>
Vector<T> cross(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    Vector<T> res(vec1.getDim());


    res(0) = (vec1(1)*vec2(2)) - (vec1(2)*vec2(1));
    res(1) = -((vec1(0)*vec2(2)) - (vec1(2)*vec2(0)));
    res(2) = (vec1(0)*vec2(1)) - (vec1(1)*vec2(0));

    return res;
}

template<typename T>
T com(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    T com;
    Vector<T> a = vec1;
    Vector<T> b = vec2;

    com = (a * b) / norm(b);

    return com;
}

template<typename T>
Vector<T> proj(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    Vector<T> res(vec1.getDim());


    res = normal(vec2) * com(vec1,vec2);

    return res;
}

template<typename T>
T area(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    T area;
    Vector<T> a = vec1;
    Vector<T> b = vec2;

    area = norm(b)*norm(a-proj(a,b))/2;

    return area;
}

template<typename T>
bool isParallel(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    if(angle(vec1 , vec2)==0)
        return true;
    else
        return false;
}

template<typename T>
bool isOrthogonal(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    T dot;

    dot = vec1 * vec2;
    if(dot==0)
        return true;
    else
        return false;
}

template<typename T>
double angle(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    double angle;
    Vector<T> a = vec1;
    Vector<T> b = vec2;

    angle = acos((a * b)/(norm(a)*norm(b))) * 180.0 / M_PI;

    return angle;
}

template<typename T>
Vector<T> pn(const Vector<T>& vec1 , const Vector<T>& vec2)
{
    Vector<T> res(vec1.getDim());

    res = cross(vec1 , vec2);

    return res;
}

template<typename T>
std::vector< Vector<T> > ob(const std::vector< Vector<T> > v)
{
    std::vector< Vector<T> > res;
    res = v;

    for(int i=1 ; i<v.size() ; ++i){
        for(int j=0 ; j<i ; ++j){
            res[i] -= proj(v[i] , res[j]);
        }
    }
    for(int i=0 ; i<res.size() ; ++i){
        res[i] = normal(res[i]);
    }
    return res;
}

template<typename T>
T max(const Vector<T>& v)
{
    int maxid = 0;

    for(int i=1; i<v.getDim(); ++i)
    {
        if(std::abs(v(maxid)) < std::abs(v(i)))
            maxid = i;
    }

    return v(maxid);
}



