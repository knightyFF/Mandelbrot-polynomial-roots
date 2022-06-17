//minimal implementation
#ifndef MY_COMPLEX_TYPE
#define MY_COMPLEX_TYPE
template<typename T>
class complex{
public:
    complex(){re=T(0.); im=T(0.);}
    complex(T rv){re=rv; im=T(0.);}
    complex(T rv, T iv){re=rv; im=iv;}
    complex(const complex<T>& cv){re=cv.re; im=cv.im;}
    complex<T>& operator=(const complex<T>& cv){re=cv.re; im=cv.im; return *this;}
    complex<T>& operator+=(const complex<T>& cv){re+=cv.re; im+=cv.im; return *this;}
    complex<T>& operator-=(const complex<T>& cv){re-=cv.re; im-=cv.im; return *this;}

    T& real(){return re;}
    const T& real()const {return re;}
    T& imag(){return im;}
    const T& imag()const {return im;}
public:
    T re,im;
};

template<typename T>
inline complex<T> operator+(const complex<T> &c1, const complex<T>& c2){
    return complex<T>(c1.re+c2.re, c1.im+c2.im);
}

template<typename T>
inline complex<T> operator-(const complex<T> &c1, const complex<T>& c2){
    return complex<T>(c1.re-c2.re, c1.im-c2.im);
}

template<typename T>
inline complex<T> operator-(const complex<T> &c){
    return complex<T>(-c.re, -c.im);
}

template<typename T>
inline complex<T> conjugate(const complex<T> &c){
    return complex<T>(c.re, -c.im);
}

template<typename T>
inline complex<T> sqr(const complex<T> &c){
    return complex<T>(c.re*c.re-c.im*c.im, T(2.)*c.re*c.im);
}

template<typename T>
inline T norm(const complex<T> &c){
    return c.re*c.re+c.im*c.im;
}

template<typename T>
inline T abs(const complex<T> &c){
    return std::sqrt(norm(c));
}

template<typename T>
inline complex<T> operator*(const complex<T>& c1, const complex<T>& c2){
    return complex<T>(c1.re*c2.re-c1.im*c2.im, c1.re*c2.im+c1.im*c2.re);
}

template<typename T>
inline complex<T> operator*(T c1, const complex<T>& c2){
    return complex<T>(c1*c2.re, c1*c2.im);
}

template<typename T>
inline complex<T> operator/(const complex<T> &c1, const complex<T>& c2){
    return T(1.)/norm(c2) * c1*conjugate(c2);
}
#endif // COMPLEX

