#ifndef MBFORMULAS
#define MBFORMULAS

//---------------------------------------------------------------------------------------------------
#define COMPLEX complex<T>

// Say F_{n}(c,z) is defined as:
// F_{n+1} = F_{n}² + c
// F_{0}(c,z) = z
// This function returns:
// - the number of iterations used: < MaxIter if bailing out and ==MaxIter otherwise
// - z = F_{MaxIter}(c,0) ; In case no bailout occure
// - dz = d/dc F_{MaxIter}(c,0)
template<typename T>
inline int MPoly_ZdZ(const COMPLEX &c, const int MaxIter, COMPLEX& z, COMPLEX& dz)
{
    z=COMPLEX(0.,0.);// Of course this line can be removed. Then, the returned z_out = F_{MaxIter}(c,z_in)
    dz=COMPLEX(0.,0.);
    int i=0;
    for(; i<MaxIter && norm(z)<T(4.); i++){
        dz = T(2.) * z * dz + COMPLEX(1.);
        z  = z * z + c;
    }
    return i;
}

// Returns true if F_{MaxIter}([c],0) contains 0. that means that we MAY have a zero.
// if F_{MaxIter}([c],0) does not contain 0 then we are SURE there is no roots
// here [c] is the complex interval which is a disc of radius r centered at c.
// [c] = c + [e]r
// [e] = [0,1] exp(i [-inf, +inf]); it is the disc of radius 1 centered at 0.
// [e] can also be seen as a function: e(t,a)= t exp(i a) where t is in [0,1] and a is in [-inf, +inf]
// Here we use a kind of Reduced Affine Arithmetic (RAA) augmented by Taylor expansion... We set:
//   [zn] = F_{n}([c],0) = zn + [e]r dzn + 1/2 [e]² r² En
//   Where:
//     zn = F_{n}(c,0)
//     dzn= d/dc F_{n}(c,0)
//     En is a real number that contains the upper order error
//  Now we have
//    [zn1] = [zn]² + [c]
//          = (zn + [e]r dzn + 1/2 [e]² r² En)² + c + [e]r
//          = zn² + 2 [e]r zn dzn + 2/2 [e]² r² zn En + [e]²r² dzn² + 2/2 [e]^3 r^3 dzn En + 1/4 [e]^4 r^4 En² + c + [e]r
//          = (zn² + c) + (2 zn dzn + 1) [e]r + 1/2 [e]² r² (2 zn En + 2 dzn² + 2 [e]r dzn En + 1/2 [e]² r² En²)
//          = zn1 + dzn1 [e]r + 1/2 [e]² r² [En1]
//  Reducing [En1]:
//    En1 = 2 |zn| En + 2 |dzn|² + 2 [e]r |dzn| En + 1/2 [e]² r² En²
//    En1 = 2 (|zn| En + |dzn|² + r |dzn| En + 1/4 r² En²)
// One can also use higher and lower taylor expansion orders. For ex:
//    zn + [e]r En ; Lower order
//    zn + [e]r dzn + 1/2 [e]² r² ddzn + 1/6 [e]^3 r^3 En
//    ... etc.
// And yes, it really looks like "perturbation theory" of superfractalthing! for which it may be useful :)
//
// The drawback of this methos is that we have to compute absolute values of complex numbers which requires sqrt()s
// but it doesn't suffer from (catastrophic) interval expansion of regular IA when doing multiplications.
// TODO: Implement an RAA-Taylor class... Anyway, for this application we won't use it too much. :/
template<typename T>
inline bool MFDiscContainsZero1(const COMPLEX &c, const T &r, const int MaxIter)
{
    COMPLEX z = COMPLEX(0.);
    COMPLEX dz= COMPLEX(0.);
    T rz = T(0.), rdz= T(0.), Ei= T(0.), rr= T(0.);
    for(int i=0; i<MaxIter &&
                 rz-rr < T(2.) && /*bailout if [z] completly outside Disc((0,0),2)*/
                 rz-rr > T(-2.) /*or if [z] contains Disc((0,0),2)*/
               ; i++){
        Ei = T(2.) * (rdz * rdz + (rz + (rdz + T(0.25) * r * Ei) * r) * Ei);
        dz = T(2.) * z * dz + COMPLEX(1.);
        z  = z * z + c;
        rdz = abs(dz);
        rz = abs(z);
        rr = r * (rdz + T(0.5) * r * Ei); // Radius of uncertainty disc around z_i
    }
    if(rz - rr <= T(0.)) return true;//if [zn] contains 0 return true
    return false;
}

//same as above but with :
// [zn] = F_{n}([c],0) = zn + [e]r dzn + [e]² r² En
template<typename T>
bool MFDiscContainsZero(const COMPLEX &c, const T &r, const int MaxIter)
{
    COMPLEX z = COMPLEX(0.);
    COMPLEX dz= COMPLEX(0.);
    T rz = T(0.), rdz= T(0.), Ei= T(0.), rr= T(0.);
    for(int i=0; i<MaxIter && rz-rr < T(2.) && rz-rr > T(-2.); i++){//bail out when [zn] is entirely outside 2... or if it covers everything
        Ei = rdz * rdz + (T(2.) * rz + r * (T(2.) * rdz + r * Ei)) * Ei;
        dz = T(2.) * z * dz + COMPLEX(1.);
        //dz *= T(2.); dz *= z; dz += COMPLEX(1.);//no need to do this. the compiler seems smart enought to optimize
        z  = z * z + c;
        //z *= z; z += c;
        rdz = abs(dz);
        rz = abs(z);
        rr = r * (rdz + r * Ei); // Radius of uncertainty disc around z_i
    }
    if(rz - rr <= T(0.)) return true;//if [zn] contains 0 return true
    return false;
}
#endif // MBFORMULAS

