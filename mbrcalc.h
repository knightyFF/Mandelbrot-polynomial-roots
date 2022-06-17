#ifndef MBRCALC
#define MBRCALC

//-------------------------------------------------------------------------------------------------------
//these are used for sorting the array containing the roots
template<typename T>
struct CompRup{
    bool operator()(const complex<T>& c0, const complex<T>& c1){
        return c0.real() < c1.real();
    }
};

template<typename T>
struct CompIup{
    bool operator()(const complex<T>& c0, const complex<T>& c1){
        return c0.imag() < c1.imag();
    }
};
//
template<typename T>
class MRootsCalc;


//The class that computes the roots
template<typename T>
class MRootsCalc
{
public:
    MRootsCalc(int MaxIter){
        m_MaxIter = MaxIter;// it's actually the order of the Mandelbrot polynomial
        // choosing the FP accuracy and setting the maximal MaxIter wrt the floating point type used
        // MaxIter should less than bits/2-1. This is due to rounding errors... it will work but will miss some roots when maxIter > bits/2-1
        m_epsilon = numeric_limits<T>::epsilon();
        m_MaxIter = min(numeric_limits<T>::digits / 2 - 1 , m_MaxIter);
        m_eps = pow(2., -2*(m_MaxIter - 1));
        m_nbrRoots = 1 << (m_MaxIter - 1);//period = m_MaxIter - 1;
        //m_tab.reserve(1024);//this is not necessary
    }

    int GetMaxIter(){ return m_MaxIter;}

    int GetExpectedNbrRoots(){ return m_nbrRoots;}

    int GetNbrRoots(){ return m_tab.size();}

    void findRoots(){
        findRoots(complex<T>(0.,0.), T(2.));//this is the stupid way :D
        /* faster way :
            findRootsRecurse(complex<T>(-1.,1.), T(1.));
            findRootsRecurse(complex<T>(1.,1.), T(1.));
           then if c is a non real root, cs = conjugate(c) is also a root.
        */
    }

    void findRoots(const complex<T> &center, const T &halfWidth){//Please use a square that can be perfectly divided by powers of two
        m_tab.empty();
        //subdivide
        findRootsRecurse(center, halfWidth);
        //filter the duplicates and sort them
        filter();
    }

    void findRootsMT(const unsigned int nThreads){
        findRootsMT(complex<T>(0.,0.), T(2.), nThreads);
    }

    void findRootsMT(const complex<T> &center, const T &halfWidth, const unsigned int nThreads){//Please use a square that can be perfectly divided by powers of two
        if (m_MaxIter < 15 ){//Don't use MT because it would be slower
            findRoots(center, halfWidth);
            return;
        }

        m_tab.empty();
        m_tabForMT.empty();

        //The roots are not distributed evenly. We need a way to balance the load between the threads
        //Collect regions that may contain roots
        m_MaxIterForMT = 11;// that would give us roughly 1000 regions to search
        m_epsForMT = pow(2., -2*(m_MaxIterForMT - 1));
        findRootsRecurseForMT(center, halfWidth);
        //Now, m_tabForMT contains the regions of interest

        //Launch working threads
        vector< thread > tThrd;
        for (unsigned int i=0; i < nThreads; i++){
            m_tMRC.push_back(MRootsCalc <T>(m_MaxIter));
        }
        for (unsigned int i=0; i < nThreads; i++){
            tThrd.push_back(thread(MTfunc, this, i, nThreads));
        }
        for (unsigned int i=0; i < nThreads; i++) tThrd[i].join();

        //Computations are done. grab results
        for (unsigned int i=0; i < nThreads; i++){
            for(unsigned int j=0; j<m_tMRC[i].m_tab.size(); j++){
                m_tab.push_back(m_tMRC[i].m_tab[j]);
            }
            m_tMRC[i].m_tab.empty();
        }
        //We have all the roots found... filter them
        filter();
    }

    void printThem(){
        cout.precision(numeric_limits<T>::max_digits10);
        cout.setf(ios::fixed);
        for (int i = 0; i < m_tab.size() /*m_cIndex*/; i++){
            cout << m_tab[i].real() << "\t" << m_tab[i].imag() << endl;
        }
    }

    void saveThem(){//we need a try catch here but I don't want to use exceptions
        string fileName = "roots_" + std::to_string(m_MaxIter) + ".txt"; //
        cout << "Saving results to: " << fileName << endl;
        ofstream SaveFile(fileName);//("save.txt");
        SaveFile.precision(numeric_limits<T>::max_digits10);
        SaveFile.setf(ios::fixed);
        for (unsigned int i = 0; i < m_tab.size() /*m_cIndex*/; i++){
            SaveFile << m_tab[i].real() << "\t" << m_tab[i].imag() << endl;
        }
    }

private:

    void MTfunc(const unsigned int id, const unsigned int nThreads){
        for (unsigned int i = id; i < m_tabForMT.size(); i += nThreads){
            m_tMRC[id].findRootsRecurse(m_tabForMT[i], m_epsForMT);
        }
    }

    void filter(){
        //sort in real direction
        std::sort(m_tab.begin(),m_tab.end(),CompRup<T>());
        //remove duplicates
        clean();
        //sort in imaginary direction
        std::sort(m_tab.begin(),m_tab.end(),CompIup<T>());
        //remove duplicates
        clean();
        //sort in real direction
        std::sort(m_tab.begin(),m_tab.end(),CompRup<T>());
    }

    void clean(){
        int i0 = 0;
        int i1 = 1;
        int cIndex = m_tab.size();
        while (i1 < cIndex){
            if(std::max(abs(m_tab[i1].real()-m_tab[i0].real()), abs(m_tab[i1].imag()-m_tab[i0].imag())) < m_eps){
                i1++;
            } else {
                i0++;
                m_tab[i0] = m_tab[i1];
                i1++;
            }
        }
        cIndex = i0+1;
        m_tab.resize(cIndex);
    }

    inline complex<T> refine(complex<T> c){
        //use Newton method
        complex<T> z, dz;
        for(int i=0; i<20; i++){
            int mi=MPoly_ZdZ(c, m_MaxIter, z, dz);
            if(mi<m_MaxIter) cout << "boum";//doesn't happen...
            complex<T> delta = z / dz;
            c -= delta;
            //if(norm(delta) < m_epsilon*m_epsilon) break;
            if((max(abs(delta.real()),abs(delta.imag()))) < T(0.5)*m_epsilon) break;
            //if(abs(delta/c) < m_epsilon) break;
            //if(abs(delta) < m_epsilon * abs(c)) break;
        }
        return c;
    }

    inline void store(complex<T> cc){
        int cIndex = m_tab.size();
        if (cIndex>0){//first time filtering
            complex<T> cc1 = m_tab[cIndex -1];
            if(std::max(abs(cc1.real()-cc.real()), abs(cc1.imag()-cc.imag())) < m_eps) return;
        }
        m_tab.push_back(cc);
    }

    void findRootsRecurse(complex<T> cc /*center of square*/, T hw /*half width of square*/){
        //verify if there is no roots in the disc (cc, hw*sqrt(2.))
        if( ! MFDiscContainsZero(cc, hw * T(1.42), m_MaxIter) ) return;
        //verify if there is only one root. if true, report it in m_tab.
        //One way of doing this is to compute an interval of the derivative of F_{n}(c,0) wrt c.
        //if this interval doesn't contain 0 then there is only one root.
        //but we use a simpler test. we know that the minimal separation between roots is m_eps.
        if(hw <= m_eps){// should be == but who knows
            //we know there is only one root. refine then save it.
            cc = refine(cc);
            store(cc);
            return;
        }
        //Many possible roots -> subdivide
        hw *= T(0.5);
        findRootsRecurse(cc + complex<T>(-hw,-hw),hw);
        findRootsRecurse(cc + complex<T>( hw,-hw),hw);
        findRootsRecurse(cc + complex<T>( hw, hw),hw);
        findRootsRecurse(cc + complex<T>(-hw, hw),hw);
    }

    void findRootsRecurseForMT(complex<T> cc /*center of square*/, T hw /*half width of square*/){
        //verify if there is no roots in the disc (cc, hw*sqrt(2.))
        if( ! MFDiscContainsZero(cc, hw * T(1.42), m_MaxIter) ) return;
        //if we attain max recursion depth report it in m_tabForMT.
        if(hw <= m_epsForMT){// should be == but who knows
            m_tabForMT.push_back(cc);
            return;
        }
        hw *= T(0.5);
        findRootsRecurseForMT(cc + complex<T>(-hw,-hw),hw);
        findRootsRecurseForMT(cc + complex<T>( hw,-hw),hw);
        findRootsRecurseForMT(cc + complex<T>( hw, hw),hw);
        findRootsRecurseForMT(cc + complex<T>(-hw, hw),hw);
    }

    int m_MaxIter, m_MaxIterForMT;
    int m_nbrRoots;
    T m_eps, m_epsilon, m_epsForMT;
    vector< complex<T> > m_tab, m_tabForMT;
    vector< MRootsCalc <T> > m_tMRC;
};

#endif // MBRCALC

