/*
  A simple Mandelbrot polynomial roots calculator by knighty. 12/13/2012
  Public domain.
  Use at your own risk.

  Adapted from an evadraw script. See attachment in: http://www.fractalforums.com/theory/the-mandelbrot-polynomial-roots-challenge/msg74387/#msg74387

  Finally using templates wasn't really a good idea. Anyway a good opportunity to discover c++11 great features.

  Many thanks to: Lycium, quazor, Adam Majewski and zebastian.

  Compilation with gcc :
    g++ -std=c++11 -Ofast -Wall -o Mandelroots.exe main.cpp
  (Using -Ofast switch is much faster when using doubles)
  It should compile fine on any compiler that supports c++11.

  Usage:
  mandelroots [order]
    default order = 10
  it saves the restults into roots_[order].txt
  (The number of roots is 2^(order-1) so don't give it silly numbers ;o) )

   Using doubles one can get all the roots up to order 25. This program could not find all
  roots for higher order polynimials because of the accuracy of doubles is insufficient. (see also line 34 in file mbrcalc.h where the order (m_MaxIter) is limited wrt the floating point type used).
   With long doubles one can go up to order 31 (> 1 billion roots). Beware: it will need lots of Gb and time. -Ofast doesn't apparently optimize for long double.

  Update:
  - using my own sloppy and minimal complex class implementation. Now when using long double, the exe is much faster. Also, using -Ofast is no longer much faster than -O3
*/
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <cmath>
#if 0
    #include <complex>
#else
    #include "complex.h"
#endif
#include <limits>
#include <thread>


using namespace std;

#include "mbformulas.h"
#include "mbrcalc.h"

int main(int argc, char* argv[])
{
    unsigned int nThreads = thread::hardware_concurrency();
    if(nThreads==0) nThreads = 1;
    cout << "Number of cores detected : " << nThreads << endl;
    int order=10;
    if(argc>1){
        unsigned int n = atoi(argv[1]);//what if argv[1] is not a number
        order = n;
    }
    //printf("Mandelbrot polynomial");
    cout << "Mandelbrot polynomial" << endl;
    cout << "Order (asked for): " << order << endl;

    MRootsCalc<double> mrc(order);

    cout << "Order : " << mrc.GetMaxIter() << endl;
    cout << "Expected number of roots : " << mrc.GetExpectedNbrRoots() << endl;
    cout << "Searchig for roots..." << endl;

    chrono::steady_clock::time_point start = chrono::steady_clock::now();
    mrc.findRootsMT(nThreads);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    //mrc.printThem();
    cout.precision(3);
    cout << "time taken : " << double(chrono::duration_cast < chrono::milliseconds >(end - start).count())/1000.0 << " s"<< endl;
    cout << "Number of roots found : " << mrc.GetNbrRoots() << endl;

    mrc.saveThem();
    //int t; cin >> t;

    return 0;
}

