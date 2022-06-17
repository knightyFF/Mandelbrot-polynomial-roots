# Mandelbrot polynomial roots
 Computing millions of Mandelbrot polynomial roots using ball arothmetics.
 See : http://www.fractalforums.com/theory/the-mandelbrot-polynomial-roots-challenge/
 It is able to compute 2^20 roots -that is, more than one million- in less than 15 seconds and 4 millions (2^22) in less than a minute whan using double type. I'm sure it is possible to do better.
 
 *Compilation*
 - Using gcc : 
	g++ -std=c++11 -Ofast -Wall -o Mandelroots.exe main.cpp
 - Using another compiler : I don't know. It is a simple command line program so... any c++11 compliant compiler should work.
 
 *Usage:*
	mandelroots [order]
  default order = 10
  it saves the restults into roots_[order].txt
  The number of roots is 2^(order-1) so don't give it silly numbers ;o). Actually it limits the order to order 25 when using double type.
