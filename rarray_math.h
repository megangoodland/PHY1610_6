// rarray_math.h
// Header for the module for doing math functions on rarrays for GW_detection.cc

#ifndef rarray_math_h
#define rarray_math_h


#include <rarray>
#include <complex>


// fft: Fast fourier transform function for 1D complex double rarray
//      Input: f rarray
//      Output: fhat rarray
rarray<std::complex<double>,1> fft(rarray<std::complex<double>,1>& f);

// sq_norm: This function returns an rarray of the square norms of the complex values in the rarray input
//      Input: 1D complex double rarray
//      Output: 1D double rarray
rarray<double,1> sq_norm(rarray<std::complex<double>,1>& fhat);

// correlation: This function calculates the correlation between two 1D rarrays via C_FG = <F,G> / sqrt(<F,F><G,G>)
//      Input: Two 1D double rarrays
//      Output: A double value
double correlation(rarray<double,1>& Fk,rarray<double,1>& Gk);


#endif /* rarray_math */
