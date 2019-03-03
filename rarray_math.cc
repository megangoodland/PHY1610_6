// rarray_math.cc
// Module for doing math functions on rarrays for GW_detection.cc

#include <iostream>
#include <rarray>
#include <complex>
#include <fftw3.h>
#include <cmath> // for pow
#include <cblas.h>
#include "rarray_math.h"

using namespace std;

// Fast fourier transform function: performs a fast fourier transform.
// Input: f rarray
// Output: fhat rarray
rarray<complex<double>,1> fft(rarray<complex<double>,1>& f){
  int f_size = f.extent(0);
  rarray<complex<double>,1> fhat(f_size); // initialize array to hold fhat
  fftw_plan p = fftw_plan_dft_1d(f_size,
                      (fftw_complex*)f.data(), (fftw_complex*)fhat.data(),
                      FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  return fhat;
}

// This function returns an rarray of the square norms of the complex values in the rarray input
rarray<double,1> sq_norm(rarray<complex<double>,1>& fhat){
  int f_size = fhat.extent(0);
  rarray<double,1> normsq(f_size);
  for (int i=0; i<f_size; i++){
    normsq[i] = pow((norm(fhat[i])),2); // calculate the norm squared for each value, put in new array
  }
  return normsq;
}
// This function calculates the correlation between two 1D rarrays via C_FG = <F,G> / sqrt(<F,F><G,G>)
double correlation(rarray<double,1>& Fk,rarray<double,1>& Gk){
  int f_size = Fk.extent(0); // they'll be the same size
  double *A = new double[f_size]; // put in a way that cblas will work with
  double *B = new double[f_size]; // put in a way that cblas will work with
  for (int i=0; i<f_size; i++) A[i] = Fk[i];
  for (int i=0; i<f_size; i++) B[i] = Gk[i];
  // calculate inner products
  double FG = cblas_ddot(f_size, A, 1, B, 1);
  double FF = cblas_ddot(f_size, A, 1, A, 1);
  double GG = cblas_ddot(f_size, B, 1, B, 1);
  // C_FG = <F,G> / sqrt(<F,F><G,G>)
  return (FG/(sqrt(FF*GG)));
}
