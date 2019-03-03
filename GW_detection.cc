// GW_detection.cc
// 
// This file will...
// 1. Read the predicted GW signal from GWprediction.nc.
// 2. Read one of the GW signal from observations detection01.nc . . . detection32.nc.
// 3. Compute the FFT of the two complex quantities, using FFTW.
// 4. Compute the power spectrum of both signals.
// 5. Compute the correlation coefficient between the power spectra using a .BLAS call for the inner product
// 6. Output the correlation coefficient
// 7. Repeat steps 2-to-6 for each of the signals in the observation set.
// 8. Finally, determine the 5 most significant candidates (those with the 5 largest values of the correlation coefficient) from the observations set.
//
// Compile with 'make'.
// On the clusters, it will be necessary to load the "gcc" module first.
//
// Run with
//
//     ./GW_detection
// 
// Megan Goodland, March 2019.
// Assignment 6 of PHY1610 at University of Toronto

#include <iostream>
#include <rarray>
#include <rarrayio>
#include <vector>
#include <netcdf>
#include <complex>
#include <fftw3.h>
#include <cmath> // for pow
#include "netCDF_reading.h"
#include <cblas.h>
using namespace std;
using namespace netCDF;

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



int main(){
  
  // Initializing rarrays and constants
  const int f_size = get_f_size("GWprediction.nc"); // Only need to do this once because f is same size in all files
  rarray<complex<double>,1> f(f_size); // initialize array to hold f
  rarray<complex<double>,1> fhat(f_size); // initialize array to hold fhat
  rarray<double,1> Fk(f_size); // initialize array to hold Fk
  rarray<double,1> Gk(f_size); // initialize array to hold Gk
  
  // Fill f with data from netCDF file
  f = get_f("GWprediction.nc");
  // Get fast fourier transform
  fhat = fft(f);
  // Get Fk
  Fk = sq_norm(fhat);

  // Get Gk with same steps, overwriting f and fhat
  f = get_f("detection01.nc");
  fhat = fft(f);
  Gk = sq_norm(fhat);
  
  double x = correlation(Fk, Gk);
  cout << x << endl;
    
  return 0;
}
