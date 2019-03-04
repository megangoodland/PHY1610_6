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
#include <string>
#include <cmath> 
#include <algorithm>
#include "netCDF_reading.h"
#include "rarray_math.h"
#include <cblas.h>
using namespace std;
using namespace netCDF;


int main(){
  
  // Initializing rarrays and constants
  const int f_size = get_f_size("GWprediction.nc"); // Only need to do this once because f is same size in all files
  int n_detections = 32; // number of detections
  rarray<complex<double>,1> f(f_size); // initialize array to hold f
  rarray<complex<double>,1> fhat(f_size); // initialize array to hold fhat
  rarray<double,1> Fk(f_size); // initialize array to hold Fk
  rarray<double,1> Gk(f_size); // initialize array to hold Gk
  rarray<double,1> C(n_detections); // array to hold correlation values
  rarray<string,1> files(n_detections); // array to hold filenames of detections
  rarray<string,1> largest_C(5); // array to hold filenames of detections that make the highest correlation values
  
  f = get_f("GWprediction.nc"); // Fill f with data from netCDF file
  fhat = fft(f); // Get fast fourier transform
  Fk = sq_norm(fhat); // Get Fk

  // This loop fills an rarray with all of the detection filenames
  for (int i=0; i<(n_detections); i++){
    int n = i+1;
    if (n<10){
        string filename = "detection0" + to_string(n) + ".nc";
        files[i]=filename; // Add to rarray that holds all of the file names
    }
    if (n>=10){
        string filename = "detection" + to_string(n) + ".nc";
        files[i]=filename; // Add to rarray that holds all of the file names    
    }
  }
  
  // This loop calculates every correlation and saves it in array C
  for (int i=0; i<(n_detections); i++){
    // Get Gk with same steps as Fk, overwriting f and fhat
    f = get_f(files[i]);
    fhat = fft(f);
    Gk = sq_norm(fhat);
    C[i]=correlation(Fk, Gk); // Add to rarray that holds all of the correlation values
  }
  
  // This identifies the detection files that make the 5 highest values in C
  double C2[] = C;
  cout << max_element(C2,C2+n_detections) << endl;
    
  
  cout << C << endl;
  return 0;
}
