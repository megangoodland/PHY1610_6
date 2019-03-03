# Makefile for GW_detection

# Define the compiler to use
CXX = g++

#all: GW_detection

GW_detection: GW_detection.o rarray_math.o netCDF_reading.o 
	${CXX} -std=c++11 -O3 -march=native GW_detection.o rarray_math.o netCDF_reading.o -o GW_detection  -I${BLAS_INC} -L${BLAS_LIB} -lopenblas -lnetcdf_c++4 -lfftw3

GW_detection.o: GW_detection.cc
	${CXX} -std=c++11 -O3 -march=native -c GW_detection.cc -I${BLAS_INC} -L${BLAS_LIB} -lopenblas
	
rarray_math.o: rarray_math.cc
	${CXX} -std=c++11 -O3 -march=native -c rarray_math.cc -I${BLAS_INC} -L${BLAS_LIB} -lopenblas
	
netCDF_reading.o: netCDF_reading.cc netCDF_reading.h
	${CXX} netCDF_reading.cc -c -o netCDF_reading.o -c

clean:
	rm *.o GW_detection
