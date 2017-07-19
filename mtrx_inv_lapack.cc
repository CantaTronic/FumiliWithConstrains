
#include <stdio.h>
// #include <stddef.h>
#include "mtrx_inv_lapack.h"

extern "C" {
  // LU decomoposition of a general matrix
  void dgetrf_(const int *, const int *, double *, const int *, int *, int *);
  // generate inverse of a matrix given its LU decomposition
  void dgetri_(const int *, double *, const int *, int *, double *, const int *, int *);
}

int mtrx_inv(const int N, double * const M) {
  int errorHandler;
  int * pivotArray = new int[N];
  dgetrf_(&N, &N, M, &N, pivotArray, &errorHandler);
  if(errorHandler != 0) {
    delete [] pivotArray;
    return errorHandler;
  }
  const int NN = N*N;
  double * lapackWorkspace = new double[NN];
  dgetri_(&N, M, &N, pivotArray, lapackWorkspace, &NN, &errorHandler);
  delete [] pivotArray;
  delete [] lapackWorkspace;
  return errorHandler;
}

void mtrx_print(const unsigned N, const double * const M) {
  for (unsigned i = 0; i < N*N; i++) {
    if ((i%N) == 0) putchar('\n');
    printf("%+12.8f ", M[i]);
  }
  putchar('\n');
}
