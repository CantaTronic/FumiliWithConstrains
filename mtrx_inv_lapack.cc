
#include <cstdio>
#include <cstring>
#include "mtrx_inv_lapack.h"

extern "C" {
  // LU decomoposition of a general matrix
  void dgetrf_(const int *, const int *, double *, const int *, int *, int *);
  // generate inverse of a matrix given its LU decomposition
  void dgetri_(const int *, double *, const int *, int *, double *, const int *, int *);
}

FSqMtrx::FSqMtrx(int n1, int n2): N(n1), M(new double[N*N]) {
  if(n1 != n2) {
    fprintf(stderr, "Warning:\n");
    fprintf(stderr, "  trying to make square matrix of size %i*%i\n", n1, n2);
    fprintf(stderr, "  creating matrix %i*%i instead\n", n1, n1);
    fflush(stderr);
  }
}

FSqMtrx::FSqMtrx(int n, double * m): N(n), M(new double[N*N]) {
  memcpy(M, m, N*N*sizeof(double));
}

FSqMtrx::FSqMtrx(FSqMtrx & other): N(other.N), M(new double[N*N]) {
  memcpy(M, other.M, N*N*sizeof(double));
}

FSqMtrx & FSqMtrx::operator= (FSqMtrx & rhs) {
  delete [] M;
  N = rhs.N;
  M = new double[N*N];
  memcpy(M, rhs.M, N*N*sizeof(double));
  return *this;
}

int FSqMtrx::Invert(const int N, double * const M) {
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

double FSqMtrx::Determinant(const int N, const double * const M) {
  int errorHandler;
  int * pivotArray = new int[N];
  double * T = new double[N*N];
  memcpy(T, M, N*N*sizeof(double));
  dgetrf_(&N, &N, T, &N, pivotArray, &errorHandler);
  if(errorHandler != 0) {
    delete [] pivotArray;
    delete [] T;
    return 0.;
  }
  double det = 1.;
  for(int i = 0; i < N; i++)
    if(pivotArray[i]-1 != i) // note: fortran indexing starts from 1!
      det *= -1;
  for(int i = 0; i < N; i++)
    det *= T[i*N+i]; // diagonal element
  delete [] pivotArray;
  delete [] T;
  return det;
}

void FSqMtrx::Print(const int N, const double * const M) {
  for (int i = 0; i < N*N; i++) {
    if ((i%N) == 0) putchar('\n');
    printf("%+12.8f ", M[i]);
  }
  putchar('\n');
  fflush(stdout);
}
