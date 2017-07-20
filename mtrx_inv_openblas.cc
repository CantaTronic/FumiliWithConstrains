
#include <cstdio>
#include <cstring>
#include <openblas/lapacke.h>
#include "mtrx_inv_openblas.h"

int mtrx_inv(const unsigned n, double * const A) {
  int ipiv[n+1];
  int ret;
  ret = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);
  if (ret !=0)
    return ret;
  ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, A, n, ipiv);
  return ret;
}

double mtrx_det(const int n, const double * const M) {
  int ret;
  int ipiv[n+1];
  double * T = new double[n*n];
  memcpy(T, M, n*n*sizeof(double));
  ret = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, T, n, ipiv);
  if(ret != 0) {
    delete [] T;
    return 0.;
  }
  double det = 1.;
  // sign
  for(int i = 0; i < n; i++)
    if(ipiv[i]-1 != i) // note: fortran indexing starts from 1!
      det *= -1;
  // value
  for(int i = 0; i < n; i++)
    det *= T[i*n+i]; // diagonal element
  delete [] T;
  return det;
}

void mtrx_print(const unsigned n, const double * const A) {
  for (unsigned i = 0; i < n*n; i++) {
    if ((i%n) == 0) putchar('\n');
    printf("%+12.8f ", A[i]);
  }
  putchar('\n');
}
