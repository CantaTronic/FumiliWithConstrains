
#include <stdio.h>
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

void mtrx_print(const unsigned n, const double * const A) {
  for (unsigned i = 0; i < n*n; i++) {
    if ((i%n) == 0) putchar('\n');
    printf("%+12.8f ", A[i]);
  }
  putchar('\n');
}
