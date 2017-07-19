
#pragma once

/**
  inplace inverse n x n matrix A.
  returns:
    ret = 0 on success
    ret < 0 illegal argument value
    ret > 0 singular matrix
*/
int mtrx_inv(const unsigned n, double * const A);

void mtrx_print(const unsigned n, const double * const A);
