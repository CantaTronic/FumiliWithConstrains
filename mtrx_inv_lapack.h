
#pragma once

/**
  inplace inverse n x n matrix A.
  returns:
    ret = 0 on success
    ret < 0 illegal argument value
    ret > 0 singular matrix
*/
int mtrx_inv(const int N, double * const M);

void mtrx_print(const unsigned N, const double * const M);
