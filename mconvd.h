
#pragma once

inline int ind(const int i, const int j) {
   if (i >= j) return (i*(i+1)/2 + j);
   else        return (j*(j+1)/2 + i);
}
void mconvd(const int n, double z[], double r[], double pl0[], int * indf);
