
#include <cmath>
//ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
//
//         mconvd
//
//ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
void mconvd(const int n, double z[], double r[], double pl0[], int * indf) {
  const double am = 3.4e138;
  const double rp = 5.0e-14;
  --pl0;
  --r;
  --z;
  if (n < 1) return;
  double aps = am / n;
  double ap = 1.0e0 / aps;
  aps = sqrt(aps);
  for (int i = 1, ir = 0; i <= n; ++i) {
    do {
      ++ir;
    } while (pl0[ir] <= 0.0e0);
    int ni = i * (i - 1) / 2;
    int ii = ni + i;
    if (z[ii] <= rp * fabs(r[ir]) || z[ii] <= ap) {
      pl0[ir] = -2.0e0;
      r[ir] = 0.0e0;
      *indf = ir - 1;
      return;
    }
    z[ii] = 1.0e0 / sqrt(z[ii]);
    for (int nl = ii - 1; nl - ni > 0; --nl) {
      z[nl] *= z[ii];
      if (fabs(z[nl]) >= aps) {
        int ir = 0;
        for (int i = 1; i <= i + nl - ii; ++i) {
          do {
            ++ir;
          } while (pl0[ir] <= 0.0e0);
        }
        pl0[ir] = -2.0e0;
        r[ir] = 0.0e0;
        *indf = ir - 1;
        return;
      }
    }
    if (i - n >= 0) break;
    int k = n + 1;
    do {
      --k;
      int nk = k * (k - 1) / 2;
      int nl = nk;
      int kk = nk + i;
      double d = z[kk] * z[ii];
      double c = d * z[ii];
      int l = k;
      do {
        int ll = nk + l;
        int li = nl + i;
        z[ll] -= z[li] * c;
        --l;
        nl -= l;
      } while (l - i > 0);
      --l;
      while (l > 0) {
        int ll = nk + l;
        int li = ni + l;
        z[ll] -= z[li] * d;
        --l;
      }
      z[kk] = -c;
    } while (k - i - 1 > 0);
  }
  for (int i = 1; i <= n; ++i) {
    for (int k = i; k <= n; ++k) {
      int nl = k * (k - 1) / 2;
      int ki = nl + i;
      double d = 0.0e0;
      for (int l = k; l <= n; ++l) {
        int li = nl + i;
        int lk = nl + k;
        d += z[li] * z[lk];
        nl += l;
      }
      ki = k * (k - 1) / 2 + i;
      z[ki] = d;
    }
  }
}
