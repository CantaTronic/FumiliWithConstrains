
#include <cmath>
//ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
//
//         mconvd
//
//ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
int ind(const int i, const  int j)
{
   if (i >= j) return (i*(i + 1) / 2 + j);
   else return (j*(j + 1) / 2 + i);
}
void mconvd(int n, double z[], double r[], double pl0[], int *indf)
{
   const double am = 3.4e138;
   const double rp = 5.0e-14;
   double  ap, aps, c, d;
   int i, k, l, ii, ki, li, kk, ni, ll, nk, nl, ir, lk;
   --pl0;
   --r;
   --z;
   if (n < 1) {
      return;
   }
   aps = am / n;
   aps = sqrt(aps);
   ap = 1.0e0 / (aps * aps);
   ir = 0;
   for (i = 1; i <= n; ++i) {
L1:
      ++ir;
      if (pl0[ir] <= 0.0e0) {
         goto L1;
      } else {
         goto L2;
      }
L2:
      ni = i * (i - 1) / 2;
      ii = ni + i;
      k = n + 1;
      if (z[ii] <= rp * fabs(r[ir]) || z[ii] <= ap) {
         goto L19;
      }
      z[ii] = 1.0e0 / sqrt(z[ii]);
      nl = ii - 1;
L3:
      if (nl - ni <= 0) {
         goto L5;
      } else {
         goto L4;
      }
L4:
      z[nl] *= z[ii];
      if (fabs(z[nl]) >= aps) {
         goto L16;
      }
      --nl;
      goto L3;
L5:
      if (i - n >= 0) {
         goto L12;
      } else {
         goto L6;
      }
L6:
      --k;
      nk = k * (k - 1) / 2;
      nl = nk;
      kk = nk + i;
      d = z[kk] * z[ii];
      c = d * z[ii];
      l = k;
L7:
      ll = nk + l;
      li = nl + i;
      z[ll] -= z[li] * c;
      --l;
      nl -= l;
      if (l - i <= 0) {
         goto L9;
      } else {
         goto L7;
      }
L8:
      ll = nk + l;
      li = ni + l;
      z[ll] -= z[li] * d;
L9:
      --l;
      if (l <= 0) {
         goto L10;
      } else {
         goto L8;
      }
L10:
      z[kk] = -c;
      if (k - i - 1 <= 0) {
         goto L11;
      } else {
         goto L6;
      }
L11:
      ;
   }
L12:
   for (i = 1; i <= n; ++i) {
      for (k = i; k <= n; ++k) {
         nl = k * (k - 1) / 2;
         ki = nl + i;
         d = 0.0e0;
         for (l = k; l <= n; ++l) {
            li = nl + i;
            lk = nl + k;
            d += z[li] * z[lk];
            nl += l;
         }
         ki = k * (k - 1) / 2 + i;
         z[ki] = d;
      }
   }
L15:
   return;
L16:
   k = i + nl - ii;
   ir = 0;
   for (i = 1; i <= k; ++i) {
L17:
      ++ir;
      if (pl0[ir] <= 0.0e0) {
         goto L17;
      }
   }
L19:
   pl0[ir] = -2.0e0;
   r[ir] = 0.0e0;
   *indf = ir - 1;
   goto L15;
}
