
#include "Fumili.h"
// #include <string>
// #include <vector>
// #include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double data[100000][2];

int SGZ(int m, double &S, double A[], double PL[], double G[], double Z[]);
void constraints(int M, double A[], double * psis, double ** dpsis);

int main(int arc, char ** argv) {
  idebug = 0;

  ifstream ifs("unif.dat");
  for(int i = 0; i < 100000; i++)
    ifs >> data[i][0] >> data[i][1];
  ifs.close();

  double AMN[2], AMX[2], PL0[2], A[2], R[2], SIGMA[2], VL[2*2], akappa, S;
  int N1, N2, N3, IT;
  double EPS;

  for(int i = 0; i < 2; i++) {
    AMN[i] = -10.;
    AMX[i] = 10.;
    PL0[i] = .1;
    A[i] = .0;
  }
  N1 = 1, N2 = 1, N3 = 30, IT = -30;
  EPS = .1;
  if(idebug)
    for(int i = 0; i < 2; i++) {
      cout << " i,A,PL0,AMX,AMN "
           << i << "  "
           << A[i] << "  "
           << PL0[i] << "  "
           << AMX[i] << "  "
           << AMN[i] << endl;
  }
  fumiliSK(2, S, N1, N2, N3,
           EPS, IT, A, PL0,
           AMX, AMN, R, SIGMA,
           SGZ, akappa, VL);
  if(akappa < EPS) {
    cout << " UnConstrained Fit: Fitted Values and their Errors " << endl;
    for(int i = 0; i < 2; i++) {
       cout << i << "  " << A[i] << " +/-  " << SIGMA[i] << endl;
    }
  }

  for(int i = 0; i < 2; i++) {
    AMN[i] = -10.;
    AMX[i] = 10.;
    PL0[i] = .1;
    A[i] = .0;
  }
  N1 = 1, N2 = 1, N3 = 30, IT = -30;
  EPS = .1;
  if(idebug)
    for(int i = 0; i < 2; i++) {
      cout << " i,A,PL0,AMX,AMN "
           << i << "  "
           << A[i] << "  "
           << PL0[i] << "  "
           << AMX[i] << "  "
           << AMN[i] << endl;
  }
  fumiliSK(2, S, N1, N2, N3,
           EPS, IT, A, PL0,
           AMX, AMN, R, SIGMA,
           SGZ, akappa, VL,
           1, constraints);
  if(akappa < EPS) {
    cout << " Constrained Fit: Fitted Values and their Errors " << endl;
    for(int i = 0;i < 2;i++) {
      cout << i << "  " << A[i] << " +/-  " << SIGMA[i] << endl;
    }
  }

  return 0;
}

int SGZ(int m, double &S, double A[], double PL[], double G[], double Z[]) {
  S = 0.;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j <= i; j++)
      Z[(i+1)*i/2 + j] = .0;
    G[i] = .0;
  }
  if(idebug) cout << " SGZ, m                   = " << m << endl;
  if(idebug) cout << " SGZ, Parameters at entry = " << endl;
  for(int i = 0; i < m;i++) {
    if(idebug)
      cout << " SGZ, i,A[i],PL[i]  = " << i << "  " << A[i] << "  " << PL[i] << endl;
  }
  for (int iev = 0; iev < 100000; iev++) {
    double CH = 1. + A[0]*data[iev][0] + A[1]*data[iev][1];
    double ZN = 1. + 0.5*A[0] + 0.5*A[1];
    double pdf = CH/ZN;//pdf
    // Derivatives
    double df[2];
    for(int i = 0; i < 2;i++) {
      double dch = data[iev][i];
      double dzn = 0.5;
      double ys = (dch - pdf*dzn)/ZN;
      df[i] = ys/pdf;
    }
    S -= log(pdf);
    int il = 0;
    for (int  ip = 0; ip < m; ip++) {
      // Calculation of gradient
      if (PL[ip] > .0) {
        G[ip] = G[ip] - df[ip];
        //if(idebug)cout << " G[ip],Z[il]: " << G[ip] << "  " << endl;
        for (int jp = 0; jp <= ip; jp++) {
          if (PL[jp] > .0) {
            Z[il] = Z[il] + df[ip]*df[jp];
            //if(idebug)cout << Z[il]  << "  " << endl;
            il = il + 1;
          }
        }
      }
    }
  }
  return 2;
}

void constraints(int M, double A[], double * psis, double ** dpsis) {
  /* First - cleaning derivatives */
  for(int j = 0;j < M;j++) {
    if(idebug) cout << " j, A[j] " <<j << "  " << A[j] << endl;
  }
  for(int i = 0; i < 1;i++) {
    for(int j = 0; j < M;j++) {
      dpsis[i][j] = .0;
    }
  }
  psis[0] = A[0] + A[1] - 1.3;
  dpsis[0][0] = 1.;
  dpsis[0][1] = 1.;
  if(idebug)
    for(int i = 0; i < 1; i++) {
      cout << " i, psis[i], derivatives.. = " << i  << "  " << psis[i] << "  " ;
      for(int j = 0; j < M; j++) {
        cout << dpsis[i][j] << "  " ;
      }
      cout << endl;
      cout << "*******************************" << endl;
    }
}
