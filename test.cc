
#include "Fumili.h"
#include <string>
// #include <vector>
// #include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int dataSize = 100000;  //just like this for that simple test example
const int parNum = 2;   //number of parameters (as far as I get)
double data[dataSize][parNum];

//TODO: В коде сплошные "магические числа".
//понять, что откуда взялось и исправить.

//TODO: Задокументировать основные действия и функции
//задокументировать имена параметров

int SGZ(int m, double &S, double A[], double PL[], double G[], double Z[]);  //внутренняя функция, используемая минимизатором Фумили. Задается пользователем
void constraints(int M, double A[], double * psis, double ** dpsis);    //Функция ограничений. Также задается пользователем по определённым правилам.
void getDebugInfo (double A[], double PL0[], double AMX[], double AMN[]);   //print out debug info
void setTestPars(double A[], double PL0[], double AMX[], double AMN[],
    int *_N1, int *_N2, int *_N3, int *_IT, double * _EPS); //set the same initial parameters for both tests we do
void fitInfo(string fitType, double akappa, double EPS, double A[], double SIGMA[]);
void runTestFit(int count, string testType, int idebug);

int main(int arc, char ** argv) {
    const int idebug = 1;   //debug is on

    //read input data
    ifstream ifs("unif.dat");
    for(int i = 0; i < dataSize; i++)
    ifs >> data[i][0] >> data[i][1];
    ifs.close();

    //run tests
    runTestFit(1, "UnConstrained Fit", idebug);
    runTestFit(2, "Constrained Fit", idebug);
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
  for (int iev = 0; iev < dataSize; iev++) {
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

void runTestFit(int count, string testType, int idebug) {
    /*  Define parameters:
    //AMN - array of minimal values
    //AMX - array of maximum values
    //A - array with start parameters values
    //SIGMA - error on second derivatives
    TODO: to be continued...
    */
   cout<<string(20,'*')<<"TEST"<<count<<": "<<testType<<string(20,'*')<<endl;
  //define parameters
  double AMN[parNum], AMX[parNum], PL0[parNum], A[parNum], R[parNum], SIGMA[parNum], VL[parNum*parNum], akappa, S;
  int N1, N2, N3, IT;
  double EPS;
  
  //set parameters (and print them out if debug)
  setTestPars (A, PL0, AMX, AMN, &N1, &N2, &N3, &IT, &EPS);
  if(idebug) {
    getDebugInfo (A, PL0, AMX, AMN);
  }
  //make fit
  if (testType == "UnConstrained Fit") {
      fumiliSK(parNum, S, N1, N2, N3, EPS, IT, A, PL0, AMX, AMN, R, SIGMA, SGZ, akappa, VL);
  } else {
      fumiliSK(parNum, S, N1, N2, N3, EPS, IT, A, PL0, AMX, AMN, R, SIGMA, SGZ, akappa, VL, 1, constraints);
  }
  //print out results
  fitInfo(testType, akappa, EPS, A, SIGMA);
  cout<<endl;
}

void setTestPars(double A[], double PL0[], double AMX[], double AMN[],
    int * N1, int * N2, int * N3, int * IT, double * EPS) {
  for(int i = 0; i < parNum; i++) {
    AMN[i] = -10.;
    AMX[i] = 10.;
    PL0[i] = .1;
    A[i] = .0;
  }
  *N1 = 1, *N2 = 1, *N3 = 30, *IT = -30;
  *EPS = .1;
}

void fitInfo(string fitType, double akappa, double EPS, double A[], double SIGMA[]){
    if(akappa < EPS) {
    cout << fitType << ": Fitted Values and their Errors " << endl;
    for(int i = 0; i < parNum; i++) {
      cout << i << "  " << A[i] << " +/-  " << SIGMA[i] << endl;
    }
  }
}

void getDebugInfo (double A[], double PL0[], double AMX[], double AMN[]) {
    cout <<"i\tA\tPL0\tAMX\tAMN "<<endl;
    for(int i = 0; i < parNum; i++) {
      cout << i << "\t"
           << A[i] << "\t"
           << PL0[i] << "\t"
           << AMX[i] << "\t"
           << AMN[i] << endl;
    }
}
