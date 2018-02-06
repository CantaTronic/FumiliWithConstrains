
#include <cstdio>
// #include <TMatrixD.h>
#include "mtrx_inv.h"
#include "test_mat.h"
#include <iostream>
#include <string>

using namespace std;

void sepString(string text = string(8,'='));

int main() {
  /*
  const int N = 5;
  double A[] = {
    0.378589, 0.971711, 0.016087, 0.037668, 0.312398,
    0.756377, 0.345708, 0.922947, 0.846671, 0.856103,
    0.732510, 0.108942, 0.476969, 0.398254, 0.507045,
    0.162608, 0.227770, 0.533074, 0.807075, 0.180335,
    0.517006, 0.315992, 0.914848, 0.460825, 0.731980
  };
*/
  const int N = 8;
  sepString("test_mat");
  double * A = boothroyd(N);    //взяли некую матрицу
  printf("\n%+12.8f\n", boothroyd_determinant(N));  //посчитали её детерминант
  sepString();
  
  /*
  printf(" ==== TMatrixD ====\n");
  TMatrixD M(N, N, A);
  FSqMtrx::Print(N, M.GetMatrixArray());
  printf("\n%+12.8f\n", M.Determinant());
  M.Invert();
  FSqMtrx::Print(N, M.GetMatrixArray());
  printf("\n%+12.8f\n", M.Determinant());
  printf(" ==================\n");
*/

  //завели некую матрицу(см. mtrx_inv.h)
#ifndef MTRX_LIB_ROOT
  FSqMtrx m(N, A);
#else
  FSqMtrx m(N, N, A);
#endif
  
  sepString("mtrx_inv");
  m.Print();  //напечатали оригинал
  printf("\n%+12.8f\n", m.Determinant());
  m.Invert();   //инвертировали
  m.Print();  //напечатали инвертированную
  printf("\n%+12.8f\n", m.Determinant()); 
  sepString();
}

void sepString(string text){
  const int len = 25;
  cout<<string(len,'=')<<text<<string(len,'=')<<endl;  
}
