
#pragma once

/** TODO const */
class FSqMtrx {
private:
  int N;
  double * M;
public:
  FSqMtrx(int n): N(n), M(new double[N*N]) { ; }
  FSqMtrx(int n1, int n2);
  FSqMtrx(int n, double * m);
  FSqMtrx(FSqMtrx & other);
  FSqMtrx & operator= (FSqMtrx & rhs);
  ~FSqMtrx() { delete [] M; }
  double * operator[] (int i) { return &M[i*N]; }
/**
  inplace inverse n x n matrix A.
  returns:
    ret = 0 on success
    ret < 0 illegal argument value
    ret > 0 singular matrix
*/
  static int Invert(const int N, double * const M);
  static double Determinant(const int N, const double * const M);
  static void Print(const int N, const double * const M);
  int Invert() { return Invert(N, M); }
  double Determinant() { return Determinant(N, M); }
  void Print() { Print(N, M); }
};
