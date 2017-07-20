
#pragma once

#include <cstdlib>

extern int idebug;

int fumiliSK( int M_User, double & S_User, int N1, int N2, int N3, double EPS,
              int IT, double A_User[], double PL0_User[],
              double AMX_User[], double AMN_User[], double R_User[],
              double SIGMA_User[],
              int (*SGZ)( int M, double & S_User,  double A_User[],
                          double PL_User[], double G_User[], double Z_User[] ),
              double & AKAPPA, double VL_User[], int nc = 0,
              void (*get_psis_and_derivatives)( int M, double A[],
                                                double * psis, double ** dpsis ) = NULL );

int fumiliSK( int M_User, double & S_User, int N1, int N2, int N3, double EPS,
              int IT, double A_User[], double PL0_User[],
              double AMX_User[], double AMN_User[], double R_User[],
              double SIGMA_User[],
              int (*SGZ)( int M, double & S_User,  double A_User[],
                          double PL_User[], double G_User[], double Z_User[], void *),
              double & AKAPPA, double VL_User[], void * data, int nc = 0,
              void (*get_psis_and_derivatives)( int M, double A[],
                                                double * psis, double ** dpsis ) = NULL );
