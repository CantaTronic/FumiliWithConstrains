
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include "mtrx_inv.h"
#include "ind.h"
#include "mconvd.h"
#include "Fumili.h"
#define PARAMETER_SIZE 10
using namespace std;

int idebug;

template<class T>
T * CopyArray(T * dst, T * src, size_t n) {
  memcpy(dst, src, n*sizeof(T));
  return dst;
}
void Calculate_Full_Matrix_For_Free_Parameters
    (int M, double PL[], double Z[], double VL[]) {
  /* First - Reset to zero VL : idea is to have zero values for fixed parameters */
  for (int ip = 0; ip < M; ip++)
    for (int jp = 0; jp <=ip; jp++)
      VL[ind(ip,jp)] = .0;
  int il = 0;
  for (int ip = 0; ip < M; ip++) {
    if (PL[ip] <= .0) continue;
    for (int jp = 0; jp <= ip; jp++) {
      if (PL[jp] <= .0) continue;
      VL[ind(ip,jp)] = Z[il];
      il++;
    }
  }
}
void Calculate_Full_Matrix_For_Substituted_Parameters
    (int nf, int nc, double ** SM, double VL[]) {
  /* Calculation of VL, SIGMA, R for Substituted_Parameters, here
  RV - diagonal terms of original Z for Substituted_Parameters,
  coming from SGZ and then from rearranged Z */
  /* Calculation of VL - full error matrix, first off diagonal terms */
  for(int k = nf; k < nf+nc; k++) { // index of row
    for(int l = 0; l < nf; l++) { // index of column
      VL[ind(k,l)] = .0;
      for(int i = 0; i < nf; i++) // index of column
        VL[ind(k,l)] += SM[k-nf][i]*VL[ind(i,l)];
    }
    for(int l = nf; l <= k; l++) { // index of column
      VL[ind(k,l)] = .0;
      for(int i = 0; i < nf; i++)
        for(int j = 0; j < nf; j++)
          VL[ind(k,l)] += SM[k-nf][i]*VL[ind(i,j)]*SM[l-nf][j];
    }
  }
}

void print_iteration( int M, double PL0[],
                      double PL[], int NN3, double S,
                      double GT, double AKAPPA, double ALAMBD,
                      double A[], double SIGMA[], double R[], double Z[]) {
  if(!idebug) return;
  cout << " Here " << endl;
  std::cout << " ITER.NO." << setw(3) << NN3
            << ", S=" << setiosflags(ios::scientific) << setw(12) << setprecision(5) << S
            << ", EC =" << setiosflags(ios::scientific) << setw(12) << setprecision(5) << GT
            << ", KAPPA=" << setiosflags(ios::scientific) << setw(9) << setprecision(4) << AKAPPA
            << ", LAMBDA=" << setiosflags(ios::scientific) << setw(9) << setprecision(4) << ALAMBD
            << endl;
  cout << "      PARAMETER      PARAMETER         STANDARD        CORRELATION "
       << endl;
  cout << "        NUMBER        VALUE           DEVIATION          FACTOR "
       << endl;
  for(int I = 0; I < M; I++) {
    if(PL0[I] > .0) {
      if(PL[I] > .0)
        cout << "        " << setw(3) << I << "    "
             << "     " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << A[I]
             << "     " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << SIGMA[I]
             << "     " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << R[I] << endl;
      else if(PL[I] == .0)
        cout << "        " << setw(3) << I << "    "
             << "     " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << A[I]
             << "     " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << SIGMA[I]
             << "     " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << R[I]
             << "  ON BOUNDARY " <<  endl;
      else if(PL[I] >= -1.)
        cout << "        " << setw(3) << I << "    "
             << "     " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << A[I]
             << "     " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << SIGMA[I]
             << "     " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << R[I]
             << "  ON BOUNDARY " <<  endl;
      else
        cout << "         " << setw(3) << I << "    "
             << "         " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << A[I]
             << " INFINITE ERROR ESTIMATED " <<  endl;
    } else {
      cout << "        " << setw(3) << I
           << "         " << setiosflags(ios::scientific) << setw(12) << setprecision(7) << A[I]
           << "     THIS PARAMETER FIXED" << endl;
    }
  }
  // Output correlation coefficients
  cout << endl;
}
void monitoVK(double S, int M, int NN3, int IT, double GT,
              double AKAPPA, double ALAMBD, double A[],
              double SIGMA[], double R[],
              double PL0[], double PL[], int ENDFLG, double Z[]) {
  if(!idebug) return;
  if(IT <= 0) return;
  if((ENDFLG == 0 && !(NN3%IT)) || (NN3 == 0 || ENDFLG == 1))
    print_iteration(M, PL0, PL, NN3, S, GT, AKAPPA, ALAMBD, A, SIGMA, R, Z);
  else if(ENDFLG == -1)
    cout <<  " MINIMIZATION TERMINATED AS NO FURTHER DECREASE IN S IS OBTAINABLE";
  else if(ENDFLG == -2)
    cout <<  " MINIMIZATION TERMINATED AS INFINITE ERRORS ESTIMATED ";
  else if(ENDFLG == -3)
    cout <<  " MINIMIZATION TERMINATED AS ITERATION LIMIT REACHED ";
  else if(ENDFLG == -4)
    cout <<  " MINIMIZATION TERMINATED AS NEGATIVE OR ZERO Y PERSISTED AS LOGARITHMIC ARGUEMENT ";
}

bool Get_RV_SM( int nf, int nc,
                double * psis, double ** dpsis,
                double * RV, double ** SM ) {
  double tolerance = 2.2204e-16;
  FSqMtrx m(nc,nc);
  /* Getting QM - last nc terms in each constraint */
  for(int i = 0; i < nc; i++) {
    for(int j = 0;j < nc; j++) {
      m[i][j] = dpsis[i][nf + j] ;
      if(i == j && fabs(m[i][j]) < tolerance)
        return 0;
    }
  }
  m.Invert();
  /* Getting RV, SM */
  for(int i = 0; i < nc; i++) {
    RV[i] = .0;
    for(int j = 0; j < nc; j++)
      RV[i] -= m[i][j]*psis[j];
      for(int l = 0; l < nf; l++) {
        SM[i][l] =.0;
        for(int j = 0;j < nc; j++)
          SM[i][l] -= m[i][j]*dpsis[j][l];
      }
    }
  return 1;
}
void Get_New_Quadratic_Form(int nf, int nc,
                            double * RV,
                            double ** SM, double & S,
                            double G[], double Z[]) {
  /*
    nf, nc - number of free parameters and constraints
    psis,dpsi - vector of constraint values and matrix of constraint derivatives  over all parameters
    RV - vector of constraint values timed by inverted matrix <-- for control
    SM - matrix for calculation of errors of substituted parameters
  */
  double G_New[PARAMETER_SIZE], Z_New[PARAMETER_SIZE*(PARAMETER_SIZE + 1)/2];
  /* Now Quadratic_Form */
  double S_New = S;
  for(int k = 0; k < nc; k++) {
    double tt = G[nf+k];
    for(int l = 0; l < nc; l++)
      tt += 0.5*Z[ind((nf+k), (nf+l))]*RV[l]; //From Old G,Z
    S_New += RV[k]*tt;
  }
  /* New Gradient */
  for(int i = 0; i < nf; i++) {
    G_New[i] = G[i] ;
    for(int k = 0; k < nc; k++) {
      G_New[i] += G[nf+k]*SM[k][i] ;
      double tt = Z[ind(nf+k, i)];
      // From Old Gradient timed by derivaties of substituted parameters overfree
      for(int l = 0; l < nc; l++)
        tt += Z[ind(nf+k, nf+l)]*SM[l][i];
        // From Old Z timed by first derivatives of the expansion
        // of substituted parameters over free and constraint values
      G_New[i] += RV[k]*tt;
    }
    for(int j = 0; j <= i; j++) { // Z_New is also symmetrical as Old Z
      Z_New[ind(i, j)] = Z[ind(i,j)];
      for(int k = 0; k < nc; k++) {
        Z_New[ind(i, j)] += Z[ind(i, nf+k)]*SM[k][j]+Z[ind(nf+k, j)]*SM[k][i] ;
        for(int l = 0; l < nc; l++)
          Z_New[ind(i, j)] += Z[ind(nf+l, nf+k)]*SM[l][j]*SM[k][i];
      }
    }
  }
  /* Resetting S,G,Z */
  S = S_New;
  for(int j = 0; j < nf; j++)
    G[j] = G_New[j];
  for(int i = 0; i < nf*(nf+1)/2; i++)
    Z[i]= Z_New[i];
}
void Get_nx(int M_User, int & nx, double PL0_User[],
            int * ind_new_to_A_User) {
  int qty_Z = 0;
  for(int i = 0; i < M_User; i++) // over derivatives in constraint
    if(PL0_User[i] > .0)
      qty_Z++;
  nx = M_User-qty_Z;
  int qty_nx = 0;
  for(int i = 0; i < M_User; i++) { // over derivatives in constraint
    if(PL0_User[i] == .0) {
      ind_new_to_A_User[M_User-nx+qty_nx] = i;
      qty_nx++;
    }
  }
}
void Fill_In_Qty_nc(int nc, int n, double PL0[], double ** dpsis, int * Qty_nc) {
  /* Setting to zero */
  for(int i = 0; i < nc; i++)
    Qty_nc[i] = 0;
  for(int i = 0; i < nc; i++)
    for(int j = 0; j < n; j++)
      if(dpsis[i][j] != .0 && PL0[j] > .0)
        Qty_nc[i]++;
  for(int i = 0; i < nc; i++)
    if(!Qty_nc[i]) {
      cout << " fumiliSK: Fill_In_Qty_nc : Constraint number "
           << i << " has zero parameter derivatives, terminated" << endl;
      terminate();
    }
}
void Fill_In_Numbers(int nc, int n, double PL0[], double ** dpsis, int ** Numbers) {
  for(int i = 0; i < nc; i++) {
    int jj = 0;
    for(int j = 0; j < n; j++) {
      if(dpsis[i][j] != .0 && PL0[j] > .0) {
        Numbers[i][jj] = j;
        jj++;
      }
    }
  }
}
void Shift_Indices(int nc, int * Qty_nc, int ** Numbers, int * Init_Index, int * Figures) {
  // Updating indices
  bool Break = true;
  int ii = nc-1;
  while(Break && ii >= 0) {
    if(!(Init_Index[ii] == Qty_nc[ii]-1)) {
      Init_Index[ii]++;
      Break = false;
    } else {
      Init_Index[ii] = 0;
      ii--;
    }
  }
  for(int i = 0; i < nc;i++)
    Figures[i] = Numbers[i][Init_Index[i]];
}
bool Check_Not_Repeatedness(int nc, int * Figures) {
  if(nc > 1) {
    for(int i = 0; i < nc - 1;i++)
      for(int j = i + 1; j < nc ;j++)
        if(Figures[i] == Figures[j])
          return false;
    return true;
  }
  if(nc == 1) return true;
  return false;
}
void Do_Indices_For_Substituted_Parameters
    ( int M, int nc, int nx,
      double PL0_User[],
      int * Figures,
      int * ind_new_to_A_User, int * ind_new_to_Z_User ) {
  for(int qty_nc = 0; qty_nc < nc; qty_nc++) { //over constraint
    int i_nc = Figures[qty_nc];
    int qty_Z = 0;
    for(int i = 0; i < M; i++) { // over parameters
      if(PL0_User[i] <= .0) continue;
      if(i == i_nc) {
        ind_new_to_A_User[M - nx - nc + qty_nc] = i_nc;
        ind_new_to_Z_User[M - nx - nc + qty_nc] = qty_Z;
      }
      qty_Z++;
    }
  }
}

void Do_Indices_For_Free_Parameters(int M_User, int & qty_nf, int nc, int nx,
                                    double PL0_User[],
                                    int * Figures,
                                    int * ind_new_to_A_User, int * ind_new_to_Z_User) {
  qty_nf = 0; // Running Qty's
  int qty_Z = 0;
  for(int i = 0; i < M_User; i++) { // over parameters
    if(PL0_User[i] <= .0) continue;
    bool Free = true;
    for(int qty_nc = 0; qty_nc < nc; qty_nc++) //over constraint
      Free = Free && (Figures[qty_nc] != i);
    if(Free) {
      ind_new_to_A_User[qty_nf] = i;
      ind_new_to_Z_User[qty_nf] = qty_Z; // change 02.09.16
      qty_nf++;
    }
    qty_Z++;
  }
}
void Sort(int nc, int * Original, int * Sorted) {
  for(int i = 0; i < nc; i++)
    Sorted[i] = Original[i];
  for(int i = 0; i < nc-1; i++)
    for(int j = i+1; j < nc ; j++)
      if(Sorted[i] > Sorted[j]) {
        int temp = Sorted[i];
        Sorted[i] = Sorted[j];
        Sorted[j] = temp;
      }
}
void SelectionDifferentCombinations
    ( int & Qty_Selected, int nc,
      int ** All_Selected_Combinations, int **& SelectedIndices ) {
  // When nc > 1 :
  // At the entry we have All_Selected_Combinations,1-st Qty_Selected are eligible
  // for further selection, because some of them are just permutations;
  // from them we should select really different combinations,put them
  // into SelectedIndices and their Qty into Qty_Selected <-- so it is changed!
  // Qty_Selected2 - Qty of entries after 2-nd(final) selection
  // When nc > 1 :
  int Sorted[Qty_Selected][nc],Qty_Selected2,IndicesOfTaken[Qty_Selected];
  for(int i = 0; i <  Qty_Selected;i++)
    Sort(nc,All_Selected_Combinations[i],Sorted[i]);
  IndicesOfTaken[0] = 0;
  Qty_Selected2 = 1;
  for(int i = Qty_Selected2; i < Qty_Selected; i++) { // Over All_Selected_Combinations
    bool Adding_All = false;
    for(int j = 0; j < Qty_Selected2; j++) { // Over temp
      bool Adding = true;
      for(int k = 0; k <  nc ;k++)
        Adding = Adding && (Sorted[i][k] == Sorted[IndicesOfTaken[j]][k]);
      Adding_All = Adding_All || Adding;
    }
    if(!Adding_All) {
      IndicesOfTaken[Qty_Selected2] = i;
      Qty_Selected2++;
    }
  }
  // Filling..
  Qty_Selected = Qty_Selected2;
  SelectedIndices = new int *[Qty_Selected];
  for(int i = 0; i < Qty_Selected; i++) { // Over SelectedIndices
    SelectedIndices[i] = new int[nc];
    for(int k = 0; k < nc; k++)
      SelectedIndices[i][k] = All_Selected_Combinations[IndicesOfTaken[i]][k];
  }
}
bool Get_Not_Repeated_Combination
    ( int nc, int * Qty_nc, int ** Numbers,
      int * Figures, int & Qty_Selected, int **& SelectedIndices ) {
  if(nc == 1) {
    Qty_Selected = Qty_nc[0];
    // Allocation of SelectedIndices
    SelectedIndices = new int*[Qty_Selected];
    for(int i = 0; i < Qty_Selected; i++) { // Over SelectedIndices
      SelectedIndices[i] = new int[nc];
      for(int k = 0; k < nc; k++)
        SelectedIndices[i][k] = Numbers[k][i];
    }
    return true;
  }
  int * Init_Index = new int[nc];
  for(int i = 0; i < nc; i++) {
    Init_Index[i] = 0;
    Figures[i] = Numbers[i][Init_Index[i]];
  }
  int Total_Number = 1;
  for(int i = 0; i < nc; i++)
    Total_Number *= Qty_nc[i];
  int ** All_Selected_Combinations = new int*[Total_Number];
  for(int i = 0; i < Total_Number; i++)
    All_Selected_Combinations[i] = new int[nc];
  Qty_Selected = 0;
  for(int i = 0; i < Total_Number; i++) {
    if(Check_Not_Repeatedness(nc,Figures)) {
      for(int j = 0;j < nc;j++)
        All_Selected_Combinations[Qty_Selected][j] = Figures[j];
      Qty_Selected++;
    }
    Shift_Indices(nc,Qty_nc,Numbers,Init_Index,Figures);
  }
  // Before this call Qty_Selected - combinations
  // after 1-st selection : list may contain the same permutated combinations
  // 2-nd selection : leave only one combination from permutated
  SelectionDifferentCombinations (Qty_Selected,nc,
                                  All_Selected_Combinations,
                                  SelectedIndices);
  // Afte this call Qty_Selected changed !
  for(int i = 0; i < Total_Number ;i++)
    delete [] All_Selected_Combinations[i];
  delete [] All_Selected_Combinations ;
  delete [] Init_Index;
  // Here we have Combinations regardless their productions, their number is Qty_Selected
  // now we should take combinations with different productions
  return true;
}
void GetDeterminant(double ** dpsis_User, int nc,
                    int * SelectedIndices, double & DeterminantValue) {
  FSqMtrx m(nc, nc);
  for(int i = 0; i < nc; i++)
    for(int j = 0; j < nc; j++)
      m[i][j] = dpsis_User[i][SelectedIndices[j]];
  DeterminantValue = m.Determinant();
}
void SelectFigures
    ( int Qty_Selected, int ** SelectedIndices,
      double * DeterminantValues, int nc, int * Figures ) {
  // Selection of a list of substituted parameters, having the biggest Detrminant
  int i_Selected = 0;
  double BiggestDeterminant = fabs(DeterminantValues[i_Selected]);
  for(int i = 1; i < Qty_Selected; i++) {
    if(fabs(DeterminantValues[i]) <= BiggestDeterminant) continue;
    i_Selected = i;
    BiggestDeterminant = fabs(DeterminantValues[i_Selected]);
  }
  for(int j = 0; j < nc;j++)
    Figures[j] = SelectedIndices[i_Selected][j];
}
void Do_Figures(int nc, int M_User, double PL0_User[], double ** dpsis_User, int * Figures) {
  int * Qty_nc;
  Qty_nc = new int[nc];
  Fill_In_Qty_nc(nc, M_User, PL0_User, dpsis_User, Qty_nc);
  // to find qty of parameters, not fixed by user and having nonzero derivatives
  int ** Numbers;
  Numbers = new int*[nc];
  for(int i = 0; i < nc;i++)
    Numbers[i] = new int[Qty_nc[i]];
  Fill_In_Numbers(nc, M_User, PL0_User, dpsis_User, Numbers); // Nomera parametrov, for which nonzero derivatives
  int Qty_Selected, ** SelectedIndices; // Qty_Selected - number
  // of combinations after
  // two selections :
  //1-st : to produce list of parameters, when each entry
  // does not contain the same parameter more than one time;
  //2 -nd : to produce list of parameters, when each entry
  // contains unique combination of parameters;
  // Array SelectedIndices is being created in SelectionDifferentCombinations,
  // called by Get_Not_Repeated_Combination and has Qty_Selected entries
  if( !Get_Not_Repeated_Combination(nc, Qty_nc, Numbers, Figures,
                                    Qty_Selected, SelectedIndices) ) {
    cout << " Do_Figures: Problem is poorly defined: constraints are linearly dependent, terminated " << endl;
    terminate();
  }
  double * DeterminantValues = new double[Qty_Selected];
  if(nc > 1)
    for(int j = 0; j < Qty_Selected;j++)
      GetDeterminant(dpsis_User,nc,SelectedIndices[j],DeterminantValues[j]);
  if(nc == 1)
    for(int j = 0; j < Qty_Selected;j++)
      DeterminantValues[j] = dpsis_User[0][SelectedIndices[j][0]];
  SelectFigures(Qty_Selected,SelectedIndices,DeterminantValues,nc,Figures);
  for(int j = 0; j < Qty_Selected;j++)
    delete SelectedIndices[j];
  delete SelectedIndices;
  delete DeterminantValues;
  for(int i = 0; i < nc;i++)
    delete Numbers[i];
  delete Qty_nc;
  delete Numbers;
}
void Get_Index_Functions (int M, int nf, int nc, int nx,
                          double A_User[], double PL0_User[],
                          double * psis, double ** dpsis_User,
                          int * ind_new_to_A_User,
                          int * ind_new_to_Z_User) {
  /*
    Meaning of indices :
    M - user's number of parameters
    nf - qty of free parameters
    nc - qty of of constrained parameters
    nx - qty of   parameters fixed by user
    PL0_User - Parameter Limiters - parameter steps if not zero, if zero, then
              parameter considered as fixed
    dpsis_User - derivatives of User constraints over parameters
    ind_new_to_A_User - gives for new A( which is rearranged A_User and has nf + nc  + nx dimensionality)
                        indices  :
                        A[i] = A_User[ind_new_to_A_User[i]], i is running from 0
                                to M - 1
    A_User - parameters given by User : they are separated in three groups :
            free parameters - i.e. nf fitted parameters they are first in rearranged parameter array A
            constrained parameters  - nc parameters which differentials are being  expressed
                                      by differentials of free parameters and
                                      as result we get new quadratic form. They are next in
                                      array A
            fixed parameters - nx parameters  which are fixed by User.They are last in array A.
    A      - rearranged A_User :
            1-st nf values are free parameters,i.e. taking part in a reduced
                  quadratic form
            next nc values are substituted parameters
            last nx values are fixed parameters
    ind_new_to_Z_User      - gives for new Z (which is rearranged Z_User and has nf + nc dimensionality)
                        indices :
                        ind_new_to_Z_User[i] - nomer parametra in matrix Z, packed by User.
                        During the packing index of Z skips fixed parameters .
                        It means that Z[ind(i,j)] =
                        Z_User[ind(ind_new_to_Z_User[i],ind_new_to_Z_User[j])],
                        i,j -
                          Remeber that User gives Z already packed
                          according to PL0_User
  */
  int * Figures;
  int M_User = M;
  Figures = new int[nc];
  Do_Figures(nc, M_User,PL0_User,dpsis_User,Figures);
  Do_Indices_For_Substituted_Parameters(M_User,nc,nx,PL0_User,Figures,
                                        ind_new_to_A_User,ind_new_to_Z_User);
  Do_Indices_For_Free_Parameters (M_User,nf,nc,nx,PL0_User,Figures,
                                  ind_new_to_A_User,ind_new_to_Z_User);
  // Index for Z - we are ready for its creation
  // exit if qty_ind_nc != nc
  if (nf + nc + nx != M_User) {
    cout << " terminated due to nf + nc + nx != M_User "
         << endl;
    terminate();
  }
  delete [] Figures;
}
void Get_SIGMA_R_VL_User (int M_User, int nf, int nc,
                          int * ind_new_to_A_User,
                          double Z0_T[], double VL[],
                          double SIGMA_User[], double R_User[], double VL_User[] ) {
  /* Cleaning VL_User */
  for(int i = 0; i < M_User; i++) {
    SIGMA_User[i] = .0;
    R_User[i] = .0;
    for(int j = 0; j <= i; j++)
      VL_User[ind(i,j)] = .0;
  }
  /* Getting R_User 24.07.14 Only for free parameters */
  for(int i = 0; i < nf+nc; i++) {
    for(int j = 0; j <= i; j++)
      VL_User[ind(ind_new_to_A_User[i],ind_new_to_A_User[j])] = VL[ind(i,j)];
    R_User[ind_new_to_A_User[i]] = Z0_T[ind(i,i)];
  }
  for(int i = 0; i < M_User; i++) {
    SIGMA_User[i] = sqrt(VL_User[ind(i,i)]);
    R_User[i] = R_User[i]*VL_User[ind(i,i)];
  }
}
void Get_PL0_PL_User (int nf,int *ind_new_to_A_User,
                      double PL0[],double PL[],
                      double PL0_User[], double PL_User[]) {
  for(int i = 0; i < nf; i++) {
    PL0_User[ind_new_to_A_User[i]] = PL0[i];
    PL_User[ind_new_to_A_User[i]] = PL[i];
  }
}
void Rearrange_User_Indices_For_A_AMNX_PL0
    ( int M, int * ind_new_to_A_User,
      double A_User[], double A[],
      double AMX_User[], double AMX[],
      double AMN_User[], double AMN[],
      double PL0_User[], double PL0[] ) {
   /* Maxim : if no constraints then everything is set to previous */
  for(int j = 0; j < M; j++) { // over constraints
    A[j] = A_User[ind_new_to_A_User[j]];
    AMX[j] = AMX_User[ind_new_to_A_User[j]];
    AMN[j] = AMN_User[ind_new_to_A_User[j]];
    PL0[j] = PL0_User[ind_new_to_A_User[j]];
  }
}
bool Rearrange_S_G_Z (double A_User[],
                      int nf, int nc, int nx,
                      int * ind_new_to_A_User, int * ind_new_to_Z_User,
                      double * psis, double ** dpsis_User, double ** dpsis,
                      double * RV, double ** SM,
                      double & S,
                      double G_User[], double G[],
                      double Z_User[], double Z[]) {
  for(int j = 0; j < nc; j++) // over constraints
    for(int i = 0; i < nf + nc; i++)
      dpsis[j][i] = dpsis_User[j][ind_new_to_A_User[i]];
  if(Get_RV_SM(nf,nc,psis,dpsis,RV,SM)) {
    /* Take out possible fixed parameter from G_User,Z_User */
    for(int i = 0; i < nf + nc +nx;i++)
      G[i] = G_User[ind_new_to_A_User[i]];
    for(int i = 0; i < nf + nc;i++)
      for(int j = 0; j <= i;j++)
        Z[ind(i,j)] = Z_User[ind(ind_new_to_Z_User[i],ind_new_to_Z_User[j])];
    Get_New_Quadratic_Form(nf, nc, RV, SM, S, G, Z); // S,G,Z - now new
    return 1;
  }
  return 0;
}
void Get_New_A_User
    ( int nf, int nc, double * RV, double ** SM, double A[], double DA[],
      int * ind_new_to_User, double A_User[] ) {
/*
  This is a change, compared with D510:(24.07.14)
  it calculates new values of parameters, based on increments DA.
  If nc = 0 then then nf is equal to user number of parameters.
  If not then we calculate first,increments for substituted parameter,then new values of
  all nonfixed parameters
*/
  if(nc) {
    for(int i = 0; i < nc; i++) {
      /* increment */
      double da = RV[i];
      for(int k = 0; k < nf;k++)
        da += SM[i][k]*DA[k];
      DA[nf + i] = da;
    }
    for(int i = 0; i < nf + nc; i++) {
      A[i] = A[i] + DA[i];
      A_User[ind_new_to_User[i]] = A[i];
    }
  } else {
    for(int i = 0; i < nf + nc;i++)
      A[i] = A[i] + DA[i];
  }
}

namespace SGZ_data {
  void * p_data = 0;
  int (*p_SGZ)( int, double &,  double [],
              double [], double [], double [], void * ) = 0;
  int SGZ(int M, double &S_User, double A_User[],
          double PL_User[], double G_User[], double Z_User[]) {
    return p_SGZ(M, S_User, A_User, PL_User, G_User, Z_User, p_data);
  }
}
int fumiliSK( int M_User, double &S_User, int N1, int N2, int N3, double EPS,
              int IT, double A_User[], double PL0_User[],
              double AMX_User[], double AMN_User[], double R_User[],
              double SIGMA_User[],
              int (*SGZ)( int M, double &S_User,  double A_User[],
                          double PL_User[], double G_User[], double Z_User[], void * ),
              double &AKAPPA, double VL_User[], void * data, int nc,
              void (*get_psis_and_derivatives)( int M, double A[],
                                                double *psis, double **dpsis ) ) {
  SGZ_data::p_data = data;
  SGZ_data::p_SGZ = SGZ;
  return fumiliSK(M_User, S_User, N1, N2, N3, EPS,
                  IT, A_User, PL0_User,
                  AMX_User, AMN_User, R_User,
                  SIGMA_User,
                  SGZ_data::SGZ,
                  AKAPPA, VL_User, nc,
                  get_psis_and_derivatives);

}

//ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
//
//         fumiliSK - version for Fimulator
//
//ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
int fumiliSK( int M_User, double & S_User, int N1, int N2, int N3, double EPS,
              int IT, double A_User[], double PL0_User[],
              double AMX_User[], double AMN_User[], double R_User[],
              double SIGMA_User[],
              int (*SGZ)( int M, double & S_User,  double A_User[],
                          double PL_User[], double G_User[], double Z_User[] ),
              double & AKAPPA, double VL_User[], int nc,
              void (*get_psis_and_derivatives)( int M, double A[],
                                                double * psis, double ** dpsis ) )
/* Main phylosophy - code is transparent to having constrains, i.e. - if no constraints all will go as before

  This is modified version of fumiliSK, supposed to be used in a case of
  a fit with constraints :
  M_User- total number of all the parameters : free, substituted  and fixed.
  nf  - qty of free parameters
  nc  - qty of constraints or or a substituted parameters
  nx  - qty of  a fixed parameters
  psis - values of constraint equations
  dpsis  - array of constraint derivatives

  07.04.14 : Get RV,SM according to formula
  \Delta X2 = RV + SM*\Delta X1(NIM A345,346,1994) <-- Keep attention !
  There were errors in there
  see also NIM A314,578,1992,

*/
{
   /*
     A_User - parameters given by User : we rearrange it
     into A which is  array of three groups parameters:
     free  - i.e. nf fitted parameters they are first in  array A
     constrained parameters  - nc parameters which differentials are being  expressed
     over by differentials of free parameters and
     as result we get new quadratic form. They are next in
     array A
     fixed parameters - nx parameters  which are fixed by User.They are last in array A.
     A      - rearranged A_User :
     1-st nf values are free parameters,i.e. taking part in a reduced
     quadratic form
     next nc values are substituted parameters
     last nx values are fixed parameters
   */
   double * RV, * ZV ;
   double S, ** dpsis_User, ** SM,** dpsis;
   double * psis; // Inserted by VK(30.05.16)
   int nf, nx;
   int * ind_new_to_A_User = 0;
/*
  ind_new_to_A_User - gives for new A(having nf  + nc  + nx dimensionality),
  which is rearranged A_User(having nf + nc  + nx dimensionality),
  i.e.  :
  A[i] = A_User[ind_new_to_A_User[i]], i is running from 0
  to M - 1
*/
  int * ind_new_to_Z_User = 0, * Ind_Z = 0;
/*
  gives for new Z (which is rearranged Z_User and has nf + nc dimensionality)
  indices :
  ind_new_to_Z_User[i] - nomer parametra in matrix Z, packed by User.
  During the packing index of Z skips fixed parameters .
  It means that Z[ind(i,j)] =
  Z_User[ind(ind_new_to_Z_User[i],ind_new_to_Z_User[j])],
  i,j -
  Remeber that User gives Z already packed
  according to PL0_User
*/
   int ENDFLG, INDFLG[5], I;
   double A[PARAMETER_SIZE], PL0[PARAMETER_SIZE], PL_User[PARAMETER_SIZE],
          G_User[PARAMETER_SIZE], Z_User[PARAMETER_SIZE*(PARAMETER_SIZE + 1)/2],
          VL[PARAMETER_SIZE*(PARAMETER_SIZE + 1)/2],
          R[PARAMETER_SIZE],SIGMA[PARAMETER_SIZE];
   int M;
   double Z[PARAMETER_SIZE*(PARAMETER_SIZE + 1)/2];
   double Z0[PARAMETER_SIZE*(PARAMETER_SIZE + 1)/2];
   double Z0_T[PARAMETER_SIZE*(PARAMETER_SIZE + 1)/2];
   double G[PARAMETER_SIZE];
   double PL[PARAMETER_SIZE];
   double DA[PARAMETER_SIZE];
   double AMX[PARAMETER_SIZE], AMN[PARAMETER_SIZE];
   double RP = 1.e-15;
   INDFLG[2] = 0;// Why?
   if(nc) { // Keep in mind that there is special logic when nc = M
    ind_new_to_A_User  = new int[M_User];
    for (int i = 0; i < M_User; i++)
      ind_new_to_A_User[i] = -1;
    Get_nx(M_User,nx,PL0_User,ind_new_to_A_User);
    // finds the number of fixed by user parameters and sets pointers to them
    if(nc > M_User-nx) {
      cout << " fumiliSK: terminated due to number of nonfixed parameters less than constraints " << endl;
      terminate();
    }
    nf = M_User - nx - nc;
    ind_new_to_Z_User = new int[M_User - nx];
    for(int i = 0; i < M_User - nx; i++)
      ind_new_to_Z_User[i] = -1 ;
    psis = new  double[nc];
    dpsis_User = new double*[nc];
    for(int i = 0; i < nc;i++)
      dpsis_User[i] = new double[M_User]; // M is not yet changed!!
    get_psis_and_derivatives (M_User,A_User,
                              psis,dpsis_User);
    Get_Index_Functions(M_User,nf,nc,nx,A_User,PL0_User,
                        psis,dpsis_User,
                        ind_new_to_A_User,ind_new_to_Z_User);
    RV = new double[nc];
    ZV = new double[nc] ;
    SM = new double*[nc];
    dpsis = new double*[nc];
    for(int i = 0; i < nc; i++) {
      SM[i] = new double[nf];
      dpsis[i] = new double[nf + nc]; // M is not yet changed!!
    }
    /* Setting PL_User to PL0_User. During iterations PL is being
      filled from dynamically changed PL and before monito PL_User should be set to PL */
    for(int i = 0; i < M_User; i++)
      PL_User[i] =  PL0_User[i];
    Rearrange_User_Indices_For_A_AMNX_PL0
        (M_User, ind_new_to_A_User, A_User, A, AMX_User, AMX, AMN_User, AMN, PL0_User, PL0);
    if(nf == 0) M = nc;
    else M = nf;
  } else {
    M = M_User;
    S = S_User;
    CopyArray(A, A_User, M);
    CopyArray(AMX, AMX_User, M);
    CopyArray(AMN, AMN_User, M);
    CopyArray(PL0, PL0_User, M);
  }
  /* Change dimensionality of problem */
  int NN2, N, FIXFLG,  IFIX1, FI, NN3, NN1, N0;
  double T1;
  NN2 = 0;
  N = M;
  FIXFLG = 0;
  ENDFLG = 0;
  INDFLG[1] = 0;
  IFIX1 = -1;
  FI = 0;
  NN3 = 0;
  for (I = 0; I < N; I++) {
    R[I] = 0.;
    if (EPS > 0.) SIGMA[I] = 0.;
    PL[I] = PL0[I];
  }

  double SP, T, OLDS, GT;
  while(true) {
    // C-----START NEW ITERATION
    NN1 = 1;
    T1 = 1.;

    while(true) {
      S = 0.;
      N0 = 0;
      for (I = 0; I < N; I++) {
        G[I] = 0.;
        if (PL0[I] > .0) {
          N0 = N0 + 1;
          if (PL[I] > .0) PL0[I] = PL[I];
        }
      }
      int NN0;
      NN0 = N0 * (N0 + 1) / 2;
      if (NN0 >= 1)
        for (I = 0; I < NN0; I++)
          Z[I] = 0.;
      INDFLG[0] = -1; // My change : 05.05.14
      int ijkl;
      if(nc) {
        ijkl = SGZ(M_User, S_User, A_User, PL0_User, G_User, Z_User);
        if(!ijkl) return -1;
        S = S_User;
        get_psis_and_derivatives(nf + nc + nx, A_User, psis, dpsis_User);
        if( Rearrange_S_G_Z(A_User, nf, nc, nx, ind_new_to_A_User, ind_new_to_Z_User,
                            psis, dpsis_User, dpsis, RV, SM, S, G_User, G, Z_User, Z) ) {
          /* Saving terms for calculation of R for Substituted_Parameters */
          for(int k = 0; k <  nc; k++) // index of row
            ZV[k] = Z[ind(nf+ k,nf + k)];
        } else {
          return 11;
        }
      } else {
        ijkl = SGZ(M_User, S, A, PL0, G, Z);
        S_User = S; // Change 16.09.14
      }
      /* Here we should have correct A_User! */
      if (!ijkl) return 10;
      SP = RP * fabs(S);
      if (NN0 >= 1)
        for (I = 0; I < NN0; I++)
          Z0[I] = Z[I];
      if (NN3 > 0 && NN1 <= N1) {
        T = 2.*(S - OLDS - GT);
        if (INDFLG[0] == -1) {
          if (fabs(S - OLDS) <= SP && -GT <= SP) goto L19;
          if (0.59*T < -GT) goto L19;
          T = -GT / T;
          if (T < 0.25) T = 0.25;
        } else {
          T = 0.25;
        }
        GT = GT * T;
        T1 = T1 * T;
        NN2 = 0;
        for (I = 0; I < N; I++) {
          if (PL[I] <= 0.) continue;
          A[I] = A[I] - DA[I];
          PL[I] = PL[I] * T;
          DA[I] = DA[I] * T;
          //  A[I] = A[I] + DA[I];
        }
        Get_New_A_User(M, nc, RV, SM, A, DA, ind_new_to_A_User, A_User);
        NN1 = NN1 + 1;
      } else {
        break;
      }
    }

L19:
    if (INDFLG[0] != -1) {
      ENDFLG = -4;
      goto L85;
    }
    int K1, K2, I1, J, L;
    K1 = 1;
    K2 = 1;
    I1 = 1;
    // Begin block 01
    for (I = 0; I < N; I++) {
      if (PL0[I] <= .0) continue;
      if (PL[I] == 0.) PL[I] = PL0[I];
      if (PL[I] > .0) {
        if ((A[I] >= AMX[I] && G[I] < 0.) ||
            (A[I] <= AMN[I] && G[I] > 0.)) {
              PL[I] = 0.;
              K1 = K1 + I1;
        } else {
          for (J = 0; J <= I; J++) {
            if (PL0[J] > .0) {
              if (PL[J] > .0) {
                Z[K2 -1] = Z0[K1 -1];
                K2 = K2 + 1;
              }
              K1 = K1 + 1;
            }
          }
        }
      } else {
        K1 = K1 + I1;
      }
      I1 = I1 + 1;
    }
    // End block 01
  /* Change as of 05.05.14
    // Below is change of block one as of (05.05.14) see upper, to simplify it in my opinion
    for (I = 0; I < N; I++) {
      if (PL0[I] <= .0) continue;
      if (PL[I] == 0.) PL[I] = PL0[I];
      if ((A[I] >= AMX[I] && G[I] < 0.) ||
          (A[I] <= AMN[I] && G[I] > 0.)) {
        PL[I] = 0.;
      }
    }
  Second -  to pack Z0 into Z by PL
    K2 = 0;
    for (I = 0; I < N; I++) {
      if (PL[I] <= .0) continue;
      for (J = 0; J <= I; J++) {
        if (PL[J] <= .0) continue;
        if(ind_Z0[I] == -1 || ind_Z0[J] == -1)
          terminate();
        Z[K2] = Z0[ind(ind_Z0[I],ind_Z0[J])];
        K2++;
      }
    }
  */
    // Begin block 02
    I1 = 1;
    L  = 1;
    for (I = 0; I < N; I++) {
      if (PL[I] <= .0) continue;
      R[I] = Z[L - 1];
      I1 = I1 + 1;
      L = L + I1;
    }
    int L1, K, IFIX;
    double BI, AIMAX, AMB;
    N0 = I1 - 1;
    for(int i = 0; i < N0*(N0 + 1)/2;i++)
      Z0_T[i] = Z[i];
    mconvd(N0, Z, R, PL, INDFLG);
    if (INDFLG[0] != -1) {
      INDFLG[0] = -1;
      INDFLG[1] = 1;
      FIXFLG = FIXFLG + 1;
      FI = 0;
      goto L19;
    }
    for(int i = 0; i < N0; i++)
      for(int j = 0; j < N0; j++) {
        double tt = .0;
        for(int k = 0 ;k < N0 ;k++)
          tt +=Z0_T[ind(i,k)]* Z[ind(k,j)];
      }
    I1 = 1;
    for (I = 0; I < N; I++) {
      DA[I] = 0.;
      if (PL[I] <= .0) continue;
      L1 = 1;
      for (L = 0; L < N; L++) {
        if (PL[L] <= .0) continue;
          if (I1 <= L1)
            K = L1 * (L1 - 1) / 2 + I1;
          else
            K = I1 * (I1 - 1) / 2 + L1;
          DA[I] = DA[I] - G[L] * Z[K - 1];
          L1 = L1 + 1;
      }
      I1 = I1 + 1;
    }
    // End block 2
    // Begin block 3
    double AFIX, SIGI, AKAP;
    AFIX = 0.;
    IFIX = -1;
    I1 = 1;
    L = I1;
    for (I = 0; I < N; I++) {
      if (PL[I] <= .0) continue;
      SIGI = sqrt(fabs(Z[L - 1]));
      R[I] = R[I] * Z[L - 1];
      if (EPS > .0) SIGMA[I] = SIGI;
      if ((A[I] >= AMX[I] && DA[I] > 0.) ||
          (A[I] <= AMN[I] && DA[I] < .0)) {
        AKAP = fabs(DA[I] / SIGI);
        if (AKAP > AFIX) {
          AFIX = AKAP;
          IFIX = I;
          IFIX1 = I;
        }
      }
      I1 = I1 + 1;
      L = L + I1;
    }
    if (IFIX != -1) {
      PL[IFIX] = -1.;
      FIXFLG = FIXFLG + 1;
      FI = 0;
      goto L19;
    }
  // End block 3
  // Begin block 4
    double ALAMBD, AL, BM, ABI, ABM;
    ALAMBD = 1.;
    AKAPPA = 0.;
    int IMAX;
    IMAX = -1;
    for (I = 0; I < N; I++) {
      if (PL[I] <= .0) continue;
      BM = AMX[I] - A[I];
      ABI = A[I] + PL[I];
      ABM = AMX[I];
      if (DA[I] <= .0) {
        BM = A[I] - AMN[I];
        ABI = A[I] - PL[I];
        ABM = AMN[I];
      }
      BI = PL[I];
      if (BI > BM) {
        BI = BM;
        ABI = ABM;
      }
      if (fabs(DA[I]) > BI) {
        AL = fabs(BI / DA[I]);
        if (ALAMBD > AL) {
          IMAX = I;
          AIMAX = ABI;
          ALAMBD = AL;
        }
      }
      AKAP = fabs(DA[I] / SIGMA[I]);
      if (AKAP > AKAPPA) AKAPPA = AKAP;
    }
    // End block 4
    // Begin block 5
    GT = 0.;
    AMB = 1.e18;
    if (ALAMBD > .0) AMB = 0.25 / ALAMBD;
    for (I = 0; I < N; I++) {
      if (PL[I] <= .0) continue;
      if (NN2 > N2 && fabs(DA[I] / PL[I]) >= AMB) {
        PL[I] = 4.*PL[I];
        T1 = 4.;
      }
      DA[I] = DA[I] * ALAMBD;
      GT = GT + DA[I] * G[I];
    }
    // End block 5
  /* This is my change/addition (24.07.14) to keep compatibility with D510 */
    if (IMAX >= 0)
      DA[IMAX] = AIMAX - A[IMAX];

    // Begin block 6
    if (-GT <= SP && T1 < 1. && ALAMBD < 1.)ENDFLG = -1;
    if (ENDFLG >= 0) {
      if (AKAPPA < fabs(EPS)) {
        if (FIXFLG == 0) {
          ENDFLG = 1;
        } else {
          if (ENDFLG == 0) {
            ENDFLG = 1;
            FIXFLG = 0;
            IFIX1 = -1;
            for (I = 0; I < M; I++) PL[I] = PL0[I];
            INDFLG[1] = 0;
            goto L19;
          } else {
            if (IFIX1 >= 0) {
              FI = FI + 1;
              ENDFLG = 0;
            }
          }
        }
      } else {
        if (FIXFLG != 0) {
          if (FI > FIXFLG) {
            ENDFLG = 1;
            FIXFLG = 0;
            IFIX1 = -1;
            for (I = 0; I < M; I++) PL[I] = PL0[I];
            INDFLG[1] = 0;
            goto L19;
          } else {
            FI = FI + 1;
            ENDFLG = 0;
          }
        } else {
          FI = FI + 1;
          ENDFLG = 0;
        }
      }
    }

L85:
    if (ENDFLG == 0 && NN3 >= N3) ENDFLG = -3;
    if (ENDFLG > 0 && INDFLG[1] > 0) ENDFLG = -2;
    //MONITO (S,M,NN3,IT,EPS,GT,AKAPPA,ALAMBD);
    /*
      07.002.14(VK):
      Change to do prints of iteration
    */
    Calculate_Full_Matrix_For_Free_Parameters(M, PL, Z, VL);
    if(nc) {
      Calculate_Full_Matrix_For_Substituted_Parameters(nf, nc, SM, VL);
      Get_SIGMA_R_VL_User(M_User,nf,nc,ind_new_to_A_User,
                          Z0_T, VL,
                          SIGMA_User, R_User, VL_User);
      /* Getting  of PL0_User and PL_User only for monitoVK */
      Get_PL0_PL_User(nf, ind_new_to_A_User, PL0, PL, PL0_User, PL_User);
    } else {
      /* Getting User values */
      CopyArray(A_User, A, M_User);
      CopyArray(SIGMA_User, SIGMA, M_User);
      CopyArray(R_User, R, M_User);
      CopyArray(PL0_User, PL0, M_User);
      CopyArray(PL_User, PL, M_User);
      CopyArray(VL_User, VL, M_User*(M_User+1)/2);
      CopyArray(Z_User, Z, M_User*(M_User+1)/2);
    }
    // Here is Rearrange new A,SIGMA,R,PL,PL0,Z to users,
    // transforming from new indices to old ones
    if (ENDFLG == 0) {
      if(nc) {
        Get_New_A_User(N,nc,RV,SM,A,DA,ind_new_to_A_User,A_User);
      } else {
        for (I = 0; I < N; I++) A[I] = A[I] + DA[I];
        if (IMAX >= 0) A[IMAX] = AIMAX;
      }
      OLDS = S;
      NN2 = NN2 + 1;
      NN3 = NN3 + 1;
    } else {
      if(nc) {
        delete ind_new_to_A_User;
        delete Ind_Z;
        delete ind_new_to_Z_User;
        for(int i = 0; i < nc; i++) {
          delete dpsis_User[i]; // M is not yet changed!!
          delete SM[i]; // M is not yet changed!!
          delete dpsis[i]; // M is not yet changed!!
        }
        delete dpsis_User;
        delete dpsis;
        delete psis;// Inserted by VK(30.05.16)
        delete SM;
        delete RV;
        delete ZV;
      }
      return ENDFLG - 1;;
    }
    // Transform fitted parameters into User's for SGZ
  }
}





















