#include "ublas.hpp"

//------------------------------------------------------------------------------
extern "C" {
  void umfpack_di_defaults(double *);
  int umfpack_di_numeric
  (
   const int Ap [ ],
   const int Ai [ ],
   const double Ax [ ],
   void *Symbolic,
   void **Numeric,
   const double Control [UMFPACK_CONTROL],
   double Info [UMFPACK_INFO]
   ) ;
  int umfpack_di_symbolic
  (
   int n_row,
   int n_col,
   const int Ap [ ],
   const int Ai [ ],
   const double Ax [ ],
   void **Symbolic,
   const double Control [UMFPACK_CONTROL],
   double Info [UMFPACK_INFO]
   ) ;
  int umfpack_di_solve
  (
   int sys,
   const int Ap [ ],
   const int Ai [ ],
   const double Ax [ ],
   double X [ ],
   const double B [ ],
   void *Numeric,
   const double Control [UMFPACK_CONTROL],
   double Info [UMFPACK_INFO]
   ) ;
  void umfpack_di_free_symbolic( void **Symbolic ) ;
  void umfpack_di_free_numeric( void **Numeric) ;
};
//------------------------------------------------------------------------------
unyque::DVector unyque::umfpackSolve(const unyque::SparseMatrix &K,
				     const unyque::DVector &rhs) {

  int size, counter, col;
  int *Kcol, *Krow;
  double *Kval, *b, *x;
  double *Info, *Control;
  void *Symbolic, *Numeric ;

  size = (int) K.size1();

  // Set UMFPACK control
  Info = new double [UMFPACK_INFO];
  Control = new double [UMFPACK_CONTROL];
  umfpack_di_defaults (Control);
  Control [UMFPACK_DENSE_ROW] = 0.4;
  Control [UMFPACK_DENSE_COL] = 0.4;

  // Get sparse matrix data in compressed column format
  Kcol = new int[size + 1];
  Krow = new int[K.nnz()];
  Kval = new double[K.nnz()];
  Kcol[0] = 0;
  counter = 0;
  col = 1;

  for (col_iter_t i2 = K.begin2(); i2 != K.end2(); ++i2, ++col) {
    for (row_iter_t i1 = i2.begin(); i1 != i2.end(); ++i1, ++counter) {
      Krow[counter] = i1.index1();
      Kval[counter] = *i1;
    }
    Kcol[col] = counter;
  }

  // Perform column pre-ordering and symbolic factorization
  (void) umfpack_di_symbolic (size, size, Kcol, Krow, Kval, &Symbolic, Control,
  			      Info);

  // Perform numeric factorization
  (void) umfpack_di_numeric (Kcol, Krow, Kval, Symbolic, &Numeric, Control,
  			     Info);

  // Initialize RHS
  b = new double[size];
  for (int i = 0; i < size; i++)
    b[i] = rhs(i);

  // Solve linear system
  x = new double[size];
  (void) umfpack_di_solve (UMFPACK_A, Kcol, Krow, Kval, x, b, Numeric, Control,
  			   Info);

  // Copy back result
  unyque::DVector result(size);
  for (int i = 0; i < size; i++)
    result(i) = x[i];

  // Check residual
  unyque::DVector residual(rhs);
  axpy_prod(K, -result, residual, false);
  if (ublas::norm_2(residual) > 1) {
    cout << "maxnorm of residual: " << ublas::norm_2(residual) << endl;
    exit(0);
  }

  delete [] Kcol;
  delete [] Krow;
  delete [] Kval;
  delete [] b;
  delete [] x;
  delete [] Info;
  delete [] Control;
  umfpack_di_free_symbolic (&Symbolic) ;
  umfpack_di_free_numeric (&Numeric) ;

  return result;

}
