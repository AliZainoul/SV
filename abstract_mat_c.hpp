#ifndef ABSTRACT_MAT_C_H
#define ABSTRACT_MAT_C_H
using namespace std;
#include "full_mat_c.hpp"


class Vector;                   // forward declaration

class AbsMtx                    // base matrix, an abstract class
{
protected:
  int nrows;                    // number of rows in the matrix
  int ncols;
  double **mx;

public:

  virtual Vector operator*(const Vector&) const = 0;  // Matrix-Vector Product

  int CG(Vector& x, const Vector& b, double& eps, unsigned int& iter, int pn=0);
  // Preconditioned Conjugate Gradient method for Ax = b.
  // Usage: A.CG = x; where Ax = b <=> x = A^(-1)b
  // A: symetric positive definite
  // x: on entry: Initial guess; on return: approximate solution
  // b: right side vector
  // eps: on entry:  stopping criterion, epsilon
  //      on return: absolute residual in two-norm for approximate solution
  // iter: on entry:  max number of iterations allowed;
  //       on return: actual number of iterations taken.
  // pn: =0 if no preconditionner, =1 if diag preconditionner, =2 if SSOR;
  // It returns 0 for sucessful return and 1 for breakdowns

  /* virtual ~AbsMtx()
  {
    for (int i = 0; i< nrows; i++) delete[]  mx[i];
    delete[] mx;
  };*/

};


#endif

//****** end abstract_mat_c.hpp  ***************/
