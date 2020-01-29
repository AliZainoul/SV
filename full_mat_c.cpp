#include <cstdlib>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "full_mat_c.hpp"
#include "error.hpp"
#include "Vector.hpp"

using namespace std;


// Definitions for members of class FullMtx

FullMtx::FullMtx(int n, int m, double** dbp) {
  nrows = n;
  ncols = m;
  mx = new double* [nrows];
  for (int i =  0; i < nrows; i++) {
    mx[i] = new double [ncols];
    for (int j = 0; j < ncols; j++) mx[i][j] = dbp[i][j];
  }
}

FullMtx::FullMtx(int n, int m, double a) {
  nrows = n;
  ncols = m;
  mx = new double* [nrows];
  for (int i =  0; i< nrows; i++) {
    mx[i] = new double [ncols];
    for (int j = 0; j < ncols; j++) mx[i][j] = a;
  }
}

FullMtx::FullMtx(initializer_list<initializer_list<double>> lst){
  nrows = lst.size();
  ncols = (*lst.begin()).size();
  mx = new double* [nrows];
  int i=0;
  for (const auto& l : lst) {
    mx[i] = new double [ncols];
	  copy(l.begin(), l.end(),mx[i] );
	  i++;
  }
}

FullMtx::FullMtx(const FullMtx& mat) {
  nrows = mat.nrows;
  ncols = mat.ncols;
  mx = new double* [nrows];
  for (int i =  0; i< nrows; i++) {
    mx[i] = new double [ncols];
    for (int j = 0; j < ncols; j++) mx[i][j] = mat[i][j];
  }
}

FullMtx & FullMtx::operator=(const FullMtx & mat) {
  if (this != &mat) {
    if (nrows !=mat.nrows || ncols !=mat.ncols)
      error("Bad matrix sizes in FullMtx::operator=()");
    for (int i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++)
    mx[i][j]  = mat.mx[i][j];
  }
  return *this;
}

// Useful for Testing Small Matrices
ostream& operator<<(ostream& s, const FullMtx& mat)
{
  for (int i =0; i< mat.nrows; i++)
  {
    s << "\n" << i << "-th row:    ";
    for (int j =0; j< mat.ncols; j++)
    {
      s << mat.mx[i][j] << "  ";
      if (j%10 == 9) s << "\n";
    }
    s << "\n";
  }
  return s;
}

Vector FullMtx::operator*(const Vector& vec) const
{
  if (ncols != vec.size())
    error("matrix and vector sizes do not match in FullMtx::operator*()");
  Vector tm(nrows);
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++) tm[i] += mx[i][j]*vec[j];
  return tm;
}

FullMtx operator*(const FullMtx& Mtx1, const FullMtx& Mtx2)
{
  if (Mtx1.ncols != Mtx2.nrows) error ("BAd matrix sizes.");
  FullMtx Mtx (Mtx1.nrows,Mtx2.ncols);
  for (int i = 0; i < Mtx1.nrows; i++)
  for (int j = 0; j < Mtx2.ncols; j++)
  for (int k = 0; k < Mtx1.ncols; k++)
  Mtx[i][j] += (Mtx1[i][k]) * (Mtx2[k][j]);
  return Mtx;
}

// ****** End full_mat_c.cpp  ****** //
