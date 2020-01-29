#ifndef FULL_MAT_C_H
#define FULL_MAT_C_H
#include "Vector.hpp"
#include "abstract_mat_c.hpp"
using namespace std;

class FullMtx: public AbsMtx
{
private:
  //int ncols; // Number of columns of a Matrix
	//int nrows; // Number of rows of a Matrix
  //double** mx; // Entries of the Matrix

public:
  FullMtx(int n, int m, double**);                          // Constructor 1
  FullMtx(int n, int m, double t = 0);                      // Constructor 2
  FullMtx(initializer_list<initializer_list<double>> lst);  // Constructor 3
  FullMtx(const FullMtx&);                                  // Copy Constructor
  ~FullMtx()                                                // Destructor
  {for (int i = 0; i < nrows; i++) delete[]  mx[i]; delete[] mx;}

  friend FullMtx operator * (const FullMtx&, const FullMtx&); // Matrix Multiplication
  FullMtx& operator=(const FullMtx&); // Overload of the Operator '='
  Vector operator*(const Vector&) const; // Matrix-Vector multiply
  double* operator[](int i) const { return mx[i]; } // Method that Returns the i-th row
  friend ostream& operator<<(ostream&, const FullMtx&);     // Overload of the Operator '<<'
  double getnrows() const { return nrows; };
  double getncols() const { return ncols; };


};


#endif
//****** End full_mat_c.hpp  ******//
