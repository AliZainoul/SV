#include <cmath>
#include <iostream>
//#include <utility>
using namespace std;
#include "error.hpp"
#include "Vector.hpp"
#include "full_mat_c.hpp"
#include "mesh.hpp"
#include "SV.hpp"
typedef double (*pfn) (double);

int gg = 10; // Gravity Constant
double f(double x) { return (x>0)?1:2;};
double fii(double x, double y) { return abs(x/y) + sqrt(gg*y); };
double maxi(const FullMtx& mat, int k, int n)
{
  double maximm=abs(mat[k][0]);
  for (int i = 1; i <= n-1; i++){
    double vi = abs(mat[k][i]);
    if (maximm <= vi) maximm = vi;
  }
  return maximm;
}

// Saint Venant Main Function
void SV(int n,int itermax)
{
  // MesH Generating [-10,10]
  Mesh m(-10,10,n);

  // Initialization of q & H;
  FullMtx H(itermax, n+2, 0.0); // Initialization of H
  for (int k=0;k<=n+1;k++) H[0][k]=f(m.getcenter()[k]);

  FullMtx q(itermax, n+2, 0.0); // Initialization of q
  FullMtx FH(itermax, n, 0.0);
  FullMtx Fq(itermax, n, 0.0);
  Vector dt(itermax);
  Vector lambda(itermax);
  H[0][1] = H[0][0];

  FullMtx M(itermax, n+2, 0.0);
  for (int k=0;k<=n+1;k++)  M[0][k] = fii(q[0][k], H[0][k]);
  lambda[0] = maxi(M,0,M.getncols());
  dt[0] = m.geth() / (3*lambda[0]);

  for (int i  = 0 ; i  < n ; i++)
  {
    Fq[0][i] =  (gg/4)*(H[0][i]*H[0][i] + H[0][i+1]*H[0][i+1]);
  }
  q[1][0] = q[0][0] - (dt[0]/ m.geth()) * (Fq[0][0]);
  cout << "*****" << endl;
  cout << q[1][0] ;
  cout << "*****" << endl;

  for (int i  = 1 ; i  <= n ; i++)
  {
    q[1][i] = q[0][i] - (dt[0] / m.geth())*(Fq[0][i]-Fq[0][i-1]);
  }

  for (int i  = 0 ; i  < n ; i++)
  {
    FH[0][i] = (1/2)*(q[0][i] + q[1][i]) - (dt[0]/2)*(H[0][i+1] - H[0][i]);
  }

  H[1][0] = H[0][0] - (dt[0] / m.geth()) * (FH[0][0]);
  for (int i  = 1 ; i  <= n ; i++)
  {
    H[1][i] = H[0][i] - (dt[0] / m.geth())*(FH[0][i]-FH[0][i-1]);
  }

  H[1][n+1] = H[1][n];
  q[1][n+1] = q[1][n];
  for (int k=0;k<=n+1;k++)  M[1][k] = fii(q[1][k], H[1][k]);

  double time=0.0;
  double Tmax = 20;

  for (int iter = 2; iter < itermax-2; iter ++)
  {
    time += dt[iter];
    if (time>=Tmax) error("Too much time.");
    else
    {
      lambda[iter-1] = maxi(M,iter-1,M.getncols());
      dt[iter-1] = m.geth() / (3*lambda[iter-1]);
      for (int i  = 0 ; i  < n ; i++)
      {
        Fq[iter-1][i] =  (1/2)*(((q[iter-1][i])*(q[iter-1][i]) / H[iter-1][i])
        + (gg/2)*((H[iter-1][i])*(H[iter-1][i]) + (H[iter-1][i+1])*(H[iter-1][i+1]))
        + (1/2)*((q[iter-1][i+1])*(q[iter-1][i+1]) / H[iter-1][i+1]))
        - (lambda[iter-1]/2)*(q[iter-1][i+1] - q[iter-1][i]);
      }
      for (int i  = 0 ; i  <= n ; i++)
      {
        q[iter][i] = q[iter-1][i] - (dt[iter-1] / m.geth())*(Fq[iter-1][i]-Fq[iter-1][i-1]);
      }

      for (int i  = 0 ; i  < n ; i++)
      {
        FH[iter-1][i] = (1/2)*(q[iter-1][i] + q[1][i]) - (dt[iter-1]/2)*(H[iter-1][i+1] - H[iter-1][i]);
      }
      for (int i  = 0 ; i  <= n ; i++)
      {
        H[iter][i] = H[iter-1][i] - (dt[iter-1] / m.geth())*(FH[iter-1][i]-FH[iter-1][i-1]);
        M[iter][i] = fii(q[iter][i], H[iter][i]);
      }
    }
    H[iter][1] = H[iter][0];
    q[iter][1] = q[iter][0];
    H[iter][n+1] = H[iter][n];
    q[iter][n+1] = q[iter][n];
  }

  cout << "Fq"  << Fq;
  cout << "q"   << q;
  cout << "FH" << FH;
  cout << "H"  << H;
  cout << "M"  << M;
}
