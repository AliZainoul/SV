#include <cmath>
#include <iostream>
#include <utility>
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
  for (int i = 1; i <= n+2; i++){
    double vi = abs(mat[k][i]);
    if (maximm < vi) maximm = vi;
  }
  return maximm;
}


// Saint Venant Main Function
bool SV(int n,int itermax)
{
  // MesH Generating [-10,10]
  Mesh m(-10,10,n);

  // Initialization of q & H;
  FullMtx H(itermax, n+2, 0.0); // Initialization of H
  for (int k=0;k<=n+1;k++) H[0][k]=f(m.getcenter()[k]);
  FullMtx q(itermax, n+2, 0.0); // Initialization of q
  FullMtx M(itermax, n+2, 0.0);
  for (int k=0;k<=n+1;k++)  M[0][k] = fii(H[0][k], q[0][k]);
  FullMtx FH(itermax, n+2, 0.0);
  FullMtx Fq(itermax, n+2, 0.0);
  Vector dt(itermax);
  Vector lambda(itermax);

  lambda[0] = maxi(M,0,M.getncols());
  dt[0] = m.geth() / (3*lambda[0]);
  for (int i  = 1 ; i  <= n+1 ; i++)
  {
    Fq[0][i] =  (1/2)*(((q[0][i])*(q[0][i]) / H[0][i])
    + (gg/2)*((H[0][i])*(H[0][i]) + (H[0][i+1])*(H[0][i+1]))
    + (1/2)*((q[0][i+1])*(q[0][i+1]) / H[0][i+1]))
    - (lambda[0]/2)*(q[0][i+1] - q[0][i]);

    q[1][i] = q[0][i] - (dt[0] / m.geth())*(Fq[0][i]-Fq[0][i-1]);
    FH[0][i] = (1/2)*(q[0][i] + q[1][i]) - (dt[0]/2)*(H[0][i+1] - H[0][i]);
    H[1][i] = H[0][i] - (dt[0] / m.geth())*(FH[0][i]-FH[0][i-1]);
  }
  double time=0.0;
  double Tmax = 200;
  for (int iter = 1; iter < itermax; iter ++)
  {
    time += dt[iter];
    if (time>=Tmax) error("Too much time.");
    else
    {
      lambda[iter] = maxi(M,iter,M.getncols());
      dt[iter] = m.geth() / (3*lambda[iter]);
      for (int i = 1; i <= n+1; i++)
      {
        Fq[iter][i] =   (1/2)*(((q[iter][i])*(q[iter][i]) / H[iter][i])
        + (gg/2)*(H[iter][i])*(H[iter][i])
        + (1/2)*((q[iter][i+1])*(q[iter][i+1]) / H[iter][i+1])
        + (gg/2)*(H[iter][i+1])*(H[iter][i+1]))
        - (lambda[iter]/2)*(q[iter][i+1] - q[iter][i]);
        q[iter+1][i] = q[iter][i] - (dt[iter] / m.geth())*(Fq[iter][i] - Fq[iter][i-1]);
        FH[iter][i] = (1/2)*(q[iter][i] + q[iter+1][i]) - (lambda[iter]/2)*(H[iter][i+1] - H[iter][i]);
        H[iter+1][i] = H[iter][i] - (dt[iter] / m.geth())*(FH[iter][i] - FH[iter][i-1]);
        M[iter][i] = fii(q[iter][i],H[iter][i]);
      }
    }
  }
  cout << "Matrix H" << H << endl;
  cout << "Matrix q" << q << endl;
  return true;
}
