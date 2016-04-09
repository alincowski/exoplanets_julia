/* Newton-Raphshon (AKA Newton's Method) Rootfinding
   by Andrew Lincowski
   for ASTR598 Exoplanets
   Comparing to ExoJulia package

   Time-stamp: <newt_meth_roots.cc on Friday, 8 April, 2016 at 13:18:13 PDT (linc)>

   compile with:
      g++ ps2.2.cc -o2 -o ps2.2

   run as:
      ./ps2.2

   sample output:

[lincowski@nimoy ps2]$ g++ ps2.2.cc -o2 -o ps2.2
[lincowski@nimoy ps2]$ ps2.2
For f(x) = x^3 - 25x^2 + 165x - 275, the roots are: 
x = 2.55406
x = 6.94707
x = 15.4989



*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>

using namespace std;

double func(double *ecc, double *E, double *M) {

  return *E-*ecc*sin(*E) - *M;

}

/*
double deriv(double ecc, double E, double M) {

    double h = 1e-12;
    return (func(ecc,E+h,M) - func(ecc,E,M) ) / h;

}
*/

double deriv(double *ecc, double *E, double *M) {

    return 1-*ecc*cos(*E);

}


double newt(double *ecc, double *E, double *M) {


  double h = -func(ecc,E,M)/deriv(ecc,E,M);
  double e = fabs(*E+h);
  *E = *E + h;

  if(fabs(h) > 1e-12) return newt(ecc,E,M);
  else return *E;
  
}

int sign(double x){

  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;

}

/*
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
*/

main() {

  int N = 100;

  double ecc[N];
  double M[N];
  double E[N][N];

  double ecc_max = 0.999;
  double M_max = 2*M_PI;
  double decc = ecc_max/N;
  double dM = M_max/N;

  for(int i=0;i<N;i++) {
    for(int j=0;j<N;j++) {
      M[i] = i*dM/M_max;
      ecc[j] = j*decc/ecc_max;
      E[i][j] = M[i] + 0.85*ecc[j]*sign(sin(M[i]));
      E[i][j] = newt(&ecc[j],&E[i][j],&M[i]);

    }
  }


}
