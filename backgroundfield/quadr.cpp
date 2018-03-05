/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*
This file is part of Vlasiator.
*/

#include <math.h>
#include <stdlib.h>

#include "quadr.hpp"

/*
  1D,2D,3D Romberg non-singular integration a'la Numerical Recipes.
  Integration bounds must be constants.
  Header file is quadr.H.
  Test program is tstquadr.C
  Same in Mathematica is tstquadr.ma
*/

static void trapez(const T1DFunction& func, double a, double b, double& S, int& it, int n)
/* Compute the nth stage of refinement of extended trapezoidal rule.
   When called with n=1, it returns S as the crudest estimate of the integral.
   Subsequent calls with n=2,3,... (in that sequential order) will improve
   the accuracy of S by adding 2^(n-2) additional interior points.
   The argument it need not be assigned a value before, but:
   *** S and it must not be modified between sequential calls! ***/
{
   int j;
   if (n == 1) {
      S = 0.5*(b-a)*(func.call(a) + func.call(b));
      it = 1;
      //recflops(4);
   } else {
      const double delta = (b-a)/it;   // the spacing of points to be added
      double x = a + 0.5*delta;
      double sum = 0;
      for (j=0; j<it; j++) {
         sum+= func.call(x);
         x+= delta;
      }
      S = 0.5*(S + (b-a)*sum/it);      // replacement of S by its refined value
      it*= 2;
      //recflops(6+2*flops_div+it*2);
   }
}

static void polint(const double xa[], const double ya[], int n, double x, double& y, double& dy)
/* Given arrays xa and ya, each of length n, and given a value x,
   return y and error estimate dy. If P(x) is the polynomial of degree n-1 such
   that P(xa[i]) == ya[i], i=0..n-1, then the returned value y is P(x).*/
{
   int i,m,ns;
   const int nmax = 20;
   double C[nmax],D[nmax];
   ns = 0;
   double dif = fabs(x-xa[0]);
   for (i=0; i<n; i++) {
      double dift = fabs(x-xa[i]);
      if (dift < dif) {
         ns = i;
         dif = dift;
      }
      C[i] = D[i] = ya[i];
   }
   //recflops(n*2);
   y = ya[ns];
   ns--;
   for (m=0; m<n-1; m++) {
      for (i=0; i<n-m-1; i++) {
         const double h0 = xa[i] - x;
         const double hp = xa[i+m+1] - x;
         const double W = C[i+1] - D[i];
         double den = h0 - hp;
         den = W/den;
         D[i] = hp*den;
         C[i] = h0*den;
      }
      if (2*(ns+1) < n-m-1)
         dy = C[ns+1];
      else {
         dy = D[ns];
         ns--;
      }
      y+= dy;
   }
   //recflops((n-1)*((6+flops_div)*(n-m-1)+1));
}

static void ratint(const double xa[], const double ya[], int n, double x, double& y, double& dy)
/* Given arrays xa and ya, each of length n, and given a value x,
   return y and error estimate dy. The value returned is that of the diagonal rational
   function, evaluated at x, which passed through the n points (xa[i],ya[i]), i=0..n-1. */
{
  int i,ns=0;
   double w,t,hh,h,dd;
   const int nmax = 20;
   const double tiny = 1e-30;      // a small number, smaller than anything...
   double C[nmax],D[nmax];
   hh = fabs(x-xa[0]);
   for (i=0; i<n; i++) {
      h = fabs(x-xa[i]);
      if (h == 0) {
         y = ya[i];
         dy = 0;
         return;
      } else if (h < hh) {
         ns = i;
         hh = h;
      }
      C[i] = ya[i];
      D[i] = ya[i] + tiny;   // avoid 0/0 situation
   }
   //recflops(2+n*3);
   y = ya[ns--];
   int m1;
   for (m1=1; m1<n; m1++) {
      for (i=0; i<n-m1; i++) {
         w = C[i+1] - D[i];
         h = xa[i+m1] - x;
         t = (xa[i] - x)*D[i]/h;   // h will never be zero since this was tested above
         dd = t - C[i+1];
         if (dd == 0) {
            cerr << "*** Error in ratint\n";
            exit(111);
         }
         // this error occurs if the interpolating function has a pole at the requested value of x
         dd = w/dd;
         D[i] = C[i+1]*dd;
         C[i] = t*dd;
      }
      //recflops((n-m1)*(7+2*flops_div));
      y+= (dy = (2*(ns+1) < (n-m1) ? C[ns+1] : D[ns--]));
   }
}

double Romberg_simple(const T1DFunction& func, double a, double b, double absacc)
{
   int j,k1, it = 0;
   const int jmax = 8/*10*//*20*/;   // maximum number of steps
   const int k = 3/*4*//*5*/;      // (maximum) number of points used in the extrapolation
   double S[jmax+1];         // successive trapezoidal approximations
   double H[jmax+1];         // and their step sizes
   H[0] = 1;
   double result = 0, dresult = 0;
   for (j=1; j<=jmax; j++) {
      trapez(func,a,b,S[j-1],it,j);
      result = S[j-1];
      if (j >= 1/*k*/) {
         k1 = (j < k) ? j : k;
         polint(&H[j-k1],&S[j-k1],k1,0.,result,dresult);
         if (fabs(dresult) < absacc) {/*cout << "dresult=" << dresult << ", absacc=" << absacc << "\n";*/ break;}
      }
      S[j] = S[j-1];
      H[j] = 0.25*H[j-1];
   }
//   if (j > jmax) clog << "warning: romberg had " << j << " steps, result=" << result << ", dresult=" << dresult << "\n";
   return result;
}

double Romberg(const T1DFunction& func, double a, double b, double absacc)
{
   int j, k1, it = 0;
   const int jmax = 8/*20*/;   // maximum number of steps
   const int k = 5/*3*//*4*//*5*/;      // Maximum number of points used in the extrapolation
   const int min_k = 2/*k*/;         // Minimum number of points used in the interpolation
   double S[jmax+1];         // successive trapezoidal approximations
   double H[jmax+1];         // and their step sizes
   H[0] = 1;
   double result=0, result_pol=0, result_rat=0;
   double dresult, dresult_pol = 0, dresult_rat, dresult_polrat;
   for (j=1; j<=jmax; j++) {
      trapez(func,a,b,S[j-1],it,j);
      result = S[j-1];
      if (j >= min_k) {
         k1 = (j < k) ? j : k;
         // Try both rational and polynomial extrapolation.
         // The error estimate is max(dresult_pol, dresult_rat, fabs(result_pol - result_rat));
         polint(&H[j-k1],&S[j-k1],k1,0.,result_pol,dresult_pol);
         ratint(&H[j-k1],&S[j-k1],k1,0.,result_rat,dresult_rat);
         dresult_pol = fabs(dresult_pol);
         dresult_rat = fabs(dresult_rat);
         dresult = (dresult_pol > dresult_rat ? dresult_pol : dresult_rat);
         dresult_polrat = fabs(result_pol - result_rat);
         if (dresult_polrat > dresult) dresult = dresult_polrat;
         if (dresult < absacc) {
                //cout << "dresult=" << dresult << ", absacc=" << absacc << "\n";
            result = 0.5*(result_pol + result_rat);
            break;
         }
      }
      S[j] = S[j-1];
      H[j] = 0.25*H[j-1];
   }
//   if (j > jmax) cout << "warning: romberg had " << j << " steps, result=" << result << ", dresult=" << dresult << "\n";
   return result;
}



class Tinty_f2D: public T1DFunction {
private:
   const T2DFunction& f;
   double ymin2D, ymax2D, the_absacc2D;
public:
   Tinty_f2D(const T2DFunction& f1, double ymin, double ymax, double absacc)
      : f(f1), ymin2D(ymin), ymax2D(ymax), the_absacc2D(absacc) {}
   virtual double call(double x) const {return Romberg(T2D_fix1(f,x),ymin2D,ymax2D,the_absacc2D);}
};

double Romberg(const T2DFunction& func, double a, double b, double c, double d, double absacc)
{
//   cout << "2d romberg a=" << a << ", b=" << b << ", c=" << c << ", d=" << d << ", absacc=" << absacc << "\n";
   return Romberg(Tinty_f2D(func,c,d,absacc/(b-a)),a,b,absacc);
}

// 3D

class Tintxy_f3D : public T1DFunction {
private:
   const T3DFunction& f;
   double xmin3D,xmax3D, ymin3D,ymax3D, the_absacc3D;
public:
   Tintxy_f3D(const T3DFunction& f1, double xmin,double xmax, double ymin,double ymax, double absacc)
      : f(f1), xmin3D(xmin),xmax3D(xmax), ymin3D(ymin),ymax3D(ymax), the_absacc3D(absacc) {}
   virtual double call(double z) const {return Romberg(T3D_fix3(f,z),xmin3D,xmax3D,ymin3D,ymax3D,the_absacc3D);}
};

double Romberg(const T3DFunction& func, double a, double b, double c, double d, double e, double f, double absacc)
{
   return Romberg(Tintxy_f3D(func,a,b,c,d,absacc/(f-e)),e,f,absacc);
}
