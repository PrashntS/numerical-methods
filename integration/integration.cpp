#include <iostream>
#include <functional>
#include <cmath>

using namespace std;

/**
 * Integrates f from a to be using trapezoidal rule.
 * @param  f Function f(x)
 * @param  a Lower Limit
 * @param  b Upper Limit
 * @return   Integral value
 */
double trapezoidal_rule(function <double (double)> f, double a, double b) {
  double h = b - a;

  return (h / 2) * (f(a) + f(b));
}

/**
 * Integrates f from a to b using simpson's rule.
 * See `trapezoidal_rule`.
 */
double simpsons_rule(function <double (double)> f, double a, double b) {
  double h = b - a,
         f1 = f(a),
         f2 = f(a + h / 2),
         f3 = f(b);

  return (h / 6) * (f1 + (4 * f2) + f3);
}

/**
 * One Point Quadrature Method in One Dimension
 * @overloaded
 */
double one_point_quadrature1D(function <double (double)> f, double a, double b) {
  // Mid Point Method, integrates 2(1) - 1 = 1 degree polynomial exactly.
  double  w = (b - a),
          x = (b + a) / 2;

  return w * f(x);
}

/**
 * One Point Quadrature Method in Two Dimensions
 * @overloaded
 */
double one_point_quadrature2D(function <double (double, double)> f,
                            double a,
                            double b,
                            double c,
                            double d) {
  double wi = (b - a),
         wj = (c - d),
         xi = (b + a) / 2,
         xj = (c + d) / 2;

  return wi * wj * f(xi, xj);
}

/**
 * Two Point Quadrature Method in One Dimension
 * @overloaded
 */
double two_point_quadrature1D(function <double (double)> f, double a, double b) {
  //: Transformation values
  double p = (b - a) / 2,
         q = (b + a) / 2,
         w = 1 / sqrt(3);
  //: Quadratures
  double x1 = q - w * p,
         x2 = q + w * p;
  //: Integral
  return p * (f(x1) + f(x2));
}

/**
 * Two Point Quadrature in Two Dimensions.
 * @overloaded Function
 */
double two_point_quadrature2D(function <double (double, double)> f,
                            double a,
                            double b,
                            double c,
                            double d) {
  //: Transformation values
  double pi = (b - a) / 2,
         qi = (b + a) / 2,
         pj = (d - c) / 2,
         qj = (d + c) / 2,
         w = 1 / sqrt(3);
  //: Quadratures
  double x1i = qi - w * pi,
         x2i = qi + w * pi,
         x1j = qj - w * pj,
         x2j = qj + w * pj;
  //: Integral
  return pi * pj * (f(x1i, x1j) + f(x1i, x2j) + f(x2i, x1j) + f(x2i, x2j));
}

/**
 * Implements the Composite Method for any methods defined here.
 * @param rule Numerical Method to be used.
 * @param f    Function f(x) to be passed on to the method.
 * @param a    Lower Limit
 * @param b    Upper Limit
 * @param subintervals  Number of Sub intervals to be made between a and b
 * @return     Integral Value
 */
double composite(function <double (function <double (double)>, double, double)> rule,
                 function <double (double)> f,
                 double a,
                 double b,
                 int subintervals) {
  double q = (b - a) / subintervals,
         integral = 0, _a = 0, _b = 0;
  for (int i = 0; i < subintervals; i += 1) {
    _a = a + i * q;
    _b = a + (i + 1) * q;
    integral += rule(f, _a, _b);
  }
  return integral;
}

/**
 * Implements Composite Method for Two Dimensional integrals.
 */
double composite(function <double (function <double (double, double)>, double, double, double, double)> rule,
                 function <double (double, double)> f,
                 double a,
                 double b,
                 double c,
                 double d,
                 int subintervals) {
  double qi = (b - a) / subintervals,
         qj = (d - c) / subintervals,
         integral = 0, _a = 0, _b = 0, _c = 0, _d = 0;
  for (int i = 0; i < subintervals; i += 1) {
    _a = a + i * qi;
    _b = a + (i + 1) * qi;
    _c = c + i * qj;
    _d = d + (i + 1) * qj;
    integral += rule(f, _a, _b, _c, _d);
  }
  return integral;
}

/**
 * One Dimesional Function
 */
double f1(double x) {
  return sin(x * x);
}

/**
 * Two Dimesional Function
 */
double f2(double x, double y) {
  return sin(x * x) + cos(y * y);
}

int main(int argc, char const *argv[]) {
  cout << "Integral of One Dimesional function f(x) = sin(x^2) from 0 to 5." << "\n\n";

  cout << "Trapezoidal Rule:      " << trapezoidal_rule(f1, 0, 5)     << "\n";
  cout << "Simpson's Rule:        " << simpsons_rule(f1, 0, 5)        << "\n";
  cout << "One Point Quadrature:  " << one_point_quadrature1D(f1, 0, 5) << "\n";
  cout << "Two Points Quadrature: " << two_point_quadrature1D(f1, 0, 5) << "\n\n";

  cout << "Composite Rule (20 subintervals)" << "\n\n";

  cout << "Trapezoidal Rule:      " << composite(trapezoidal_rule, f1, 0, 5, 20)     << "\n";
  cout << "Simpson's Rule:        " << composite(simpsons_rule, f1, 0, 5, 20)        << "\n";
  cout << "One Point Quadrature:  " << composite(one_point_quadrature1D, f1, 0, 5, 20) << "\n";
  cout << "Two Points Quadrature: " << composite(two_point_quadrature1D, f1, 0, 5, 20) << "\n\n\n";

  cout << "Integral of Two Dimesional function f(x, y) = sin(x^2) + cos(y^2) from x, y 0 to 5." << "\n\n";

  cout << "One Point Quadrature:  " << one_point_quadrature2D(f2, 0, 5, 0, 5) << "\n";
  cout << "Two Points Quadrature: " << two_point_quadrature2D(f2, 0, 5, 0, 5) << "\n\n";

  cout << "Composite Rule (20 subintervals)" << "\n\n";

  cout << "One Point Quadrature:  " << composite(one_point_quadrature2D, f2, 0, 5, 0, 5, 20) << "\n";
  cout << "Two Points Quadrature: " << composite(two_point_quadrature2D, f2, 0, 5, 0, 5, 20) << "\n";

  return 0;
}
