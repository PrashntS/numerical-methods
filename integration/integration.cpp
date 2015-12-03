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

double f(double x) {
  return sin(x * x);
}

int main(int argc, char const *argv[]) {
  cout << composite(trapezoidal_rule, f, 0, 5, 20) << '\n';
  cout << composite(simpsons_rule, f, 0, 5, 20) << '\n';
  return 0;
}
