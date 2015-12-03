#include <iostream>

/**
 * Used Internally to build polynomial from roots.
 */
double multiply_combinations(double *arr, int n, int r, double *buffer, int i, int j) {
    if (j == r) {
        double mul = 1;
        for(int p = 0; p < r; p++)
            mul *= buffer[p];
        return mul;
    }

    for (int q = i; q <= n - r + j; q++) {
        buffer[j] = arr[q];
        return multiply_combinations(arr, n, r, buffer, q + 1, j + 1);
    }
}

/**
 * Used Internally to build polynomial from roots.
 */
double sum_combinations(double *arr, int n, int r) {
    double buffer[r], sum = 0;
    for (int p = 0; p <= n - r; p++) {
        buffer[0] = arr[p];
        sum += multiply_combinations(arr, n, r, buffer, p + 1, 1);
    }
    return sum;
}

/**
 * Lagrange Polynomial
 * @param  n Number of interpolating points.
 * @param  x X values
 * @param  y f(x) Values
 * @return   String representation of the polynomial.
 */
std::string lagrange(int n, double *x, double *y) {
  int m = n - 1; // n-1 points
  // Calculate "Lm,0, Lm,1, Lm,2 .. Lm,n"
  double L[n];
  double P[n][m];
  double Coeffs[n];

  for (int i = 0; i < n; i++) {
    double denominator = 1;
    double polynomial_roots[m];
    for (int j = 0; j < n; j++) {
      if (i == j) {
        continue;
      }
      polynomial_roots[j] = x[j];
      denominator *= x[i] - x[j];
    }
    L[i] = y[i] / denominator;
    for (int j = 1; j <= n; j++) {
      P[i][j] = L[i] * sum_combinations(polynomial_roots, m, j);
    }
  }

  for (int i = 0; i < n; ++i) {
    Coeffs[i] = 0;
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      Coeffs[i] += P[i][j];
    }
  }

  std::string polynomial = "";

  for (int i = 0; i < n; ++i) {
    polynomial += std::to_string(Coeffs[i]) + "x^" + std::to_string(i) + " + ";
  }

  return polynomial;
}

/**
 * Newton's Divided Difference Method.
 * @param  n Number of interpolating points.
 * @param  x X values
 * @param  y f(x) Values
 * @param  k Estimation point.
 * @return   Estimated value
 */
double divided_difference(int n, double *x, double *y, double k) {
  double p[100],     // To store intermideate values.
        f1 = 0,     // Temporary Variable
        f2 = 0,     // Stores approximation for f(k)
        f  = y[1];  // Keeps initial value of y[1]
  int i, j = 1;

  do {
    for (i = 1; i <= n - 1; i++) {
      p[i] = ((y[i + 1] - y[i]) / (x[i + j] - x[i]));
      y[i] = p[i];
    }

    f1 = 1;
    for(i=1; i <= j; i++) {
      f1 *= (k - x[i]);
    }

    f2 += (y[1]*f1);

    n--;
    j++;
  } while(n != 1);

  return f + f2;
}

/**
 * Cubic Spline, Not a Knot condition.
 * @param  n Number of interpolating points.
 * @param  x X values
 * @param  f f(x) Values
 * @param  t Estimation point.
 * @return   Estimated value
 */
double cubic_not_a_knot(int n, double * x, double * f, double t) {
  double * h, * dl, * dd, * du, * b, * c, * d;
  int i, found = 0;

  h = new double[n];
  dl = new double[n];
  dd = new double[n];
  du = new double[n];
  b = new double[n];
  c = new double[n];
  d = new double[n];

  for (i = 0; i < n - 1; i++) {
    h[i] = x[i + 1] - x[i];
  }

  for (i = 0; i < n - 3; i++) {
    dl[i] = du[i] = h[i + 1];
  }

  for (i = 0; i < n - 2; i++) {
    dd[i] = 2.0 * (h[i] + h[i + 1]);
    c[i] = (3.0 / h[i + 1]) * (f[i + 2] - f[i + 1]) -
      (3.0 / h[i]) * (f[i + 1] - f[i]);
  }

  dd[0] += (h[0] + h[0] * h[0] / h[1]);
  dd[n - 3] += (h[n - 2] + h[n - 2] * h[n - 2] / h[n - 3]);
  du[0] -= (h[0] * h[0] / h[1]);
  dl[n - 4] -= (h[n - 2] * h[n - 2] / h[n - 3]);

  int m = n - 2;


  for (i = 0; i < m - 1; i++) {
    du[i] /= dd[i];
    dd[i + 1] -= dl[i] * du[i];
  }

  c[0] /= dd[0];
  for (i = 1; i < m; i++)
    c[i] = (c[i] - dl[i - 1] * c[i - 1]) / dd[i];

  for (i = m - 2; i >= 0; i--)
    c[i] -= c[i + 1] * du[i];


  for (i = n - 3; i >= 0; i--)
    c[i + 1] = c[i];
  c[0] = (1.0 + h[0] / h[1]) * c[1] - h[0] / h[1] * c[2];
  c[n - 1] = (1.0 + h[n - 2] / h[n - 3]) * c[n - 2] - h[n - 2] / h[n - 3] * c[n - 3];
  for (i = 0; i < n - 1; i++) {
    d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
    b[i] = (f[i + 1] - f[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
  }


  i = 1;

  while (!found && (i < n - 1)) {
    if (t < x[i])
      found = 1;
    else
      i++;
  }
  t = f[i - 1] + (t - x[i - 1]) * (b[i - 1] + (t - x[i - 1]) * (c[i - 1] +
    (t - x[i - 1]) * d[i - 1]));
  return (t);
}

/**
 * Cubic Spline, Clamped condition.
 * @param  n Number of interpolating points.
 * @param  x X values
 * @param  f f(x) Values
 * @param  t Estimation point.
 * @param  fpa Value of derivative at a
 * @param  fpb Value of derivative at b
 * @return   Estimated value
 */
double cubic_clamped(int n, double * x, double * f, double t, double fpa, double fpb) {
  double * h, * dl, * dd, * du, * b, * c, * d;
  int i, found = 0;

  h = new double[n];
  dl = new double[n];
  dd = new double[n];
  du = new double[n];
  b = new double[n];
  c = new double[n];
  d = new double[n];

  for (i = 0; i < n - 1; i++) {
    h[i] = x[i + 1] - x[i];
    dl[i] = du[i] = h[i];
  }

  dd[0] = 2.0 * h[0];
  dd[n - 1] = 2.0 * h[n - 2];
  c[0] = (3.0 / h[0]) * (f[1] - f[0]) - 3.0 * fpa;
  c[n - 1] = 3.0 * fpb - (3.0 / h[n - 2]) * (f[n - 1] - f[n - 2]);
  for (i = 0; i < n - 2; i++) {
    dd[i + 1] = 2.0 * (h[i] + h[i + 1]);
    c[i + 1] = (3.0 / h[i + 1]) * (f[i + 2] - f[i + 1]) -
      (3.0 / h[i]) * (f[i + 1] - f[i]);
  }

  int m = n - 2;


  for (i = 0; i < m - 1; i++) {
    du[i] /= dd[i];
    dd[i + 1] -= dl[i] * du[i];
  }

  c[0] /= dd[0];
  for (i = 1; i < m; i++)
    c[i] = (c[i] - dl[i - 1] * c[i - 1]) / dd[i];

  for (i = m - 2; i >= 0; i--)
    c[i] -= c[i + 1] * du[i];


  for (i = 0; i < n - 1; i++) {
    d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
    b[i] = (f[i + 1] - f[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
  }


  i = 1;

  while (!found && (i < n - 1)) {
    if (t < x[i])
      found = 1;
    else
      i++;
  }
  t = f[i - 1] + (t - x[i - 1]) * (b[i - 1] + (t - x[i - 1]) * (c[i - 1] +
    (t - x[i - 1]) * d[i - 1]));
  return (t);
}

/**
 * Main Function
 */
int main() {
    double x[10] = {1, 2, 3, 4, 5};
    double y[10] = {1, 4, 9, 16, 25};
    // Function value at 3.4
    std::cout << cubic_not_a_knot(5, x, y, 3.4);
    std::cout << '\n';
    // Using normal conditions.
    std::cout << cubic_clamped(5, x, y, 3.4, 0, 0);
    std::cout << '\n';
    std::cout << lagrange(5, x, y);
    std::cout << '\n';
    std::cout << divided_difference(5, x, y, 3.4); // Function value at 3.4
    std::cout << '\n';
    return 0;
}
