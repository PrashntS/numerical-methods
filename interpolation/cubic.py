#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

dat = [
    [300, 0.04],
    [400, 0.035],
    [500, 0.046],
    [600, 0.058],
    [700, 0.067],
    [800, 0.083],
    [900, 0.097],
    [1000, 0.111],
    [1100, 0.125],
]

class CubicSpline(object):

    """
    """

    def __init__(self, dat):
        self.dat = dat
        self._make_h()
        self._make_x()
        self._make_a()
        self._make_c()
        self._make_d()
        self._make_b()

    def _make_h(self):
        self.h = []
        for i in range(0, len(self.dat) - 1):
            #: x[i+1] - x[i]
            self.h.append(self.dat[i + 1][0] - self.dat[i][0])

    def _make_a(self):
        #: ai = yi
        self.a = [_[1] for _ in self.dat]

    def _make_x(self):
        #: ai = yi
        self.x = [_[0] for _ in self.dat]

    def _make_c(self):
        A, C = self._make_tridiagonal()
        self.c = list(np.linalg.solve(A, C))

    def _make_d(self):
        self.d = []

        for i in range(0, len(self.dat) - 3):
            _d = (self.c[i + 1] - self.c[i]) / 3 * self.h[i]
            self.d.append(_d)

    def _make_b(self):
        self.b = []

        for i in range(0, len(self.dat) - 3):
            _b = ((self.a[i + 1] - self.a[i]) / self.h[i]) - ((2 * self.c[i] + self.c[i + 1]) * self.h[i] / 3)
            self.b.append(_b)

    def _make_tridiagonal(self):
        d = []
        lim = len(self.dat) - 2
        A = []
        B = []
        C = []

        for i in range(lim):
            row = []
            if i == 0:
                s = self.Eq9()
                row.append(s[0])
                row.append(s[1])
                C.append(s[2])
                row += [0] * (lim - 2)
            elif i == (lim - 1):
                s = self.Eq10()
                row += [0] * (lim - 2)
                row.append(s[0])
                row.append(s[1])
                C.append(s[2])
            else:
                s = self.Astar(i)
                row += [0] * (i - 1)
                row.append(s[0])
                row.append(s[1])
                row.append(s[2])
                C.append(s[3])
                row += [0] * (lim - (i - 1) - 3)
            A.append(row)
        return A, C

    def Astar(self, i):
        a = self.h[i - 1]
        b = 2 * (self.h[i - 1] + self.h[i])
        c = self.h[i]
        d = 3 * ((self.a[i + 1] - self.a[i]) / self.h[i]) - 3 * ((self.a[i] - self.a[i - 1]) / self.h[i - 1])
        return a, b, c, d

    def Eq9(self):
        a = (3 * self.h[0]) + (2 * self.h[1]) + (self.h[0] ** 2 / self.h[1])
        b = self.h[1] - (self.h[0] ** 2 / self.h[1])
        c = (3 * (self.a[2] - self.a[1]) / self.h[1]) - (3 * (self.a[1] - self.a[0]) / self.h[0])
        return a, b, c

    def Eq10(self):
        n = len(self.dat) - 1
        a = self.h[n - 2] - (self.h[n - 1] ** 2 / self.h[n - 2])
        b = 3 * self.h[n - 1] + 2 * self.h[n - 1] + (self.h[n - 1] ** 2 / self.h[n - 2])
        c = (3 * (self.a[n] - self.a[n - 1]) / self.h[n - 1]) - (3 * (self.a[n - 1] - self.a[n - 2]) / self.h[n - 2])
        return a, b, c

    def print_eq(self):
        eq = "<{0}> + <{1}>(x - <{2}>) + <{3}>(x - <{4}>)^2 + <{5}>(x - <{6}>)^3"
        e = []
        for i in range(len(self.dat) - 3):
            e.append(eq.format(self.a[i], self.b[i], self.x[i], self.c[i], self.x[i], self.d[i], self.x[i]))
        return e

if __name__ == "__main__":
    sol = CubicSpline(dat)
    A,C = sol._make_tridiagonal()
    print(np.array(A))
    print(np.array(C))
    print(np.array(sol.print_eq()))
