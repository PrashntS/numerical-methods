#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math

def f0(x):
    """Function 1"""
    return (math.tan(math.pi * x) - x - 6)

def f1(x):
    """Function 2"""
    return (x**3 + (2 * x**2) - (3 * x) - 1)

def fluid(x):
    """Function 3"""
    rho_f = 0.890
    rho_o = 0.120
    r     = 5
    return (((rho_f / 3) * x**3) - (r * rho_f * x**2) + ((4 * r**3 * rho_o) / 3))

class num_method(object):
    """
    This is a base class incorporating some common functions of all the three methods.
    Other methods inherit this class.
    """

    def __init__(self, f, r = None):
        """
        Constructor Method.
        """

        self.f = f
        self.MAX_ITER = 99
        self.root = None
        self.round = r
        self.temp_iter = 0

    def choose_interval(self, left, right):
        """
        Returns correct interval.
        """

        left_multiple  = self.f(left[0])  * self.f(left[1])
        right_multiple = self.f(right[0]) * self.f(right[1])

        if left_multiple < 0:
            return left

        elif right_multiple < 0:
            return right

        else:
            return None

    def routine(self, interval):
        """
        Defines iteration scheme.
        """

        self.temp_iter += 1
        approx_root  = self.root_approx(interval[0], interval[1])

        if self.round:
            approx_root = round(approx_root, self.round)

        left         = [interval[0], approx_root]
        right        = [approx_root, interval[1]]
        new_interval = self.choose_interval(left, right)

        if self.root == approx_root:
            iter_count = self.temp_iter
            self.temp_iter = 0
            return self.root, iter_count
        elif self.temp_iter >= self.MAX_ITER:
            iter_count = self.temp_iter
            self.temp_iter = 0
            return self.root, iter_count
        elif new_interval is None:
            iter_count = self.temp_iter
            self.temp_iter = 0
            return self.root, iter_count
        else:
            self.root = approx_root
            return self.routine(new_interval)

class bisection(num_method):
    """
    Implements Previous Class, to approximate the Root using Bisection Method.
    """

    def root_approx(self, lower, upper):
        root = (lower + upper) / 2
        print ("Method: Bisection -> Count: {0}, Root: {1}".format(self.temp_iter, root))
        return root

class false_position(num_method):
    """
    Implements False Position Method.
    """

    def __init__(self, f, r = None):
        self.pn = []
        self.er = []
        super(false_position, self).__init__(f, r)

    def root_approx(self, a, b):
        root = (b - self.f(b) * ((b - a) / (self.f(b) - self.f(a))))
        self.pn.append(root)
        print ("Method: False Posn -> Count: {0}, Root: {1} -> Error: {2}".format(self.temp_iter, root, self.error(self.temp_iter - 1)))
        return root

    def error(self, n):
        if n < 2:
            return None

        try:
            l = (self.pn[n] - self.pn[n - 1]) / (self.pn[n - 1] - self.pn[n - 2])

            return math.fabs((l / (l - 1)) * (self.pn[n] - self.pn[n - 1]))
        except ZeroDivisionError:
            return None

class newton_method(object):
    """
    Implements Newtons Method.
    """

    def __init__(self, func, func_, r = 4):
        self.f  = func
        self.f_ = func_
        self.r  = r
        self.count = 0

    def iterate(self, xo):
        self.count += 1
        val = round(self.root(xo), 10)
        print ("Method -> Newton's, Count -> {0}, Root -> {1}".format(self.count, val))
        if val == xo:
            "Converged"
            return val, self.count
        elif self.count == 200:
            return val, self.count
        else:
            return self.iterate(val)

    def root(self, xn):
        return (xn - self.f(xn) / self.f_(xn))

def func(x):
    return x * (1 - math.cos(x))

def func_(x):
    return (1 - math.cos(x)) + x * math.sin(x)

if __name__ == "__main__":
    """
    Main Function and Example Usage.
    """

    ob = bisection(func, 5)
    ob.routine([-2, 1])

    ob = newton_method(func, func_)
    ob.iterate(1)
