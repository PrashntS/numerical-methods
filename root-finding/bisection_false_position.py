#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math

def f0(x):
    return (math.tan(math.pi * x) - x - 6) 

def f1(x):
    return (x**3 + (2 * x**2) - (3 * x) - 1)

class num_method(object):
    def __init__(self, f, r = None):
        self.f = f
        self.MAX_ITER = 99
        self.root = None
        self.round = r
        self.temp_iter = 0

    def choose_interval(self, left, right):
        left_multiple  = self.f(left[0])  * self.f(left[1])
        right_multiple = self.f(right[0]) * self.f(right[1])
        if left_multiple < 0:
            return left

        elif right_multiple < 0:
            return right

        else:
            raise Exception

    def routine(self, interval):
        self.temp_iter   += 1
        approx_root  = self.root_approx(interval[0], interval[1])

        if self.round:
            approx_root = round(approx_root, self.round)

        left         = [interval[0], approx_root]
        right        = [approx_root, interval[1]]
        new_interval = self.choose_interval(left, right)

        if self.temp_iter >= self.MAX_ITER:
            iter_count = self.temp_iter
            self.temp_iter = 0
            return self.root, iter_count

        else:
            self.root = approx_root
            return self.routine(new_interval)

class bisection(num_method):
    def root_approx(self, lower, upper):
        return (lower + upper) / 2

class false_position(num_method):
    def __init__(self, f, r = None):
        self.pn = []
        self.er = []
        self.count = 0
        super(false_position, self).__init__(f, r)

    def root_approx(self, a, b):
        root = (b - self.f(b) * ((b - a) / (self.f(b) - self.f(a))))
        self.count += 1
        self.pn.append(root)
        print ("Count: {0}, Root: {1} -> Error: {2}".format(self.count, root, self.error(self.count - 1)))
        return root

    def error(self, n):
        if n < 2:
            return None

        try:
            l = (self.pn[n] - self.pn[n - 1]) / (self.pn[n - 1] - self.pn[n - 2])

            return math.fabs((l / (l - 1)) * (self.pn[n] - self.pn[n - 1]))
        except ZeroDivisionError:
            return None

if __name__ == "__main__":
    obj2 = false_position(f1, 20)
    print (obj2.routine([1, 2]))
