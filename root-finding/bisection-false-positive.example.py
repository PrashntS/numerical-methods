#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math

def f(x):
    return (math.tan(math.pi * x) - x - 6) 

class num_method(object):
    def __init__(self, f, r = None):
        self.f = f
        self.MAX_ITER = 990
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
        elif self.root == approx_root:
            # CONVERGED!
            iter_count = self.temp_iter
            self.temp_iter = 0
            return self.root, iter_count
        else:
            self.root = approx_root
            return self.routine(new_interval)

class bisection(num_method):
    def root_approx(self, lower, upper):
        return (lower + upper) / 2

class false_positive(num_method):
    def root_approx(self, a, b):
        return (b - f(b) * ((b - a) / (f(b) - f(a))))

if __name__ == "__main__":
    obj1 = bisection(f, 9)
    obj2 = false_positive(f, 9)
    print (obj1.routine([0, 0.49]))
    print (obj2.routine([0, 0.49]))
