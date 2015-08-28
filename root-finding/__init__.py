#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from bisection_false_position import *

def f(x):
    return (x**3 + (2 * x**2) - (3 * x) - 1)

if __name__ == "__main__":
    ob = false_position(f, 9)
    ob.routine([1, ])
