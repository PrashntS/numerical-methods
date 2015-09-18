#!/usr/bin/env python3
# -*- coding: utf-8 -*-

A = [
    [1,  1,  1,  1],
    [1,  1,  2,  3],
    [-1, 0,  2,  1],
    [3,  2, -1,  0]
]

B = [1, 2, 1, 1]


class Matrix(object):

    def __init__(self, A, B):
        self.A = A
        self.B = B
        self.AUG = self.augment()
        self.scale_vector = self.scale_vector_calc()

    def augment(self):
        self.AUG = [self.A[_] + [self.B[_]] for _ in range(len(self.A))]
        return self.AUG

    def interchange(self, _from, _to, target = "aug", push = True):
        mat = self.AUG if target is "aug" else self.A

        #: Swap the rows
        mat[_to], mat[_from] = mat[_from], mat[_to]

        if push:
            if target is "aug":
                self.AUG = mat
            else:
                self.A = A

        return mat

    def const_product(self, row_index, constant, target = "aug", push = True):
        mat = self.AUG if target is "aug" else self.A

        #: Swap the rows
        mat[row_index] = [_ * constant for _ in mat[row_index]]

        if push:
            if target is "aug":
                self.AUG = mat
            else:
                self.A = A

        return mat

    def add_rows(self, _from, _to, constant = 1, target = "aug", push = True):
        #: mat[_to] = mat[_to] + constant * mat[_from]

        mat = self.AUG if target is "aug" else self.A

        #: Swap the rows
        mat[_to] = [self._mul(mat[_to][_], round(constant, 3) * mat[_from][_]) for _ in range(len(mat[_from]))]

        if push:
            if target is "aug":
                self.AUG = mat
            else:
                self.A = A

        return mat

    def scale_vector_calc(self):
        self.scale_vector = [max([abs(__) for __ in _]) for _ in A]
        return self.scale_vector

    def upper_triangle(self):
        """
        Without Scaled Pivot - Partial Pivot
        """
        for offset in range(len(self.AUG)):
            """
            We have the Row
            """
            row = self.AUG[offset][offset:]
            for i in range(1, len(self.AUG) - offset):
                const = -1 * self.AUG[offset + i][offset:][0] / row[0]
                self.add_rows(offset, offset + i, const)

            self.scaled_pivot()

    def scaled_pivot(self):
        """
        """
        for offset in range(len(self.AUG)):
            column = [(_, self.AUG[_][offset:][0]) for _ in range(offset, len(self.AUG))]
            col_lis = sorted(column, key = lambda x: abs(x[1] / self.scale_vector[offset]), reverse=True)

            self.AUG = self.AUG[0:offset] + [self.AUG[_] for (_, __) in col_lis]

    def _mul(self, a, b):
        res = a + b

        if abs(res) < 0.0000000005:
            return 0
        else:
            return round(res, 10)

a = Matrix([
    [0.7, 1725],
    [0.4352, -5.433]
    ], [1739, 3.271])
a.upper_triangle()
print(a.AUG)
print(a.AUG[1][2] / a.AUG[1][1])