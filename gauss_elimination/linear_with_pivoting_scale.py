#!/usr/bin/env python3
# -*- coding: utf-8 -*-

A = [
    [5,  7, 14, -8],
    [2, -2, -1,  2],
    [3,  1,  4, -1],
    [1,  3,  2,  4]
]

B = [20, 1, 7, -4]


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
        mat[_to] = [mat[_to][_] + constant * mat[_from][_] for _ in range(len(mat[_from]))]

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
            print(self.AUG)
            #print(row[0])
            for i in range(1, len(self.AUG) - offset):
                const = -1 * self.AUG[offset + i][offset:][0] / row[0]
                self.add_rows(offset, offset + i, const)

            # const = -1 * self.AUG[offset + 2][offset:][0] / row[0]
            # self.add_rows(offset, offset + 2, const)
            # const = -1 * self.AUG[offset + 3][offset:][0] / row[0]
            # self.add_rows(offset, offset + 3, const)

    def zero_transpose(self):
        """
        """
        pass


a = Matrix(A, B)

a.upper_triangle()
print(a.AUG)