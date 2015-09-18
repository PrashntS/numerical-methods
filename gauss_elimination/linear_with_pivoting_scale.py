#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Class Implementing Gauss Elimination with Pivots.
Uses in built list operators.
"""

class GaussElimination(object):
    """
    Implements the Gaussian Elimination Method to solve the System of Linear Equations.
    Uses the Scaled Pivoting technique to find the Upper Triangular matrix of the given Augmented Matrix.
    """

    def __init__(self, A, B):
        """
        Inititalised and Solves the Equation of form AX = B.
        Args:
            A (list): Matrix A
            B (list): Matrix B
        Example: eqn = GaussElimination(A, B)
                 solution = eqn.solution

        Corner Cases:
            This class DOES NOT check the conformity of the given Matrices. Expect Exceptions on every non - applicable cases.
        """

        self.A = A
        self.B = B
        self.AUG = self.augment()
        self.scale_vector = self.scale_vector_calc()
        self.solution = self.solve()

    def augment(self):
        """
        Creates the Augmented Matrix: AUG = [A | B]
        """
        self.AUG = [self.A[_] + [self.B[_]] for _ in range(len(self.A))]
        return self.AUG

    def interchange(self, _from, _to, target = "aug", push = True):
        """
        Performs the "Row Interchange" operation in a Matrix.
        Args:
            _from (int): Index of the "R1" row.
            _to (int): Index of the "R2" row.
            target (str): Defaults at `aug`. Flag to determine if the A matrix is manipulated or AUG.
            push (bool): If set True, pushes the changes to the original matrix. Defaults at True.

        Returns:
            (list): Operated Matrix
        """
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
        """
        Performs the "Product with a Constant" operation in a Matrix.
        Args:
            row_index (int): Index of the row.
            constant (float): The multiple.
            target (str): Defaults at `aug`. Flag to determine if the A matrix is manipulated or AUG.
            push (bool): If set True, pushes the changes to the original matrix. Defaults at True.

        Returns:
            (list): Operated Matrix
        """
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
        """
        Performs the "Row Addition" operation in a Matrix.
        Args:
            _from (int): Index of the "R1" row.
            _to (int): Index of the "R2" row.
            constant (float): Multiply the "R1" row with the constant. Defaults at 1.0.
            target (str): Defaults at `aug`. Flag to determine if the A matrix is manipulated or AUG.
            push (bool): If set True, pushes the changes to the original matrix. Defaults at True.

        Returns:
            (list): Operated Matrix
        """
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
        """
        Calculates the Scale Vector to be used for pivoting.

        Returns:
            (list): Scale Vector
        """
        self.scale_vector = [max([abs(__) for __ in _]) for _ in A]
        return self.scale_vector

    def upper_triangle(self):
        """
        Finds the Upper Triangular form of the Augmented Matrix.

        Returns:
            (list): Upper Triangular form of the Matrix
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
        Performs the Scaled Pivoting and Row transformation of the Matrix.
        Acts upon the common object attribute `AUG`.
        """
        for offset in range(len(self.AUG)):
            column = [(_, self.AUG[_][offset:][0]) for _ in range(offset, len(self.AUG))]
            col_lis = sorted(column, key = lambda x: abs(x[1] / self.scale_vector[offset]), reverse=True)

            self.AUG = self.AUG[0:offset] + [self.AUG[_] for (_, __) in col_lis]

    def _mul(self, a, b):
        """
        Multiplies with a threshold on zero accuracy.

        Args:
            a, b (float): Operands

        Returns:
            (float)
        """
        res = a + b

        if abs(res) < 0.0000000005:
            return 0
        else:
            return round(res, 10)

    def solve(self):
        """
        Solves the Augmented matrix for Linear Equation solution.
        """
        s = []

        #: Find the Upper Triangular Representation (Scaled Pivoting applied)
        self.upper_triangle()

        for i in range(0, len(self.AUG)):
            #: Reverse the Iteration
            j = len(self.AUG) - i - 1

            #: Find the remaining blocks
            rem = self.AUG[j][0:i]

            #: Multiply Solution and `rem` vectors
            v_mul = sum([_ * __ for (_, __) in zip(rem, s)])

            #: Find Current Element
            cur = self.AUG[j][j]

            #: Solution is This:
            sol = (self.AUG[j][len(self.AUG)] - v_mul) / cur
            
            #: Add solution
            s.append(sol)

        # Report Solution
        return s

    def mapped_sol(self):
        sol = self.solve()
        return [("x{0}".format(_), sol[_]) for _ in range(len(sol))]

if __name__ == "__main__":
    # AX = B
    A = [
        [1,  1,  1,  1, 2],
        [1,  1,  2,  3, 3],
        [-1, 0,  2,  1, 4],
        [3,  2, -1,  0, 5],
        [2,  3, -1,  0, 5]
    ]

    B = [1, 2, 1, 1, 3]

    print(GaussElimination(A, B).mapped_sol())