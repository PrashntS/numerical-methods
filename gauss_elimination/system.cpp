#include <iostream>
#include <armadillo>
#include <stdlib.h>

#define CLOSE_TO_ZERO 0.00001

using namespace std;
using namespace arma;

/**
 * Gauss Elimination Without Pivoting
 * @param  aug_mat Augmented Matrix
 * @return         Solution Vector
 */
vec gauss_elimination(mat aug_mat) {
    //: Check for Singularity
    mat coeff_matrix = aug_mat.submat(0, 0, aug_mat.n_rows - 1, aug_mat.n_cols - 2);
    double determinant_coeff_mat = det(coeff_matrix);

    if (abs(determinant_coeff_mat) < CLOSE_TO_ZERO) {
        cout << "[ERR] Singular Coeffecient Matrix.";
        exit (EXIT_FAILURE);
    }

    //: Row Eliminations
    for (int i = 0; i < aug_mat.n_rows - 1; ++i) {
        for (int j = i+1; j < aug_mat.n_rows; ++j) {
            double ratio = aug_mat(j, i) / aug_mat(i, i);
            aug_mat.row(j) = aug_mat.row(j) - ratio * aug_mat.row(i);
        }
    }

    for (int i = aug_mat.n_rows - 1; i > 0; --i) {
        for (int j = i-1; j >= 0; --j) {
            double ratio = aug_mat(j, i) / aug_mat(i, i);
            aug_mat.row(j) = aug_mat.row(j) - ratio * aug_mat.row(i);
        }
    }

    // Calculating the X. Divide the last column (B) of the Augmented Matrix by the
    // Diagonal of the Augmenterd matrix.
    vec x = aug_mat.col(aug_mat.n_cols - 1) / aug_mat.diag();

    return x;
}

/**
 * Gauss Elimination With Pivoting
 * @param  aug_mat Augmented Matrix
 * @return         Solution Vector
 */
vec gauss_elimination_pivoted(mat aug_mat) {
    //: Check for Singularity
    mat coeff_matrix = aug_mat.submat(0, 0, aug_mat.n_rows - 1, aug_mat.n_cols - 2);
    double determinant_coeff_mat = det(coeff_matrix);

    if (abs(determinant_coeff_mat) < CLOSE_TO_ZERO) {
        cout << "[ERR] Singular Coeffecient Matrix.";
        exit (EXIT_FAILURE);
    }

    //: Pivoting
    for (int i = 0; i < aug_mat.n_cols-1; ++i) {
        uword r;
        aug_mat.col(i).max(r);
        aug_mat.swap_rows(r, i);
    }

    //: Row Eliminations
    for (int i = 0; i < aug_mat.n_rows - 1; ++i) {
        for (int j = i+1; j < aug_mat.n_rows; ++j) {
            double ratio = aug_mat(j, i) / aug_mat(i, i);
            aug_mat.row(j) = aug_mat.row(j) - ratio * aug_mat.row(i);
        }
    }

    for (int i = aug_mat.n_rows - 1; i > 0; --i) {
        for (int j = i-1; j >= 0; --j) {
            double ratio = aug_mat(j, i) / aug_mat(i, i);
            aug_mat.row(j) = aug_mat.row(j) - ratio * aug_mat.row(i);
        }
    }

    // Calculating the X. Divide the last column (B) of the Augmented Matrix by the
    // Diagonal of the Augmenterd matrix.
    vec x = aug_mat.col(aug_mat.n_cols - 1) / aug_mat.diag();

    return x;
}

/**
 * L U Decomposition
 * @param  aug_mat Augmented Matrix
 * @return         Solution Vector
 */
vec l_u_decomposition(mat aug_mat) {
    mat coeff_matrix = aug_mat.submat(0, 0, aug_mat.n_rows - 1, aug_mat.n_cols - 2);

    //: Declare l and u matrices.
    mat l(size(coeff_matrix)), u(size(coeff_matrix));

    //: Declare column vectors.
    colvec b, x, y;

    //: Decompose B from Augmented Matrix.
    b = aug_mat.col(aug_mat.n_cols - 1);

    //: Check for Singularity
    double determinant_coeff_mat = det(coeff_matrix);

    if (abs(determinant_coeff_mat) < CLOSE_TO_ZERO) {
        cout << "[ERR] Singular Coeffecient Matrix.";
        exit (EXIT_FAILURE);
    }

    //: Identity matrix
    l.eye();
    //: Prepare U
    u = coeff_matrix;

    //: Calculate L and U
    for (int i = 0; i < coeff_matrix.n_rows-1; ++i) {
        for (int j = i+1; j < coeff_matrix.n_rows; ++j) {
            double ratio = u(j, i) / u(i, i);

            u.row(j) = u.row(j) - ratio*u.row(i);
            l(j, i) = ratio;
        }
    }

    //: b is B vector for y
    y = (l.i()) * b;

    //: y is B vector for x
    x = (u.i()) * y;

    return x;
}

/**
 * Gauss Jacobi Method
 * NOTE: This method does NOT guarantee a solution.
 * @param  aug_mat Augmented Matrix
 * @return         Solution Vector
 */
vec gauss_jacobi(mat aug_mat) {
    mat coeff_matrix = aug_mat.submat(0, 0, aug_mat.n_rows - 1, aug_mat.n_cols - 2);

    //: Declare column vectors.
    colvec x, xi, b;

    x.resize(coeff_matrix.n_cols);
    x.zeros();
    xi.resize(coeff_matrix.n_cols);
    xi.zeros();

    //: Decompose B from Augmented Matrix.
    b = aug_mat.col(aug_mat.n_cols - 1);

    //: Check for Singularity
    double determinant_coeff_mat = det(coeff_matrix);

    if (abs(determinant_coeff_mat) < CLOSE_TO_ZERO) {
        cout << "[ERR] Singular Coeffecient Matrix.";
        exit (EXIT_FAILURE);
    }

    int count=0;
    while(++count < 50) {
        for (int i = 0 ; i < coeff_matrix.n_cols ; ++i) {
            double sigma = 0;
            for (int j = 0; j < coeff_matrix.n_cols ; ++j) {
                if (i != j) {
                    sigma = sigma + coeff_matrix(i, j) * x(j);
                }
            }
            xi(i) = (b(i) - sigma) / coeff_matrix(i, i);
        }
        //: Update X
        x = xi;
    }

    return x;
}

/**
 * Gauss Seidel Method
 * NOTE: This method does NOT guarantee a solution.
 * @param  aug_mat Augmented Matrix
 * @return         Solution Vector
 */
vec gauss_seidel(mat aug_mat) {
    mat coeff_matrix = aug_mat.submat(0, 0, aug_mat.n_rows - 1, aug_mat.n_cols - 2);

    //: Declare column vectors.
    colvec x, xi, b;

    x.resize(coeff_matrix.n_cols);
    x.zeros();
    xi.resize(coeff_matrix.n_cols);
    xi.zeros();

    //: Decompose B from Augmented Matrix.
    b = aug_mat.col(aug_mat.n_cols - 1);

    //: Check for Singularity
    double determinant_coeff_mat = det(coeff_matrix);

    if (abs(determinant_coeff_mat) < CLOSE_TO_ZERO) {
        cout << "[ERR] Singular Coeffecient Matrix.";
        exit (EXIT_FAILURE);
    }

    int count=0;
    while(++count < 50) {
        for (int i = 0 ; i < coeff_matrix.n_cols ; ++i) {
            double sigma = 0;
            for (int j = 0; j < coeff_matrix.n_cols ; ++j) {
                if (i != j) {
                    sigma = sigma + coeff_matrix(i, j) * xi(j);
                }
            }
            xi(i) = (b(i) - sigma) / coeff_matrix(i, i);
        }
        //: Update X
        x = xi;
    }

    return x;
}

/**
 * Main Function
 * @param  argc Commandline Argument Counts
 * @param  argv Commandline Argument Vectors
 * @return      Exit Status
 */
int main(int argc, char** argv) {

    //: The Augmented Matrix needs to be put here.
    mat aug_mat = {
        {2, 1, 11},
        {5, 7, 13}
    };

    gauss_elimination(aug_mat).print("Gauss Elimination Without Pivoting:");
    gauss_elimination_pivoted(aug_mat).print("Pivoted Gaussian Elimination");
    l_u_decomposition(aug_mat).print("L U Decomposition");
    gauss_jacobi(aug_mat).print("Gauss Jacobi");
    gauss_seidel(aug_mat).print("Gauss Seidel");

    return 0;
}
