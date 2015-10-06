#include <iostream>
#include <iomanip>
#include <math.h>

#define MAXITER 50

using namespace std;

/**
 * PLEASE NOTE:
 * This script CANNOT calculate the derivative of the Function, hence, the
 * function and Derivatve MUST be hard coded by the user.
 *
 * This program accepts Command Line Arguments.
 */

// Program Interfaces declaration. The rest of the script implemets these.
void usage();
float f(float);
float f_p(float);
float newton(float, float, float, float);

int main(int argc, char * argv[]) {
    if(argc != 5) usage();

    float val1, val2, precision, root, est;
    val1 = atof(argv[1]);
    val2 = atof(argv[2]);
    est = atof(argv[3]);
    precision = atof(argv[4]);

    root = newton(precision, val1, val2, est); // approximation for root using netwon-raphson method

    return 0;
}

/**
 * Implementation of Newton's Method.
 * Arguments:
 *     dx_accuracy: Desired Accuracy
 *     x1, x2: Initial Range
 *     est: Estimated Root
 * Returns:
 *     float (root)
 */
float newton(float dx_accuracy, float x1, float x2, float est) {
    float fx, df, point_x, dx;
    point_x = est;
    for(unsigned int i = 1; i <= MAXITER; i++) {
        fx = f(point_x);
        df = f_p(point_x);

        dx = fx/df;
        point_x -= dx;

        cout << setw(20) << left << i << left << setw(20) << point_x << endl;

        // Accuracy Check
        if(fabs(dx) < dx_accuracy) {
            return point_x;
        }
        if((point_x < x1) || (point_x > x2)) {
            cerr << "Out of Bounds. ERROR." << endl;
            exit(1);
        }

    }

    return point_x;
}

float f(float x) {
    return (x*x*x) - 2*x - 2;
}

float f_p(float x) {
    return (3*x*x) -2;
}

void usage() {
    cerr << "Usage: <x_a> <x_b> <estimate of root between x_a and x_b> <tolerance>" << endl;
    exit(1);
}
