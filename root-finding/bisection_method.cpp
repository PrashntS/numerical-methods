#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;

// Program Interfaces declaration. The rest of the script implemets this.
float f(float x);
float bisection_method(float (*f)(float), float, float, float);
void usage();

int main(int argc, char *argv[]) {
    if(argc != 4) usage();

    float val1, val2, precision, root;
    val1 = atof(argv[1]);
    val2 = atof(argv[2]);
    precision = atof(argv[3]);

    root = bisection_method(f, precision, val1, val2);
    return 0;
}

/**
 * Implementation of Bisection Method.
 * Arguments:
 *     f: The Function
 *     dx_accuracy: Desired Accuracy
 *     x1, x2: Initial Range
 * Returns:
 *     float (root)
 */
float bisection_method (float (*f)(float), float dx_accuracy, float x1, float x2) {
    if((f(x1)*f(x2) > 0)) {
        cerr << "bad range" << endl;
        exit(1);
    }

    float Eo = x2 - x1;

    unsigned int NumOfIter = 100; // Maximum 100 iterations.
    float dx, point1, x_mid, f_mid;

    dx = x2 - x1;

    for (unsigned int iter = 1; iter <= NumOfIter; iter++) {
        x_mid = x1 + (dx*=0.5); // determine the midpoint
        f_mid = f(x_mid);   // deterine f(midpoint)

        cout << left << setw(20) << iter << left << setw(20) << x_mid << endl;

        //root reaches the precision indicated, or root has been found

        if((dx < dx_accuracy) || (f_mid == 0.0)) {
            return x_mid;
        }

        point1 = f(x1);

        // Left or Right Interval
        if((f_mid * point1) < 0)
            x2 = x_mid;
        else
            x1 = x_mid;
    }

    return x_mid;
}


float f(float x) {
    return (x*x*x) - (2*x) - 2;
}

void usage() {
    cerr << "Usage: <x_a> <x_b> <tolerance>" << endl;
    exit(1);
}
