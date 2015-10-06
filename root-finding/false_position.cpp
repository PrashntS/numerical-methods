#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;

// Program Interfaces declaration. The rest of the script implemets this.
float f(float x);
float false_position_method(float (*f)(float), float, float, float);
void usage();

int main(int argc, char *argv[]) {
    if(argc != 4) usage();

    float val1, val2, precision, root;
    val1 = atof(argv[1]);
    val2 = atof(argv[2]);
    precision = atof(argv[3]);

    root = false_position_method(f, precision, val1, val2);
    return 0;
}

/**
 * Implementation of False Positin Method.
 * Arguments:
 *     f: The Function
 *     dx_accuracy: Desired Accuracy
 *     x1, x2: Initial Range
 * Returns:
 *     float (root)
 */
float false_position_method (float (*f)(float), float dx_accuracy, float x1, float x2) {
    if((f(x1)*f(x2) > 0)) {
        cerr << "bad range" << endl;
        exit(1);
    }

    float Eo = x2 - x1;

    unsigned int NumOfIter = 100; // Maximum 100 iterations.
    float dx, point1, x_new, f_rt;

    dx = x2 - x1;

    for (unsigned int iter = 1; iter <= NumOfIter; iter++) {
        x_new = x2 - f(x2) * ((x2 - x1) / (f(x2) - f(x1)));
        f_rt = f(x_new);   // deterine f(midpoint)

        cout << left << setw(20) << iter << left << setw(20) << x_new << endl;

        //root reaches the precision indicated, or root has been found

        if((dx < dx_accuracy) || (f_rt == 0.0)) {
            return x_new;
        }

        point1 = f(x1);

        // Left or Right Interval
        // Have to do this only once.
        // So, we do it for first iteration only.
        if (iter == 1);
            if((f_rt * point1) < 0)
                x2 = x_new;
            else
                x1 = x_new;

    }

    return x_new;
}


float f(float x) {
    return (x*x*x) - (2*x) - 2;
}

void usage() {
    cerr << "Usage: <x_a> <x_b> <tolerance>" << endl;
    exit(1);
}
