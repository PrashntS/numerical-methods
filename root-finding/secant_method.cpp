#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;

// Program Interfaces declaration. The rest of the script implemets this.
float f(float x);
float secant_method(float (*f)(float), float, float);
void usage();

int main(int argc, char *argv[]) {
    if(argc != 3) usage();

    float val1, val2, root;
    val1 = atof(argv[1]);
    val2 = atof(argv[2]);

    root = secant_method(f, val1, val2);
    return 0;
}

/**
 * Implementation of Secant Method
 * Arguments:
 *     f: The Function
 *     xA, xB: Initial Range
 * Returns:
 *     float (root)
 */
float secant_method (float (*f)(float), float xA, float xB) {
    float e = 1.0e-12;
    float fA, fB;
    float d;
    int i;
    int limit = 50;

    fA=(*f)(xA);
    for (i=0; i<limit; i++) {
        fB=(*f)(xB);
        d = (xB - xA) / (fB - fA) * fB;
        if (fabs(d) < e)
            break;
        xA = xB;
        fA = fB;
        xB -= d;
        cout << left << setw(20) << i << left << setw(20) << xB << endl;
    }
    if (i==limit) {
        printf("Function is not converging near (%7.4f,%7.4f).\n", xA,xB);
        return -99.0;
    }
    return xB;
}


float f(float x) {
    return (x*x*x) - (2*x) - 2;
}

void usage() {
    cerr << "Usage: <x_a> <x_b> <tolerance>" << endl;
    exit(1);
}
