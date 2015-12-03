#include<iostream>
using namespace std;

float multiply_combinations(float *arr, int n, int r, float *buffer, int i, int j) {
    if (j == r) {
        float mul = 1;
        for(int p = 0; p < r; p++)
            mul *= buffer[p];
        return mul;
    }

    for (int q = i; q <= n - r + j; q++) {
        buffer[j] = arr[q];
        return multiply_combinations(arr, n, r, buffer, q + 1, j + 1);
    }
}

float sum_combinations(float *arr, int n, int r) {
    float buffer[r], sum = 0;
    for (int p = 0; p <= n - r; p++) {
        buffer[0] = arr[p];
        sum += multiply_combinations(arr, n, r, buffer, p + 1, 1);
    }
    return sum;
}

int main() {
    float arr[] = {1,2,3, 4,5};
    std::cout << sum_combinations(arr, 5, 1) << '\n';
    std::cout << sum_combinations(arr, 5, 3) << '\n';
    std::cout << sum_combinations(arr, 5, 5) << '\n';
    return 0;
}
