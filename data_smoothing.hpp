#pragma once

#include <cmath>
#include <vector>
#include <assert.h>

inline double gauss(double sigma, double x) {
    double expVal = -1 * (pow(x, 2) / pow(2 * sigma, 2));
    double divider = sqrt(2 * M_PI * pow(sigma, 2));
    return (1 / divider) * exp(expVal);
}

inline std::vector<double> gaussKernel(int samples, double sigma) {
    std::vector<double> v;
    
    bool doubleCenter = false;
    if (samples % 2 == 0) {
        doubleCenter = true;
        samples--;
    }
    int steps = (samples - 1) / 2;
    double stepSize = (3 * sigma) / steps;
    
    for (int i = steps; i >= 1; i--) {
        v.push_back(gauss(sigma, i * stepSize * -1));
    }
    
    v.push_back(gauss(sigma, 0));
    if (doubleCenter) {
        v.push_back(gauss(sigma, 0));
    }
    
    for (int i = 1; i <= steps; i++) {
        v.push_back(gauss(sigma, i * stepSize));
    }

    assert(v.size() == samples);
    
    return v;
}

inline std::vector<double> gaussSmoothen(std::vector<double> values, double sigma, int samples) {
    std::vector<double> out;
    auto kernel = gaussKernel(samples, sigma);
    int sampleSide = samples / 2;
    unsigned long ubound = values.size();
    for (unsigned long i = 0; i < ubound; i++) {
        double sample = 0;
        int sampleCtr = 0;
        for (long j = i - sampleSide; j <= i + sampleSide; j++) {
            if (j > 0 && j < ubound) {
                int sampleWeightIndex = int(sampleSide + (j - i));
                sample += kernel[sampleWeightIndex] * values[j];
                sampleCtr++;
            }
        }
        double smoothed = sample / (double)sampleCtr;
        out.push_back(smoothed);
    }
    return out;
}
