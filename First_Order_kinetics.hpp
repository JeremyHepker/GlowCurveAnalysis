//
//  First_Order_kinetics.hpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#ifndef First_Order_kinetics_hpp
#define First_Order_kinetics_hpp

#include <string>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>
class First_Order_Kinetics{
private:
    std::vector<double> temp_data;
    std::vector<double> count_data;
    std::vector<double> curve;
    std::vector<std::vector<double>> glow_curves;
    const double DERIV_STEP = 1e-2;
    const int MAX_ITER = 1000;
    double k = .000086173303;
public:
    First_Order_Kinetics(std::pair<std::vector<int>,std::vector<double>>);
    double activation_energy(int TL_index,int TM_index,int TR_index);
    double frequency_factor();
    void glow_curve();
    //double gauss_newton(std::vector<double> const &temp,int max_index,double E,double Tm,double Im);
    
    double Func(const std::vector<double> input, const std::vector<double> params);
    
    std::vector<std::vector<double>> jacobian(int max_index,int TL_index,int TR_index,double E,double Tm,double Im);
    void transpose(std::vector<std::vector<double>> const &A,std::vector<std::vector<double>> &B, int n, int m);
    std::vector<std::vector<double>> multiply(std::vector<std::vector<double>> const &A,std::vector<std::vector<double>> const &B);
    void invert(std::vector<std::vector<double>> &A, bool neg);
    double Deriv(const std::vector<double> input, const std::vector<double> params, int n);
    double dotProduct(std::vector<double> A, std::vector<double> B);
    std::vector<std::vector<double>> Identity(int num, double lambda);
    std::vector<double> vec_matrix_multi(std::vector<std::vector<double>> const &A,std::vector<double> const &B);
    void LevenbergMarquardt(const std::vector<std::vector<double>> &inputs, const std::vector<std::vector<double>> &outputs, std::vector<double> &params);
    std::vector<std::vector<double>> jacobian(const std::vector<std::vector<double>> &inputs, std::vector<double> params);
    std::vector<std::vector<double>> return_glow_curve(){
        return glow_curves;
    }
};

#endif /* First_Order_kinetics_hpp */
