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
#include <cstdlib>
class First_Order_Kinetics{
private:
    std::vector<int> temp_data;
    std::vector<double> count_data;
    double k = .000086173303;
public:
    First_Order_Kinetics(std::pair<std::vector<int>,std::vector<double>>);
    double activation_energy(std::vector<double> &data, bool first);
    double frequency_factor();
    std::vector<std::vector<double>> glow_curve();
    double gauss_newton(std::vector<double> const &temp,int max_index,double E,double Tm,double Im);
    std::vector<double> kernal(double E,double Tm,double Im,double dm);
    std::vector<std::vector<double>> jacobian(int max_index,int TL_index,int TR_index,double E,double Tm,double Im);
    void transpose(std::vector<std::vector<double>> const &A,std::vector<std::vector<double>> &B, int n, int m);
    std::vector<std::vector<double>> multiply(std::vector<std::vector<double>> const &A,std::vector<std::vector<double>> const &B);
    void invert(std::vector<std::vector<double>> &A);

};

#endif /* First_Order_kinetics_hpp */
