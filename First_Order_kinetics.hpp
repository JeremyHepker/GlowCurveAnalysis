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
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>

class First_Order_Kinetics{
private:
    std::vector<int> temp_data;
    std::vector<double> count_data;
    std::vector<std::pair<int, int>> temp_range = {{450,500},{420,470},{390,450},{300,450}};
public:
    First_Order_Kinetics(std::pair<std::vector<int>,std::vector<double>> &data);
    double activation_energy(std::vector<double> &data);
    double frequency_factor();
    std::vector<std::vector<double>> glow_curve();
    std::vector<double> kernal(int max_index,double E,std::vector<int> &range);
};

#endif /* First_Order_kinetics_hpp */
