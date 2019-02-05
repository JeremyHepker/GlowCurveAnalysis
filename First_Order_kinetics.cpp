//
//  First_Order_kinetics.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//
#include "File_Manager.hpp"
#include "First_Order_kinetics.hpp"
using namespace std;
First_Order_Kinetics::First_Order_Kinetics(std::pair<std::vector<int>,std::vector<double>> &data):count_data(data.second),temp_data(&data.first[1300],&data.first[1800]){
};
double First_Order_Kinetics::activation_energy(vector<double> &data){
    double m_g = 0.0,E = 0.0,C = 0.0, K = .00008617,b = 0.0,t = 0.0,d = 0.0,w = 0.0;
    
    vector<double>::iterator TM = max_element(data.begin(), data.end());
    int half_intensity = *TM/2;
    vector<double>::iterator TR = upper_bound(TM,data.end(),half_intensity,greater<double>());
    
    int TM_index  = int(TM-data.begin());
    int TR_index = int(TR-data.begin());
    int TL_index = int(TM_index - (TR-TM));
    
    double T_M = double(temp_data[TM_index])+273.15;
    double T_1 = double(temp_data[TL_index])+273.15;
    double T_2 = double(temp_data[TR_index])+273.15;
    
    t = T_M-T_1;
    if(t == 0){
        t = 1;
    }
    w = T_2-T_1;
    d = T_2-T_M;
    m_g = d/w;
    b = 1.58+4.2*(m_g -0.42);
    C = 1.51+(3*(m_g -0.42));
    E = ((C*K*(T_M*T_M))/t) - (b*(2*K*T_M));
    return E;
};

vector<vector<double>> First_Order_Kinetics::glow_curve(){
    File_Manager manage("/Users/jeremyhepker/Documents/NERS499/GlowCurveAnalsys/GlowCurveAnalsys/test_1.csv");
    double threshold = 0.9;
    vector<vector<double>> sum;
    vector<double> curve = count_data;
    double E = activation_energy(curve);
    auto a = max_element(count_data.begin(), count_data.end());
    int max_index = int(a - count_data.begin());
    threshold = (*a) * threshold;
    int count = 0;
    while (*a > threshold && *a > 1){
        vector<double> temp =kernal(max_index,E,temp_data);
        sum.push_back(temp);
        transform(curve.begin(), curve.end(), temp.begin(), curve.begin(), minus<double>());
        a = max_element(curve.begin(), curve.end());
        max_index = int(a - curve.begin());
        E = activation_energy(curve);
        count += 1;
    }
    return sum;
};

vector<double> First_Order_Kinetics::kernal(int max_index,double E,vector<int> &range){
    double k = pow(8.617,-5);
    vector<double> peak;
    peak.reserve(range.size());
    double Tm = double(temp_data[max_index])+273.15;
    double Im = double(count_data[max_index]);
    double I_t = 0.0, T = 0.0;
    double dm = (2.0*k*Tm)/E;
    for(int i = 0;i <= range.size();++i){
        T = double(range[i])+273.15;
        I_t = Im*exp(1.0 +(E/(k*T))*((T-Tm)/Tm)-((T*T)/(Tm*Tm))*exp((E/(k*T))*((T-Tm)/Tm))*(1.0-((2.0*k*T)/E))-dm);
        peak.push_back(I_t);
    };
    return peak;
}
