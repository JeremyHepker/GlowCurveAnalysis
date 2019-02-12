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
First_Order_Kinetics::First_Order_Kinetics(std::pair<std::vector<int>,std::vector<double>> data):count_data(data.second),temp_data(data.first){
};
double First_Order_Kinetics::activation_energy(vector<double> &data, bool first){
    vector<double>::iterator TR = data.end();
    double m_g = 0.0,E = 0.0,C = 0.0, K = .000086173303,b = 0.0,t = 0.0,d = 0.0,w = 0.0;
    vector<double>::iterator TM = max_element(data.begin(), data.end());
    double half_intensity = *TM/2.0;
    if(first){
        TR = upper_bound(TM,data.end(),half_intensity,greater<double>());
    }else{
        TR = upper_bound(data.begin(),TM,half_intensity,less<double>());
    }
    int TM_index = int(TM-data.begin());
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
    vector<vector<double>> residual;
    vector<double> curve = count_data;
    bool first =true;
    int c = 0;
    while (c< 4){
        double E = activation_energy(curve,first);
        int max = int(max_element(curve.begin(), curve.end()) - curve.begin());
        double Tm = double(temp_data[max])+273.15;
        double Im = double(count_data[max]);
        double min_E=E,min_T = Tm, min_t_error = INFINITY,min_E_error =INFINITY;
        double min_Im = Im, dm=0.0;
        vector<double> temp;
        int half_intensity = Im/2;
        vector<double>::iterator TR = lower_bound(count_data.begin()+max,count_data.end(),half_intensity,greater<double>());
        int TR_index = int(TR-count_data.begin());
        int TL_index = int(max - (TR_index-max));
        
        for(double inc = double(temp_data[TL_index]); inc < temp_data[TR_index]; inc += .1){
            vector<int>::iterator t_Im = lower_bound(temp_data.begin()+TL_index,temp_data.begin()+TR_index,int(inc),less<int>());
            double temp_Im = double(count_data[int(t_Im - temp_data.begin())]);
            dm = (2.0*k*(inc+273.15))/min_E;
            temp =kernal(min_E,inc+273.15,temp_Im,dm);
            vector<double> r(temp.size(), 0);
            for(int i = TL_index; i < TR_index; ++i){
                r[i] += temp[i] - count_data[i];
            }
            double sum = 0.0;
            for(auto x:r){
                sum += abs(pow(x,2.0));
            }
            if(sqrt(sum) < min_t_error){
                min_t_error =sqrt(sum);
                min_Im = temp_Im;
                min_T = inc+273.15;
            }
            //residual.push_back(temp);
        }
        dm = (2.0*k*min_T)/min_E;
        temp =kernal(min_E,min_T,min_Im,dm);
        //residual.push_back(temp);
         
        for(double inc = E*(.20); inc < E*(1.80); inc += .01){
            //Calculate the residual function r(x)
            double dm = (2.0*k*Tm)/min_E;
            temp =kernal(inc,Tm,Im,dm);
            vector<double> r(temp.size(), 0);
            for(int i = TL_index; i < TR_index; ++i){
                r[i] += temp[i] - count_data[i];
            }
            double sum = 0.0;
            for(auto x:r){
                sum += pow(x,2.0);
            }
            if(sqrt(sum) < min_E_error){
                min_E_error = sqrt(sum);
                min_E = inc;
            }
            if(c ==1) residual.push_back(temp);
        }
        //E = gauss_newton(temp,max_index,E,Tm,Im);
        //push back the iteration
        dm = (2.0*k*min_T)/min_E;
        temp =kernal(min_E,min_T,min_Im,dm);
        if(c ==1)residual.push_back(temp);
        transform(curve.begin()+TL_index,curve.end(),temp.begin()+TL_index,curve.begin()+TL_index,minus<double>());
        //if(c ==1)residual.push_back(curve);
        first = false;
        ++c;
    }
    return residual;
};


/*transform(curve.begin(),curve.end(),temp.begin(),curve.begin(),minus<double>());
 E = activation_energy(curve);
 a = max_element(curve.begin(), curve.end());
 max_index = int(a - curve.begin());*/
vector<double> First_Order_Kinetics::kernal(double E,double Tm,double Im,double dm){
    vector<double> peak;
    double T=0.0;
    double I_t = 0.0;
    peak.reserve(temp_data.size());
    for(int i = 0;i <= temp_data.size();++i){
        T = double(temp_data[i])+273.15;
        I_t = Im*exp(1.0 +(E/(k*T))*((T-Tm)/Tm)-((T*T)/(Tm*Tm))*exp((E/(k*T))*((T-Tm)/Tm))*(1.0-((2.0*k*T)/E))-dm);
        peak.push_back(I_t);
    };
    return peak;
}





//-------------------------GAUSS NEWTON METHOD----------------------------------------//

double First_Order_Kinetics::gauss_newton(vector<double> const &temp,int max_index,double E,double Tm,double Im){
    double E_n = E;
    
    int half_intensity = Im/2;
    vector<double>::iterator TR = lower_bound(count_data.begin()+max_index,count_data.end(),half_intensity,greater<double>());
    int TR_index = int(TR-count_data.begin());
    int TL_index = int(max_index - (TR_index-max_index));
    vector<double> r(temp.size(), 0);
    transform(temp.begin()+TL_index,temp.end(),count_data.begin()+TL_index,r.begin()+TL_index,minus<double>());

    //Compute the Jacobian matrix which is a 2x1800 matrix
    /*
     {dr/dE_1, dr/dE_1, ...,dr/dE_1800
     dr/dT_1, dr/dT_1, ...,dr/dT_1800]
     */
    vector<vector<double>> jacobian_1 = jacobian(max_index,TL_index,TR_index,E,Tm,Im);

    //Reserve space and take the transpose of the jacobian.
    vector<vector<double>> jacobian_2(jacobian_1[0].size(), vector<double>(jacobian_1.size(),0));
    transpose(jacobian_1, jacobian_2, int(jacobian_1.size()), int(jacobian_1[0].size()));

    //calculate the gradient, for comparing to the tolerance
    // if ||g(x)|| = ||J(x)^T*r(x)|| <= epsilon : stop;
    vector<double> gradient(jacobian_2.size(),0);
    double ssq = 0.0;
    for (int i=0;i<int(jacobian_2.size());i++){
        for (int j=0;j<int(jacobian_2[0].size());j++){
            gradient[i] += (jacobian_2[i][j]*r[j+TL_index]);
            ssq+=pow(gradient[i],2.0);
        }
    };
    //Calculate the transpose of the jacobian times the jacobian
    //Resulting in 1800x1800 matrix
    jacobian_1 = multiply(jacobian_2, jacobian_1);

    //Invert the matrix just calculated
    //  ---HERE IS THE PROBLEM---
    invert(jacobian_1);

    //calculate s_k for the next iteration.
    // x_k+1 = x_k + s_k
    double sk = 0.0;
    vector<double> s(jacobian_1.size(),0);
    for (int i=0;i<int(jacobian_1.size());i++){
        for (int j=0;j<int(jacobian_1[0].size());j++){
            s[i] += -(jacobian_1[i][j]*gradient[j]);
            sk+=pow(s[i],2.0);
        }
    };
    E_n += sqrt(sk);
    for(auto x: jacobian_1){
        x.clear();
    };
    for(auto x: jacobian_2){
        x.clear();
    };
    return E_n;
}
//---------------------GAUSS NEWTON METHOD HELPER FUNCTIONS----------------------//

vector<vector<double>> First_Order_Kinetics::jacobian(int max_index,int TL_index,int TR_index,double E,double Tm,double Im){
    double T=0.0;
    vector<double> dT, dE;
    for(int i = TL_index;i <= TR_index;++i){
        T = double(temp_data[i])+273.15;
        double d_dT = exp(1.0-(2.0*k*Tm)/E+(E*(-Tm+T))/(k*Tm*T)-(exp((E*(-Tm+T))/(k*Tm*T))*(T*T)*(1.0-(2.0*k*T)/E))/(Tm*Tm))* ((2.0*k*Tm)/(E*E)-(2.0*exp((E*(-Tm+T))/(k*Tm*T))*k*(T*T*T))/((E*E)*(Tm*Tm))+(-Tm+T)/(k*Tm*T)-(exp((E*(-Tm+T))/(k*Tm*T))*T*(-Tm+T)*(1.0-(2.0*k*T)/E))/(k*(Tm*Tm*Tm)))*Im;
    
        if (d_dT  == +1.0/0.0 || d_dT  == -1.0/0.0){
            d_dT = 0.0;
        }
        double d_dE = E*(-(2*k*(T*T*T)*exp((E*(T - Tm))/(k*Tm*T)))/((E*E)*(Tm*Tm)) + (2*k*Tm)/(E*E) - (T*(T - Tm)*(1 - (2*k*T)/E)*exp((E*(T - Tm))/(k*Tm*T)))/(k*(Tm*Tm*Tm)) + (T - Tm)/(k*Tm*T))*exp(-((T*T)*(1 - (2*k*T)/E)*exp((E*(T - Tm))/(k*Tm*T)))/(Tm*Tm) + (E*(T - Tm))/(k*Tm*T) - (2*k*Tm)/E + 1);
        if (d_dE  == +1.0/0.0 || d_dE  == -1.0/0.0){
            d_dE = 0.0;
        }
        dE.push_back(d_dE);
        dT.push_back(d_dT);
     };
    vector<vector<double>> jacob;
    jacob.push_back(dT);
    jacob.push_back(dE);
    return jacob;
}

void First_Order_Kinetics::transpose(vector<vector<double>> const &A,vector<vector<double>> &B, int n, int m){
    for(auto i = B.begin(); i != B.end();++i) i->resize(A.size());
    for (int i=0; i < n; i++) {
        for (int j=0; j < m; j++) {
            B[j][i] = A[i][j];
        }
    }
}
vector<vector<double>> First_Order_Kinetics::multiply(vector<vector<double>> const &A,vector<vector<double>> const &B){
    size_t row_A = A.size();
    size_t col_A = A[0].size();
    size_t col_B = B[0].size();
    double stop = 1* pow(double(10),double(-35));
    vector<vector<double>> C(row_A, vector<double>(col_B,0));
    for (size_t i = 0; i <row_A; i++){
        for (size_t j = 0; j < col_B; j++){
            for (size_t k = 0; k < col_A; k++){
                C[i][j] += A[i][k] * B[k][j];
                ((C[i][j] < stop) ? C[i][j] = 0.0 :false);
            }
        }
    }
    return C;
};

void First_Order_Kinetics::invert(vector<vector<double>> &A){
    int n = int(A.size());
    for(int i = 0; i < n; i++){
        for(int j = n; j < 2*n; j++){
            if(i==(j-n))
                A[i][j] = 1.0;
            else
                A[i][j] = 0.0;
        }
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i!=j){
                double ratio = 0.0;
                if(A[i][i] != 0)ratio = A[j][i]/A[i][i];
                for(int k = 0; k < 2*n; k++){
                    A[j][k] -= ratio * A[i][k];
                }
            }
        }
    }
    for(int i = 0; i < n; i++){
        double a = A[i][i];
        for(int j = 0; j < 2*n; j++){
            if(a != 0) A[i][j] /= a;
        }
    }
}
