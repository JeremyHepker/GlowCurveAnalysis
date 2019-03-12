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
First_Order_Kinetics::First_Order_Kinetics(std::pair<std::vector<int>,std::vector<double>> data):count_data(data.second),temp_data(data.first.begin(),data.first.end()){

};
/*---------------------------------Activation Energy---------------------------------------*/
double First_Order_Kinetics::activation_energy(int TL_index,int TM_index,int TR_index){
    double m_g = 0.0,E = 0.0,C = 0.0, K = .000086173303,b = 0.0,t = 0.0,d = 0.0,w = 0.0;
    
    //Calculate Max temperature and half widths
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

/*---------------------------------Main Deconvolution---------------------------------------*/
void First_Order_Kinetics::glow_curve(){
    curve = count_data;
    vector<double>::iterator TR = curve.end(),TL = curve.begin(),TM = curve.end();
    int TM_index = 0, TR_index = 0, TL_index = 0, c = 0;
    double curve_sum = std::accumulate(curve.begin(), curve.end(), 0);
    double curve_sum_start=curve_sum;
    vector<double> sum(curve.size(), 0.0);
    vector<vector<double>> updated_params;
    vector<vector<double>> inputs(int(temp_data.size()), vector<double>(1,1.0));
    for(int i=0; i < int(temp_data.size()); i++) {
        inputs[i][0] = double(temp_data[i]);
    }
    double FOM = 1;
    
    while(FOM == 1){
        while (curve_sum > curve_sum_start*.05){
            //Find the range of the curve
            TL = curve.begin();
            TR = curve.end();
            TM = max_element(TL, TR);
            double half_intensity = *TM/2.0;
            TL = lower_bound(TL,TM,half_intensity,less<double>());
            TR = lower_bound(TM,TR,half_intensity,greater<double>());
            int diff1 = int(TM - TL);
            int diff2 = int(TR - TM);
            if(TR == curve.end()) diff2 += diff1;
            if(diff1 <= diff2){
                TM_index = int(TM - curve.begin());
                TL_index = int(TL - curve.begin());
                TR_index = int(TM_index + (TM_index - TL_index));
            }else{
                TM_index = int(TM - curve.begin());
                TR_index = int(TR - curve.begin());
                TL_index = int(TM_index - (TR_index - TM_index));
            }
            if(TR_index > int(curve.size())){
                TR_index = int(curve.size())-1;
            }
            //Caluclate Initial Guess for both variable
            double Tm = double(temp_data[TM_index]);
            double E = activation_energy(TL_index,TM_index,TR_index);
            
            vector<vector<double>> outputs(int(temp_data.size()),vector<double>(1,1.0));
            
            for(int i=0; i < int(temp_data.size()); i++) {
                outputs[i][0] = curve[i];
            }
            
            // Guess the parameters, it should be close to the true value, else it can fail
            vector<double> params = {E,Tm,double(TL_index),double(TR_index)};
            cout<<"----------------Curve: "<<c+1<<"----------------"<<endl;
            LevenbergMarquardt(inputs,outputs, params);
            
            glow_curves.push_back(vector<double>(temp_data.size(),0));
            for(int i = 0; i < int(temp_data.size()); ++i ){
                vector<double> input(1,0.0);
                input[0]= double(temp_data[i]);
                glow_curves[c][i] = Func(input,params);
            }
            transform(glow_curves[c].begin(),glow_curves[c].end(),sum.begin(),sum.begin(),plus<double>());
            transform(curve.begin(),curve.end(),sum.begin(),curve.begin(),minus<double>());
            
            for(int i = 0; i < int(temp_data.size());++i){
                if(curve[i]<0.0) curve[i] = 0.0;
                if(i > TM_index && c == 0) curve[i] = 0.0;
                if(i > TM_index && i < TR_index) curve[i] = 0.0;
            }
            ++c;
            curve_sum = accumulate(curve.begin(), curve.end(), 0);
        }
        double integral = accumulate(sum.begin(), sum.end(), 0);
        for(int i = 0; i < int(curve.size()); ++i){
            FOM += abs(count_data[i] - sum[i])/integral;
        }
        cout<<FOM<<endl;
        glow_curves.push_back(sum);
    }
};

/*-------------------------First Order Kinetics Function--------------------------------*/
double First_Order_Kinetics::Func(const vector<double> input, const vector<double> params){
    double T=0.0;
    double I_t = 0.0;
    double E = params[0];
    double Tm = params[1]+273.15;
    double dm = (2.0*k*(Tm))/E;
    vector<double>::iterator t_Im = upper_bound(temp_data.begin(),temp_data.end(),Tm-273.15,less<int>());
    double Im = double(curve[int(t_Im - temp_data.begin())]);
    T = double(input[0]+273.15);
    I_t = Im*exp(1.0 +(E/(k*T))*((T-Tm)/Tm)-((T*T)/(Tm*Tm))*exp((E/(k*T))*((T-Tm)/Tm))*(1.0-((2.0*k*T)/E))-dm);
    return I_t;
}
//-------------------------Levenburg Marquardt METHOD----------------------------------------//
void First_Order_Kinetics::LevenbergMarquardt(const vector<vector<double>> &inputs, const vector<vector<double>> &outputs, vector<double> &params){
    int input_size = int(inputs.size());
    int num_params = 2;
    vector<vector<double>> Jf_T(num_params, vector<double>(input_size,0.0)); // Jacobian of Func()
    vector<double> error(input_size,0.0);
    vector<vector<double>> input(1, vector<double>(1,1)); // single row input
    vector<vector<double>> H;
    double lambda = 0.01;
    double updateJ = 1;
    double e = 0.0;
    double E_est = params[0];
    double Tm_est = params[1];
    for(int i = 0; i < MAX_ITER; ++i){
        if(updateJ == 1){
            //Evaluate the jacobian matrix at the current paramater.
            vector<double> t_params = {E_est,Tm_est};
            //Jf_T = jacobian(inputs, params);
            for(int j=0; j < input_size; j++) {
                input[0][0] = inputs[j][0];
                for(int k = 0; k < 2; ++k){
                    Jf_T[k][j] = Deriv(input[0],t_params , k);
                }
                //Evaluate the distance error at the current paramters.
                if(j > params[2] && j < params[3]){
                    error[j] = outputs[j][0] - Func(input[0], t_params);
                }else{
                    error[j] = 0.0;
                }
            }
            //Calculate the hessian mat
            vector<vector<double>> Jf(Jf_T[0].size(), vector<double>(Jf_T.size(),0.0));
            transpose(Jf_T,Jf,int(Jf_T.size()), int(Jf_T[0].size()));
            H = multiply(Jf_T,Jf);
            if(i == 0){
                e = dotProduct(error, error);
            }
        }
        
        //apply the damping factor to the hessian matrix
        vector<vector<double>> I = Identity(num_params, lambda);
        vector<vector<double>> H_lm(num_params, vector<double>(num_params,0.0));
        for(int j = 0; j < int(H.size()); ++j){
            for(int s = 0; s < int(H.size()); ++s){
                H_lm[j][s] = H[j][s] + I[j][s];
            }
        }
        
        //Update parameters
        invert(H_lm, true);
        vector<double> Jf_error = vec_matrix_multi(Jf_T,error);
        vector<double> delta = vec_matrix_multi(H_lm,Jf_error);
        E_est += delta[1];
        Tm_est += delta[0];
        
        //Evaluate the total distance error at the updated paramaters.
        vector<double> temp_params;
        vector<double> temp_error(input_size,0.0);
        temp_params.push_back(E_est);
        temp_params.push_back(Tm_est);
        bool bad = false;
        for(int j=params[2]; j < params[3]; j++) {
            input[0][0] = inputs[j][0];
            double err = outputs[j][0] - Func(input[0], temp_params);
            temp_error[j] = err;
        }
        
        double temp_e = bad? e:dotProduct(temp_error, temp_error);
        
        if(temp_e < e){
            lambda /= 10;
            params[0] = E_est;
            params[1] = Tm_est;
            e = temp_e;
            updateJ = 1;
            vector<double> temp_outputs;
            for(int j = 0; j < int(temp_data.size()); ++j ){
                vector<double> input;
                input.push_back(double(temp_data[j]));
                temp_outputs.push_back(Func(input,params));
            }
            //glow_curves.push_back(temp_outputs);
        }else{
            updateJ = 0;
            lambda *= 10;
        }
    }
}

//---------------------Levenburg Marquardt Method Helper Functions----------------------//
double First_Order_Kinetics::Deriv(const vector<double> input, const vector<double> params, int n)
{
    // Assumes input is a single row matrix
    
    // Returns the derivative of the nth parameter
    vector<double> params1 = params;
    vector<double> params2 = params;
    double d = 0;
    // Use central difference  to get derivative
    params1[n] -= DERIV_STEP;
    params2[n] += DERIV_STEP;
    double p1 = Func(input, params1);
    double p2 = Func(input, params2);
    d = (p2 - p1) / (2*DERIV_STEP);
    
    return d;
}

//Function to transpose a matrix
void First_Order_Kinetics::transpose(vector<vector<double>> const &A,vector<vector<double>> &B, int n, int m){
    for(auto i = B.begin(); i != B.end();++i) i->resize(A.size());
    for (int i=0; i < n; i++) {
        for (int j=0; j < m; j++) {
            B[j][i] = A[i][j];
        }
    }
}
//Functtion to multiply two matricies together
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
//Function to muliply an array with a matrix
vector<double> First_Order_Kinetics::vec_matrix_multi(vector<vector<double>> const &A,vector<double> const &B){
    vector<double> output;
    for(auto i = A.begin(); i != A.end();++i){
        double sum = 0.0;
        for(int j = 0; j != int(B.size()); ++j){
            sum += i->at(j) * B[j];
        }
        output.push_back(sum);
    }
    return output;
}
//Function to invert a matrix, with option for negtive inverse.
void First_Order_Kinetics::invert(vector<vector<double>> &A, bool neg){
    double det = A[0][0] *A[1][1] - A[1][0] * A[0][1];
    vector<vector<double>> temp(2, vector<double>(2,0));
    if(neg){
        det *= (-1.0);
    }
    temp[0][0] = (1/det)*A[1][1];
    temp[0][1] = -(1/det)*A[0][1];
    temp[1][0] = -(1/det)*A[1][0];
    temp[1][1] = (1/det)*A[0][0];
    A = temp;
}
//Functin to preform the dot product between two vectors
double First_Order_Kinetics::dotProduct(vector<double> A, vector<double> B){
    double product = 0;
    // Loop for calculate cot product
    for (int i = 0; i < int(A.size()); i++)
        product = product + A[i] * B[i];
    return product;
}
//Function to create an identity matrix multiplied by give lambda
vector<vector<double>> First_Order_Kinetics::Identity(int num, double lambda){
    int row, col;
    vector<vector<double>> I(num, vector<double>(num,0));
    for (row = 0; row < num; row++){
        for (col = 0; col < num; col++){
            // Checking if row is equal to column
            if (row == col)
                I[row][col] = lambda;
        }
    }
    return I;
}
//Function for finding the Jacobian
vector<vector<double>> First_Order_Kinetics::jacobian(vector<vector<double>> const &inputs, vector<double> params){
    double T=0.0, E = params[0],Tm = params[1] + 273.15;
    vector<double> dT, dE;
    vector<double>::iterator t_Im = upper_bound(temp_data.begin(),temp_data.end(),Tm- 273.15,less<double>());
    double Im = double(curve[t_Im - temp_data.begin()]);
    
    for(int i = 0;i < inputs.size();++i){
        T = double(inputs[i][0])+273.15;
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
    jacob.push_back(dE);
    jacob.push_back(dT);
    return jacob;
}

