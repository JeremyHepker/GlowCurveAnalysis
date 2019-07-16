//
//  First_Order_kinetics.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "Levenberg–Marquardt.hpp"
//#include "smartPeakDetect.hpp"

using namespace std;

First_Order_Kinetics::First_Order_Kinetics(std::pair<std::vector<double>,std::vector<double>> data, std::vector<std::vector<double>> peakParams):count_data(data.second),temp_data(data.first), peakParams(peakParams){};

/*---------------------------------Main Deconvolution---------------------------------------*/
double First_Order_Kinetics::glow_curve(){
    //int peak = 0;
    double FOM = 1.0;
    //curve = count_data;
    //const double curve_sum_start = std::accumulate(curve.begin(), curve.end(), 0);
    //double curve_sum = curve_sum_start*2;
    //vector<vector<double>> param_list;
    //findPeaks(temp_data,count_data,param_list);
    //vector<vector<double>> __tempPeaks__;
    vector<double> sum(curve.size(), 0.0);
    //int main_peak_limit = 0;
    double integral = 0.0;
    cout<<".";
    cout.flush();
//    while(curve_sum > curve_sum_start *.03 && int(param_list.size())<4){
//        param_list.push_back(initial_guess(curve,peak));
//        vector<double> params = {param_list[peak][0],param_list[peak][1],double(param_list[peak][2]),double(param_list[peak][3])};
//        //LevenbergMarquardt(curve, params);
//        auto t_Im = upper_bound(temp_data.begin(),temp_data.end(),params[1],less<int>());
//        double Im = double(curve[int(t_Im - temp_data.begin())]);
//        param_list[peak] = params;
//        param_list[peak].insert(param_list[peak].begin()+2,Im);
//        __tempPeaks__.push_back(vector<double>(temp_data.size(),0));
//        for(int i = 0; i < int(temp_data.size()); ++i ){
//            double output = Func(temp_data[i],param_list[peak]);
//            __tempPeaks__[peak][i] = output;
//        }
//        auto begin = __tempPeaks__[peak].begin();
//        auto end = __tempPeaks__[peak].end();
//        transform(begin,end,sum.begin(),sum.begin(),plus<double>());
//        transform(curve.begin(),curve.end(),sum.begin(),curve.begin(),minus<double>());
//        int TM_index = (param_list[peak][4] + param_list[peak][3])/2;
//        for(int i = 0; i < int(temp_data.size());++i){
//            if(curve[i]<0.0) curve[i] = 0.0;
//            if(i > TM_index && peak == 0) curve[i] = 0.0;
//            if(i > TM_index && i < param_list[peak][4]) curve[i] = 0.0;
//        }
//        if(peak == 0) main_peak_limit = int(param_list[peak][4]);
//        peak++;
//        curve_sum = accumulate(curve.begin(), curve.end(), 0);
//    }
//    for(int i = main_peak_limit; i < int(curve.size()); ++i){
//        if(__tempPeaks__[0][i] < count_data[i]) count_data[i] = __tempPeaks__[0][i];
//    }
    curve = count_data;
    //vector<vector<double>> peakPrams;
    //findPeaks(count_data, temp_data, peakPrams);
    cout<<".";
    cout.flush();
    LevenbergMarquardt(curve, peakParams, FOM);
    if(FOM > 1.0){
        cout<<"."<<endl<<"----- Levenberg-Marquardt failed to converged -----"<<endl;
        cout<<"data most likely contains to much noise, or is improperly formatted"<<endl;
        return -1;
    }
    cout<<"."<<endl<<"----- Levenberg-Marquardt converged to a FOM of "<<(FOM*100)<<"% -----"<<endl;
    integral = 0.0;
    sum = vector<double>(temp_data.size(),0.0);
    vector<double> peak_areas = vector<double>(4, 0.0);
    for(int i = 0;i < int(peakParams.size()); ++i){
        glow_curves.push_back(vector<double>(temp_data.size(),0.0));
    }
    for(int i = 0; i < int(temp_data.size()); ++i ){
        double output = 0.0;
        for(int x = 0; x < int(peakParams.size()); ++x){
            double out = Func2(temp_data[i],peakParams[x]);
            peak_areas[x] += out;
            output += out;
            glow_curves[x][i] = out;
        }
        sum[i] = output;
        integral += output;
    }
    glow_curves.push_back(sum);
    glow_curves.push_back(count_data);
    for(int i = 0; i < int(peakParams.size()); ++i){
        cout<<"----- Area Under Curve #"<<i+1<<" :"<<peak_areas[i]<<" -----"<<endl;
    }
    curve_areas = peak_areas; 
    return FOM;
};

//vector<double> First_Order_Kinetics::initial_guess(vector<double> &curve, int j){
//    vector<double>::iterator TR = curve.end(),TL = curve.begin(),TM = curve.end();
//    int TM_index = 0, TR_index = 0, TL_index = 0;
//    double half_intensity = 0;
//    //Find the range of the curve
//    TL = curve.begin();
//    TR = curve.end();
//    TM = max_element(TL, TR);
//    half_intensity = *TM/2.0;
//    TL = lower_bound(TL,TM,half_intensity,less<double>());
//    TR = lower_bound(TM,TR,half_intensity,greater<double>());
//    int diff1 = int(TM - TL);
//    int diff2 = int(TR - TM);
//    if(TR == curve.end()) diff2 += diff1;
//    TM_index = int(TM - curve.begin());
//    if(diff1 <= diff2){
//        TL_index = int(TL - curve.begin());
//        TR_index = int(TM_index + (TM_index - TL_index));
//    }else{
//        TR_index = int(TR - curve.begin());
//        TL_index = int(TM_index - (TR_index - TM_index));
//    }
//    if(TR_index > int(curve.size())){
//        TR_index = int(curve.size())-1;
//    }
//    if(TL_index < 0){
//        TL_index = 0;
//    }
//    double Tm = double(temp_data[TM_index]);
//    double E = activation_energy(TL_index,TM_index,TR_index);
//    vector<double> params = {E,Tm,double(TL_index),double(TR_index)};
//    return params;
//}

/*-------------------------First Order Kinetics Function--------------------------------*/

//double First_Order_Kinetics::Func(const double input, const vector<double> params){
//    double T=0.0;
//    double I_t = 0.0;
//    double E = params[0];
//    double Tm = params[1]+273.15;
//    double dm = (2.0*k*(Tm))/E;
//    auto t_Im = upper_bound(temp_data.begin(),temp_data.end(),Tm-273.15,less<int>());
//    double Im = double(curve[int(t_Im - temp_data.begin())]);
//    T = double(input+273.15);
//    I_t = Im*exp(1.0 +(E/(k*T))*((T-Tm)/Tm)-((T*T)/(Tm*Tm))*exp((E/(k*T))*((T-Tm)/Tm))*(1.0-((2.0*k*T)/E))-dm);
//    return I_t;
//}
double First_Order_Kinetics::Func2(const double input, const vector<double> params){
    double T=0.0;
    double I_t = 0.0;
    double E = params[0];
    double Tm = params[1]+273.15;
    double dm = (2.0*k*(Tm))/E;
    double Im = params[2];
    T = double(input+273.15);
    I_t = Im*exp(1.0 +(E/(k*T))*((T-Tm)/Tm)-((T*T)/(Tm*Tm))*exp((E/(k*T))*((T-Tm)/Tm))*(1.0-((2.0*k*T)/E))-dm);
    return I_t;
}
//-------------------------Levenburg Marquardt METHOD----------------------------------------//
//void First_Order_Kinetics::LevenbergMarquardt(const vector<double> &outputs, vector<double> &params){
//    int input_size = int(temp_data.size());
//    int num_params = 2;
//    vector<vector<double>> Jf_T(num_params, vector<double>(input_size,0.0));
//    vector<double> error(input_size,0.0);
//    vector<vector<double>> input(1, vector<double>(1,1));
//    vector<vector<double>> H;
//    double lambda = 0.0001;
//    double updateJ = 1;
//    double e = 0.0;
//    double E_est = params[0];
//    double Tm_est = params[1];
//    for(int i = 0; i < MAX_ITER; ++i){
//        if(updateJ == 1){
//            //Evaluate the jacobian matrix at the current paramater.
//            vector<double> t_params = {E_est,Tm_est};
//            for(int j=0; j < input_size; j++) {
//                for(int k = 0; k < 2; ++k){
//                    Jf_T[k][j] = Deriv(temp_data[j],t_params , k);
//                }
//                //Evaluate the distance error at the current paramters.
//                if( j>params[2] && j < params[3]) {
//                    error[j] = outputs[j] - Func(temp_data[j], t_params);
//                }
//            }
//            //Calculate the hessian mat
//            vector<vector<double>> Jf(Jf_T[0].size(), vector<double>(Jf_T.size(),0.0));
//            transpose(Jf_T,Jf,int(Jf_T.size()), int(Jf_T[0].size()));
//            H = multiply(Jf_T,Jf);
//            e = dotProduct(error, error);
//            
//        }
//        
//        //apply the damping factor to the hessian matrix
//        vector<vector<double>> I = Identity(num_params, lambda);
//        vector<vector<double>> H_lm(num_params, vector<double>(num_params,0.0));
//        for(int j = 0; j < int(H.size()); ++j){
//            for(int s = 0; s < int(H.size()); ++s){
//                H_lm[j][s] = H[j][s] + I[j][s];
//            }
//        }
//        
//        //Update parameters
//        invert(H_lm, true);
//        vector<double> Jf_error = vec_matrix_multi(Jf_T,error);
//        vector<double> delta = vec_matrix_multi(H_lm,Jf_error);
//        E_est += delta[1];
//        Tm_est += delta[0];
//        
//        //Evaluate the total distance error at the updated paramaters.
//        vector<double> temp_params;
//        vector<double> temp_error(input_size,0.0),temp_temp(input_size,0.0);
//        temp_params.push_back(E_est);
//        temp_params.push_back(Tm_est);
//        bool bad = false;
//        for(int j=params[2]; j < params[3]; j++){
//            temp_error[j] = outputs[j] - Func(temp_data[j], temp_params);
//        }
//        double temp_e = bad? e:dotProduct(temp_error, temp_error);
//        if(temp_e < e){
//            lambda /= 10;
//            params[0] = E_est;
//            params[1] = Tm_est;
//            e = temp_e;
//            updateJ = 1;
//        }else{
//            updateJ = 0;
//            lambda *= 10;
//        }
//    }
//}

//-------------------------Levenburg Marquardt 3 iteration METHOD----------------------------------------//
void First_Order_Kinetics::LevenbergMarquardt(const vector<double> &outputs, vector<vector<double>> &params, double &FOM){
    int input_size = int(temp_data.size());
    int num_params = int(params.size());
    int d = 1;
    int main_hold = 0;
    double main_FOM =FOM;
    while(FOM > .02){
        main_FOM =FOM;
        if(main_hold > 3){
            break;
        }
        for(int param_num = 0; param_num < 3 ; ++param_num){
            vector<double> temp_params;
            vector<double> temp_output(curve.size(), 0.0);
            for(int i = 0; i < params.size();++i){
                temp_params.push_back(params[i][param_num]);
            }
            int other_param1 = 0;
            int other_param2 = 0;
            if(param_num == 0){
                other_param1 = 1;
                other_param2 = 2;
            }else if(param_num == 1){
                other_param1 = 0;
                other_param2 = 2;
            }else{
                other_param1 = 0;
                other_param2 = 1;
            }
            vector<vector<double>> Jf_T(num_params, vector<double>(input_size,0.0));
            vector<double> error(input_size,0.0);
            vector<vector<double>> H;
            double lambda = 0.01;
            double updateJ = 1;
            double e = 0.0;
            int i = 0;
            int inner_hold = 0;
            while(FOM > .03 && i < 500){
                if(updateJ == 1){
                    //Evaluate the jacobian matrix at the current paramater.
                    vector<double>t_parms(3, 0.0);
                    double integral = 0.0;
                    for(int j=0; j < input_size; j++) {
                        double output = 0.0;
                        for(int k = 0; k < num_params;++k){
                            t_parms[param_num] = temp_params[k];
                            t_parms[other_param1] = params[k][other_param1];
                            t_parms[other_param2] = params[k][other_param2];
                            for(int x = 0; x < 1; x++){
                                Jf_T[k - x][j] = Deriv(temp_data[j],t_parms, x);
                            }
                            output += Func2(temp_data[j],t_parms);
                        }
                        temp_output[j] = output;
                        error[j] = outputs[j] - output;
                        integral += output;
                    }
                    FOM = 0.0;
                    for(int z = 0; z < int(curve.size()); ++z){
                        FOM += abs(curve[z] - temp_output[z])/integral;
                    }
                    //Calculate the hessian mat
                    vector<vector<double>> Jf(Jf_T[0].size(), vector<double>(Jf_T.size(),0.0));
                    transpose(Jf_T,Jf,int(Jf_T.size()), int(Jf_T[0].size()));
                    H = multiply(Jf_T,Jf);
                    e = dotProduct(error, error);
                }
                
                //apply the damping factor to the hessian matrix
                vector<vector<double>> I = Identity(num_params, lambda);
                vector<vector<double>> H_lm(num_params, vector<double>(num_params,0.0));
                for(int j = 0; j < int(H.size()); ++j){
                    for(int s = 0; s < int(H.size()); ++s){
                        H_lm[j][s] = H[j][s] + I[j][s];
                    }
                }
                invert(H_lm, true);
                vector<double> Jf_error = vec_matrix_multi(Jf_T,error);
                vector<double> delta = vec_matrix_multi(H_lm,Jf_error);
                vector<double> t_params = temp_params;
                for(int x = 0; x < delta.size(); ++x){
                    t_params[x] += delta[x];
                }
                double integral = 0.0;
                //Evaluate the total distance error at the updated paramaters.
                vector<double> temp_error(input_size,0.0);
                vector<double> t_param(3,0.0);
                for(int j=0; j < input_size; j++){
                    double output = 0.0;
                    for(int k = 0; k < num_params;++k){
                        t_param[param_num] = t_params[k];
                        t_param[other_param1] = params[k][other_param1];
                        t_param[other_param2] = params[k][other_param2];
                        double peak = Func2(temp_data[j],t_param);
                        output += peak;
                    }
                    temp_output[j] = output;
                    integral += output;
                    temp_error[j] = outputs[j] - output;
                }
                double temp_FOM = 0.0;
                for(int z = 0; z < int(curve.size()); ++z){
                    temp_FOM += abs(curve[z] - temp_output[z])/integral;
                }
                double temp_e = dotProduct(temp_error, temp_error);
                if(temp_FOM < FOM){
                    lambda /= 10;
                    temp_params = t_params;
                    e = temp_e;
                    updateJ = 1;
                    inner_hold = 0;
                }else{
                    inner_hold += 1;
                    updateJ = 0;
                    lambda *= 10;
                }
                if(inner_hold > 25) i = 500;
                ++i;
            }
            for(int i = 0; i < int(temp_params.size());++i){
                params[i][param_num]= temp_params[i];
            }
        }
        if(abs(main_FOM - FOM) < (1e-4)){
            main_hold += 1;
        }else{
            main_hold = 0;
        }
        ++d;
        cout<<".";
        cout.flush();
    }
}
