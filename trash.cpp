//
//  trash.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 2/11/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include <stdio.h>
//Compute the Jacobian matrix which is a 2x1800 matrix
/*
 {dr/dE_1, dr/dE_1, ...,dr/dE_1800
 dr/dT_1, dr/dT_1, ...,dr/dT_1800]
 */
/*vector<vector<double>> jacobian_1 = jacobian(max_index,E,Tm,Im,I_t);
 
 //Reserve space and take the transpose of the jacobian.
 vector<vector<double>> jacobian_2(jacobian_1[0].size(), vector<double>(jacobian_1.size(),0));
 transpose(jacobian_1, jacobian_2, int(jacobian_1.size()), int(jacobian_1[0].size()));
 
 //calculate the gradient, for comparing to the tolerance
 // if ||g(x)|| = ||J(x)^T*r(x)|| <= epsilon : stop;
 vector<double> gradient(jacobian_2.size(),0);
 cout << setprecision(25);
 for (int i=0;i<int(jacobian_2.size());i++){
 for (int j=0;j<int(jacobian_2[0].size());j++){
 gradient[i] += (jacobian_2[i][j]*temp[j]);
 }
 };
 
 //Calculate the transpose of the jacobian times the jacobian
 //Resulting in 1800x1800 matrix
 jacobian_2 = multiply(jacobian_2, jacobian_1);
 
 //Invert the matrix just calculated
 //  ---HERE IS THE PROBLEM---
 vector<vector<double>> jacobian_3;
 jacobian_3 = jacobian_2;
 //invert(jacobian_3);
 
 //calculate s_k for the next iteration.
 // x_k+1 = x_k + s_k
 vector<double> s(jacobian_3.size(),0);
 for (int i=0;i<int(jacobian_3.size());i++){
 for (int j=0;j<int(jacobian_3[0].size());j++){
 s[i] += (jacobian_3[i][j]*gradient[j]);
 }
 };*/


/*void GaussNewton(double(*Func)(const vector<double> input, const vector<double> params), const vector<vector<double>> inputs, const vector<vector<double>> outputs, vector<double> &params)
 {
 int m = int(inputs.size());
 int n = int(inputs[0].size());
 int num_params = int(params.size());
 
 vector<vector<double>> r(m, vector<double>(1,1.0)); // residual matrix
 vector<vector<double>> Jf(num_params, vector<double>(m,1)); // Jacobian of Func()
 vector<vector<double>> input(1, vector<double>(n,1)); // single row input
 
 double last_mse = 0;
 
 for(int i=0; i < MAX_ITER; i++) {
 double mse = 0;
 
 for(int j=0; j < m; j++) {
 for(int k=0; k < n; k++) {
 input[0][k] = inputs[j][k];
 }
 
 r[j][0] = outputs[j][0] - Func(input[0], params);
 
 mse += r[j][0]*r[j][0];
 
 for(int k=0; k < num_params; k++) {
 Jf[k][j] = Deriv(Func, input[0], params, k);
 }
 }
 
 mse /= m;
 
 // The difference in mse is very small, so quit
 if(fabs(mse - last_mse) < 1e-8) {
 break;
 }
 
 //vector<double> delta = ((Jf.t()*Jf)).inv() * Jf.t()*r;
 vector<vector<double>> trans(Jf[0].size(), vector<double>(Jf.size(),0.0));
 transpose(Jf,trans,int(Jf.size()), int(Jf[0].size()));
 vector<vector<double>> trans_jf = multiply(Jf, trans);
 invert(trans_jf,false);
 vector<vector<double>> trans_r = multiply(Jf,r);
 vector<vector<double>> delta = multiply(trans_jf,trans_r);
 for(int z = 0; z < params.size(); ++z){
 params[z] += delta[z][0];
 }
 //printf("%d: mse=%f\n", i, mse);
 //printf("%d %f\n", i, mse);
 
 last_mse = mse;
 }
 }*/

/*double First_Order_Kinetics::gauss_newton(vector<double> const &temp,int max_index,double E,double Tm,double Im){
    double E_n = E;
    
    int half_intensity = Im/2;
    vector<double>::iterator TR = lower_bound(count_data.begin()+max_index,count_data.end(),half_intensity,greater<double>());
    int TR_index = int(TR-count_data.begin());
    int TL_index = int(max_index - (TR_index-max_index));
    vector<double> r(temp.size(), 0);
    transform(temp.begin()+TL_index,temp.end(),count_data.begin()+TL_index,r.begin()+TL_index,minus<double>());
    
    //Compute the Jacobian matrix which is a 2x1800 matrix
    
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
}*/
/*
 vector<vector<double>> residual;
 vector<double> curve = count_data;
 bool first =true;
 int c = 0;
 while (c< 4){
 vector<double>::iterator TR = curve.end(),TL = curve.end();
 vector<double>::iterator TM = max_element(curve.begin(), curve.end());
 int TM_index, TR_index, TL_index;
 double half_intensity = *TM/2.0;
 if(first){
 TR = upper_bound(TM,curve.end(),half_intensity,greater<double>());
 TM_index = int(TM-curve.begin());
 TR_index = int(TR-curve.begin());
 TL_index = int(TM_index - (TR-TM));
 }else{
 TL = upper_bound(curve.begin(),TM,half_intensity,less<double>());
 TM_index = int(TM-curve.begin());
 TL_index = int(TL-curve.begin());
 TR_index = int(TM_index + (TM-TL));
 }
 double Tm = double(temp_data[TM_index])+273.15;
 double Im = double(count_data[TM_index]);
 double E = activation_energy(curve,TL_index,TM_index,TR_index);
 double min_t_error = INFINITY,min_E_error =INFINITY;
 double min_Im = Im, dm=0.0,min_E=E,min_T = Tm;
 vector<double> temp;
 
 for(double inc = double(temp_data[TL_index]); inc < temp_data[TR_index]; inc += .1){
 vector<int>::iterator t_Im = lower_bound(temp_data.begin()+TL_index,temp_data.begin()+TR_index,int(inc),less<int>());
 double temp_Im = double(count_data[int(t_Im - temp_data.begin())]);
 dm = (2.0*k*(inc+273.15))/min_E;
 temp =kernal(min_E,inc+273.15,temp_Im,dm);
 vector<double> r(temp.size(), 0);
 for(int i = TL_index; i < TR_index; ++i){
 r[i] += temp[i] - curve[i];
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
 }
 
 for(double inc = E*(.20); inc < E*(1.80); inc += .01){
 //Calculate the residual function r(x)
 double dm = (2.0*k*min_T)/min_E;
 temp =kernal(inc,Tm,Im,dm);
 vector<double> r(temp.size(), 0);
 for(int i = TL_index; i < TR_index; ++i){
 r[i] += temp[i] - curve[i];
 }
 double sum = 0.0;
 for(auto x:r){
 sum += pow(x,2.0);
 }
 if(sqrt(sum) < min_E_error){
 min_E_error = sqrt(sum);
 min_E = inc;
 }
 }
 //push back the iteration
 dm = (2.0*k*min_T)/min_E;
 temp =kernal(min_E,min_T,min_Im,dm);
 double damper = double(temp_data.back()/double(temp_data.size()));
 if(damper>0.4){
 damper = .01;
 }else if(damper<0.2){
 damper = .1;
 }else{
 damper = 0.1;
 }
 int plus = int(ceil(double(TL_index) * damper));
 transform(curve.begin()+TL_index-plus,curve.end(),temp.begin()+TL_index-plus,curve.begin()+TL_index-plus,minus<double>());
 //residual.push_back(temp);
 if( c==0)residual.push_back(temp);
 if(c == 4 )residual.push_back(curve);
 
 /*
 int plus = int(ceil(double(temp_data[TL_index]) * 0.01));
 plus = int(distance(temp_data.begin(), find(temp_data.begin(), temp_data.end(), temp_data[TL_index]-plus)));
 transform(curve.begin()+plus,curve.end(),temp.begin()+plus,curve.begin()+plus,minus<double>());
 
first = false;
++c;
}

return residual;
*/







for(int j=0; j < input_size; j++) {
    input[0][0] = inputs[j][0];
    for(int k=0; k < num_params; k++) {
        Jf[k][j] = Deriv(input[0], params, k);
    }
    //Evaluate the distance error at the current paramters.
    error[j] = outputs[j][0] - Func(input[0], params);
    }




void First_Order_Kinetics::glow_curve(){
    curve = count_data;
    double FOM = 1;
    bool first_iter = true;
    vector<vector<double>> temp_curves;
    vector<vector<double>> updated_params;
    vector<double> params(4,0.0);
    vector<vector<double>> inputs(int(temp_data.size()), vector<double>(1,1.0));
    for(int i=0; i < int(temp_data.size()); i++) {
        inputs[i][0] = double(temp_data[i]);
    }
    while(FOM > .13){
        int c = 0;
        FOM = 1;
        curve = count_data;
        double curve_sum = std::accumulate(curve.begin(), curve.end(), 0);
        double curve_sum_start=curve_sum;
        vector<double> sum(curve.size(), 0.0);
        while (curve_sum > curve_sum_start*.05){
            vector<vector<double>> outputs(int(temp_data.size()),vector<double>(1,1.0));
            for(int i=0; i < int(temp_data.size()); i++) {
                outputs[i][0] = curve[i];
            }
            
            // Guess the parameters, it should be close to the true value, else it can fail
            if(first_iter) params = initial_guess(curve);
            else params = updated_params[c];
            cout<<"----------------Curve: "<<c+1<<"----------------"<<endl;
            LevenbergMarquardt(inputs,outputs, params);
            
            if(first_iter)temp_curves.push_back(vector<double>(temp_data.size(),0));
            for(int i = 0; i < int(temp_data.size()); ++i ){
                vector<double> input(1,0.0);
                input[0]= double(temp_data[i]);
                temp_curves[c][i] = Func(input,params);
            }
            transform(temp_curves[c].begin(),temp_curves[c].end(),sum.begin(),sum.begin(),plus<double>());
            transform(curve.begin(),curve.end(),sum.begin(),curve.begin(),minus<double>());
            
            for(int i = 0; i < int(temp_data.size());++i){
                if(curve[i]<0.0) curve[i] = 0.0;
                int TM_index = params[2] + (params[3]-params[2]);
                if(i > TM_index && c == 0) curve[i] = 0.0;
                if(i > TM_index && i < params[3]) curve[i] = 0.0;
            }
            if(first_iter) updated_params.push_back(params);
            else updated_params[c] = params;
            ++c;
            curve_sum = accumulate(curve.begin(), curve.end(), 0);
        }
        double integral = accumulate(sum.begin(), sum.end(), 0);
        for(int i = 0; i < int(curve.size()); ++i){
            FOM += abs(count_data[i] - sum[i])/integral;
        }
        first_iter = false;
        cout<<FOM<<endl;
    }
    glow_curves = temp_curves;
    cout<<"Final FOM " << FOM<<endl;
};

/*-------------------------First Order Kinetics Function--------------------------------*/
vector<double> First_Order_Kinetics::initial_guess(vector<double> &curve){
    vector<double>::iterator TR = curve.end(),TL = curve.begin(),TM = curve.end();
    int TM_index = 0, TR_index = 0, TL_index = 0;
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
    // Guess the parameters, it should be close to the true value, else it can fail
    vector<double> params = {E,Tm,double(TL_index),double(TR_index)};
    return params;
}




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
    double FOM = 0.0;
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
                if( j>params[2] && j < params[3]) {
                    error[j] = Func(input[0], t_params);
                }else{
                    error[j] = 0.0;
                }
                double integral = accumulate(error.begin(), error.end(), 0);
                if(j>params[2] && j < params[3]) {
                    FOM += abs(outputs[j][0] - error[j])/integral;
                }
            }
            //Calculate the hessian mat
            vector<vector<double>> Jf(Jf_T[0].size(), vector<double>(Jf_T.size(),0.0));
            transpose(Jf_T,Jf,int(Jf_T.size()), int(Jf_T[0].size()));
            H = multiply(Jf_T,Jf);
            //e = dotProduct(error, error);
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
        //bool bad = false;
        for(int j=params[2]; j < params[3]; j++) {
            input[0][0] = inputs[j][0];
            temp_error[j] = Func(input[0], temp_params);
        }
        double integral = accumulate(temp_error.begin(), temp_error.end(), 0);
        double FOM_temp = 0.0;
        for(int j=params[2]; j < params[3]; j++) {
            input[0][0] = inputs[j][0];
            FOM_temp += abs(outputs[j][0] - temp_error[j])/integral;
        }
        //double temp_e = bad? e:dotProduct(temp_error, temp_error);
        
        if(FOM_temp <  FOM){
            cout<< FOM<<" "<< FOM_temp<<endl;
            lambda /= 10;
            params[0] = E_est;
            params[1] = Tm_est;
            //e = temp_e;
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




void First_Order_Kinetics::glow_curve(){
    curve = count_data;
    int c = 0;
    double curve_sum = std::accumulate(curve.begin(), curve.end(), 0);
    double curve_sum_start=curve_sum;
    vector<double> sum(curve.size(), 0.0);
    vector<vector<double>> updated_params;
    vector<vector<double>> temp_curves;
    bool first = true;
    vector<double> inputs(int(temp_data.size()), 1.0);
    for(int i=0; i < int(temp_data.size()); i++) {
        inputs[i] = double(temp_data[i]);
    }
    while (curve_sum > curve_sum_start*.05){
        if(first){
            updated_params.push_back(initial_guess(curve));
        }
        
        vector<double> params = {updated_params[c][0],updated_params[c][1],double(updated_params[c][2]),double(updated_params[c][3])};
        cout<<"----------------Curve: "<<c+1<<"----------------"<<endl;
        LevenbergMarquardt(inputs,curve, params);
        updated_params[c] = params;
        temp_curves.push_back(vector<double>(temp_data.size(),0));
        for(int i = 0; i < int(temp_data.size()); ++i ){
            vector<double> input(1,0.0);
            input[0]= double(temp_data[i]);
            temp_curves[c][i] = Func(input,params);
        }
        transform(temp_curves[c].begin(),temp_curves[c].end(),sum.begin(),sum.begin(),plus<double>());
        transform(curve.begin(),curve.end(),temp_curves[c].begin(),curve.begin(),minus<double>());
        int TM_index =  (updated_params[c][3] + updated_params[c][2])/2;
        for(int i = 0; i < int(temp_data.size());++i){
            if(curve[i]<0.0) curve[i] = 0.0;
            if(i > TM_index && c == 0) curve[i] = 0.0;
            if(i > TM_index && i < params[3]) curve[i] = 0.0;
        }
        ++c;
        curve_sum = accumulate(curve.begin(), curve.end(), 0);
    }
    glow_curves = temp_curves;
    double FOM = 0;
    double integral = accumulate(sum.begin(), sum.end(), 0);
    for(int i = 0; i < int(curve.size()); ++i){
        FOM += abs(count_data[i] - sum[i])/integral;
    }
    cout<<FOM<<endl;
    glow_curves.push_back(sum);
};


vector<double> First_Order_Kinetics::initial_guess(vector<double> &curve){
    vector<double>::iterator TR = curve.end(),TL = curve.begin(),TM = curve.end();
    int TM_index = 0, TR_index = 0, TL_index = 0;
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
    // Guess the parameters, it should be close to the true value, else it can fail
    vector<double> params = {E,Tm,double(TL_index),double(TR_index)};
    return params;
}

for(int i =0; i < int(glow_curves.size());++i){
    double delta = .01;
    double temp_FOM = FOM;
    int iter = 0;
    transform(curve.begin(),curve.end(),glow_curves[i].begin(),curve.begin(),plus<double>());
    while(iter < 10){
        transform(sum.begin(),sum.end(),glow_curves[i].begin(),sum.begin(),minus<double>());
        if(updated_params[i][2]-delta > 0)updated_params[i][2]-= delta;
            if(updated_params[i][3]+delta < temp_data.size())updated_params[i][3]+= delta;
                LevenbergMarquardt(inputs,curve, updated_params[i]);
                vector<double> temp_curves(temp_data.size(),0);
                for(int j = 0; j < int(temp_data.size()); ++j ){
                    vector<double> input(1,0.0);
                    input[0] = double(temp_data[j]);
                    double output = Func(input,updated_params[i]);
                    temp_curves[i] = output;
                }
        transform(temp_curves.begin(),temp_curves.end(),sum.begin(),sum.begin(),plus<double>());
        temp_FOM = 0.0;
        double integral = accumulate(sum.begin(), sum.end(), 0);
        for(int i = 0; i < int(curve.size()); ++i){
            temp_FOM += abs(count_data[i] - sum[i])/integral;
        }
        if(temp_FOM < FOM){
            delta *= 10;
            glow_curves[i] = temp_curves;
        }
        else delta /= 10;
            ++iter;
    }
    FOM = temp_FOM;
    cout<<FOM<<endl;
    }



    if(first){
        bool main = c==0? true:false;
        updated_params.push_back(initial_guess(curve,main));
        vector<double> params = {updated_params[c][0],updated_params[c][1],double(updated_params[c][2]),double(updated_params[c][3])};
        cout<<"----------------Curve: "<<c+1<<"----------------"<<endl;
        LevenbergMarquardt(inputs,curve, params,false);
        updated_params[c] = params;
    }else if(c != 0){
        if(c == glow_curves.size())break;
        int TM_index = (updated_params[c][2] + updated_params[c][3])/2;
        updated_params[c][0] = activation_energy(updated_params[c][2],TM_index,updated_params[c][3]);
        vector<double> temp_curve(curve.size(), 0.0);
        
        for(int x = updated_params[c][2]; x <= updated_params[c][3]; ++x){
            temp_curve[x] = count_data[x];
        }
        int range = updated_params[c][3] -updated_params[c][2];
        TM = max_element(temp_curve.begin(),temp_curve.end());
        /*if(*TM == temp_curve[updated_params[c][3]]){
         range = updated_params[c][3];
         for(int x = updated_params[c][3]; x > updated_params[c][2]; --x){
         temp_curve[range] = count_data[x];
         range++;
         }
         range = range - updated_params[c][2];
         updated_params[c][3] = updated_params[c][2] + range;
         TM_index = (updated_params[c][3] + updated_params[c][2])/2;
         updated_params[c][0] = activation_energy(updated_params[c][2],TM_index,updated_params[c][3]);
         updated_params[c][1] = double(temp_data[TM_index]);
         }*/
        //subtractions.push_back(temp_curve);
        vector<double> params = {updated_params[c][0],updated_params[c][1],double(updated_params[c][2]),double(updated_params[c][3])};
        cout<<"----------------Curve: "<<c+1<<"----------------"<<endl;
        temp_curve = vector<double>(temp_curve.begin()+updated_params[c][2], temp_curve.begin()+updated_params[c][3]);
        inputs = vector<double>(range, 1.0);
        for(int i=0; i < range; i++) {
            inputs[i] = double(temp_data[updated_params[c][2]+i]);
        }
        LevenbergMarquardt(inputs,temp_curve, params,first);
        updated_params[c] = params;
        }
    
    
    /*-------------------------Levenburg Marquardt second iteration METHOD----------------------------------------//
     vector<double> First_Order_Kinetics::LevenbergMarquardt2(const vector<double> &inputs, const vector<double> &outputs, vector<vector<double>> &params){
     int input_size = int(inputs.size());
     int num_params = int(params.size())*2;
     vector<double> temp_params;
     for(int i = 0; i < params.size();++i){
     for(int j = 0; j < 2; ++j){
     temp_params.push_back(params[i][j]);
     }
     }
     vector<vector<double>> Jf_T(num_params, vector<double>(input_size,0.0)); // Jacobian of Func()
     vector<double> error(input_size,0.0);
     vector<vector<double>> input(1, vector<double>(1,1)); // single row input
     vector<vector<double>> H;
     double lambda = 0.01;
     double updateJ = 1;
     double e = 0.0;
     for(int i = 0; i < 100; ++i){
     if(updateJ == 1){
     //Evaluate the jacobian matrix at the current paramater.
     vector<double>t_parms(2, 0.0);
     for(int j=0; j < input_size; j++) {
     input[0][0] = inputs[j];
     double output = 0.0;
     for(int k = 0; k < num_params;++k){
     t_parms[0] = temp_params[k];
     ++k;
     t_parms[1] = temp_params[k];
     for(int x = 0; x < 2; x++){
     Jf_T[k - x][j] = Deriv(input[0],t_parms , x);
     }
     output += Func(input[0],t_parms);
     }
     error[j] = outputs[j] - output;
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
     //Update parameters
     double deter = 0;
     
     //finding determinant
     deter = determinant(H_lm, int(H_lm.size()));
     
     cout<<"\n\ndeterminant: "<<deter;
     
     cout<<"\n\nInverse of matrix is: \n";
     for(int t = 0; t < 8; t++){
     for(int z = 0; z < 8; z++){
     cout<<((H_lm[(z+1)%8][(t+1)%8] * H_lm[(z+2)%8][(t+2)%8]) - (H_lm[(z+1)%8][(t+2)%8] * H_lm[(z+2)%8][(t+1)%8]))/ deter<<"\t";
     }
     }
     //invert(H_lm, true);
     vector<double> Jf_error = vec_matrix_multi(Jf_T,error);
     vector<double> delta = vec_matrix_multi(H_lm,Jf_error);
     vector<double> t_params = temp_params;
     for(int x = 0; x < delta.size(); ++x){
     t_params[x] += delta[x];
     }
     
     //Evaluate the total distance error at the updated paramaters.
     vector<double> temp_error(input_size,0.0);
     bool bad = false;
     vector<double> t_param(2,0.0);
     for(int j=0; j < input_size; j++) {
     input[0][0] = inputs[j];
     double output = 0.0;
     for(int k = 0; k < num_params;++k){
     t_param[0] = t_params[k];
     ++k;
     t_param[1] = t_params[k];
     output += Func(input[0],t_param);
     }
     double err = outputs[j] - output;
     temp_error[j] = err;
     }
     double temp_e = bad? e:dotProduct(temp_error, temp_error);
     vector<double> temp_outputs;
     if(temp_e < e){
     cout<<"HERE"<<endl;
     lambda /= 10;
     temp_params = t_params;
     e = temp_e;
     updateJ = 1;
     }else{
     updateJ = 0;
     lambda *= 10;
     }
     }
     return temp_params;
     }
    
     /*for(int k = 0; k < num_params;++k){
     vector<double> ttemp_output(input_size,0.0);
     for(int j=0; j < input_size; j++) {
     t_param[param_num] = t_params[k];
     t_param[other_param] = params[k][other_param];
     ttemp_output[j] = Func(temp_data[j],t_param);
     }
     auto begin = ttemp_output.begin();
     auto end = ttemp_output.end();
     transform(begin,end,temp_output.begin(),ttemp_output.begin(),plus<double>());
     auto TM = max_element(ttemp_output.begin(), ttemp_output.end());
     int TM_index = int(TM - ttemp_output.begin());
     maximums.push_back(TM_index);
     }
     //glow_curves.push_back(temp_output);
     for(int max = 0; max < int(maximums.size()); max++){
     if(temp_output[maximums[max]] > curve[maximums[max]]*(1.05)){
     cout<<"too high"<<endl;
     //glow_curves.push_back(temp_output);
     cout<< "Old max: "<<temp_params[max]<<endl;
     temp_params[max] = temp_params[max]*0.9;
     cout<< "New max: "<<temp_params[max]<<endl;
     }
     }*/
    //double err = outputs[j] - output;
    //temp_error[j] = err;
