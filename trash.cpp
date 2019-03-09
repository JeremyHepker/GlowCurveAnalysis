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
