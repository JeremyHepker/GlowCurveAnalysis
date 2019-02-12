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
