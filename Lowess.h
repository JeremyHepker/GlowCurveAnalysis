//
//  Lowess.h
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 3/7/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#ifndef CPPLOWESS_LOWESS_H
#define CPPLOWESS_LOWESS_H

//Polynomial Fit
#include <iostream>
#include <vector>

using namespace std;

class boxFIR
{
    int numCoeffs; //MUST be > 0
    vector<double> b; //Filter coefficients
    vector<double> m; //Filter memories
    
public:
    boxFIR(int _numCoeffs) :
    numCoeffs(_numCoeffs)
    {
        if (numCoeffs<1)
            numCoeffs = 1; //Must be > 0 or bad stuff happens
        
        double val = 1./numCoeffs;
        for (int ii=0; ii<numCoeffs; ++ii) {
            b.push_back(val);
            m.push_back(0.);
        }
    }
    
    void filter(vector<double> &a)
    {
        double output;
        
        for (int nn=0; nn<a.size(); ++nn)
        {
            //Apply smoothing filter to signal
            output = 0;
            m[0] = a[nn];
            for (int ii=0; ii<numCoeffs; ++ii) {
                output+=b[ii]*m[ii];
            }
            
            //Reshuffle memories
            for (int ii = numCoeffs-1; ii!=0; --ii) {
                m[ii] = m[ii-1];
            }
            
            a[nn] = output;
        }
    }
    
    
};

#endif //CPPLOWESS_LOWESS_H
