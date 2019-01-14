//
//  main.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include <iostream>
#include "Standard Mode.hpp"
using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    standard stan;
    stan.read("/Users/jeremyhepker/Documents/NERS499/GlowCurveAnalsys/GlowCurveAnalsys/test_1.csv");
    stan.write("outputfile.csv");
    cout<<"we did it";
    return 0;
    
}

