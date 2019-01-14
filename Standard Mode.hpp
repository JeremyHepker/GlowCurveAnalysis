//
//  Standard Mode.hpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#ifndef Standard_Mode_hpp
#define Standard_Mode_hpp

#include <stdio.h>
#include <vector>
#include <string>

using namespace std;
class standard{
private:
    string filename, out_filename, date_time, filename_str;
    vector<int> raw_temp_data, raw_count_data, raw_time_data;
    int counts_total;
    float heat_rate; 
    bool time;
public:
    standard();
    //This function displays the software logo, prints the usage, and pompts the user
    //to re-run the program with proper arguments.
    void greeting();
    
    //This function reads in the .xls file and parses the data into vector of coordinate pairs.
    void read(string filename);
    
    //This is a function to write the output to a new CSV file.
    void write(string out_filename);
    
    //This function calculates the rate at which the TLD was heating. 
    void heat_rate_calc();
    
    //Deconvolve the input data using first order kinetics.
    void deconvolve_first_order_kinetics();
    
    //deconvolve the input data using the OTOR model. 
    void deconvolve_OTOR();
    
    //Global abort function for error managment. 
    void stop_program();
};
#endif /* Standard_Mode_hpp */
