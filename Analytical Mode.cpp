//
//  Analytical Mode.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/21/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "Analytical Mode.hpp"
#include "CSV_iterator.cpp"
#include <stdio.h>
#include <iostream>
#include <getopt.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <sstream>
using namespace std;
analytic::analytic(string &filename, string &out_filename, string &model_type, string &verbose_mode,string &material){
    if(verbose_mode == "true") verbose = true;
    read(filename);
    if(model_type == "FOK"){
        //deconvolve_first_order_kinetics();
    }else{
        //deconvolve_OTOR();
    }
    write(out_filename);
};


void analytic::greeting(){
    //cout<<usage<<endl;
};

//This function reads in the .xls file and parses the data into vector of coordinate pairs.
void analytic::read(string &filename){
    //Open and test the user input file.
    ifstream file(filename);
    if(!file.is_open()){
        cout<<"error reading file"<<endl;
        stop_program();
    }
    //Parse the input file header data.
    for(int i = 0; i < 6; ++i){
        string line;
        getline(file, line);
        if(i == 2){
            line.erase(0, 18);
            auto comma = line.find_first_of(",");
            line.erase(comma,3);
            filename_str = line;
        }else if(i == 3){
            line.erase(0,17);
            auto comma = line.find_first_of(",");
            line.erase(comma,3);
            date_time = line;
        }else if(i == 5){
            if(line.find("Time") != string::npos){
                time = true;
            }
        }
    }
    //read in the csv file and parse date into two raw input vectors.
    auto i = csv_iterator<int>( file );
    ++i;
    raw_temp_data.push_back(*i);
    ++i;
    raw_count_data.push_back(*i);
    ++i;
    if(time){
        raw_time_data.push_back(*i);
        ++i;
    }
    counts_total = *i;
    while(file){
        ++i;
        raw_temp_data.push_back(*i);
        ++i;
        raw_count_data.push_back(*i);
        ++i;
        if(time){
            raw_time_data.push_back(*i);
            ++i;
        }
    }
    file.close();
};

//This is a function to write the output to a new CSV file.
void analytic::write(string &out_filename){
    ofstream file;
    file.open (out_filename+".csv");
    if(!file.is_open()){
        cout<<"error reading file"<<endl;
        stop_program();
    }
    for(int i = 0; i<raw_temp_data.size();++i){
        file << raw_temp_data[i]<<","<<raw_count_data[i]<<",\n";
    }
    file.close();
};

void analytic::heat_rate_calc(){
    auto size = raw_temp_data.size();
    size--;
    float measure1 = (raw_temp_data[size]-raw_temp_data[0])/((raw_time_data[size]-raw_time_data[0]));
    float measure2 = (raw_temp_data[floor(size/2)]-raw_temp_data[0])/((raw_time_data[floor(size/2)]-raw_time_data[0]));
    float measure3 = (raw_temp_data[size]-raw_temp_data[floor(size/2)])/((raw_time_data[size]-raw_time_data[floor(size/2)]));
    heat_rate = (measure1 + measure2 + measure3)/3;
}

//Deconvolve the input data using first order kinetics.
void deconvolve_first_order_kinetics();

//deconvolve the input data using the OTOR model.
void deconvolve_OTOR();

void analytic::stop_program(){
    exit(1);
};
