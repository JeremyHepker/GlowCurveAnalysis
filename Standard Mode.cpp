//
//  Standard Mode.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "Standard Mode.hpp"
#include "CSV_iterator.cpp"
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <sstream>
using namespace std;

/*const char* usage = R"foo(
This tool is designed to automatically determine the TLD material glow curve which is being analyized, deconvolves, and outputs the glow peaks. TLD glow curve data provided by the user in a properly fromatted .xls file.

Results: When in non-verbose mode, no output will be displayed while running analysis. A message will dispaly successful analysis along with the goodness of fit value for the analysis. Output file containing glow peaks will be output to a unless specified default named zip file. If you require a full output of all analysis progress, please use the verbose setting.

Output:  [\x1B[35;40m+\x1B[0m] Added Headers, [\x1B[35;40m-\x1B[0m] Removed Headers, [\x1B[35;40m!\x1B[0m] Altered Headers, [ ] No Change
                                                                                           
   Usage .:
   -u / --url Complete URL
   -f / --file <Path to User Agent file> / If no file is provided, -d options must be present
   -s / --single provide single user-agent string (may need to be contained within quotes)
   -o / --output <Path to output file> CSV formated output (FILE WILL BE OVERWRITTEN[\x1B[31;40m!\x1B[0m])
   -v / --verbose results (Displays full headers for each check) >> Recommended
   --debug See debug messages (This isnt the switch youre looking for)\n
   Example .:)foo";*/

//This function displays the software logo, prints the usage, and pompts the user
//to re-run the program with proper arguments.

standard::standard():time(false){
    greeting();
};
                                                                                          
                                                                                          
void standard::greeting(){
    //cout<<usage<<endl;
};

//This function reads in the .xls file and parses the data into vector of coordinate pairs.
void standard::read(string filename){
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
            line.find_first_of("Time");
            time = true;
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
void standard::write(string out_filename){
    ofstream file;
    file.open (out_filename);
    if(!file.is_open()){
        cout<<"error reading file"<<endl;
        stop_program();
    }
    for(int i = 0; i<raw_temp_data.size();++i){
        file << raw_temp_data[i]<<","<<raw_count_data[i]<<",\n";
    }
    file.close();
};

void standard::heat_rate_calc(){
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

void standard::stop_program(){
    exit(1);
};
