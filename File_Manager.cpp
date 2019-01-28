//
//  File_Manager.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "File_Manager.hpp"
#include "CSV_iterator.cpp"
using namespace std;
File_Manager::File_Manager(std::string given_filename):filename(given_filename){};

//This function reads in the .csv file and parses the data into vector of coordinate pairs.
pair<vector<int>,vector<double>> File_Manager::read(){
    //Open and test the user input file.
    ifstream file(filename);
    if(!file.is_open()){
        cerr<<"Error opening file: "<<filename<<endl;
        exit(1);
    }
    //Parse the input file header data.
    string line;
    getline(file, line);
    while(line.find("Counts") == string::npos){
        header += (line += "\n");
        getline(file, line);
    }
    line += "\n";
    header += line;
    //read in the csv file and parse date into two raw input vectors.
    auto i = csv_iterator<int>( file );
    ++i;
    raw_temp_data.push_back(*i);
    ++i;
    raw_count_data.push_back(*i);
    ++i;
    total_counts = *i;
    while(file){
        ++i;
        raw_temp_data.push_back(*i);
        ++i;
        raw_count_data.push_back(*i);
        ++i;
    }
    file.close();
    return make_pair(raw_temp_data, raw_count_data);
};

//This is a function to write the output to a new CSV file.
void File_Manager::write(vector<vector<double>> glow_curves, string output_name){
    ofstream file;
    file.open("/Users/jeremyhepker/Documents/NERS499/GlowCurveAnalsys/GlowCurveAnalsys/"+output_name+".csv");
    if(!file.is_open()){
        cerr<<"Could not open output file : "<<output_name<<endl;
        exit(1);
    }
    for(int i = 0; i<raw_temp_data.size();++i){
        file << raw_temp_data[i]<<","<<raw_count_data[i];
        for(int j = 0; j<glow_curves.size();++j){
            file<<","<<glow_curves[j][i];
        }
        file<<",\n";
    }
    file.close();
};
double File_Manager::temp_rate(){
    return time;
}
