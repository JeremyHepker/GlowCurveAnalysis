//
//  File_Manager.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "File_Manager.hpp"
#include "CSV_iterator.cpp"
#include <iomanip>

using namespace std;
File_Manager::File_Manager(std::string given_filename):filename(given_filename){};

//This function reads in the .csv file and parses the data into vector of coordinate pairs.
pair<vector<double>,vector<double>> File_Manager::read(){
    //Open and test the user input file.
    ifstream file(filename);
    if(!file.is_open()){
        cerr<<"Error opening file: "<<filename<<endl;
        exit(1);
    }
    //Parse the input file header data.
    string line;
    getline(file, line,',');
    while(line.find("Counts") == string::npos){
        header += (line += "\n");
        getline(file, line,',');
    }
    header += line;
    //read in the csv file and parse date into two raw input vectors.
    getline(file, line,'\r');
    getline(file, line,'\r');
    auto i = csv_iterator<int>( file );
    ++i;
    while(file){
        raw_temp_data.push_back(*i);
        ++i;
        raw_count_data.push_back(*i);
        ++i;
        if(!file.eof()) ++i;
    }
    file.close();
    //raw_count_data = gaussSmoothen(raw_count_data,10.0, int(raw_count_data.size())-1);
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
    file<<"temp,count_original";
    for(int j = 0; j<glow_curves.size();++j){
        string ster = "count_" + to_string(j);
        file<<","<<ster;
    }
    file<<",\n";
    file.setf(ios_base::fixed);
    file<<setprecision(5);
    for(int i = 0; i<raw_temp_data.size();++i){
        file << raw_temp_data[i]<<",";
        for(int j = 0; j<glow_curves.size();++j){
            file<<","<<double(glow_curves[j][i]);
        }
        file<<",\n";
    }
    file.close();
};
double File_Manager::temp_rate(){
    return time;
}

File_Manager::~File_Manager()
{
    raw_temp_data.clear();
    raw_count_data.clear();
}
double File_Manager::gauss(double sigma, double x) {
    double expVal = -1 * (pow(x, 2) / pow(2 * sigma, 2));
    double divider = sqrt(2 * M_PI * pow(sigma, 2));
    return (1 / divider) * exp(expVal);
}

vector<double> File_Manager::gaussKernel(int samples, double sigma) {
    vector<double> v;
    
    bool doubleCenter = false;
    if (samples % 2 == 0) {
        doubleCenter = true;
        samples--;
    }
    int steps = (samples - 1) / 2;
    double stepSize = (3 * sigma) / steps;
    
    for (int i = steps; i >= 1; i--) {
        v.push_back(gauss(sigma, i * stepSize * -1));
    }
    
    v.push_back(gauss(sigma, 0));
    if (doubleCenter) {
        v.push_back(gauss(sigma, 0));
    }
    
    for (int i = 1; i <= steps; i++) {
        v.push_back(gauss(sigma, i * stepSize));
    }
    
    assert(v.size() == samples);
    
    return v;
}

vector<double> File_Manager::gaussSmoothen(vector<double> values, double sigma, int samples) {
    vector<double> out;
    auto kernel = gaussKernel(samples, sigma);
    int sampleSide = samples / 2;
    //int valueIdx = samples / 2 + 1;
    unsigned long ubound = values.size();
    for (unsigned long i = 0; i < ubound; i++) {
        double sample = 0;
        int sampleCtr = 0;
        long tempj =(i - sampleSide);
        long tempi =(i + sampleSide);
        for (long j = tempj; j <= tempi; j++) {
            if (j > 0 && j < ubound) {
                int sampleWeightIndex = int(sampleSide + (j - i));
                sample += kernel[sampleWeightIndex] * values[j];
                sampleCtr++;
            }
        }
        double smoothed = sample / (double)sampleCtr;
        out.push_back(smoothed);
    }
    return out;
}
