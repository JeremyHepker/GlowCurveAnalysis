//
//  main.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include <iostream>
#include <getopt.h>
#include <string>
#include <locale>
#include "polyfit.hpp"
#include "Usage.h"
#include <fstream>
#include <iomanip>
#include "File_Manager.hpp"
#include "data_smoothing.hpp"
#include "batch_handler.hpp"
#include "First_Order_kinetics.hpp"
using namespace std;
vector<string> getArguments(int argc,char * argv[]) {
    string input,mode,filename,out_filename = "",verbose = "false", model = "";
    opterr = true;
    int choice,option_index = 0;
    option long_options[] = {
        { "mode",    required_argument,      nullptr, 'm' },
        { "filename", required_argument,      nullptr, 'f' },
        { "model",   required_argument,      nullptr, 'l' },
        { "output",   required_argument,      nullptr, 'o' },
        { "verbose", no_argument,            nullptr, 'v' },
        { "help",   no_argument,            nullptr, 'h' },
        { nullptr, 0,                     nullptr, '\0'}
    };
    while ((choice = getopt_long(argc, argv, "m:f:l:o:vh",long_options, &option_index)) != -1) {
        switch (choice) {
            case 'm':
                mode = optarg;
                if(mode == "standard" || mode == "smart" || mode == "analytical"){
                    break;
                }else{
                    cerr<<"Invalid Argument Provided : -m <mode>"<<endl;
                    cerr<<usage;
                    exit(0);
                }
            case 'f':
                if(string(optarg).find(".csv") == string::npos){
                    cerr<<"Invalid Input Provided : -f <filename.csv>"<<endl;
                    cerr<<usage;
                    exit(0);
                }
                filename = optarg;
                break;
            case 'l':
                if(strcmp(optarg,"FOK")!=0 && strcmp(optarg,"OTOR")!=0){
                    cerr<<"Invalid Input Provided : -l <model>"<<endl;
                    cerr<<usage;
                    exit(0);
                }
                model = optarg;
                break;
            case 'o':
                out_filename = optarg;
                break;
            case 'v':
                verbose = "true";
                break;
            case 'h':
                cout<<usage<<endl;
                break;
            default:
                cerr << "Invalid Input Provided"<< endl;
                cerr<<usage<<endl;
                exit(1);
        }
    }
    vector<string> args = {mode,filename,out_filename,model,verbose};
    return args;
} // getScheme()

int material_intake(){
    int material = -1;
    for(size_t i = 0; i <= materials.size()-1; ++i){
        cout<<"["<<i<<"] "<<materials[i]<< endl;
    }
    cout<<"Enter Material Index: ";
    cin>>material;
    while(material < 0 || material > 8){
        cerr<<"Enter a valid material."<<endl;
        material = material_intake();
    }
    return material;
}
vector<pair<int, int>> peak_intake(){
    vector<pair<int, int>> thresholds;
    int peaks =0;
    cout<<"Enter Number of Peaks: ";
    cin>> peaks;
    return thresholds;
}
void heating_FOM(vector<double> &rates, vector<double> &FOMs, string dir){
    std::ofstream myfile;
    dir += "/temp_v_FOM.csv";
    myfile.open(dir);
    if(!myfile.is_open()){
        cerr<<"Could not open output file : Heating_rate_vs_FOM.csv"<<endl;
        exit(1);
    }
    for(int i = 0; i < int(rates.size()); ++i){
        myfile<<rates[i]<<","<<FOMs[i]<<"\n";
    }
    myfile.close();
}

void barcodesOutput(vector<double> &codes, vector<double> &FOMs, string dir){
    std::ofstream myfile;
    dir += "/barcode_v_FOM.csv";
    myfile.open(dir);
    if(!myfile.is_open()){
        cerr<<"Could not open output file : Heating_rate_vs_FOM.csv"<<endl;
        exit(1);
    }
    for(int i = 0; i < int(codes.size()); ++i){
        myfile<<codes[i]<<","<<FOMs[i]<<"\n";
    }
    myfile.close();
}

double count_sum(vector<double> &counts){
    double sum = std::accumulate(counts.begin(), counts.end(), 0);
    return sum;
}
void countOutput(vector<double> &counts, vector<string> &names, string dir){
    std::ofstream myfile;
    dir += "/file_counts.csv";
    myfile.open(dir);
    if(!myfile.is_open()){
        cerr<<"Could not open output file : Heating_rate_vs_FOM.csv"<<endl;
        exit(1);
    }
    for(int i = 0; i < int(counts.size()); ++i){
        myfile<<names[i+1]<<","<<counts[i]<<"\n";
    }
    myfile.close();
}

int main(int argc, char * argv[]) {
    string dir = argv[1];
    if(dir.back() == '/') dir.pop_back();
    vector<string> files = batch_handler(dir);
    vector<double> FOMs, counts, heating_rates,barcodes;
    int count = 0;
    auto i = files.begin();
    string output_dir = *i++;
    for(; i != files.end(); ++i){
        if(i->find("temp.csv") !=string::npos){
            files.erase(i);
            continue;
        }
        cout<<"----------------------------"<<endl;
        string filename = i->substr((i->find_last_of("/\\"))+1);
        cout<<"Processing: "<<filename<<"("<<count<<" of "<<files.size()<<")"<<endl<<"Reading in File  .";
        cout.flush();
        File_Manager manager = *new File_Manager(*i);
        cout<<".";
        cout.flush();
        pair<vector<double>, vector<double>> data = manager.read();
        remove( dir + "/temp.csv" );
        vector<double > idk = polyfit(data.first, data.second, 5);
        //gaussSmoothen(data.second, .01, int(data.second.size()));
        cout<<"."<<endl;
        cout.flush();
        cout<<"Deconvoluting Glow Peak  .";
        cout.flush();
        First_Order_Kinetics FOK_Model = *new First_Order_Kinetics(data);
        FOMs.push_back(FOK_Model.glow_curve());
        if(FOMs.back() < 0 || FOMs.back() > 1.0 || isnan(FOMs.back())){
            FOMs.pop_back();
            count++;
            continue;
        }
        barcodes.push_back(manager.barcode());
        heating_rates.push_back(manager.temp_rate(filename));
        vector<vector<double>> peaks = FOK_Model.return_glow_curve();
        counts.push_back(count_sum(data.second));
        filename = output_dir +"/"+ filename;
        manager.write(peaks, filename, count);
        cout<<"----------------------------"<<endl;
        count++;
    }
    barcodesOutput(barcodes, FOMs,output_dir);
    if(FOMs.size() > 0)heating_FOM(heating_rates, FOMs,output_dir);
    return 0;
}
