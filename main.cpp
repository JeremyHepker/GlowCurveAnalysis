//
//  main.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include <iostream>
#include <getopt.h>
#include "Usage.h"
#include <string>
#include "Standard Mode.hpp"
using namespace std;
vector<string> getArguments(int argc,char * argv[]) {
    if(argc < 2){
        cout<<usage<<endl; 
    }
    string input;
    string mode,filename,out_filename = "GCA_output",verbose = "false", model = "FOK";
    // These are used with getopt_long()
    opterr = true; // Give us help with errors
    int choice;
    int option_index = 0;
    option long_options[] = {
        { "mode",    required_argument,      nullptr, 'm' },
        { "filename", required_argument,      nullptr, 'f' },
        { "model",   required_argument,      nullptr, 'l' },
        { "material", required_argument,      nullptr, 's' },
        { "output",   required_argument,      nullptr, 'o' },
        { "verbose", no_argument,            nullptr, 'v' },
        { "help",   no_argument,            nullptr, 'h' },
        { nullptr, 0,                     nullptr, '\0'}
    };
    
    // TODO: Fill in the double quotes, to match the mode and help options.
    while ((choice = getopt_long(argc, argv, "m:f:l:s:o:vh",long_options, &option_index)) != -1) {
        switch (choice) {
            case 'm':
                if(strcmp(optarg,"standard")==0){
                    mode = "standard";
                }else if(strcmp(optarg,"analytical")==0){
                    mode = "analytical";
                }else if(strcmp(optarg,"smart")==0){
                    mode = "smart";
                }else{
                    cerr<<"The mode provided was invalide, please select one of the folling: standard, analytical, or smart."<<endl<< "Example : -m standard"<<endl<<"Example : -m analytical"<<endl<<"Example : -m smart"<<endl;
                    exit(0);
                }
                break;
            case 'f':
                if(string(optarg).find(".csv") == string::npos){
                    cerr<<"Invalide input file, please ensure file is of the extension .csv"<<"Example : -f <filename>.csv"<<endl<<"If error persists include full file path"<<endl<<"Example : -f ~/<Filename>.csv"<<endl;
                    exit(0);
                }
                filename = optarg;
                break;
            case 'l':
                if(strcmp(optarg,"FOK")!=0 && strcmp(optarg,"OTOR")!=0){
                    cerr<<"Invalide model input, please select either FOK (first order kinetics) or OTOR (OTOR level)"<<endl<<"Example : -l FOK"<<endl<<"Example : -l OTOR"<<endl;
                    exit(0);
                }
                model = optarg;
                break;
            case 'o':
                if(string(optarg).find(".csv") != string::npos){
                    string output = optarg;
                    output.erase(string(optarg).find(".csv"),4);
                    out_filename = output;
                }else{
                    out_filename = optarg;
                }
                break;
            case 'v':
                verbose = "true";
                break;
            case 'h':
                cout<<usage<<endl;
                break;
            default:
                cerr << "Invalid user input, see -h or --help for usage." << endl;
                exit(1);
        } // switch
    } // while
    vector<string> args = {mode,filename,out_filename,model,verbose};
    return args;
} // getScheme()

string analytic_input(){
    string material;
    cout<<"You have selected to use the analytical mode which deconvolves and outputs the glow peaks based on user provided information. If this is not desired enter \"exit\". Otherwise either the material being used, or the number of peaks and their approximate threshold temperatures are required"<<endl<<"Enter the number of peaks : ";
    getline(cin,material);
    return material;
}

int main(int argc, char * argv[]) {
    //standard stan;
    vector<string> args = getArguments(argc, argv);
    if(args[0]=="standard"){
        standard stand(args[1],args[2],args[3],args[4]);
    }else if(args[0]=="analytical"){
        cout<<"analytical"<<endl;
        string material = analytic_input();
        cout<<material;
    }else{
        cout<<"smart"<<endl;
    }
    cout<<"we did it";
    return 0;
    
}

