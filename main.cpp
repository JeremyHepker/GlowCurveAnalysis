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
#include "Usage.h"
#include "File_Manager.hpp"
#include "Lowess.h"
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


int main(int argc, char * argv[]) {
    vector<string> args = getArguments(argc, argv);
    int material = -1;
    vector<pair<int, int>> thresholds; 
    if(args[0]=="standard"){
        //standard stand(args[1],args[2],args[3],args[4]);
    }else if(args[0]=="analytical"){
        string choice;
        cout<<"You have selected Analytical mode, if this was incorrect enter 'quit'"<<endl;
        cout<<"Enter either:"<<endl<<"[1] To select the material which was analysized"<<endl;
        cout<<"[2] To enter number of peaks and threashhold temperatures"<<endl;
        cout<<"Enter option here: ";
        cin>>choice;
        material = material_intake();
        cout<< material;
        if(choice == "quit"){
            cout<< usage;
            exit(1);
        }else if(choice == "2"){
            cout<<"2";
        }else{
            //material = material_intake();
        }
        
    }else{
        for(int i = 0; i < 1; ++i){
            File_Manager manager = *new File_Manager(args[1]);
            pair<vector<int>, vector<double>> data = manager.read();
            boxFIR box(10);
            //box.filter(data.second);
            //vector<double> smoothed = data.second;
            First_Order_Kinetics FOK_Model = *new First_Order_Kinetics(data);
            FOK_Model.glow_curve();
            vector<vector<double>> peaks = FOK_Model.return_glow_curve();
            //vector<vector<double>> peaks;
            //peaks.push_back(smoothed);
            manager.write(peaks, args[2]);
        }
    }
    return 0;
}
