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
#include "Standard Mode.hpp"
#include "Analytical Mode.hpp"
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
    
    if(args[0]=="standard"){
        //standard stand(args[1],args[2],args[3],args[4]);
    }else if(args[0]=="analytical"){
        string choice;
        cout<<"You have selected Analytical mode, if this was incorrect enter 'quit'"<<endl;
        cout<<"Enter either:"<<endl<<"[1] To select the material which was analysized"<<endl;
        cout<<"[2] To enter number of peaks and threashhold temperatures"<<endl;
        cout<<"Enter option here: ";
        cin>>choice;
        if(choice == "quit"){
            cout<< usage;
            exit(1);
        }else if(choice == "2"){
            cout<<"2";
        }else{
            int material = material_intake();
        }
        //analytic stand(args[1],args[2],args[3],args[4],args[5]);
    }else{
        cout<<"smart"<<endl;
    }
    return 0;
    
}
