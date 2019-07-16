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
#include <fstream>
#include <iomanip>
#include "src/File_Manager.hpp"
#include "src/batch_handler.hpp"
#include "src/smartPeakDetect.hpp"
#include "src/DataSmoothing.hpp"
#include "src/Levenberg–Marquardt.hpp"
//#include "First_Order_kinetics.hpp"


using namespace std;

int main(int argc, char * argv[]) {
    string dir = argv[1];
    if(dir.back() == '/') dir.pop_back();
    vector<string> filenames,files = batch_handler(dir);
    vector<vector<double>> stats(files.size(), vector<double>(0,0.0)), peakParams;
    int count = 0;
    auto i = files.begin();
    string output_dir = *i++;
    for(; i != files.end(); ++i){
        double integral= 0;
        if(i->find("temp.csv") !=string::npos){
            files.erase(i);
            continue;
        }
        cout<<"----------------------------"<<endl<<"Processing: ";
        string filename = i->substr((i->find_last_of("/\\"))+1);
        cout<<filename<<" ("<<count+1<<" of "<<files.size()-1<<")"<<endl<<"Reading in File  .";
        cout.flush();
        File_Manager fileManager = *new File_Manager(*i);
        cout<<".";
        cout.flush();
        pair<vector<double>, vector<double>> data = fileManager.read();
        remove( dir + "/temp.csv" );
        findPeaks(data.first,data.second, peakParams);
        stats[count].push_back(fileManager.barcode());
        //for(auto i = data.second.begin(); i != data.second.end();i++) integral += *i;
        remove( dir + "/temp.csv" );
        cout<<"."<<endl;
        cout.flush();
        cout<<"Deconvoluting Glow Peak  .";
        cout.flush();
        First_Order_Kinetics FOK_Model = *new First_Order_Kinetics(data,peakParams);
        stats[count].push_back(FOK_Model.glow_curve());
        stats[count].push_back(integral);
        stats[count].push_back(fileManager.temp_rate(filename));
        vector<double> peak_integral = FOK_Model.return_curve_areas();
        for(auto i = peak_integral.begin(); i != peak_integral.end();++i){
            stats[count].push_back(*i);
        }
        vector<vector<double>> returnedPeaks = FOK_Model.return_glow_curve();
        filenames.push_back(filename);
        filename = output_dir +"/"+ filename;
        fileManager.write(returnedPeaks, filename, count);
        fileManager.statistics(stats,filenames,output_dir);
        cout<<"----------------------------"<<endl;
        count++;
    }
    return 0;
}
