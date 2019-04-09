//
//  File_Manager.hpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#ifndef File_Manager_hpp
#define File_Manager_hpp

#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <stdio.h> 
#include <iostream>
#include <cstdlib>

class File_Manager{
private:
    std::vector<double> raw_temp_data;
    std::vector<double> raw_count_data;
    int total_counts = 0;
    double time = 0.1;
    std::string filename,header;
public:
    
    File_Manager(std::string filename);
    File_Manager();
    
    //This function reads in the .csv file and parses the data into vector of coordinate pairs.
    std::pair<std::vector<double>,std::vector<double>> read();
    
    //This is a function to write the output to a new CSV file.
    void write(std::vector<std::vector<double>> glow_curves, std::string output_name);

    double temp_rate();
    
    ~File_Manager();
    
};


#endif /* File_Manager_hpp */
