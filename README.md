# Glow Curve Analysis
Information about the signal read out from a thermoluminescent (TL) material is often obtained through the fitting of glow curves to the emitted spectra. This fitting is done using a general nth order kinetics formula.

Most GCA software is dependent on human input for the initial fitting parameters, thus consistent convergence becomes more of an art than a science. This work consists of developing a GCA software to automate the process of identifying and fitting the individual glow peaks of the TL spectrum to eliminate human participation.

Most GCA software is written using Matlab or other similar tools or libraries for ease of the programmer. This project uses only C++ and its standard library for high portability and accessibility for the user of the code. The intention of this project is to create an automated, easy-to-use, highly portable, open-source GCA software.

## Files

- CSV_iterator.cpp
	- Optimized .csv file handler
- Levenberg–Marquardt.hpp
	- Iterative Levenberg Marquardt algorithm 
- dataSmoothing.cpp
	- Average smoothing Algorithm & Data spike smoothing
- OTORModel.cpp
	- Experimental OTOR (not complete, more research and design required ) 
- dataSmoothing.hpp
	-  Average smoothing Algorithm & Data spike smoothing
- OTORModel.hpp
	- Experimental OTOR (not complete, more research and design required ) 
- firstOrderKineticsModel.cpp
	- Fitting model utilized to fit the glow peaks
- batch_handler.hpp
	- Batch file handling which employs the Boost Filesystem Library
- firstOrderKineticsModel.hpp
	- Fitting model utilized to fit the glow peaks
- main.cpp
	- Program Driver
- fileManager.cpp
- matrix_arithmetic.cpp
	- Optimized Matrix algebra, combining C and C++
- fileManager.hpp
- smartPeakDetect.cpp
- Levenberg–Marquardt.cpp
- smartPeakDetect.hpp

## Boost Filesystem 
- Download the boost libraries that relates to your system 
	- unix : MacOS & Linux 
	- windows : Well Windows..
	-  (automated process is under construction) 
	-	[https://www.boost.org/users/history/version_1_62_0.html](https://www.boost.org/users/history/version_1_62_0.html)
	-	Almost any errors that occur at compile time are a result of improper boost download. Follow the instructions. 

## Run the Program 
- Clone the GitHub Repo onto your computer
-- Open your terminal
-- Navegate to your desired directory
-- Clone the repository
```
git clone https://github.com/JeremyHepker/GlowCurveAnalsys.git
```
- Once the repository is successfully cloned navigate into the folder
```
cd GlowCurveAnalsys
``` 
- To build to project run...
--Any compile issues will most likely be improper installation of boost. It can be tricky. Feel free to reach out. The only Boost used is Boost FileSystem, google installing this specifically if problems occur.
```
make all
```
- To run the project now just run the following command
```
./GCA
```
- The program will prompt you for the path to the directory of spectra files you wish to process. Put in the path and it will notify you if it has found the directory or not, how many .CSV files it found in the directory, and ask if the information looks correct before deconvolution. 
- Enter y or n and hit enter to begin. 
- That is all folks. 
- To delete the build and rebuild, if changes to the files are made run...
```
make clean
make all
```
- This will clean the old build and its components out, rebuild with the new code. 
## Adding New Files
- Add your files to the src directory ONLY
- Edit the makefile, add any .cpp files (src/< filename >) to BOTH lines 14 and 16, with the other two list of files. 

## Test Data
- There is test data provided in the testSpectra directory, it is a bit of a mess but there is alot. To run the data in this folder AFTER BUILD..  
```
./GCA testSpectra/< desired test directory > 
```
- Example
```
./GCA testSpectra/TLD_100
```
## Getting Better Results
- To yeild better results quick 20< line python programs can be used to scan over directories of csv files and smooth data. 

## Plotting Results
- I used a wide variety of short 20< line python programs to plot entire directories of output data. Don't get wrapped up in c++ potting. 

## FAQ
- No one has used this but me there are not FAQs. 
- Reach out at jhepker@umich.edu with questions and concerns 
- I run this on a Mac but have run it in Linux and Windows successfully. The Boost installation is a bit more tricky but manageable.
## Running without the Makefile
- If for some reason you wish to build the program without the makefile here, run this from inside the top level GlowCurveAnalysis directory, NOT THE src folder.
```
g++ -std=c++0x -g -Wall -lboost_system -lboost_filesystem -o GCA src/smartPeakDetect.cpp src/OTORModel.cpp src/matrix_arithmetic.cpp src/Levenberg–Marquardt.cpp src/DataSmoothing.cpp src/CSV_iterator.cpp src/main.cpp src/FOKModel.cpp src/File_Manager.cpp  
```
