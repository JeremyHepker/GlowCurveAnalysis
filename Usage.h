const char* usage = R"(
This tool is designed to automatically determine the TLD material glow curve which is being analyized, deconvolves, and outputs the glow peaks. TLD glow curve data provided by the user in a properly fromatted .csv file.

Results: When in non-verbose mode, no output will be displayed while running analysis. A message will display successful analysis along with the goodness of fit value for the analysis. Output file containing glow peaks will be output to a unless specified default named zip file. If you require a full output of all analysis progress, please use the verbose setting.

Usage .:
-m / --mode : Mode to be run, standard, analytical, or smart. Unless specified defult Standard.
-f / --filename : <Path to .csv file> >> Required.
-l / --model : to be used for analysis, First Order Kinetics Model(FOK) or OTOR Level(OTOR) model. Unless specified defult First Frder Kinetics.
-s / --material : Enter the material being used >> Required unless -n tag provided.
-o / --output : <Path to output file> CSV formated output, defult glow_curve_report.csv (FILE WILL BE OVERWRITTEN).
-v / --verbose : results (Displays full progress of analysis) >> Recommended.
-h / --help : Display usage.
Example .:
>> GCA -m analytical -f </path/filename> -l FOK -s TLD100 -v
)";
