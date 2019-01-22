const char* usage = R"(
This tool is designed to automatically determine the TLD material glow curve which is being analyized, deconvolves, and outputs the glow peaks. TLD glow curve data provided by the user in a properly fromatted .xls file.

Results: When in non-verbose mode, no output will be displayed while running analysis. A message will dispaly successful analysis along with the goodness of fit value for the analysis. Output file containing glow peaks will be output to a unless specified default named zip file. If you require a full output of all analysis progress, please use the verbose setting.

Output:  [\x1B[35;40m+\x1B[0m] Added Headers, [\x1B[35;40m-\x1B[0m] Removed Headers, [\x1B[35;40m!\x1B[0m] Altered Headers, [ ] No Change

Usage .:
-u / --url Complete URL
-f / --file <Path to User Agent file> / If no file is provided, -d options must be present
-s / --single provide single user-agent string (may need to be contained within quotes)
-o / --output <Path to output file> CSV formated output (FILE WILL BE OVERWRITTEN[\x1B[31;40m!\x1B[0m])
-v / --verbose results (Displays full headers for each check) >> Recommended
--debug See debug messages (This isnt the switch youre looking for)\n
Example .:)";
