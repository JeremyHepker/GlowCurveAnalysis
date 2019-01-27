#include <vector>
#include <string>
std::vector<std::string> materials =
    {
    "TLD-100  : Lithium Fluoride Magnesium Titanium(LiF: Mg,Ti)",
    "TLD-200  : Calcium Fluoride Dyprosium (CaF2: Dy)",
    "TLD-300  : Calcium Fluoride Thulium (CaF2:Tm)",
    "TLD-400  : Calcium Fluoride Manganese (CaF2:Mn)",
    "TLD-600  : Lithium Fluoride Magnesium Titanium (Lithium isotope) (LiF:Mg,Ti(Li-6))",
    "TLD-700  : Lithium Fluoride Magnesium Titanium (Lithium isotope) (LiF:Mg,Ti(Li-7))",
    "TLD-100H : High-Sensitivity Lithium Fluoride Magnesium Copper (Lithium isotope) (LiF:Mg,Cu,P)",
    "TLD-700H : High-Sensitivity Lithium Fluoride Magnesium Copper (Lithium isotope) (LiF:Mg,Cu,P(Li-7))",
    "TLD-900  : Calcium Sulfate Dysprosium (CaSO4(Dy))",
    };
const char* usage = R"(
This tool is designed to automatically deconvolves the Thermoluminescent Dosimetry Materials glow curve which is being analyized and outputs the glow peaks. TLD glow curve data provided by the user in a properly fromatted .csv file(See provided sample file "GCA_Sample.csv")

Results: When in non-verbose mode, no output will be displayed while running analysis. A message will display successful analysis along with the goodness of fit value for the analysis. Output file containing glow peaks will be output to a unless specified default named zip file. If you desire a full output of the analysis progress, please use the verbose setting.

Usage .:
>> GCA -m <standard/analytical/smart> -f <filename.csv> -l <FOK/OTOR> -o <output_filename.csv>
Required-
-m / --mode : Mode to be run, standard, analytical, or smart. Unless specified defult Standard.
-f / --filename : <Path to .csv file> >> Required.
Optional-
-l / --model : to be used for analysis, First Order Kinetics Model(FOK) or OTOR Level(OTOR) model. Unless specified defult First Frder Kinetics.
-o / --output : <Path to output file> CSV formated output, defult glow_curve_report.csv (FILE WILL BE OVERWRITTEN).
-v / --verbose : results (Displays full progress of analysis) >> Recommended.
-h / --help : Display usage.
Example .:
>> GCA -m analytical -f </path/filename> -l FOK -v
)";

