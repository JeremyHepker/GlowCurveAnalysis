# Run a regression test
test: main.exe
	@./GCA /Users/jeremyhepker/Documents/NERS499/GlowCurveAnalsys/GlowCurveAnalsys/PeakDetection/
# Compile the main executable
main.exe: main.cpp First_Order_Kinetics.cpp File_Manager.cpp CSV_iterator.cpp
	@g++ -std=c++11 -o GCA main.cpp First_Order_Kinetics.cpp File_Manager.cpp matrix_arithmetic.cpp CSV_iterator.cpp -lboost_system -lboost_filesystem
# Remove automatically generated files
clean :
	rm -rvf *.exe *~ *.out *.dSYM *.stackdump output_*.csv

# Disable built-in rules
.SUFFIXES:
