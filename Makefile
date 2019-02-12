# Run a regression test
test: main.exe
	./GCA -f test_1.csv -m smart

# Compile the main executable
main.exe: main.cpp First_Order_Kinetics.cpp File_Manager.cpp CSV_iterator.cpp
	g++ -std=c++0x -o GCA main.cpp First_Order_Kinetics.cpp File_Manager.cpp CSV_iterator.cpp

# Remove automatically generated files
clean :
	rm -rvf *.exe *~ *.out *.dSYM *.stackdump

# Disable built-in rules
.SUFFIXES:
