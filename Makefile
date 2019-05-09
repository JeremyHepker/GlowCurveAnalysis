# Run a regression test
test: main.exe
	@./GCA -f Tests/test1.csv -m smart -o output_1
	# #@./GCA -f Tests/test2.csv -m smart -o output_2
	@./GCA -f Tests/test3.csv -m smart -o output_3
	@./GCA -f Tests/test4.csv -m smart -o output_4
	@./GCA -f Tests/test5.csv -m smart -o output_5
	@./GCA -f Tests/test6.csv -m smart -o output_6
	@./GCA -f Tests/test7.csv -m smart -o output_7
	@./GCA -f Tests/test8.csv -m smart -o output_8
	@./GCA -f Tests/test9.csv -m smart -o output_9
	@./GCA -f Tests/test10.csv -m smart -o output_10
	@./GCA -f Tests/test11.csv -m smart -o output_11
	@./GCA -f Tests/test12.csv -m smart -o output_12
	@./GCA -f Tests/test13.csv -m smart -o output_13
	# @./GCA -f Tests/test15.csv -m smart -o output_14
	# @./GCA -f Tests/test16.csv -m smart -o output_15
# Compile the main executable
main.exe: main.cpp First_Order_Kinetics.cpp File_Manager.cpp CSV_iterator.cpp
	@g++ -std=c++0x -o GCA main.cpp First_Order_Kinetics.cpp File_Manager.cpp matrix_arithmetic.cpp CSV_iterator.cpp

# Remove automatically generated files
clean :
	rm -rvf *.exe *~ *.out *.dSYM *.stackdump output_*.csv

# Disable built-in rules
.SUFFIXES:
