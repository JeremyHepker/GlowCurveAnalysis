# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -std=c++0x -g -Wall -lboost_system -lboost_filesystem

# the build target executable:
TARGET = GCA

all: $(TARGET)

$(TARGET): src/smartPeakDetect.cpp src/OTORModel.cpp src/matrix_arithmetic.cpp src/Levenberg–Marquardt.cpp src/DataSmoothing.cpp src/CSV_iterator.cpp src/main.cpp src/FOKModel.cpp src/File_Manager.cpp
	$(CC) $(CFLAGS) -o $(TARGET) src/smartPeakDetect.cpp src/OTORModel.cpp src/matrix_arithmetic.cpp src/Levenberg–Marquardt.cpp src/DataSmoothing.cpp src/CSV_iterator.cpp src/main.cpp src/FOKModel.cpp src/File_Manager.cpp

clean:
	$(RM) $(TARGET)


