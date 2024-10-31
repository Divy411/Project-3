# Include the general make settings
include ../Make.include

# Common compiler flags
CXXFLAGS = -g -Wall
OPTFLAGS = -O3
BLASFLAGS = -framework Accelerate  # For MacOS. Use -lopenblas on Linux

# Targets
all: naive.exe blas.exe

# Compile the naive version
naive.exe: benchmark.o dgemm-naive.o
	$(CXX) $(CXXFLAGS) $(BLASFLAGS) -o naive.exe benchmark.o dgemm-naive.o

# Compile the BLAS-optimized version with optimization flags
blas.exe: benchmark.o dgemm-blas.o
	$(CXX) $(OPTFLAGS) $(BLASFLAGS) -o blas.exe benchmark.o dgemm-blas.o

# Compile individual source files
benchmark.o: benchmark.cpp
	$(CXX) $(CXXFLAGS) -c benchmark.cpp

dgemm-naive.o: dgemm-naive.cpp
	$(CXX) $(CXXFLAGS) -c dgemm-naive.cpp

dgemm-blas.o: dgemm-blas.cpp
	$(CXX) $(OPTFLAGS) -c dgemm-blas.cpp $(BLASFLAGS)

# Clean up
clean:
	rm -rf *.o *.exe *.d *.dSYM
