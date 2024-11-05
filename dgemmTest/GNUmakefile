include ../Make.include

CXXFLAGS = -Wall
OPTFLAGS = -O3
BLASFLAGS = -framework Accelerate

all: naive.exe blas.exe

naive.exe: benchmark.o dgemm-naive.o
	$(CXX) $(CXXFLAGS) $(BLASFLAGS) -o naive.exe benchmark.o dgemm-naive.o

blas.exe: benchmark.o dgemm-blas.o
	$(CXX) $(OPTFLAGS) $(BLASFLAGS) -o blas.exe benchmark.o dgemm-blas.o

benchmark.o: benchmark.cpp
	$(CXX) $(CXXFLAGS) -c benchmark.cpp

dgemm-naive.o: dgemm-naive.cpp
	$(CXX) $(CXXFLAGS) -c dgemm-naive.cpp

dgemm-blas.o: dgemm-blas.cpp
	$(CXX) $(OPTFLAGS) -c dgemm-blas.cpp $(BLASFLAGS)

naive.out: naive.exe
	./naive.exe > naive.out

naive_opt.out: dgemm-naive.o benchmark.o
	$(CXX) $(OPTFLAGS) $(BLASFLAGS) -o naive.exe benchmark.o dgemm-naive.o
	./naive.exe > naive_opt.out

blas.out: blas.exe
	./blas.exe > blas.out

clean:
	rm -rf *.o *.exe *.d *.dSYM
