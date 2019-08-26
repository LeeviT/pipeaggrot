CC = mpiicpc
# Default
#COPTS = -pg -DBOOST_DISABLE_ASSERTS
## Optimize 
## COPTS = -std=c++11 -fopenmp -pg -O3 -mssse3 -finline-functions -DBOOST_DISABLE_ASSERTS
COPTS = -std=c++11 -qopenmp -fopenmp -pg -O2 -DBOOST_DISABLE_ASSERTS
## Madman
# COPTS = -O3 -march=native -funroll-loops -ftree-loop-linear -fprefetch-loop-arrays -DBOOST_DISABLE_ASSERTS
# Debug options
# COPTS = -DDEBUG -g -pg -O5 -m3dnow -mssse3 -finline-functions -fprofile-arcs -ftest-coverage -ffast-math -mrecip
INCL = -I. -I/u/83/tuikkal1/unix/csm/pipeaggrot/cppsrc/include/
# Linking needs -lm to enable certain mathematical functions
LD_LIBRARY_PATH = /u/83/tuikkal1/unix/csm/pipeaggrot/cppsrc/lib/
LINK = -L/u/83/tuikkal1/unix/csm/pipeaggrot/cppsrc/lib -lsundials_cvode -lsundials_sunlinsolspgmr -lsundials_nvecserial -lsundials_nvecopenmp -lm -llapack -lblas -lgcov
SOURCES = main.cpp combo.cpp orientation.cpp velo.cpp io.cpp point.cpp
OUTFILE = pipeaggrot.exe

$(OUTFILE): $(SOURCES)
	$(CC) -o $(OUTFILE) $(COPTS) $(INCL) $(SOURCES) $(LINK) 
	mv -f $(OUTFILE) ../

clean:
	rm -rf *.exe *.o

pack:
	tar -cf ao_src.tar *.cpp *.hpp Doxyfile makefile 

documentation:
	doxygen

profile:
	gprof pipeaggrot.exe | ./gprof2dot.py | dot -Tpng -o profile.png
