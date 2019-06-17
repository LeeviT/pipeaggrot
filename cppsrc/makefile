CC = g++
#CC = icc
# Default
#COPTS = -pg -DBOOST_DISABLE_ASSERTS
## Optimize 
COPTS = -std=c++11 -pg -O3 -m3dnow -mssse3 -finline-functions -DBOOST_DISABLE_ASSERTS
## Madman
# COPTS = -O3 -march=native -funroll-loops -ftree-loop-linear -fprefetch-loop-arrays -DBOOST_DISABLE_ASSERTS
# Debug options
# COPTS = -DDEBUG -g -pg -O5 -m3dnow -mssse3 -finline-functions -fprofile-arcs -ftest-coverage -ffast-math -mrecip
INCL = -I. -I/u/83/tuikkal1/unix/csm/pipeaggrot/cppsrc
# Linking needs -lm to enable certain mathematical functions
LINK = -L/u/83/tuikkal1/unix/csm/pipeaggrot/cppsrc/include/lib/ -lsundials_nvecserial -lsundials_cvode -lm -llapack -lblas -lgcov
# LINK = -L/home/mmo/lib/ -lm -llapack -lblas -lgcov
SOURCES = main.cpp Combo_class.cpp orientation_control_vol_acorrected.cpp
#SOURCES = main.cpp Combo_class.cpp orientation_control_vol.cpp
OUTFILE = run_control_vol_acorrected.exe
#OUTFILE = run_control_vol.exe

$(OUTFILE): $(SOURCES)
	$(CC) -o $(OUTFILE) $(COPTS) $(INCL) $(SOURCES) $(LINK)
clean:
	rm -rf *.exe *.o
pack:
	tar -cf ao_src.tar *.cpp *.hpp Doxyfile makefile 
documentation:
	doxygen
profile:
	gprof run_control_vol_acorrected.exe | ./gprof2dot.py | dot -Tpng -o profile.png
