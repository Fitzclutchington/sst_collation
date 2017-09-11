PROG = ahil2c
CXX = g++
CXXFLAGS = -g -Wall -Wpedantic -Wfatal-errors -Wunused-parameter -march=native -O2 -fopenmp -std=c++11 
LD = g++
LDFLAGS =\
	-fopenmp\
	-lm\
	-lnetcdf\
	-lopencv_core\
	-lopencv_imgproc\
	-lopencv_highgui\
	-Wl,--no-as-needed\
	
OFILES = \
	main.o\
	hermite.o\

HFILES = \
	hermite.h

LOC =\
	eigen-eigen-26667be4f70b\

all: $(PROG)

%.o: %.cc $(HFILES)
	$(CXX) -I $(LOC) -c $(CXXFLAGS) $<

$(PROG): $(OFILES)
	$(LD)  -o $(PROG) $(OFILES) $(LDFLAGS)

clean:
	rm -f $(PROG) $(OFILES)
