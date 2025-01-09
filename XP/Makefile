CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3 -Wno-comment

IGSL = ${shell gsl-config --cflags}
LGSL = ${shell gsl-config --libs}
# IGSL = -I/Users/duanpanpan/Desktop/Work/Model/GSL/gsl-install/include
# LGSL = -L/Users/duanpanpan/Desktop/Work/Model/GSL/gsl-install/lib -lgsl -lgslcblas -lm


TARGET = test.exe
all: $(TARGET)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(IGSL) -c main.cpp

$(TARGET): main.o
	$(CXX) $(CXXFLAGS) $^ $(LGSL) -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f *.o *.*~ *.exe