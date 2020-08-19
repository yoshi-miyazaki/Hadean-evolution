# Makefile

CXX = /usr/local/bin/g++-9
CXXFLAGS = -c -Wall -std=c++11 -fopenmp

SRC = ../GibbsFE_minimization

early: main.o system.o atmos1d.o mantle.o element.o gibbse.o	massb.o solution.o CGgibbsmin.o
	$(CXX) -Wall -fopenmp -o early main.o system.o atmos1d.o mantle.o element.o gibbse.o	massb.o solution.o CGgibbsmin.o
main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp
system.o: system.cpp
	$(CXX) $(CXXFLAGS) system.cpp
atmos1d.o: atmos1d.cpp
	$(CXX) $(CXXFLAGS) atmos1d.cpp
mantle.o: mantle.cpp
	$(CXX) $(CXXFLAGS) mantle.cpp
element.o: $(SRC)/element.cpp
	$(CXX) $(CXXFLAGS) $(SRC)/element.cpp
gibbse.o: $(SRC)/gibbse.cpp
	$(CXX) $(CXXFLAGS) $(SRC)/gibbse.cpp
massb.o: $(SRC)/massb.cpp
	$(CXX) $(CXXFLAGS) $(SRC)/massb.cpp
solution.o: $(SRC)/solution.cpp
	$(CXX) $(CXXFLAGS) $(SRC)/solution.cpp
CGgibbsmin.o: $(SRC)/CGgibbsmin.cpp
	$(CXX) $(CXXFLAGS) $(SRC)/CGgibbsmin.cpp
clean:
	rm -f *.o early
