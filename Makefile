# Makefile

CXX = /usr/local/bin/g++-9
CXXFLAGS = -c -Wall -std=c++11 -fopenmp

early: main.o system.o mantle.o element.o gibbse.o	massb.o solution.o CGgibbsmin.o
	$(CXX) -Wall -fopenmp -o early main.o system.o mantle.o element.o gibbse.o	massb.o solution.o CGgibbsmin.o
main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp
system.o: system.cpp
	$(CXX) $(CXXFLAGS) system.cpp
mantle.o: mantle.cpp
	$(CXX) $(CXXFLAGS) mantle.cpp
element.o: element.cpp
	$(CXX) $(CXXFLAGS) element.cpp
gibbse.o: gibbse.cpp
	$(CXX) $(CXXFLAGS) gibbse.cpp
massb.o: massb.cpp
	$(CXX) $(CXXFLAGS) massb.cpp
solution.o: solution.cpp
	$(CXX) $(CXXFLAGS) solution.cpp
CGgibbsmin.o: CGgibbsmin.cpp
	$(CXX) $(CXXFLAGS) CGgibbsmin.cpp
clean:
	rm -f *.o early
