all: fdm2d

#fdm2d: domain.o field.o main.o cntrl.o wave1d.o source.o 
#	g++ -o fdm2d main.o domain.o field.o cntrl.o wave1d.o source.o  
fdm2d: domain.o main.o cntrl.o field.o source.o wave1d.o
	g++ -o fdm2d main.o domain.o cntrl.o  field.o  source.o wave1d.o

main.o: main.cpp fdm2d.h
	g++ -c main.cpp
domain.o: domain.cpp fdm2d.h
	g++ -c domain.cpp
cntrl.o: cntrl.cpp fdm2d.h
	g++ -c cntrl.cpp

field.o: field.cpp fdm2d.h
	g++ -c field.cpp
source.o: source.cpp fdm2d.h 
	g++ -c source.cpp
wave1d.o: wave1d.cpp  wave1d.h 
	g++ -c wave1d.cpp

