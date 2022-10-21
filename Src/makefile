all: ../fdm2d ../plot_inwv ../bnd ../bnd_tight

../fdm2d: fdm2d.h domain.o main.o indx_cod.o field.o\
	PML.o InWv.o cmplx.o source.o recs.o crack.o COD.o icon.o
	g++ -o ../fdm2d domain.o main.o indx_cod.o field.o PML.o InWv.o cmplx.o source.o recs.o crack.o COD.o icon.o

../plot_inwv: cmplx.o plot_inwv.o InWv.o
	g++ -o ../plot_inwv cmplx.o plot_inwv.o InWv.o

../bnd: bnd.o
	g++ -o ../bnd bnd.o
../bnd_tight: bnd_tight.o
	g++ -o ../bnd_tight bnd_tight.o

plot_inwv.o: plot_inwv.cpp
	g++ -c plot_inwv.cpp
domain.o: domain.cpp fdm2d.h
	g++ -c domain.cpp
main.o: main.cpp fdm2d.h crack.h
	g++ -c main.cpp
indx_cod.o: indx_cod.cpp fdm2d.h
	g++ -c indx_cod.cpp
field.o: field.cpp fdm2d.h
	g++ -c field.cpp
PML.o: PML.cpp fdm2d.h
	g++ -c PML.cpp
InWv.o: InWv.cpp fdm2d.h
	g++ -c InWv.cpp
source.o: source.cpp fdm2d.h
	g++ -c source.cpp
recs.o: recs.cpp fdm2d.h
	g++ -c recs.cpp
crack.o: crack.cpp fdm2d.h crack.h
	g++ -c crack.cpp
icno.o: icon.cpp fdm2d.h
	g++ -c icon.cpp
COD.o: COD.cpp fdm2d.h crack.h
	g++ -c COD.cpp
cmplx.o: cmplx.cpp fdm2d.h
	g++ -c cmplx.cpp
rayleigh.o: rayleigh.cpp
	g++ -c rayleigh.cpp

bnd.o:bnd.cpp
	g++ -c bnd.cpp
bnd_tight.o:bnd_tight.cpp
	g++ -c bnd_tight.cpp
