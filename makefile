all: fdm2d 
#../plot_inwv

fdm2d: domain.o main.o indx_cod.o field.o\
	PML.o InWv.o cmplx.o source.o recs.o crack.o COD.o icon.o
	g++ -o ../fdm2d domain.o main.o indx_cod.o field.o PML.o InWv.o cmplx.o source.o recs.o crack.o COD.o icon.o

#../plot_inwv: cmplx.o plot_inwv.o InWv.o
#	g++ -o ../plot_inwv cmplx.o plot_inwv.o InWv.o

plot_inwv.o: plot_inwv.cpp
	g++ -c plot_inwv.cpp
domain.o: domain.cpp
	g++ -c domain.cpp
main.o: main.cpp
	g++ -c main.cpp
indx_cod.o: indx_cod.cpp
	g++ -c indx_cod.cpp
field.o: field.cpp
	g++ -c field.cpp
PML.o: PML.cpp
	g++ -c PML.cpp
InWv.o: InWv.cpp
	g++ -c InWv.cpp
source.o: source.cpp
	g++ -c source.cpp
recs.o: recs.cpp
	g++ -c recs.cpp
crack.o: crack.cpp
	g++ -c crack.cpp
icno.o: icon.cpp
	g++ -c icon.cpp
COD.o: COD.cpp
	g++ -c COD.cpp
cmplx.o: cmplx.cpp
	g++ -c cmplx.cpp
#rayleigh.o: rayleigh.cpp
	#g++ -c rayleigh.cpp
