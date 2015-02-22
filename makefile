FC = gfortran
LFLAGS = -O3 -fbounds-check 

SUBFILES = global.f90 ./ADImods/modA.f90 ./ADImods/modB.f90 ./ADImods/modC.f90 ./ADImods/modE.f90 ./ADImods/modD.f90 mod1.f90 modfreeE.f90 modphi.f90 modhomfE.f90 modallocate.f90 TDMAmod.f90 secantmod.f90

main: main.f90 
	$(FC) $(LFLAGS) -o $@ $(MOBLIB) $(SUBFILES) main.f90 $(LIBS)

clean:
	rm *~ -o -rf
