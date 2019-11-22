# Start of the makefile
# Defining variables
objects = global.o mpiown.o init.o diff.o diffall.o eqeuler.o derivs.o rhs.o  cc.o outt.o filtering.o rk.o coefficients.o euler.o 
#
##f90comp = gfortran -O3
#f90comp =  gfortran -O3 -fcheck= -Wline-truncation
f90comp =  mpif90 -mcmodel=medium -O3 -fcheck= -Wline-truncation
#f90comp =  mpif90 -mcmodel=medium -Ofast
#f90comp =  mpif90  -Ofast

#
## Makefile
#
euler: $(objects)
	$(f90comp) -o euler $(objects)
	
global.o:  global.f90
	$(f90comp) -c global.f90

mpiown.o:  global.f90 mpiown.f90
	$(f90comp) -c mpiown.f90

coefficients.o: global.o coefficients.f90
	$(f90comp) -c coefficients.f90	

diff.o: global.o   diff.f90
	$(f90comp) -c diff.f90

diffall.o: global.o   diffall.f90
	$(f90comp) -c diffall.f90

init.o: global.o coefficients.o  diff.o outt.o init.f90
	$(f90comp) -c init.f90

eqeuler.o: global.o   eqeuler.f90
	$(f90comp) -c eqeuler.f90

derivs.o: global.o  diff.o derivs.f90
	$(f90comp) -c derivs.f90

rhs.o:  global.o coefficients.o derivs.o   eqeuler.o rhs.f90
	$(f90comp) -c rhs.f90	

filtering.o:  global.o  filtering.f90
	$(f90comp) -c filtering.f90

cc.o:  global.o  cc.f90
	$(f90comp) -c cc.f90

outt.o: global.o diffall.o  outt.f90
	$(f90comp) -c outt.f90

rk.o: global.o  rhs.o  diff.o cc.o rk.f90
	$(f90comp) -c rk.f90

euler.o: global.o mpiown.o init.o rk.o  outt.o filtering.o euler.f90
	$(f90comp) -c euler.f90	


## Cleaning everything
clean:
#	rm MOD1.mod FLAT
	rm -r -f *.o 
	rm -r -f *.mod 
	rm -r -f *.dat 
	rm -r -f fort.* 
	rm -r $(objects)
#
## End of the makefile

#objects = global.o  globalq.o global2.o init.o cc.o diff.o \
#          derivs.o  eqpml.o eqeuler.o rkpml.o \
#          rhs.o rhsq.o outt.o filtering.o euler.o 
#
##f90comp = gfortran -O3
#f90comp = gfortran -mcmodel=medium -O3
#
## Makefile
#
#euler: $(objects)
#	$(f90comp) -o euler $(objects)
#	
#global.o:  global.f90
#	$(f90comp) -c global.f90
#
#outt.o: global.o diff.o  outt.f90
#	$(f90comp) -c outt.f90
#
#globalq.o: global.o  globalq.f90	
#	$(f90comp) -c globalq.f90
#
#global2.o: global.o  global2.f90	
#	$(f90comp) -c global2.f90
#
#init.o: global.o  globalq.o  derivs.o init.f90
#	$(f90comp) -c init.f90
#
#eqeuler.o: global.o  global2.o eqeuler.f90
#	$(f90comp) -c eqeuler.f90
#
#diff.o: global.o   diff.f90
#	$(f90comp) -c diff.f90
#
#derivs.o: global.o  global2.o  diff.o derivs.f90
#	$(f90comp) -c derivs.f90
#
#eqpml.o: global.o  global2.o eqpml.f90
#	$(f90comp) -c eqpml.f90
#
#rhs.o:  global.o  derivs.o   eqeuler.o  eqpml.o rhs.f90
#	$(f90comp) -c rhs.f90	
#
#rhsq.o:  global.o   globalq.o global2.o rhsq.f90
#	$(f90comp) -c rhsq.f90	
#
#rkpml.o: global.o  globalq.o rhs.o rhsq.o cc.o rkpml.f90
#	$(f90comp) -c rkpml.f90	
#
#cc.o:  global.o  cc.f90
#	$(f90comp) -c cc.f90
#
#filtering.o:  global.o  globalq.o filtering.f90
#	$(f90comp) -c filtering.f90
#
#euler.o: global.o globalq.o init.o filtering.o outt.o  rkpml.o   euler.f90
#	$(f90comp) -c euler.f90	
#
## Cleaning everything
#clean:
#	#rm MOD1.mod FLAT
#	rm -r *.dat 
#	rm -r $(objects)
#
## End of the makefile
