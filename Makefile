FC = gfortran
FFLAGS = -lstdc++ -g

all: test_gas_phase run

test_gas_phase: test_gas_phase.f90 /home/tapas/MBX_20200325_v0.2.2a/install/lib/libmbx.so
	$(FC) test_gas_phase.f90 /home/tapas/MBX_20200325_v0.2.2a/install/lib/libmbx.so $(FFLAGS) -o test_gas_phase

run: caleng_mbx.f90 /home/tapas/MBX_20200325_v0.2.2a/install/lib/libmbx.so
	$(FC) caleng_mbx.f90 /home/tapas/MBX_20200325_v0.2.2a/install/lib/libmbx.so $(FFLAGS) -o run

clean: 
	rm -f test_gas_phase run

.PHONY: all clean
