FC=gfortran
TARGETS=ising-2d Tc_ising matrix_ising

all: ${TARGETS}

ising-2d: ising-2d.f90
	${FC} -o $@ $^ -O3

Tc_ising: Tc_ising.f90
	${FC} -o $@ $^ -O3

matrix_ising: matrix_ising.f90
	${FC} -o $@ $^ -O3

	@echo "--------------------------------------------------------------------------"
	@echo "If you now wish to obtain useful results from the compiled codes please run shell script ./run_ising.sh"
	@echo "--------------------------------------------------------------------------"

clean:
	rm -f ${TARGETS}

