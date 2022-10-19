FC = mpiifort
FFLAGS = -xHost -g -O0

default: one_e_slater two_e_slater

%.o: %.f90
	$(FC) $(FFLAGS) -c $^

ONE_FILES= \
 build_c_g_gauss_STO.o \
 chi2.o \
 one_e_slater.o \
 u_gauss.o \
 write_STO_in_file_info_basis.o

TWO_FILES= \
 build_c_g_gauss_STO.o \
 chi2.o \
 u_gauss.o \
 up_gauss.o \
 two_e_slater.o \
 write_STO_in_file_info_basis.o

one_e_slater: $(ONE_FILES)
	$(FC) $(FFLAGS) $^ -o $@

two_e_slater: $(TWO_FILES)
	$(FC) $(FFLAGS) $^ -o $@

clean:
	rm *.o

.PHONY: clean default
