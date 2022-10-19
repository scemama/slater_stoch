FC = mpiifort
FFLAGS = -xHost -g -O0

default: one_e_slater two_e_slater

%.o: %.f90
	$(FC) $(FFLAGS) -c $^

ONE_FILES= \
 build_c_g_gauss_STO.o \
 chi2.o \
 number_atom.o \
 one_e_slater.o \
 read_fit_SMILES.o \
 read_basis.o \
 u_orb.o \
 u_gauss.o \
 write_STO_in_file_info_basis.o

TWO_FILES= \
 build_c_g_gauss_STO.o \
 build_mapping_ijkl.o \
 build_gaussian_expansion_of_orbitals.o \
 count_multi_center_integrals.o \
 chi2.o \
 gauss_ijkl.o \
 number_atom.o \
 one_elect.o \
 read_basis.o \
 read_fit_SMILES.o \
 u_orb.o \
 up_orb.o \
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
