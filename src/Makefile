FC = mpiifort
FFLAGS = -xHost -g -O0 -mkl=sequential

default: one_e_slater two_e_slater

%.o: %.f90
	$(FC) $(FFLAGS) -c $^

ONE_FILES= \
 build_gaussian_expansion_of_orbitals.o \
 fact.o \
 fact_inv.o \
 bigF.o \
 build_c_g_gauss_STO.o \
 chi2.o \
 eigsrt.o \
 gaunt.o \
 jacobi.o \
 norm_sto.o \
 number_atom.o \
 one_e_slater.o \
 one_elect.o \
 read_fit_SMILES.o \
 read_basis.o \
 read_geometry.o \
 smallg.o \
 u_orb.o \
 u_gauss.o \
 write_STO_in_file_info_basis.o

TWO_FILES= \
 compare_nuclei.o \
 dierfc.o \
 fact.o \
 fact_inv.o \
 bigF.o \
 build_c_g_gauss_STO.o \
 build_mapping_ijkl.o \
 build_gaussian_expansion_of_orbitals.o \
 compute_lmin_lmax.o \
 count_multi_center_integrals.o \
 chi2.o \
 eigsrt.o \
 gauss_ijkl.o \
 gaunt.o \
 ijkl_slater.o \
 ijkl_slater_xnynzn.o \
 jacobi.o \
 kinetic_slater.o \
 norm_sto.o \
 number_atom.o \
 one_elect.o \
 read_basis.o \
 read_geometry.o \
 read_fit_SMILES.o \
 smallg.o \
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
