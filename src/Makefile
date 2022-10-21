FC = mpiifort
#FFLAGS = -xHost -g -O0 -mkl=sequential
FFLAGS = -xHost -g -O2 -mkl=sequential 

default: one_e_slater two_e_slater

%.o: %.f90
	$(FC) $(FFLAGS) -c $^

ONE_FILES= \
 add_poly_multiply.o \
 bielec_integral.o \
 bigF.o \
 binom_func.o \
 binom.o \
 binom_transp.o \
 build_c_g_gauss_STO.o \
 build_gaussian_expansion_of_orbitals.o \
 chi2.o \
 doublefact.o \
 eigsrt.o \
 fact_inv.o \
 fact.o \
 gaunt.o \
 gaussian_product.o \
 general_primitive_integral.o \
 give_explicit_poly_and_gaussian.o \
 give_polynom_mult_center_mono_elec.o \
 give_polynom_mult_center_x.o \
 hermite.o \
 I_x1_new.o \
 I_x2_new.o \
 I_x1_pol_mult_a1.o \
 I_x1_pol_mult_a2.o \
 I_x1_pol_mult_mono_elec.o \
 I_x1_pol_mult.o \
 I_x1_pol_mult_recurs.o \
 I_x2_pol_mult_mono_elec.o \
 I_x2_pol_mult.o \
 jacobi.o \
 logfact.o \
 multiply_poly_2.o \
 multiply_poly.o \
 NAI_pol_mult.o \
 norm_sto.o \
 number_atom.o \
 one_elect.o \
 one_electron_I.o \
 one_electron_Kij.o \
 one_electron_sij.o \
 one_electron_Vij.o \
 one_e_slater.o \
 read_basis.o \
 read_fit_SMILES.o \
 read_geometry.o \
 recentered_poly2.o \
 rint1.o \
 rint_large_n.o \
 rint.o \
 rint_sum.o \
 smallg.o \
 u_gauss.o \
 u_orb.o \
 write_STO_in_file_info_basis.o


TWO_FILES= \
 add_poly_multiply.o \
 bielec_integral.o \
 bigF.o \
 binom_func.o \
 binom.o \
 binom_transp.o \
 build_c_g_gauss_STO.o \
 build_gaussian_expansion_of_orbitals.o \
 build_mapping_ijkl.o \
 chi2.o \
 compare_nuclei.o \
 compute_densities.o \
 compute_int_slater.o \
 compute_jacobian_and_pi0.o \
 compute_lmin_lmax.o \
 compute_r_tilde.o \
 compute_u_tilde.o \
 count_multi_center_integrals.o \
 dierfc.o \
 doublefact.o \
 draw_configuration.o \
 eigsrt.o \
 fact_inv.o \
 fact.o \
 gaunt.o \
 gaussian_product.o \
 gauss_ijkl.o \
 general_primitive_integral.o \
 give_explicit_poly_and_gaussian.o \
 give_polynom_mult_center_mono_elec.o \
 give_polynom_mult_center_x.o \
 hermite.o \
 ijkl_slater.o \
 ijkl_slater_xnynzn.o \
 I_x1_pol_mult_a1.o \
 I_x1_pol_mult_a2.o \
 I_x1_pol_mult_mono_elec.o \
 I_x1_pol_mult.o \
 I_x1_pol_mult_recurs.o \
 I_x2_pol_mult_mono_elec.o \
 I_x2_pol_mult.o \
 jacobi.o \
 kinetic_slater.o \
 logfact.o \
 multiply_poly_2.o \
 multiply_poly.o \
 NAI_pol_mult.o \
 norm_sto.o \
 number_atom.o \
 one_elect.o \
 one_electron_I.o \
 one_electron_Kij.o \
 one_electron_sij.o \
 one_electron_Vij.o \
 read_basis.o \
 read_fit_SMILES.o \
 read_geometry.o \
 recentered_poly2.o \
 rint1.o \
 rint_large_n.o \
 rint.o \
 rint_sum.o \
 smallg.o \
 two_e_slater.o \
 u_gauss.o \
 u_orb.o \
 up_gauss.o \
 up_orb.o \
 write_STO_in_file_info_basis.o

one_e_slater: $(ONE_FILES)
	$(FC) $(FFLAGS) $^ -o $@

two_e_slater: $(TWO_FILES)
	$(FC) $(FFLAGS) $^ -o $@

clean:
	rm *.o

.PHONY: clean default
