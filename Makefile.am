ACLOCAL_AMFLAGS = -I m4

MODULES = trexio.mod common_data.mod

bin_PROGRAMS = one_e_slater two_e_slater
BUILT_SOURCES = $(MODULES)

trexio.mod: src/trexio.o  src/common_data.o
common_data.mod: src/common_data.o

CLEANFILES= $(MODULES)

COMMON_SOURCE = \
 src/trexio.F90 \
 src/common_data.F90 \
 src/add_poly_multiply.F90 \
 src/allocate_basis.F90 \
 src/bielec_integral.F90 \
 src/bigF.F90 \
 src/binom_func.F90 \
 src/binom.F90 \
 src/build_c_g_gauss_STO.F90 \
 src/build_gaussian_expansion_of_orbitals.F90 \
 src/build_mapping_ijkl.F90 \
 src/davidson.F90 \
 src/doublefact.F90 \
 src/fact.F90 \
 src/gaunt.F90 \
 src/gaussian_product.F90 \
 src/gauss.F90 \
 src/general_primitive_integral.F90 \
 src/give_explicit_poly_and_gaussian.F90 \
 src/give_polynom_mult_center_mono_elec.F90 \
 src/give_polynom_mult_center_x.F90 \
 src/hermite.F90 \
 src/I_x1_pol_mult_a1.F90 \
 src/I_x1_pol_mult_a2.F90 \
 src/I_x1_pol_mult_mono_elec.F90 \
 src/I_x1_pol_mult.F90 \
 src/I_x1_pol_mult_recurs.F90 \
 src/I_x2_pol_mult_mono_elec.F90 \
 src/I_x2_pol_mult.F90 \
 src/logfact.F90 \
 src/multiply_poly_2.F90 \
 src/multiply_poly.F90 \
 src/NAI_pol_mult.F90 \
 src/norm_sto.F90 \
 src/number_atom.F90 \
 src/one_elect.F90 \
 src/one_electron_I.F90 \
 src/one_electron_Kij.F90 \
 src/one_electron_sij.F90 \
 src/one_electron_Vij.F90 \
 src/read_basis.F90 \
 src/read_fit_SMILES.F90 \
 src/read_geometry.F90 \
 src/rint1.F90 \
 src/rint_large_n.F90 \
 src/rint.F90 \
 src/rint_sum.F90 \
 src/smallg.F90 \
 src/u_gauss.F90 \
 src/u_orb.F90

src/trexio.$(OBJEXT): src/common_data.$(OBJEXT)

MODULES_O = src/common_data.$(OBJEXT) src/trexio.$(OBJEXT)

one_e_slater_SOURCES = $(MODULES_O) $(COMMON_SOURCE) \
 src/one_e_slater.F90

if HAVE_MMAP
MMAP_FILES = src/fortran_mmap.c src/mmap.F90
else
MMAP_FILES =
endif

two_e_slater_SOURCES = $(MODULES_O) $(COMMON_SOURCE) $(MMAP_FILES) \
 src/compare_nuclei.F90 \
 src/compute_densities.F90 \
 src/compute_int_slater.F90 \
 src/compute_jacobian_and_pi0.F90 \
 src/compute_lmin_lmax.F90 \
 src/count_multi_center_integrals.F90 \
 src/davidson_clean.F90 \
 src/dierfc.F90 \
 src/draw_configuration.F90 \
 src/gauss_ijkl.F90 \
 src/ijkl_slater.F90 \
 src/ijkl_slater_xnynzn.F90 \
 src/kinetic_slater.F90 \
 src/ortho_qr.F90 \
 src/ortho_svd.F90 \
 src/randomized_svd.F90 \
 src/svd.F90 \
 src/svd_clean.F90 \
 src/two_e_slater.F90


