implicit double precision (a-h,o-z)
parameter(nbasis_max=100)
common/basis/nbasis

character*3 basis_type
common/bas/basis_type

parameter(ng_max=30,is_type_max=3)
parameter(nw=40)
parameter(nint_max=13000000)
parameter(n_contract_max=20)

logical finite_range,finite_range_save,one_elec_ZV,two_elec_ZV,one_elec_exact,one_elec_STO_exact
common/log/finite_range,one_elec_ZV,two_elec_ZV
common/fr/eps_rc,alpha(0:16),beta(0:16)
common/fr_int/n_eps

!***************************************************************************************
character*80 ATOM
parameter(number_nuclei_max=100)
common/nuclei_1/charge(number_nuclei_max),centers_nuclei(3,number_nuclei_max),enucl
common/nuclei_2/number_nuclei,nucleus_number(nbasis_max)
common/nuclei_3/ATOM(number_nuclei_max)
common/nuclei_4/G_center(3,nbasis_max,nbasis_max)
common/nuclei_5/dist_ij(nbasis_max,nbasis_max)
!****************************************************************************************

common/kcp/kcp_ijkl(nbasis_max,nbasis_max,nbasis_max,nbasis_max)
common/opt_integers_2/nint,is(nint_max),ks(nint_max),js(nint_max),ls(nint_max)

character*80 orb_b
parameter(n_atoms_max=92)
common/integer_orb/n_b(n_atoms_max),  &
                   n_cont_b(nbasis_max,n_atoms_max)
common/real_orb/gamma_b(nbasis_max,n_contract_max,n_atoms_max), &
                coef_b(nbasis_max,n_contract_max,n_atoms_max)
common/character_orb/orb_b(nbasis_max,n_atoms_max)

character*80 orb_name,orb_name_full
common/name/orb_name(nbasis_max),orb_name_full(nbasis_max)

!*****************************************************************************************************************
common/orbitals_real/center(3,nbasis_max),g_contract(n_contract_max,nbasis_max),c_contract(n_contract_max,nbasis_max)
common/orbitals_integer/npower(3,nbasis_max),i_type(nbasis_max),n_contract(nbasis_max)
parameter(n_gauss_max=ng_max*n_contract_max)
common/orbitals_real2/n_gauss(nbasis_max),c_gauss(2,n_gauss_max*ng_max,nbasis_max), &
g_gauss(n_gauss_max*ng_max,nbasis_max),g_min(nbasis_max),r_c(nbasis_max),g_slater(nbasis_max),r_infty(nbasis_max),&
r_c_prim(n_contract_max,nbasis_max),gauss_min(nbasis_max)
common/seuil_g/g_thr
common/i_seuil_g/i_add_large_g,i_add_2p,i_add_3d

common/star/n_star(n_contract_max,nbasis_max),ng_star(16)
common/act2/n_kt(n_contract_max,nbasis_max),ng_kt,ng(n_contract_max,nbasis_max),ng0,level
common/fit1/c_fit_GTO(ng_max,ng_max,0:16),g_fit_GTO(ng_max,ng_max,0:16)
common/fit2/c_fit_STO(ng_max,ng_max,0:16),g_fit_STO(ng_max,ng_max,0:16)
common/fit3/c_fit_exp(ng_max,ng_max),g_fit_exp(ng_max,ng_max)
common/fit4/n_STO(nbasis_max)

common/ZV/alpha_ZV,a_ZV(nbasis_max,nbasis_max)
!common/c_SCF/niter_SCF,nocc





