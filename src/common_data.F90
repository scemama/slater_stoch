module common_data
  implicit none

  integer, parameter             :: nbasis_max=100
  integer, parameter             :: ng_max=30, is_type_max=3
  integer, parameter             :: nw=40
  integer, parameter             :: nint_max=13000000
  integer, parameter             :: n_contract_max=20
  integer, parameter             :: number_nuclei_max=100
  integer, parameter             :: n_atoms_max=92
  integer, parameter             :: n_gauss_max=ng_max*n_contract_max

  character(3) basis_type
  character(80)                  :: ATOM(number_nuclei_max)
  character(80)                  :: orb_b(nbasis_max,n_atoms_max)
  character(80)                  :: orb_name(nbasis_max), orb_name_full(nbasis_max)

  logical                        :: finite_range
  logical                        :: finite_range_save
  logical                        :: one_elec_ZV
  logical                        :: two_elec_ZV
  logical                        :: one_elec_exact
  logical                        :: one_elec_STO_exact

  integer                        :: nbasis
  integer                        :: n_eps
  integer                        :: n_b(n_atoms_max), n_cont_b(nbasis_max,n_atoms_max)
  integer                        :: number_nuclei
  integer                        :: nucleus_number(nbasis_max)
  integer                        :: kcp_ijkl(nbasis_max,nbasis_max,nbasis_max,nbasis_max)
  integer                        :: nint
  integer                        :: is(nint_max)
  integer                        :: ks(nint_max)
  integer                        :: js(nint_max)
  integer                        :: ls(nint_max)
  integer                        :: n_STO(nbasis_max)
  integer                        :: n_star(n_contract_max,nbasis_max),ng_star(16)
  integer                        :: n_kt(n_contract_max,nbasis_max),ng_kt,ng(n_contract_max,nbasis_max),ng0,level
  integer                        :: i_add_large_g,i_add_2p,i_add_3d
  integer                        :: npower(3,nbasis_max),i_type(nbasis_max),n_contract(nbasis_max)

  double precision               :: charge(number_nuclei_max)
  double precision               :: centers_nuclei(3,number_nuclei_max)
  double precision               :: enucl
  double precision               :: G_center(3,nbasis_max,nbasis_max)
  double precision               :: dist_ij(nbasis_max,nbasis_max)
  double precision               :: gamma_b(nbasis_max,n_contract_max,n_atoms_max)
  double precision               :: coef_b(nbasis_max,n_contract_max,n_atoms_max)
  double precision               :: eps_rc,alpha(0:16),beta(0:16)
  double precision               :: alpha_ZV,a_ZV(nbasis_max,nbasis_max)
  double precision               :: c_fit_GTO(ng_max,ng_max,0:16),g_fit_GTO(ng_max,ng_max,0:16)
  double precision               :: c_fit_STO(ng_max,ng_max,0:16),g_fit_STO(ng_max,ng_max,0:16)
  double precision               :: c_fit_exp(ng_max,ng_max),g_fit_exp(ng_max,ng_max)
  double precision               :: g_thr
  double precision               :: n_gauss(nbasis_max),c_gauss(2,n_gauss_max*ng_max,nbasis_max)
  double precision               :: g_gauss(n_gauss_max*ng_max,nbasis_max),g_min(nbasis_max)
  double precision               :: r_c(nbasis_max),g_slater(nbasis_max),r_infty(nbasis_max)
  double precision               :: r_c_prim(n_contract_max,nbasis_max),gauss_min(nbasis_max)
  double precision               :: center(3,nbasis_max),g_contract(n_contract_max,nbasis_max)
  double precision               :: c_contract(n_contract_max,nbasis_max)

end module
