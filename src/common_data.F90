module common_data
  implicit none

  integer, parameter             :: nbasis_max=100     ! Max number of basis functions

  integer, parameter             :: nw=40              ! Number of independent walkers for the stochastic calculations

  ! Gaussian fits
  integer, parameter             :: ng_max=30                              ! Max number of Gaussians for fits
  integer                        :: ng0                                    ! Number of Gaussians for the STO fit
  integer                        :: i_add_2p                               ! How many functions to add for 2p's
  integer                        :: i_add_3d                               ! How many functions to add for 3d's
  integer                        :: i_add_large_g                          ! How many functions to add for large exponents

  ! Two-electron integrals. allocated in build_mapping_ijkl 
  integer                        :: nint                                   ! Number of <ij|kl> integrals
  integer, allocatable           :: is(:), ks(:), js(:), ls(:)             ! Indices i,j,k,l for each integral

  ! Basis per atom type
  integer, parameter             :: n_atoms_max=92                         ! Max atomic number
  integer, parameter             :: n_contract_max=20                      ! Max number of contractions
  integer, parameter             :: n_gauss_max=ng_max*n_contract_max
  integer, parameter             :: n_slater_max=100                       ! Max number of STO per atom type
  integer                        :: n_b(n_atoms_max)                       ! number of basis functions per atom type (S,P,D,F, or G) 
  integer                        :: n_cont_b(n_contract_max, n_atoms_max)  ! number of contracted primitives for the radial part
  double precision               :: gamma_b(n_slater_max,n_contract_max,n_atoms_max) ! Exponent in basis
  double precision               :: coef_b(n_slater_max,n_contract_max,n_atoms_max)  ! Contraction coefficient in basis
  character(80)                  :: orb_b(n_slater_max,n_atoms_max)                  ! name of the basis function, that is S,P,D,F or G


  ! Allocated in read_geometry
  integer                        :: number_nuclei                            ! Number of nuclei in molecule
  character(80), allocatable     :: ATOM(:)                                  ! Atom labels
  double precision, allocatable  :: charge(:)                                ! Nuclear charges
  double precision, allocatable  :: centers_nuclei(:,:)                      ! XYZ coordinates of nuclei
  double precision               :: enucl                                    ! Nuclear repulsion energy

  ! Allocated in allocate_basis
  integer                        :: nbasis                            ! Number of basis functions
  character(80), allocatable     :: orb_name(:)                       ! Name of the orbital (1S, 2S, 3P) 
  character(80), allocatable     :: orb_name_full(:)                  ! Full name of the orbital (5G_XYZZ)
  integer, allocatable           :: nucleus_number(:)                 ! Nucleus on which the function is centered
  integer, allocatable           :: n_STO(:)                          ! Value of n in STO, in r^n.exp(-a.r)
  integer, allocatable           :: npower(:,:)                       ! Polynomial part: powers of x,y,z
  integer, allocatable           :: i_type(:)                         ! 1: c.exp(-a.r)  2: c.r^2.exp(-a.r)
  integer, allocatable           :: n_contract(:)                     ! Number of contracted STOs
  integer, allocatable           :: ng(:,:)                           ! Number of Gaussians used to fit a contracted STO
  integer, allocatable           :: n_gauss(:)                        ! Number of Gaussians used to fit an STO
  double precision, allocatable  :: center(:,:)                       ! Where the basis functions are centered

  double precision               :: a_ZV(nbasis_max,nbasis_max)              ! Zero-variance parameter
  double precision               :: c_fit_exp(ng_max,ng_max)                 ! SMILES fit
  double precision               :: g_fit_exp(ng_max,ng_max)                 ! SMILES fit
  double precision               :: g_thr                                    ! Threshold to identify large exponents
  double precision               :: c_gauss(n_gauss_max*ng_max,nbasis_max,2) ! Coefficients of contracted Gaussians
  double precision               :: g_gauss(n_gauss_max*ng_max,nbasis_max)   ! Exponents of contracted Gaussians
  double precision               :: g_slater(nbasis_max)                     ! Exponent of STO
  double precision               :: g_contract(n_contract_max,nbasis_max)    ! Exponents of contracted STOs
  double precision               :: c_contract(n_contract_max,nbasis_max)    ! Contraction coefficients of STO

end module
