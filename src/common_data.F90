module common_data
  implicit none

  integer, parameter             :: nw=40             ! Number of independent walkers for the stochastic calculations

  ! Gaussian fits
  integer, parameter             :: ng_max=30                              ! Max number of Gaussians for fits
  integer                        :: ng0                                    ! Number of Gaussians for the STO fit
  double precision               :: g_thr                                  ! Threshold to identify large exponents
  integer                        :: i_add_2p                               ! How many functions to add for 2p's
  integer                        :: i_add_3d                               ! How many functions to add for 3d's
  integer                        :: i_add_large_g                          ! How many functions to add for large exponents

  ! Two-electron integrals. allocated in build_mapping_ijkl 
  integer*8                      :: nint                                   ! Number of <ij|kl> integrals
  integer, allocatable           :: is(:), ks(:), js(:), ls(:)             ! Indices i,j,k,l for each integral

  ! Basis per atom type
  integer, parameter             :: n_atoms_max=92                         ! Max atomic number
  integer, parameter             :: n_contract_max=20                      ! Max number of contractions
  integer, parameter             :: n_gauss_max=ng_max*n_contract_max
  integer, parameter             :: n_slater_max=32                        ! Max number of STO per atom type
  integer, allocatable           :: n_b(:)                                 ! number of basis functions per atom type (S,P,D,F, or G) 
  integer, allocatable           :: n_cont_b(:,:)                          ! number of contracted primitives for the radial part
  double precision, allocatable  :: gamma_b(:,:,:)                         ! Exponent in basis
  double precision, allocatable  :: coef_b(:,:,:)                          ! Contraction coefficient in basis
  character(80), allocatable     :: orb_b(:,:)                             ! Name of the basis function, that is S,P,D,F or G

  ! Allocated in read_geometry
  integer                        :: number_nuclei                            ! Number of nuclei in molecule
  character(2), allocatable      :: ATOM(:)                                  ! Atom labels
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
  double precision, allocatable  :: c_gauss(:,:,:)                           ! Coefficients of contracted Gaussians
  double precision, allocatable  :: g_gauss(:,:)                             ! Exponents of contracted Gaussians
  double precision, allocatable  :: g_slater(:)                              ! Exponent of STO
  double precision, allocatable  :: g_contract(:,:)                          ! Exponents of contracted STOs
  double precision, allocatable  :: c_contract(:,:)                          ! Contraction coefficients of STO

  double precision               :: c_fit_exp(ng_max,ng_max)                 ! SMILES fit
  double precision               :: g_fit_exp(ng_max,ng_max)                 ! SMILES fit
  double precision               :: a_ZV

end module
