program integrals
include 'j.inc'
!!!!!!!!!!!COMMENT IF NOT MPI
!  include 'mpif.h'
!!!!!!!!!!!COMMENT IF NOT MPI

character*80,charabid
character*80,MOLECULE

logical normalization_OA

! Arrays for one-electron integrals
double precision S_ij_gaus(nbasis_max,nbasis_max), &
                 V_ij_gaus(nbasis_max,nbasis_max),K_ij_gaus(nbasis_max,nbasis_max),Kij
double precision S_ij_ZV(nbasis_max,nbasis_max), &
                 V_ij_ZV(nbasis_max,nbasis_max),K_ij_ZV(nbasis_max,nbasis_max)
double precision S_ij_ZV2(nbasis_max,nbasis_max), &
                 V_ij_ZV2(nbasis_max,nbasis_max),K_ij_ZV2(nbasis_max,nbasis_max)
double precision S_ij_ex(nbasis_max,nbasis_max), &
                 V_ij_ex(nbasis_max,nbasis_max),K_ij_ex(nbasis_max,nbasis_max)
double precision S_ij_STO_ex(nbasis_max,nbasis_max), &
                 V_ij_STO_ex(nbasis_max,nbasis_max),K_ij_STO_ex(nbasis_max,nbasis_max)
double precision precond(nbasis_max,nbasis_max)
double precision sqrtSii_gaus  (nbasis_max)
double precision sqrtSii_STO_ex(nbasis_max)

!! BEGIN BIG ARRAYS nint_max 
double precision ijkl_gaus(nint_max)
double precision ijkl(nint_max),ijkl2(nint_max)

dimension e_S_ijkl(nint_max)
dimension e_G_ijkl(nint_max)
dimension e_test_ijkl(nint_max)
dimension error_ijkl(nint_max)
dimension mono_center(nint_max)
double precision moy(nint_max), moy2(nint_max), moy2t(nint_max)
!! END BIG ARRAYS nint_max 

dimension e_Sij(nbasis_max,nbasis_max)
dimension e_Vij(nbasis_max,nbasis_max)
dimension e_Kij(nbasis_max,nbasis_max)
dimension e_G_Sij(nbasis_max,nbasis_max)
dimension e_G_Vij(nbasis_max,nbasis_max)
dimension e_G_Kij(nbasis_max,nbasis_max)
dimension error_Sij(nbasis_max,nbasis_max)
dimension error_Vij(nbasis_max,nbasis_max)
dimension error_Kij(nbasis_max,nbasis_max)
double precision moy_S(nbasis_max,nbasis_max), moy_S2(nbasis_max,nbasis_max), moy_S2t(nbasis_max,nbasis_max)
double precision moy_V(nbasis_max,nbasis_max), moy_V2(nbasis_max,nbasis_max), moy_V2t(nbasis_max,nbasis_max)
double precision moy_K(nbasis_max,nbasis_max), moy_K2(nbasis_max,nbasis_max), moy_K2t(nbasis_max,nbasis_max)

!! MONTE CARLO PART
integer i_tab_mc(nbasis_max,nbasis_max)
double precision :: d_x(4)
logical clean,mono_center_Slater

dimension r12_inv(nw),r1(nw,3),r2(nw,3)
dimension pi_0(nw),pot(nw),rkin(nw),rkin_G(nw)
dimension rt1(nw,3),rt2(nw,3)
dimension ut1(3,nw,nbasis_max*nbasis_max),ut2(3,nw,nbasis_max*nbasis_max)
dimension rho(nw,nbasis_max*nbasis_max,2),poly(nw,nbasis_max*nbasis_max,2)
dimension weight(nw),weight_kin(nw) 
dimension rho_G(nw,nbasis_max*nbasis_max,2)
dimension weight_G(nw),weight_kin_G(nw)
dimension rjacob(nw)

character*(128) :: filename_in
character*(128) :: filename_out_gaus_s
character*(128) :: filename_out_gaus_v
character*(128) :: filename_out_gaus_k
character*(128) :: filename_out_gaus_ijkl
character*(128) :: filename_out_s
character*(128) :: filename_out_v
character*(128) :: filename_out_k
character*(128) :: filename_out_s_ex
character*(128) :: filename_out_v_ex
character*(128) :: filename_out_k_ex
character*(128) :: filename_out_ijkl
character*(128) :: filename_basis

integer mpi_rank
logical MPI
integer ierr, mpi_size
double precision mpi_size_inv

integer*4                      :: seed(33)
integer*4                      :: put(33)
!!!!!!!!!!!COMMENT IF NOT MPI
!integer*4                      :: seed(12)
!!!!!!!!!!!COMMENT IF NOT MPI

! call getarg(1,filename_in)

MPI=.false.
mpi_rank=0

!!!!!!!!!!!COMMENT IF NOT MPI
! call mpi_init(ierr)
! call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_rank, ierr)
! call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_size, ierr)
! call sleep(mpi_rank/20)
! MPI=.true.
! write(*,*)'mpi_rank=',mpi_rank
!if(mpi_rank.eq.0)write(*,*)'mpi_size=',mpi_size
!!!!!!!!!!!COMMENT IF NOT MPI

write(filename_in,'(A4)') 'j_in'

if(mpi_rank.eq.0)write(*,*)'INPUT FILE USED=',filename_in
open(unit=5,file=filename_in)

if(mpi_rank.eq.0)print*,'Simulation number?'
read(5,'(a80)')charabid
read(5,*)num_simulation
if(mpi_rank.eq.0)write(*,*)'num_simulation=',num_simulation

if(mpi_rank.eq.0)then
 write(filename_out_gaus_s,'(A15)')'overlap_ao_gaus'
 write(filename_out_gaus_v,'(A15)')'nuclear_ao_gaus'
 write(filename_out_gaus_k,'(A15)')'kinetic_ao_gaus'
 write(filename_out_gaus_ijkl,'(A14)')'bielec_ao_gaus'
 write(filename_out_s,'(A10)')'overlap_ao'
 write(filename_out_v,'(A10)')'nuclear_ao'
 write(filename_out_k,'(A10)')'kinetic_ao'
 write(filename_out_s_ex,'(A13)')'overlap_ao_ex'
 write(filename_out_v_ex,'(A13)')'nuclear_ao_ex'
 write(filename_out_k_ex,'(A13)')'kinetic_ao_ex'
 write(filename_out_ijkl,'(A9)')'bielec_ao'
endif
 
!!!!!!!!!!!COMMENT IF NOT MPI
seed(1) = mpi_rank+1+(num_simulation-1)*10000
!!!!!!!!!!!COMMENT IF NOT MPI
do i_rand=2,12
 seed(i_rand) = i_rand
enddo
call random_seed(put=seed)

!*****************
! BEGIN READ INPUT
!*****************
if(mpi_rank.eq.0)print*,'MOLECULE?'
read(5,'(a80)')charabid
read(5,*)MOLECULE
if(mpi_rank.eq.0)write(*,*)MOLECULE
do i=1,3
 read(5,'(a80)')charabid
 if(mpi_rank.eq.0)write(6,'(a80)')charabid
enddo
read(5,*)basis_type
if(mpi_rank.eq.0)write(*,*)'basis_type=',basis_type
read(5,'(a80)')charabid
if(mpi_rank.eq.0)write(6,'(a80)')charabid
if(mpi_rank.eq.0)print*,'file name for basis set?'
read(5,*)filename_basis
if(mpi_rank.eq.0)print*,'*******'
if(mpi_rank.eq.0)write(*,*)trim(filename_basis)
if(mpi_rank.eq.0)print*,'*******'
do i=1,3
 read(5,'(a80)')charabid
 if(mpi_rank.eq.0)write(6,'(a80)')charabid
enddo
read(5,*)finite_range
if(mpi_rank.eq.0)write(*,*)'finite_range=',finite_range
read(5,'(a80)')charabid
if(mpi_rank.eq.0)write(6,'(a80)')charabid
read(5,*)n_eps
if(mpi_rank.eq.0)write(6,*)n_eps
if(mpi_rank.eq.0)then
 if(finite_range)write(*,*)'n determining place where orbitals vanish=',n_eps
endif
do i=1,5
 read(5,'(a80)')charabid
 if(mpi_rank.eq.0)write(6,'(a80)')charabid
enddo
read(5,*)ng0,i_add_2p,i_add_3d,g_thr,i_add_large_g
if(mpi_rank.eq.0)write(6,*)ng0,i_add_2p,i_add_3d,g_thr,i_add_large_g
do i=1,3
 read(5,'(a80)')charabid
 if(mpi_rank.eq.0)write(6,'(a80)')charabid
enddo
read(5,*)level
if(mpi_rank.eq.0)then
 if(finite_range)write(*,*)'level=',level
endif
do i=1,2
 read(5,'(a80)')charabid
 if(mpi_rank.eq.0)write(6,'(a80)')charabid
enddo
read(5,*)i_print_orb
if(mpi_rank.eq.0)write(6,*)i_print_orb
if(i_print_orb.ne.0)write(*,*)'ORBITALS ARE PRINTED'
do i=1,2
 read(5,'(a80)')charabid
 if(mpi_rank.eq.0)write(6,'(a80)')charabid
enddo
read(5,*)npts_one_elec
if(mpi_rank.eq.0)write(*,*)'Number of MC steps for one-electron integrals=',npts_one_elec, ' times number of blocks'
if(mod(npts_one_elec,nw).ne.0)stop 'npts_one_elec must be a multiple of nw'
read(5,'(a80)')charabid 
if(mpi_rank.eq.0)write(6,'(a80)')charabid 
read(5,*)npts_two_elec
if(mpi_rank.eq.0)write(*,*)'Number of MC steps for two-electron integrals=',npts_two_elec, ' times number of blocks'
if(mod(npts_two_elec,nw).ne.0)stop 'npts_two_elec must be a multiple of nw'
do i=1,2
 read(5,'(a80)')charabid
 if(mpi_rank.eq.0)write(6,'(a80)')charabid
enddo
read(5,*)aread,bread
if(mpi_rank.eq.0)print*,'aread=',aread
if(mpi_rank.eq.0)print*,'bread=',bread
if(mpi_rank.eq.0)print*,'Number of occupied orbitals?'
read(5,'(a80)')charabid
read(5,*)nocc
if(mpi_rank.eq.0)print*,'nocc=',nocc
if(mpi_rank.eq.0)print*,'Number iter for SCF calculations?'
read(5,'(a80)')charabid
read(5,*)niter_SCF
if(mpi_rank.eq.0)print*,'niter=',niter_SCF
!! Determination of one-center bielectronic
if(mpi_rank.eq.0)print*,'mono_center_Slater?'
read(5,'(a80)')charabid
read(5,*)mono_center_Slater
if(mpi_rank.eq.0)print*,'compute momo-center Slater=',mono_center_Slater
read(5,'(a80)')charabid
read(5,*)normalization_OA
if(mpi_rank.eq.0)print*,'orbital are normalized ?',normalization_OA

!*****************
! END READ INPUT
!*****************

call check_basis_type(filename_basis)
call read_basis(filename_basis)
if(finite_range.eqv..true.)call read_ng_star
if(mpi_rank.eq.0)print*,'READING geometry in angstrom'
call read_geometry(MOLECULE)
if(mpi_rank.eq.0)print*,'ENUCL=',enucl

call build_gaussian_expansion_of_orbitals(mpi_rank)
if(i_print_orb.eq.1)call print_orb

if(mpi_rank.eq.0)write(*,*)'**********************************'
if(mpi_rank.eq.0)write(*,*)'Number of basis functions ',nbasis

!! build arrays is(kcp),js(kcp),ks(kcp),ls(kcp) and kcp_ijkl(i,j,k,l) ;  kcp=1,nint and  i=1,nbasis (idem for j,k,l)
!! nint= total number of two-electron integrals
call build_mapping_ijkl

if(mpi_rank.eq.0)print*,'**************************************************************************'
if(mpi_rank.eq.0)print*,'Analytic calculation of APPROXIMATE GAUSSIAN one-electron integrals for ZV'
if(mpi_rank.eq.0)print*,'**************************************************************************'

call cpu_time(t0)
seuil_cauchy=1.d-15
nint_zero=0
i_value=1
n_zero_cauchy=0
k_sort=0
call cpu_time(t0)
kkk=0

num=nbasis**2
do i=1,nbasis
 do k=1,nbasis
  call one_elect(i,k,Sij,Vij,Kij)
  S_ij_gaus(i,k)=Sij
  V_ij_gaus(i,k)=Vij
  K_ij_gaus(i,k)=Kij
  if(finite_range)then
   i_value=int_ij_zero_rc(i,k)
   if(i_value.eq.0)nint_zero=nint_zero+1
  endif
  k_sort=k_sort+1
  if(num.ge.10.and.mod(k_sort,num/10).eq.0)then
   kkk=kkk+1
   call cpu_time(t1)
   if(mpi_rank.eq.0)write(*,'(a,i3,a,e22.15)')' CPU TIME block',kkk,' of GAUSSIAN ONE-electron',t1-t0
   t0=t1
   k_sort=0
  endif
 enddo
enddo

do i=1,nbasis
 sqrtSii_gaus(i)=dsqrt(dabs(S_ij_gaus(i,i)))
enddo

do i=1,nbasis
 do k=1,nbasis
  rnorm=sqrtSii_gaus(i)*sqrtSii_gaus(k)
  if(.not.normalization_OA)rnorm=1.d0
  S_ij_gaus(i,k)=S_ij_gaus(i,k)/rnorm
  V_ij_gaus(i,k)=V_ij_gaus(i,k)/rnorm
  K_ij_gaus(i,k)=K_ij_gaus(i,k)/rnorm
 enddo
enddo


if(mpi_rank.eq.0)then
 open(unit=11,file=filename_out_gaus_s)
 open(unit=12,file=filename_out_gaus_v)
 open(unit=13,file=filename_out_gaus_k)
 rewind 11
 rewind 12
 rewind 13
 do i=1,nbasis
  do k=1,nbasis
   write(11,'(2(I5,X),2D22.15)')i,k,S_ij_gaus(i,k)
   write(12,'(2(I5,X),D22.15)') i,k,V_ij_gaus(i,k)
   write(13,'(2(I5,X),D22.15)') i,k,K_ij_gaus(i,k)
  enddo
 enddo
 close(11)
 close(12)
 close(13)
endif

if(mpi_rank.eq.0)print*,'TOTAL NUMBER OF ONE-ELECTRON INTEGRALS =',nbasis**2
if(mpi_rank.eq.0)print*,'NUMBER OF ONE-ELECTRON INTEGRALS REMOVED BY RC=',nint_zero
call cpu_time(t1)
if(mpi_rank.eq.0)print*,'DONE!',' TIME=',t1-t0

if(mpi_rank.eq.0)print*,'********************************************************************'
if(mpi_rank.eq.0)print*,'ANALYTIC CALCULATION OF APPROXIMATE GAUSSIAN TWO-ELECTRON INTREGRALS'
if(mpi_rank.eq.0)print*,'********************************************************************'

if(mpi_rank.eq.0)call count_multi_center_integrals

do i=1,nbasis
 do k=1,nbasis
  precond(i,k)=gauss_ijkl(i,k,i,k)
 enddo
enddo

k_sort=0
k_sort2=0
nint_zero=0
i_value=1

n_zero_cauchy=0
nint_zero=0
seuil_cauchy=1.d-15

call cpu_time(t0)

if(mpi_rank.eq.0)then
 open(unit=14,file=filename_out_gaus_ijkl)
 rewind 14
endif

kkk=0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,l,i_value, k_sort, t1, rnorm) &
!$OMP  REDUCTION(+:n_zero_cauchy,nint_zero) SCHEDULE(dynamic,1000)
do kcp=nint,1,-1

 i=is(kcp)
 k=ks(kcp)
 j=js(kcp)
 l=ls(kcp)

 if(dsqrt(precond(i,k)*precond(j,l)).lt.seuil_cauchy)then
  n_zero_cauchy=n_zero_cauchy+1
  ijkl_gaus(kcp)=0.d0
 else
  if(finite_range)then
   i_value=int_ijkl_zero_rc(i,k,j,l)
   if(i_value.eq.0)nint_zero=nint_zero+1
  endif
  ijkl_gaus(kcp)=gauss_ijkl(i,k,j,l)
 endif

 k_sort=k_sort+1
 if( nint.ge.10.and.mod(k_sort,nint/10).eq.0)then
  !$OMP CRITICAL
  kkk=kkk+1
  call cpu_time(t1)
  if(mpi_rank.eq.0)write(*,'(a,i3,a,e22.15)')' CPU TIME block',kkk,' of GAUSSIAN TWO-electron',t1-t0
  t0=t1
  k_sort=0
  !$OMP END CRITICAL
 endif

! if(dabs(ijkl_gaus(kcp)).gt.-1.d-15)write(10,'(4(I5,X),2D22.15)') i,j,k,l, ijkl_gaus(kcp)

  rnorm=sqrtSii_gaus(i)*sqrtSii_gaus(k)*sqrtSii_gaus(j)*sqrtSii_gaus(l)
  if(.not.normalization_OA)rnorm=1.d0
  ijkl_gaus(kcp)=ijkl_gaus(kcp)/rnorm
  if(mpi_rank.eq.0)write(14,'(4(I5,X),D22.15)') i,j,k,l, ijkl_gaus(kcp)

enddo !kcp
!$OMP END PARALLEL DO

close(14)

if(mpi_rank.eq.0)print*,'*********************************************************************'
if(mpi_rank.eq.0)print*,'IJKL REMOVED DUE TO CAUCHY-SCHWARZ=',n_zero_cauchy,' OVER',nint,'IJKL'
if(finite_range)then
if(mpi_rank.eq.0)print*,'IJKL REMOVED IN ADDITION DUE TO RC=',nint_zero,' OVER',nint,'IJKL'
endif
if(mpi_rank.eq.0)print*,'TOTAL NUMBER OF IJKL REMOVED    =',nint_zero+n_zero_cauchy,' LEFT',nint-(nint_zero+n_zero_cauchy)
if(mpi_rank.eq.0)print*,'*******************************************************************************'


if(mpi_rank.eq.0)then
 print*
 print*,'***************************************************'
 print*,'SCF CALCULATION WITH APPROXIMATE GAUSSIAN INTEGRALS'
 print*,'***************************************************'

 call cpu_time(t0)
 iread_c0=0
 call SCF(nocc,S_ij_gaus,V_ij_gaus,K_ij_gaus,ijkl_gaus,iread_c0,ehf,niter_SCF)
 call cpu_time(t1)

!*********************************************'
if(mpi_rank.eq.0) print*,'BEGIN WRITING DATA_RED FOR FUTURE ZVMC'
 open(100,file='data_for_zvmc')
 rewind 100
 write(100,*)nint
 do i=1,nbasis
  do k=1,nbasis
   write(100,'(3e22.15)')S_ij_gaus(i,k),V_ij_gaus(i,k),K_ij_gaus(i,k)
  enddo
  enddo
 do kcp=1,nint
  write(100,'(e22.15)')ijkl_gaus(kcp)
 enddo
 close(100)
 print*,' DONE!'
if(mpi_rank.eq.0)then
 if( (basis_type.eq.'GTO').and.(finite_range.eqv..false.))then
  stop 'exact gaussian GTO calculation DONE'
 endif
endif

endif !mpi_rank=0

!!******************
!! BEGIN MONTE CARLO 
!!******************
!!
call cpu_time(t0)
if(mpi_rank.eq.0)then
 print*,'BEGIN READING OF DATA_FOR_ZVMC'
 open(100,file='data_for_zvmc')
 rewind 100
 read(100,*)nint
 do i=1,nbasis
  do k=1,nbasis
   read(100,'(3e22.15)')S_ij_gaus(i,k),V_ij_gaus(i,k),K_ij_gaus(i,k)
  enddo
  enddo
 do kcp=1,nint
  read(100,'(e22.15)')ijkl_gaus(kcp)
 enddo
 close(100)
 print*,'READING OF DATA_FOR_ZVMC DONE!'
endif
!**************************************
call build_mapping_ijkl

do i=1,nbasis
 do k=1,nbasis
  a_ZV(i,k)=aread
 enddo
enddo

alpha_ZV=1.d0
nbl=10

if(mpi_rank.eq.0)then
 print*,'*********************************'
 print*,'ZVMC EXACT ONE-ELECTRON INTEGRALS'
 print*,'*********************************'
 print*,' Number of Monte Carlo steps= ',nbl,' X ',npts_one_elec
 print*,'*******************************************************'
 if(mod(npts_one_elec,nw).ne.0)stop 'npts_one_elec must be a multiple of nw'
endif

one_elec_STO_exact=.true.
if( (basis_type.eq.'STO').and.(finite_range.eqv..false.).and.(one_elec_STO_exact.eqv..true.))then
 ng0_save=ng0
 ng0=20
 call build_gaussian_expansion_of_orbitals(mpi_rank)
 do i=1,nbasis
  do k=1,nbasis
   call one_elect(i,k,Sij,Vij,Kij)
   S_ij_STO_ex(i,k)=Sij
   V_ij_STO_ex(i,k)=Vij
   K_ij_STO_ex(i,k)=Kij
  enddo
  sqrtSii_STO_ex(i)=dsqrt(dabs(S_ij_STO_ex(i,i)))
 enddo

 do i=1,nbasis
  do k=1,nbasis
   rnorm=sqrtSii_STO_ex(i)*sqrtSii_STO_ex(k)
   if(.not.normalization_OA)rnorm=1.d0
   S_ij_STO_ex(i,k)=S_ij_STO_ex(i,k)/rnorm
   V_ij_STO_ex(i,k)=V_ij_STO_ex(i,k)/rnorm
   K_ij_STO_ex(i,k)=K_ij_STO_ex(i,k)/rnorm
  enddo
 enddo

 ng0=ng0_save
 call build_gaussian_expansion_of_orbitals(mpi_rank)

endif !one_elec_STO_exact

! if normalization_OA = true we de-normalize for the ZV step
if(normalization_OA)then
 do i=1,nbasis
  do k=i,nbasis
   rnorm=sqrtSii_gaus(i)*sqrtSii_gaus(k)
   S_ij_gaus(i,k)=S_ij_gaus(i,k)*rnorm
   V_ij_gaus(i,k)=V_ij_gaus(i,k)*rnorm
   K_ij_gaus(i,k)=K_ij_gaus(i,k)*rnorm
  enddo
 enddo
endif

do i=1,nbasis
 do k=1,nbasis
 S_ij_ZV(i,k)=0.d0
 S_ij_ZV2(i,k)=0.d0
 V_ij_ZV(i,k)=0.d0
 V_ij_ZV2(i,k)=0.d0
 K_ij_ZV(i,k)=0.d0
 K_ij_ZV2(i,k)=0.d0
 enddo
enddo

do kkk=1,nbl

 call cpu_time(t0)

 do i=1,nbasis
  do k=1,nbasis
   e_Sij(i,k)=0.d0
   e_Vij(i,k)=0.d0
   e_Kij(i,k)=0.d0
   e_G_Sij(i,k)=0.d0
   e_G_Vij(i,k)=0.d0
   e_G_Kij(i,k)=0.d0
  enddo
 enddo

 do kk=1,npts_one_elec/nw

  call draw_configuration(1,r1,r2)
  call compute_pi0(1,r1,r2,pi_0)

  do i=1,nbasis
   do k=i,nbasis
    i_value=1
    if(finite_range)i_value=int_ij_zero_rc(i,k)
    if(i_value.ne.0)then
     call compute_r_tilde(1,i,k,r1,r2,rt1,rt2)
     call compute_u_tilde(1,i,k,rt1,rt2,ut1,ut2)
     call compute_densities(1,i,k,ut1,ut2,rho,rho_G,poly)

     ik = (k-1)*nbasis_max+i
     call compute_jacobian(1,i,k,j,l,r1,r2,rjacob)
     do kw=1,nw
      factor=rjacob(kw)/pi_0(kw)/(g_min(i)+g_min(k))**1.5d0
      weight  (kw)    =factor* rho  (kw,ik,1)
      weight_kin(kw)  =factor* poly (kw,ik,1)
      weight_G(kw)    =factor* rho_G(kw,ik,1)
      weight_kin_G(kw)=factor* poly (kw,ik,1)
     enddo
     call compute_pot_nuclear(i,k,ut1,pot)
     call compute_kinetic(i,k,ut1,rkin,rkin_G)
     e_Sij(i,k)=e_Sij(i,k)     +     sum(weight  (1:nw))
     e_Vij(i,k)=e_Vij(i,k)     +     sum(pot(1:nw)*weight  (1:nw))
     e_Kij(i,k)=e_Kij(i,k)     +     sum(rkin(1:nw)  *weight_kin  (1:nw))
     e_G_Sij(i,k)=e_G_Sij(i,k) +     sum(weight_G(1:nw))
     e_G_Vij(i,k)=e_G_Vij(i,k) +     sum(pot(1:nw)*weight_G(1:nw))
     e_G_Kij(i,k)=e_G_Kij(i,k) +     sum(rkin_G(1:nw)*weight_kin_G(1:nw))
    endif !i_value
   enddo !k
  enddo !i

 enddo !npts_one_elec

 do i=1,nbasis
  do k=i,nbasis
   e_Sij(i,k)=e_Sij(i,k)/npts_one_elec
   e_Vij(i,k)=e_Vij(i,k)/npts_one_elec
   e_Kij(i,k)=e_Kij(i,k)/npts_one_elec
   e_G_Sij(i,k)=e_G_Sij(i,k)/npts_one_elec
   e_G_Vij(i,k)=e_G_Vij(i,k)/npts_one_elec
   e_G_Kij(i,k)=e_G_Kij(i,k)/npts_one_elec

   if(finite_range)i_value=int_ij_zero_rc(i,k)
    if(i_value.ne.0)then
     e_tot_Sij=e_Sij(i,k)-e_G_Sij(i,k)+S_ij_gaus(i,k)
     e_tot_Vij=e_Vij(i,k)-e_G_Vij(i,k)+V_ij_gaus(i,k)
     e_tot_Kij=e_Kij(i,k)-e_G_Kij(i,k)+K_ij_gaus(i,k)
    else
     e_tot_Sij=0.d0
     e_tot_Vij=0.d0
     e_tot_Kij=0.d0
    endif
    S_ij_ZV(i,k) =S_ij_ZV(i,k) +e_tot_Sij
    S_ij_ZV2(i,k)=S_ij_ZV2(i,k)+e_tot_Sij**2
    V_ij_ZV(i,k) =V_ij_ZV(i,k) +e_tot_Vij
    V_ij_ZV2(i,k)=V_ij_ZV2(i,k)+e_tot_Vij**2
    K_ij_ZV(i,k) =K_ij_ZV(i,k) +e_tot_Kij
    K_ij_ZV2(i,k)=K_ij_ZV2(i,k)+e_tot_Kij**2
  enddo
 enddo

 do i=1,nbasis
  do k=i,nbasis
   e_Sij(k,i)=e_Sij(i,k)
   e_Vij(k,i)=e_Vij(i,k)
   e_Kij(k,i)=e_Kij(i,k)
   e_G_Sij(k,i)=e_G_Sij(i,k)
   e_G_Vij(k,i)=e_G_Vij(i,k)
   e_G_Kij(k,i)=e_G_Kij(i,k)
   S_ij_ZV(k,i) =S_ij_ZV(i,k) 
   S_ij_ZV2(k,i)=S_ij_ZV2(i,k)
   V_ij_ZV(k,i) =V_ij_ZV(i,k) 
   V_ij_ZV2(k,i)=V_ij_ZV2(i,k)
   K_ij_ZV(k,i) =K_ij_ZV(i,k) 
   K_ij_ZV2(k,i)=K_ij_ZV2(i,k)
  enddo
 enddo

 call cpu_time(t1)
 if(mpi_rank.eq.0)write(*,'(a,i3,a,e22.15)')' CPU TIME block',kkk,' of EXACT ZV ONE-electron',t1-t0

enddo !kkk=1,nbl

do i=1,nbasis
 do k=1,nbasis
  S_ij_ZV (i,k) =S_ij_ZV(i,k)/nbl
  S_ij_ZV2(i,k)= S_ij_ZV2(i,k)/nbl
  V_ij_ZV (i,k) =V_ij_ZV(i,k)/nbl
  V_ij_ZV2(i,k)= V_ij_ZV2(i,k)/nbl
  K_ij_ZV (i,k) =K_ij_ZV(i,k)/nbl
  K_ij_ZV2(i,k)= K_ij_ZV2(i,k)/nbl
  error_Sij(i,k)=dsqrt( dabs(S_ij_ZV2(i,k)-S_ij_ZV(i,k)**2))/dsqrt(dfloat(nbl))
  error_Vij(i,k)=dsqrt( dabs(V_ij_ZV2(i,k)-V_ij_ZV(i,k)**2))/dsqrt(dfloat(nbl))
  error_Kij(i,k)=dsqrt( dabs(K_ij_ZV2(i,k)-K_ij_ZV(i,k)**2))/dsqrt(dfloat(nbl))
 enddo
enddo

!!!!!!!!!!!COMMENT IF NOT MPI
! mpi_size_inv = 1.d0/dble(mpi_size)
! call MPI_AllReduce(S_ij_ZV, moy_S, nbasis_max*nbasis, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
! moy_S(:,:)=moy_S(:,:) * mpi_size_inv
! call MPI_AllReduce(V_ij_ZV, moy_V, nbasis_max*nbasis, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
! moy_V(:,:)=moy_V(:,:) * mpi_size_inv
! call MPI_AllReduce(K_ij_ZV, moy_K, nbasis_max*nbasis, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
! moy_K(:,:)=moy_K(:,:) * mpi_size_inv

! moy_S2t(:,:) = (S_ij_ZV(:,:) -  moy_S(:,:))**2
! call MPI_reduce(moy_S2t, moy_S2, nbasis_max*nbasis, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
! moy_S2(:,:)=dsqrt(moy_S2(:,:) * mpi_size_inv/ (dble(mpi_size-1)) )

! moy_V2t(:,:) = (V_ij_ZV(:,:) -  moy_V(:,:))**2
! call MPI_reduce(moy_V2t, moy_V2, nbasis_max*nbasis, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
! moy_V2(:,:)=dsqrt(moy_V2(:,:) * mpi_size_inv/ (dble(mpi_size-1)) )

! moy_K2t(:,:) = (K_ij_ZV(:,:) -  moy_K(:,:))**2
! call MPI_reduce(moy_K2t, moy_K2, nbasis_max*nbasis, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
! moy_K2(:,:)=dsqrt(moy_K2(:,:) * mpi_size_inv/ (dble(mpi_size-1)) )
!!!!!!!!!!!COMMENT IF NOT MPI

if(MPI.eqv..false.)then
 moy_S(:,:)=S_ij_ZV(:,:)
 moy_S2(:,:)=error_Sij(:,:)
 moy_V(:,:)=V_ij_ZV(:,:)
 moy_V2(:,:)=error_Vij(:,:)
 moy_K(:,:)=K_ij_ZV(:,:)
 moy_K2(:,:)=error_Kij(:,:)
endif

if (mpi_rank == 0) then

 open(unit=101,file=filename_out_s)
 open(unit=102,file=filename_out_v)
 open(unit=103,file=filename_out_k)

 open(unit=104,file=filename_out_s_ex)
 open(unit=105,file=filename_out_v_ex)
 open(unit=106,file=filename_out_k_ex)

 error_mean_S=0.d0
 error_max_S=-1.d0
 error_mean_V=0.d0
 error_max_V=-1.d0
 error_mean_K=0.d0
 error_max_K=-1.d0

 if(normalization_OA)then
  do i=1,nbasis
   do k=1,nbasis
    rnorm=sqrtSii_STO_ex(i)*sqrtSii_STO_ex(k)
    moy_S(i,k)= moy_S(i,k)/rnorm
    moy_S2(i,k)= moy_S2(i,k)/rnorm
    moy_V(i,k)= moy_V(i,k)/rnorm
    moy_V2(i,k)= moy_V2(i,k)/rnorm
    moy_K(i,k)= moy_K(i,k)/rnorm
    moy_K2(i,k)= moy_K2(i,k)/rnorm
   enddo
  enddo
 endif

 do i=1,nbasis
  do k=1,nbasis
   write(101,'(2(I5,X),2D22.15)') i,k,moy_S(i,k), moy_S2(i,k)
   write(104,'(2(I5,X),D22.15)') i,k,S_ij_STO_ex(i,k)
   write(102,'(2(I5,X),2D22.15)') i,k,moy_V(i,k), moy_V2(i,k)
   write(105,'(2(I5,X),D22.15)') i,k,V_ij_STO_ex(i,k)
   write(103,'(2(I5,X),2D22.15)') i,k,moy_K(i,k), moy_K2(i,k)
   write(106,'(2(I5,X),D22.15)') i,k,K_ij_STO_ex(i,k)
   error_mean_S=error_mean_S+dabs(moy_S2(i,k))
   error_mean_V=error_mean_V+dabs(moy_V2(i,k))
   error_mean_K=error_mean_K+dabs(moy_K2(i,k))
   if(moy_S2(i,k).gt.error_max_S)error_max_S=moy_S2(i,k)
   if(moy_V2(i,k).gt.error_max_V)error_max_V=moy_V2(i,k)
   if(moy_K2(i,k).gt.error_max_K)error_max_K=moy_K2(i,k)
  enddo
 enddo
 close(101)
 close(102)
 close(103)
 close(104)
 close(105)
 close(106)

 print*,'Error_mean S_ij ZV=',error_mean_S/(nbasis**2),' error_max ',error_max_S
 print*,'Error_mean V_ij ZV=',error_mean_V/(nbasis**2),' error_max ',error_max_V
 print*,'Error_mean K_ij ZV=',error_mean_K/(nbasis**2),' error_max ',error_max_K
endif !mpi_rank

!! Determination of one-center bielectronic
do kcp=1,nint
 i=is(kcp)
 k=ks(kcp)
 j=js(kcp)
 l=ls(kcp)
 call compare_nuclei(nucleus_number(i),nucleus_number(j),nucleus_number(k),nucleus_number(l),ndiff)
 mono_center(kcp)=0
 if(mono_center_Slater)then 
  if(ndiff.eq.1)mono_center(kcp)=1
 endif
enddo

if(mpi_rank.eq.0)then
 print*,'*********************************'
 print*,'ZVMC EXACT TWO-ELECTRON INTEGRALS'
 print*,'*********************************'
 print*,' Number of Monte Carlo steps= ',nbl,' X ',npts_two_elec
 print*,'*******************************************************'
 if(mod(npts_two_elec,nw).ne.0)stop 'npts_two_elec must be a multiple of nw'
endif

do i=1,nbasis
 do j=1,nbasis
  do ll=1,3
   G_center(ll,i,j)= (g_min(i)*center(ll,i)+g_min(j)*center(ll,j))/(g_min(i)+g_min(j))
  enddo
 enddo
enddo

!**************************************************
!! Exact calculation of monocenter Slater integrals
do kcp=1,nint
 i=is(kcp)
 k=ks(kcp)
 j=js(kcp)
 l=ls(kcp)
 if(mono_center(kcp).eq.1)then
  call compute_int_slater(i,k,j,l,ijkl(kcp))
  if(normalization_OA)then
   rnorm=sqrtSii_STO_ex(i)*sqrtSii_STO_ex(k)*sqrtSii_STO_ex(j)*sqrtSii_STO_ex(l)
   ijkl(kcp)=ijkl(kcp)/rnorm
  endif 
 endif
enddo
!! end calculation ***********************************

! if normalization_OA = true we de-normalize for the ZV step
 if(normalization_OA)then
  do kcp=1,nint
   i=is(kcp)
   k=ks(kcp)
   j=js(kcp)
   l=ls(kcp)
   rnorm=sqrtSii_gaus(i)*sqrtSii_gaus(k)*sqrtSii_gaus(j)*sqrtSii_gaus(l)
   ijkl_gaus(kcp)=ijkl_gaus(kcp)*rnorm
  enddo
 endif

do kcp=1,nint
 if(mono_center(kcp).eq.0)then
  ijkl(kcp)=0.d0
  ijkl2(kcp)=0.d0
 endif
enddo

!! Determinate which i j are used in computation of W_S and W_G
!! i_tab_mc(i,k)=1  W_S and W_G can be computed
do k=1,nbasis
 do i=k,nbasis
  i_tab_mc(i,k)=0
 enddo
enddo
do kcp=1,nint
 i=is(kcp)
 k=ks(kcp)
 j=js(kcp)
 l=ls(kcp)
 i_tab_mc(i,k)=1
 i_tab_mc(j,l)=1
enddo
!! end 

do kkk=1,nbl

 call cpu_time(t0)

 do kcp=1,nint
  e_S_ijkl(kcp)=0.d0
  e_G_ijkl(kcp)=0.d0
  e_test_ijkl(kcp)=0.d0
 enddo

 k_sort=0
 k_sort2=0

 do kk=1,npts_two_elec/nw

  call draw_configuration(2,r1,r2)
  call compute_pi0(2,r1,r2,pi_0)

  do k=1,nbasis
   do i=k,nbasis
    if(i_tab_mc(i,k).eq.1)then
     call compute_r_tilde(2,i,k,r1,r2,rt1,rt2)
     call compute_u_tilde(2,i,k,rt1,rt2,ut1,ut2)
     call compute_densities(2,i,k,ut1,ut2,rho,rho_G,poly)
    endif
   enddo !k
  enddo !i

  do kcp=1,nint

   i=is(kcp)
   k=ks(kcp)
   ik = (k-1)*nbasis_max+i
   j=js(kcp)
   l=ls(kcp)
   jl = (l-1)*nbasis_max+j

   if(dsqrt(precond(i,k)*precond(j,l)).lt.seuil_cauchy)then

    n_zero_cauchy=n_zero_cauchy+1

   else

    if(finite_range)then
     i_value=int_ijkl_zero_rc(i,k,j,l)
     if(i_value.eq.0)then
      nint_zero=nint_zero+1
     else
      call compute_jacobian(2,i,k,j,l,r1,r2,rjacob)
      do kw=1,nw
       d_x(1) = ut1(1,kw,ik)-ut2(1,kw,jl)
       d_x(2) = ut1(2,kw,ik)-ut2(2,kw,jl)
       d_x(3) = ut1(3,kw,ik)-ut2(3,kw,jl)
       r12_2 = d_x(1)*d_x(1) + d_x(2)*d_x(2) + d_x(3)*d_x(3)
       factor=rjacob(kw)/pi_0(kw)/( (g_min(i)+g_min(k))*(g_min(j)+g_min(l)) )**1.5d0
       weight  (kw)=factor* rho  (kw,ik,1)*rho  (kw,jl,2)
       weight_G(kw)=factor* rho_G(kw,ik,1)*rho_G(kw,jl,2)

       r12_inv(kw) = real( 1./sqrt(real(r12_2,4)), 8) !simple precision
!      r12_inv(kw) = 1.d0/dsqrt(r12_2)                !double precision
      enddo
      e_S_ijkl(kcp)=e_S_ijkl(kcp) + sum(weight  (1:nw)*r12_inv(1:nw))
      e_G_ijkl(kcp)=e_G_ijkl(kcp) + sum(weight_G(1:nw)*r12_inv(1:nw))
      e_test_ijkl(kcp)=e_test_ijkl(kcp) + sum(r12_inv(1:nw))

     endif
    else
     call compute_jacobian(2,i,k,j,l,r1,r2,rjacob)
     do kw=1,nw
      d_x(1) = ut1(1,kw,ik)-ut2(1,kw,jl)
      d_x(2) = ut1(2,kw,ik)-ut2(2,kw,jl)
      d_x(3) = ut1(3,kw,ik)-ut2(3,kw,jl)
      r12_2 = d_x(1)*d_x(1) + d_x(2)*d_x(2) + d_x(3)*d_x(3)
      factor=rjacob(kw)/pi_0(kw)/( (g_min(i)+g_min(k))*(g_min(j)+g_min(l)) )**1.5d0
      weight  (kw)=factor* rho  (kw,ik,1)*rho  (kw,jl,2)
      weight_G(kw)=factor* rho_G(kw,ik,1)*rho_G(kw,jl,2)
      r12_inv(kw) = real( 1./sqrt(real(r12_2,4)), 8) !simple precision
!      r12_inv(kw) = 1.d0/dsqrt(r12_2)                !double precision
     enddo
     e_S_ijkl(kcp)=e_S_ijkl(kcp) + sum(weight  (1:nw)*r12_inv(1:nw))
     e_G_ijkl(kcp)=e_G_ijkl(kcp) + sum(weight_G(1:nw)*r12_inv(1:nw))
     e_test_ijkl(kcp)=e_test_ijkl(kcp) + sum(r12_inv(1:nw))
    endif
   endif
  enddo !kcp

  k_sort=k_sort+1
  k_sort2=k_sort2+nw
  if( npts_two_elec.ge.10.and.mod(k_sort,npts_two_elec/10).eq.0)then
   write(*,*)'mpi_rank= ',mpi_rank,' nsteps=',k_sort2
   k_sort=0
  endif
 enddo !npts_two_elec

 do kcp=1,nint
  if(mono_center(kcp).eq.0)then
   e_S_ijkl(kcp)=e_S_ijkl(kcp)/npts_two_elec
   e_G_ijkl(kcp)=e_G_ijkl(kcp)/npts_two_elec
   e_test_ijkl(kcp)=e_test_ijkl(kcp)/npts_two_elec
  endif
 enddo

 do kcp=1,nint
  i_value=1
  i=is(kcp)
  k=ks(kcp)
  j=js(kcp)
  l=ls(kcp)
  if(finite_range)i_value=int_ijkl_zero_rc(i,k,j,l)
  if(i_value.ne.0)then
   if(mono_center(kcp).eq.0)then
    e_tot_ijkl=e_S_ijkl(kcp)-e_G_ijkl(kcp)+ijkl_gaus(kcp)
    ijkl(kcp)=ijkl(kcp)+e_tot_ijkl
    ijkl2(kcp)=ijkl2(kcp)+e_tot_ijkl**2
   endif 
  endif 
 enddo !kcp

 call cpu_time(t1)
 if(mpi_rank.eq.0)write(*,'(a,i3,a,e22.15)')' CPU TIME block',kkk,' of EXACT ZV TWO-electron',t1-t0

enddo !kkk=1,nbl

do kcp=1,nint
 if(mono_center(kcp).eq.0)then
  ijkl(kcp)=ijkl(kcp)/nbl
  ijkl2(kcp)=ijkl2(kcp)/nbl
  error_ijkl(kcp)=dsqrt( dabs(ijkl2(kcp)-ijkl(kcp)**2))/dsqrt(dfloat(nbl))
 else
  error_ijkl(kcp)=0.d0
 endif
enddo

!!!!!!!!!!!COMMENT IF NOT MPI
! mpi_size_inv = 1.d0/dble(mpi_size)
! call MPI_AllReduce(ijkl, moy, nint, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
! moy(:)=moy(:) * mpi_size_inv

! moy2t(:) = (ijkl(:) -  moy(:))**2
! call MPI_reduce(moy2t, moy2, nint, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
! moy2(:)=dsqrt(moy2(:) * mpi_size_inv/ (dble(mpi_size-1)) )
!!!!!!!!!!!COMMENT IF NOT MPI

 if(MPI.eqv..false.)then
  moy(:)=ijkl(:)
  moy2(:)=error_ijkl(:)
  if(normalization_OA)then
   do kcp=1,nint
    if(mono_center(kcp).eq.0)then
     i=is(kcp)
     k=ks(kcp)
     j=js(kcp)
     l=ls(kcp)
     rnorm=sqrtSii_STO_ex(i)*sqrtSii_STO_ex(k)*sqrtSii_STO_ex(j)*sqrtSii_STO_ex(l)
     moy(kcp)=moy(kcp)/rnorm
     moy2(kcp)=moy2(kcp)/rnorm
    endif
   enddo
  endif
 endif

 if (mpi_rank == 0) then

  open(unit=104,file=filename_out_ijkl)

  error_max=-1.d0
  errmoy_ijkl=0.d0
  dmoy_ijkl=0.d0
  diff_ijklmax=-1.d0
  a_ijkln=0.d0
  n_ijkl=0

  do kcp=1,nint
   i=is(kcp)
   k=ks(kcp)
   j=js(kcp)
   l=ls(kcp)
   write(104,'(4(I5,X),2D22.15)') i,j,k,l, moy(kcp), moy2(kcp)
   if(moy2(kcp).gt.error_max)error_max=moy2(kcp)

   if( dabs(moy(kcp)).gt.1.d-6.and.(moy2(kcp).ne.0.d0))then
    errmoy_ijkl=errmoy_ijkl+moy2(kcp)
    dmoy_ijkl=dmoy_ijkl+1.d0
    diff_ijkl=dabs(moy(kcp)-ijkl_gaus(kcp))/moy2(kcp)
    if(diff_ijkl.gt.diff_ijklmax)diff_ijklmax=diff_ijkl
    if(diff_ijkl.gt.3.d0)then
     a_ijkln=a_ijkln+diff_ijkl
     n_ijkl=n_ijkl+1
    endif
   endif

  enddo

  close(104)

  write(*,*)'Error mean ijkl=',errmoy_ijkl/dmoy_ijkl,' error_max ',error_max
  if(idebug.eq.1)then
   write(*,*)'IJKL'
   write(*,*)'****'
   write(*,*)'Max DELTA=[ |ijkl_ZV - ijkl_gaus|/error] =',diff_ijklmax
   write(*,'(a,i8,a,i6,a)')' Numb of integrals such that DELTA > 3 ',n_ijkl,' over',nint,' two-elec integrals'
  endif

 endif

 print*,'ZVMC DONE for mpi_rank=',mpi_rank

if(mpi_rank.eq.0)then

 print*
 print*,'***************************************'
 print*,'SCF CALCULATION WITH EXACT ZV INTEGRALS'
 print*,'***************************************'
!! EXACT one-electron integrals are needed for the SCF calculation 

 call cpu_time(t0)
 iread_c0=0
 call SCF(nocc,moy_S,moy_V,moy_K,moy,iread_c0,ehf,niter_SCF)
 call cpu_time(t1)

 print*,' SCF Done! TIME=',t1-t0

endif

!!!!!!!!!!!COMMENT IF NOT MPI
! call MPI_BARRIER (MPI_COMM_WORLD, ierr)
! call mpi_finalize(ierr)
!!!!!!!!!!!COMMENT IF NOT MPI

end

integer function int_ijkl_zero_rc(i,k,j,l)
include 'j.inc'
int_ijkl_zero_rc=1
dik=r_c(i)+r_c(k)
djl=r_c(j)+r_c(l)

if(dik.lt.dist_ij(i,k))then
 int_ijkl_zero_rc=0
 return
endif

if(djl.lt.dist_ij(j,l))then
 int_ijkl_zero_rc=0
 return
endif

end

integer function int_ij_zero_rc(i,k)
include 'j.inc'
int_ij_zero_rc=1
dik=r_c(i)+r_c(k)
if(dik.lt.dist_ij(i,k))then
 int_ij_zero_rc=0
 return
endif
end

double precision function phi_orb(i_orb,x)
include 'j.inc'
dimension x(3)
dx=x(1)-center(1,i_orb)
dy=x(2)-center(2,i_orb)
dz=x(3)-center(3,i_orb)
r=dsqrt(dx**2+dy**2+dz**2)
phi_orb=dx**npower(1,i_orb)*dy**npower(2,i_orb)*dz**npower(3,i_orb)* u_orb(i_orb,r)
end

double precision function phi_gaus(i_orb,x)
include 'j.inc'
dimension x(3)
dx=x(1)-center(1,i_orb)
dy=x(2)-center(2,i_orb)
dz=x(3)-center(3,i_orb)
r=dsqrt(dx**2+dy**2+dz**2)
phi_gaus=dx**npower(1,i_orb)*dy**npower(2,i_orb)*dz**npower(3,i_orb) * u_gauss(i_orb,r)
end

subroutine read_ng_star
include 'j.inc'
if(level.gt.4)stop 'level greater than 4 not coded'

if(level.eq.1)then
 ng_star(1)=5
 ng_star(2)=3
 ng_star(3)=3
 ng_star(4)=2
 do l=5,16
  ng_star(l)=1
 enddo
endif
if(level.eq.2)then
 ng_star(1)=5
 ng_star(2)=4
 ng_star(3)=4
 ng_star(4)=2
 do l=5,16
  ng_star(l)=1
 enddo
endif
if(level.eq.3)then
 ng_star(1)=5
 ng_star(2)=4
 ng_star(3)=4
 ng_star(4)=3
 ng_star(5)=2
 do l=6,16
  ng_star(l)=1
 enddo
endif
if(level.eq.4)then
 ng_star(1)=5
 ng_star(2)=4
 ng_star(3)=4
 ng_star(4)=3
 ng_star(5)=2
 ng_star(6)=2
 do l=7,16
  ng_star(l)=1
 enddo
endif
end

subroutine compute_int_slater(i,k,j,l,ijkl_gaus)
include 'j.inc'
dimension gam(4)
integer n_i(3),n_k(3),n_j(3),n_l(3),nix(4),niy(4),niz(4),power(4)
double precision ijkl_slater_xnynzn
double precision ijkl_gaus
character*80 orb

gam(1)=g_slater(i)
gam(2)=g_slater(k)
gam(3)=g_slater(j)
gam(4)=g_slater(l)

nix(1)=npower(1,i)
nix(2)=npower(1,k)
nix(3)=npower(1,j)
nix(4)=npower(1,l)

niy(1)=npower(2,i)
niy(2)=npower(2,k)
niy(3)=npower(2,j)
niy(4)=npower(2,l)

niz(1)=npower(3,i)
niz(2)=npower(3,k)
niz(3)=npower(3,j)
niz(4)=npower(3,l)

orb=orb_name(i)
if(orb.ne.'1S'.and.orb.ne.'2S'.and.orb.ne.'3S'.and.orb.ne.'2P'.and.orb.ne.'3P'.and.orb.ne.'3D')then
print*,'orb=',orb
stop 'pb in compute_in_slater'
endif

power(1)=0
if(orb.eq.'2S')power(1)=1
if(orb.eq.'3S')power(1)=2
if(orb.eq.'3P')power(1)=1

orb=orb_name(k)
power(2)=0
if(orb.eq.'2S')power(2)=1
if(orb.eq.'3S')power(2)=2
if(orb.eq.'3P')power(2)=1

orb=orb_name(j)
power(3)=0
if(orb.eq.'2S')power(3)=1
if(orb.eq.'3S')power(3)=2
if(orb.eq.'3P')power(3)=1

orb=orb_name(l)
power(4)=0
if(orb.eq.'2S')power(4)=1
if(orb.eq.'3S')power(4)=2
if(orb.eq.'3P')power(4)=1

ijkl_gaus=ijkl_slater_xnynzn(nix,niy,niz,power,gam)
end

subroutine draw_configuration(i_el,r1,r2)
include 'j.inc'
dimension r1(nw,3),r2(nw,3)
do kw=1,nw
 do ll=1,3
  call random_number(ranf1)
  r1(kw,ll)=dsqrt(2.d0)*dierfc(2.d0-2.d0*ranf1)
 enddo
enddo
if(i_el.eq.1)return
do kw=1,nw
 do ll=1,3
  call random_number(ranf2)
  r2(kw,ll)=dsqrt(2.d0)*dierfc(2.d0-2.d0*ranf2)
 enddo
enddo
end

subroutine compute_pi0(i_el,r1,r2,pi_0)
include 'j.inc'
dimension r1(nw,3),r2(nw,3),pi_0(nw)
double precision, parameter :: pi=dacos(-1.d0)
double precision, parameter :: f = 1.d0/dsqrt(2.d0*pi)**3
do kw=1,nw
 r1_mod_2=r1(kw,1)**2+r1(kw,2)**2+r1(kw,3)**2
 pi_0(kw)=dexp(-0.5d0*r1_mod_2)*f
enddo
if(i_el.eq.1)return
do kw=1,nw
 r2_mod_2=r2(kw,1)**2+r2(kw,2)**2+r2(kw,3)**2
 pi_0(kw)=pi_0(kw)*dexp(-0.5d0*r2_mod_2)*f
enddo
end

subroutine compute_densities(i_el,i,k,ut1,ut2,rho,rho_G,poly)
include 'j.inc'
dimension ut1(3,nw,nbasis_max*nbasis_max)
dimension ut2(3,nw,nbasis_max*nbasis_max)
dimension rho  (nw,nbasis_max*nbasis_max,2)
dimension poly  (nw,nbasis_max*nbasis_max,2)
dimension rho_G(nw,nbasis_max*nbasis_max,2)
ik = (k-1)*nbasis_max+i
do kw=1,nw
 dx=ut1(1,kw,ik)-center(1,i)
 dy=ut1(2,kw,ik)-center(2,i)
 dz=ut1(3,kw,ik)-center(3,i)
 poly_i=dx**npower(1,i)*dy**npower(2,i)*dz**npower(3,i)
 r_i=dsqrt(dx**2+dy**2+dz**2)

 dx=ut1(1,kw,ik)-center(1,k)
 dy=ut1(2,kw,ik)-center(2,k)
 dz=ut1(3,kw,ik)-center(3,k)
 poly_k=dx**npower(1,k)*dy**npower(2,k)*dz**npower(3,k)
 r_k=dsqrt(dx**2+dy**2+dz**2)

 rho  (kw,ik,1)=poly_i*u_orb(i,r_i)*poly_k*u_orb(k,r_k)
 poly (kw,ik,1)=poly_i*poly_k
 rho_G(kw,ik,1)=poly_i*u_gauss(i,r_i)*poly_k*u_gauss(k,r_k)

enddo
if(i_el.eq.1)return
do kw=1,nw
 rho  (kw,ik,2)=phi_orb (i,ut2(1,kw,ik))*phi_orb (k,ut2(1,kw,ik))
 rho_G(kw,ik,2)=phi_gaus(i,ut2(1,kw,ik))*phi_gaus(k,ut2(1,kw,ik))
enddo
end

subroutine compute_r_tilde(i_el,i,k,r1,r2,rt1,rt2)
include 'j.inc'
dimension r1(nw,3),r2(nw,3)
dimension rt1(nw,3),rt2(nw,3)
do kw=1,nw
 r1_mod=dsqrt(r1(kw,1)**2+r1(kw,2)**2+r1(kw,3)**2)
 call compute_f_fp_jacob(r1_mod,i,k,0,f1,f1p)
 do ll=1,3
  rt1(kw,ll)=f1*r1(kw,ll)
 enddo !!ll
enddo !!kw
if(i_el.eq.1)return
do kw=1,nw
 r2_mod=dsqrt(r2(kw,1)**2+r2(kw,2)**2+r2(kw,3)**2)
 call compute_f_fp_jacob(r2_mod,i,k,0,f2,f2p)
 do ll=1,3
  rt2(kw,ll)=f2*r2(kw,ll)
 enddo !!ll
enddo !!kw
end

subroutine compute_pot_nuclear(i,k,ut1,pot)
include 'j.inc'
dimension ut1(3,nw,nbasis_max*nbasis_max),pot(nw)
ik = (k-1)*nbasis_max+i
do kw=1,nw
 pot(kw)=0.d0
 do kn=1,number_nuclei
  dR=dsqrt( (ut1(1,kw,ik)-centers_nuclei(1,kn))**2+ &
  (ut1(2,kw,ik)-centers_nuclei(2,kn))**2+ &
  (ut1(3,kw,ik)-centers_nuclei(3,kn))**2 )
  pot(kw)=pot(kw)-charge(kn)/dR
enddo
enddo
end

subroutine compute_kinetic(i,k,ut1,rkin,rkin_G)
include 'j.inc'
dimension ut1(3,nw,nbasis_max*nbasis_max),rkin(nw),rkin_G(nw)
dimension dr_i(3),dr_k(3)
ik = (k-1)*nbasis_max+i
do kw=1,nw

 do l=1,3
  dr_i(l)=ut1(l,kw,ik)-center(l,i)
  dr_k(l)=ut1(l,kw,ik)-center(l,k)
 enddo
  r_i=dsqrt(dr_i(1)**2+ dr_i(2)**2+dr_i(3)**2)
  r_k=dsqrt(dr_k(1)**2+ dr_k(2)**2+dr_k(3)**2)

 u_i   = u_orb   (i,r_i)
 u_G_i = u_gauss (i,r_i)

 up_i  = up_orb  (i,r_i)
 up_G_i= up_gauss(i,r_i)

 u_k   = u_orb   (k,r_k)
 u_G_k = u_gauss (k,r_k)

 up_k  = up_orb  (k,r_k)
 up_G_k= up_gauss(k,r_k)

 rkin(kw)=0.d0
 rkin_G(kw)=0.d0
 do l=1,3
  rkin(kw)=rkin(kw)    +0.5d0*(npower(l,i)/dr_i(l)*u_i  +dr_i(l)*up_i/r_i)  *(npower(l,k)/dr_k(l)*u_k  +dr_k(l)*up_k/r_k)
  rkin_G(kw)=rkin_G(kw)+0.5d0*(npower(l,i)/dr_i(l)*u_G_i+dr_i(l)*up_G_i/r_i)*(npower(l,k)/dr_k(l)*u_G_k+dr_k(l)*up_G_k/r_k)
 enddo !l
enddo !kw
end

subroutine compute_u_tilde(i_el,i,k,rt1,rt2,ut1,ut2)
include 'j.inc'
dimension rt1(nw,3),rt2(nw,3)
dimension ut1(3,nw,nbasis_max*nbasis_max),ut2(3,nw,nbasis_max*nbasis_max)
ik = (k-1)*nbasis_max+i
do kw=1,nw
 do ll=1,3
  ut1(ll,kw,ik)= 1.d0/dsqrt(g_min(i)+g_min(k))*rt1(kw,ll) + G_center(ll,i,k)
 enddo !!ll
enddo !!kw
if(i_el.eq.1)return
do kw=1,nw
 do ll=1,3
  ut2(ll,kw,ik)= 1.d0/dsqrt(g_min(i)+g_min(k))*rt2(kw,ll) + G_center(ll,i,k)
 enddo !!ll
enddo !!kw
end

subroutine compute_jacobian(i_el,i,k,j,l,r1,r2,rjacob)
include 'j.inc'
dimension r1(nw,3),r2(nw,3)
dimension rjacob(nw)
do kw=1,nw
 r1_mod=dsqrt(r1(kw,1)**2+r1(kw,2)**2+r1(kw,3)**2)
 call compute_f_fp_jacob(r1_mod,i,k,1,f1,f1p)
 rjacob(kw)=f1**2*dabs(f1+r1_mod*f1p)
enddo
if(i_el.eq.1)return
do kw=1,nw
 r2_mod=dsqrt(r2(kw,1)**2+r2(kw,2)**2+r2(kw,3)**2)
 call compute_f_fp_jacob(r2_mod,j,l,1,f2,f2p)
 rjacob(kw)=rjacob(kw)*f2**2*dabs(f2+r2_mod*f2p)
enddo
end

subroutine compute_f_fp_jacob(r,i,k,ic,f,fp)
include 'j.inc'

 if(basis_type.eq.'STO')then
 f=a_ZV(i,k)*r**alpha_ZV/(2.d0*dsqrt(g_min(i)+g_min(k)))
 if(ic.eq.0)return
 fp=alpha_ZV*a_ZV(i,k)*r**(alpha_ZV-1.d0)/(2.d0*dsqrt(g_min(i)+g_min(k)))
 endif

!return

if(basis_type.eq.'GTO'.and.(finite_range.eqv..true.))then
 f=a_ZV(i,k)/dsqrt(2.d0)*dsqrt((g_min(i)+g_min(k))/(gauss_min(i)+gauss_min(k)))
 if(ic.eq.0)return
 fp=0.d0
endif

end


subroutine compare_nuclei(ni,nj,nk,nl,ndiff)
implicit double precision (a-h,o-z)

! 1111
if( (ni.eq.nj).and.(nj.eq.nk).and.(nk.eq.nl) ) then
 ndiff=1
 return
endif

! 2111
if( (ni.ne.nj).and.(nj.eq.nk).and.(nk.eq.nl) ) then
 ndiff=2
 return
endif

! 1211
if( (ni.ne.nj).and.(nk.eq.ni).and.(nk.eq.nl) ) then
 ndiff=2
 return
endif

! 1121
if( (ni.eq.nj).and.(nj.ne.nk).and.(nl.eq.ni) ) then
 ndiff=2
 return
endif

! 1112
if( (ni.eq.nj).and.(nj.eq.nk).and.(nk.ne.nl) ) then
 ndiff=2
 return
endif

! 1212
if((ni.ne.nj).and.(nk.ne.nl).and.(ni.eq.nk).and.(nj.eq.nl) ) then
 ndiff=2
 return
endif

! 2211
if((ni.eq.nj).and.(nk.eq.nl).and.(ni.ne.nk) ) then
 ndiff=2
 return
endif

! 1221
if((ni.ne.nj).and.(nj.eq.nk).and.(nk.ne.nl).and.(ni.eq.nl) ) then
 ndiff=2
 return
endif


!11xy
if( (ni.eq.nj).and.(ni.ne.nk).and.(ni.ne.nl).and.(nk.ne.nl) ) then
 ndiff=3
 return
endif

!1x1y
if( (ni.eq.nk).and.(ni.ne.nj).and.(ni.ne.nl).and.(nj.ne.nl) ) then
 ndiff=3
 return
endif

!1xy1
if( (ni.eq.nl).and.(ni.ne.nj).and.(ni.ne.nk).and.(nj.eq.nk) ) then
 ndiff=3
 return
endif

!x1y1
if( (nj.eq.nl).and.(nj.ne.ni).and.(nj.ne.nk).and.(ni.ne.nk) ) then
 ndiff=3
 return
endif

!x11y
if( (nj.eq.nk).and.(nj.ne.ni).and.(nj.ne.nl).and.(ni.ne.nl) ) then
 ndiff=3
 return
endif

!xy11
if( (nk.eq.nl).and.(nk.ne.ni).and.(nk.ne.nj).and.(ni.ne.nj) ) then
 ndiff=3
 return
endif

!1234
if( (ni.ne.nj).and.(nj.ne.nk).and.(nk.ne.nl).and.(nl.ne.ni).and.(nl.ne.nj).and.(ni.ne.nk) ) then
 ndiff=4
 return
endif

end

      double precision function halton_vec(i)
      implicit double precision(a-h,o-z)
!*** Halton generator
integer*8 nra1,nrb1,n1h
common/genhaltv/n1h(8),nra1(8),nrb1(8)
!***********! propre au generateur de nombre aleatoire 
      integer*8 nxa1,nxb1

      nxa1=nrb1(i)-nra1(i)
      if(nxa1.eq.1)then
        nra1(i)=1
        nrb1(i)=nrb1(i)*n1h(i)
        goto 19
      endif
      nxb1=nrb1(i)/n1h(i)
      kcp=0
10    if(nxa1.le.nxb1)then
        nxb1=nxb1/n1h(i)
        kcp=kcp+1
        if(mod(kcp,100000).eq.0)then
          write(*,*)'kcp=',kcp
          write(*,*)'nxa1=',nxa1
          write(*,*)'nxb1=',nxa1
          write(*,*)'i=',i
          write(*,*)'n1h(i)=',n1h(i)
          stop
        endif 
        goto 10
      else
        nra1(i)=(1+n1h(i))*nxb1-nxa1
      endif
19    halton_vec=dfloat(nra1(i))/dfloat(nrb1(i))
      end

! inverse of error function in double precision
!
      double precision function dierfc(y)
      implicit real*8 (a - h, o - z)
      parameter (&
          qa = 9.16461398268964d-01, &
          qb = 2.31729200323405d-01, &
          qc = 4.88826640273108d-01, &
          qd = 1.24610454613712d-01, &
          q0 = 4.99999303439796d-01, &
          q1 = 1.16065025341614d-01, &
          q2 = 1.50689047360223d-01, &
          q3 = 2.69999308670029d-01, &
          q4 = -7.28846765585675d-02)
      parameter (&
          pa = 3.97886080735226000d+00, &
          pb = 1.20782237635245222d-01, &
          p0 = 2.44044510593190935d-01, &
          p1 = 4.34397492331430115d-01, &
          p2 = 6.86265948274097816d-01, &
          p3 = 9.56464974744799006d-01, &
          p4 = 1.16374581931560831d+00, &
          p5 = 1.21448730779995237d+00, &
          p6 = 1.05375024970847138d+00, &
          p7 = 7.13657635868730364d-01, &
          p8 = 3.16847638520135944d-01, &
          p9 = 1.47297938331485121d-02, &
          p10 = -1.05872177941595488d-01, &
          p11 = -7.43424357241784861d-02)
      parameter (&
          p12 = 2.20995927012179067d-03, &
          p13 = 3.46494207789099922d-02, &
          p14 = 1.42961988697898018d-02, &
          p15 = -1.18598117047771104d-02, &
          p16 = -1.12749169332504870d-02, &
          p17 = 3.39721910367775861d-03, &
          p18 = 6.85649426074558612d-03, &
          p19 = -7.71708358954120939d-04, &
          p20 = -3.51287146129100025d-03, &
          p21 = 1.05739299623423047d-04, &
          p22 = 1.12648096188977922d-03)
      z = y
      if (y .gt. 1.d0) z = 2.d0 - y
      w = qa - log(z)
      u = sqrt(w)
      s = (qc + log(u)) / w
      t = 1.d0 / (u + qb)
      x = u * (1 - s * (0.5d0 + s * qd)) - &
          ((((q4 * t + q3) * t + q2) * t + q1) * t + q0) * t
      t = pa / (pa + x)
      u = t - 0.5d0
      s = (((((((((p22 * u + p21) * u + p20) * u + &
          p19) * u + p18) * u + p17) * u + p16) * u + &
          p15) * u + p14) * u + p13) * u + p12
      s = ((((((((((((s * u + p11) * u + p10) * u +  &
          p9) * u + p8) * u + p7) * u + p6) * u + p5) * u + &
          p4) * u + p3) * u + p2) * u + p1) * u + p0) * t - &
          z * exp(x * x - pb)
      x = x + s * (1.d0 + x * s)
      if (y .gt. 1.d0) x = -x
      dierfc = x
      end

double precision function I_n_special(epsilo,gama,n,lmax)
! function that calculates the following integral involved in the bielectronic integrals :
!      int {t,[0,epsilo]} of exp(-gama t^2)*t^n
 implicit none
 integer :: n,lmax
 double precision :: epsilo,gama
 integer :: i,j,k,l
 double precision :: accu,fact,power,sq_gama,pouet
 sq_gama = dsqrt(gama)
 pouet = sq_gama * epsilo
 accu = 0.d0
 do l = 0,lmax
  accu = accu + power(l,-1.d0)*power(2*l,pouet) / (fact(l) * (dble(n+2*l+1)))
 enddo
!print*,'accu = ',accu
 I_n_special = 1.d0/dsqrt(gama**(n+1)) * accu * power(n+1,pouet)
end


double precision function I_n_bibi(n,lmax,alpha)
 implicit none
! integral of exp(-alpha * t^2) * t^n between [0:1]
 integer :: lmax,n
 double precision :: alpha
 integer :: l
 double precision :: fact,power
 I_n_bibi = 0.d0
 do l = 0, lmax,1
   I_n_bibi = I_n_bibi + (-1.d0)**l * power(2*l,alpha) / (fact(l) * (n+2 * l +1))
 enddo

end

double precision function I_n_michel(epsilo,n,lmax)
! function that calculates the following integral involved in the bielectronic integrals :
!      int {t,[0,epsilo]} of exp(-t^2)*t^n
!      this integral is estimated as a taylor expension of length "lmax"
!      lmax = 20 is suffisent to get the 8 first decimals
 implicit none
 integer :: n,lmax
 double precision :: epsilo
 integer :: i,j,k,l
 double precision :: accu,fact,power
 accu = 0.d0
 do l = 0,lmax
  accu = accu + power(l,-1.d0)*power(2*l,epsilo) / (fact(l) * (dble(n+2*l+1)))
 enddo
 I_n_michel = power(n+1,epsilo) * accu

end

 double precision function power(n,x)
 implicit none
 integer :: i,n
 double precision :: x,accu
 power= 1.d0

 do i = n,1,-1
  power = power* x
 enddo
 end

 double precision function I_n_special_exact(n,epsilo,gama)
 implicit none
 integer :: n
 double precision :: epsilo_prim,epsilo,gama,I_n
 epsilo_prim = dsqrt(gama)*epsilo
!print*,'epsilo_prim = ',epsilo_prim
!print*,'I_0 = ',I_0
!print*,'I_n(n,epsilo_prim,I_0) = ',I_n(n,epsilo_prim)
!print*,'gama**(-dble(n+1)/2.d0) = ',gama**(-dble(n+1)/2.d0)
 I_n_special_exact = gama**(-dble(n+1)/2.d0) * I_n(n,epsilo_prim)

 end

 double precision function I_n (n,epsilo)
! exact value of 
!  integral [0:epsilo] dt { exp(-t^2)  t^n }
 implicit none
 integer :: n
 double precision :: epsilo
 integer  :: i,j,k,l
 double precision :: accu,delta,b,prod
 double precision :: power,pouet
 double precision :: I_0
 double precision :: sqpi
 sqpi = dsqrt(dacos(-1.d0))
 I_0 =  0.5d0 * sqpi*erf(epsilo)
 accu = 0.d0
 delta = -0.5d0*exp(-epsilo*epsilo)
 b = 0.5d0
 do i = 1, n/2
  l =  i+i-1
  pouet = epsilo**l
  prod = 1.d0
  do j = 1, n/2 - i
   prod = prod * b * (2.d0*dble(i+j)-1.d0)
  enddo
  accu = accu + prod * pouet
 enddo
 accu = accu * delta

 prod = 1.d0
 do i = 1, n/2
  prod = prod  * dble(2*i-1)
 enddo
 prod = prod * I_0 * b **(n/2)
 I_n = accu + prod
 end

 double precision function integrate_bourrin_standard(n,gama,nx)
 implicit none
 integer, intent(in) :: n,nx
 double precision, intent(in) :: gama
 integer :: i,j
 double precision :: dx,x,power
 dx = 1.d0/dble(nx)
 integrate_bourrin_standard = 0.d0
 x = 0.d0
 do i = 1,nx
  x = x + dx
  integrate_bourrin_standard = integrate_bourrin_standard + dexp(-gama * x * x) * power(n,x)
 enddo
 integrate_bourrin_standard = integrate_bourrin_standard * dx
 end


 double precision function integrate_bourrin(n,epsilo,gama,nx)
 implicit none
 integer, intent(in) :: n,nx
 double precision, intent(in) :: epsilo,gama
 integer :: i,j
 double precision :: dx,x,power
 dx = epsilo/dble(nx)
 integrate_bourrin = 0.d0
 x = 0.d0
 do i = 1,nx
  x = x + dx
  integrate_bourrin = integrate_bourrin + dexp(-gama * x * x) * x ** n
 enddo
 integrate_bourrin = integrate_bourrin * dx
 end

recursive subroutine I_x2_new(c,B_10,B_01,B_00,res,n_pt) 
  implicit none
  !EGIN_DOC
  !  recursive function involved in the bielectronic integral
  !ND_DOC
  integer, intent(in)            :: c, n_pt
  double precision, intent(in)   :: B_10(1000),B_01(1000),B_00(1000)
  double precision, intent(out)  :: res(1000)
  integer                        :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter             :: max_dim = 511
  double precision, parameter    :: pi =  dacos(-1.d0)
  double precision, parameter    :: sqpi =  dsqrt(dacos(-1.d0))
  double precision, parameter    :: pi_5_2 =  34.9868366552d0
  double precision, parameter    :: dfour_pi =  4.d0*dacos(-1.d0)
  double precision, parameter    :: dtwo_pi =  2.d0*dacos(-1.d0)
  double precision, parameter    :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: thresh = 1.d-15
  double precision, parameter    :: cx_lda = -0.73855876638202234d0
  double precision, parameter    :: c_2_4_3 = 2.5198420997897464d0
  double precision, parameter    :: cst_lda = -0.93052573634909996d0
  double precision, parameter    :: c_4_3 = 1.3333333333333333d0
  double precision, parameter    :: c_1_3 = 0.3333333333333333d0
  
  if(c==1)then
    do i=1,n_pt
      res(i) = 0.d0
    enddo
  elseif(c==0) then
    do i=1,n_pt
      res(i) = 1.d0
    enddo
  else
    call I_x1_new(0,c-2,B_10,B_01,B_00,res,n_pt)
    do i=1,n_pt
      res(i) =  (c-1) * B_01(i) * res(i)
    enddo
  endif
end


integer function n_pt_sup(a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z)
  implicit none
  !EGIN_DOC
  ! Returns the upper boundary of the degree of the polynomial involved in the
  ! bielctronic integral :
  !       Ix(a_x,b_x,c_x,d_x) * Iy(a_y,b_y,c_y,d_y) * Iz(a_z,b_z,c_z,d_z)
  !ND_DOC
  integer                        :: a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z
  n_pt_sup =  ishft( a_x+b_x+c_x+d_x + a_y+b_y+c_y+d_y + a_z+b_z+c_z+d_z,1 )
end

!  Y_l^m(theta,phi) = i^(m+|m|) ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2  P_l^|m|(cos(theta))  exp(i m phi)
!  l=0,1,2,....
!  m=0,1,...,l
! Here:
!  n=l (n=0,1,...)  
!  m=0,1,...,n 
!  x=cos(theta) 0 < x < 1
!
!  
!  This routine computes:   PM(m,n) for n=0,...,N (number N in input) and m=0,..,n

!   Exemples (see 'Associated Legendre Polynomilas wikipedia')
!    P_{0}^{0}(x)=1
!    P_{1}^{-1}(x)=-1/2 P_{1}^{1}(x)
!    P_{1}^{0}(x)=x
!    P_{1}^{1}(x)=-(1-x^2)^{1/2}
!    P_{2}^{-2}(x)=1/24 P_{2}^{2}(x)
!    P_{2}^{-1}(x)=-1/6 P_{2}^{1}(x)
!    P_{2}^{0}(x)=1/2 (3x^{2}-1)
!    P_{2}^{1}(x)=-3x(1-x^2)^{1/2}
!    P_{2}^{2}(x)=3(1-x^2)


        SUBROUTINE LPMN(MM,M,N,X,PM)
!
! Here N = LMAX
! Here M= MMAX (we take M=LMAX in input)
!
!       =====================================================
!       Purpose: Compute the associated Legendre functions Pmn(x)
!       Input :  x  --- Argument of Pmn(x)
!                m  --- Order of Pmn(x),  m = 0,1,2,...,n
!                n  --- Degree of Pmn(x), n = 0,1,2,...,N
!                mm --- Physical dimension of PM 
!       Output:  PM(m,n) --- Pmn(x)
!       =====================================================
!
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PM(0:MM,0:(N+1))
        DO 10 I=0,N
        DO 10 J=0,M
10         PM(J,I)=0.0D0
        PM(0,0)=1.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
15            PM(0,I)=X**I
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1)
        DO 35 I=0,M
35         PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)

        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)- (I+J-1.0D0)*PM(I,J-2))/(J-I)
40      CONTINUE

        END

!! alternative way of computing associated legendre polynomial
!! CALL LPMN(100,l,l,X,PM)  ---> d0_asso_legendre_poly(l,m,x)=pm(m,l)
recursive double precision function d0_asso_legendre_poly(l,m,x) result(res)
implicit none
integer l,m
double precision x
double precision f,expo,doublefact
res=0.d0
if(m.gt.l)return
if(m.lt.0)return
if(l.eq.0.and.m.eq.0)then
res=1.d0
return
endif
if(m.eq.l)then
 f=(-1.d0)**l*doublefact(2*l-1)
 expo=dfloat(l)/2.d0
 res=f*(1.d0-x**2)**expo
 return
endif
res = ((2*l-1)*x*d0_asso_legendre_poly(l-1,m,x)  &
       -(l-1+m)*d0_asso_legendre_poly(l-2,m,x))  &
      /(l-m)
end

!! Computing the first derivative wrt x of associated legendre polynomial
!! d1_asso_legendre_poly(l,m,x)
recursive double precision function d1_asso_legendre_poly(l,m,x) result(res)
implicit none
integer l,m
double precision x
double precision f,expo,doublefact,d0_asso_legendre_poly
res=0.d0
if(m.gt.l)return
if(m.lt.0)return
if(l.eq.0.and.m.eq.0)then
res=0.d0
return
endif
if(m.eq.l)then
 f=(-1.d0)**l*doublefact(2*l-1)
 expo=dfloat(l)/2.d0-1.d0
 res=-f*l*x*(1.d0-x**2)**expo
 return
endif
res =( (2*l-1)*x*d1_asso_legendre_poly(l-1,m,x) &
      +(2*l-1)*d0_asso_legendre_poly(l-1,m,x)-(l-1+m)*d1_asso_legendre_poly(l-2,m,x)) &
    /(l-m)
end

double precision function dblefact(n)
implicit none
integer :: n,k
double precision prod
dblefact=1.d0
if(n.lt.0)return
if(mod(n,2).eq.1)then
  prod=1.d0
  do k=1,n,2
   prod=prod*dfloat(k)
  enddo
  dblefact=prod
  return
 endif
 if(mod(n,2).eq.0)then
  prod=1.d0
  do k=2,n,2
   prod=prod*dfloat(k)
  enddo
  dblefact=prod
  return
 endif
end

!c Explicit representation of Legendre polynomials P_n(x) 
!!
!! P_n0(x) = P_n(x)= \sum_{k=0}^n a_k x^k
!!
!!   with  a_k= 2^n  binom(n,k) binom( (n+k-1)/2, n )
!! coef_pm(n,k) is the k_th coefficient of P_n(x)
double precision function coef_pm(n,k)
implicit none
integer n,k
double precision arg,binom,binom_gen
if(n.eq.0.and.k.ne.0)stop 'coef_pm not defined'
if(n.eq.0.and.k.eq.0)then
coef_pm=1.d0
return
endif
arg=0.5d0*dfloat(n+k-1)
coef_pm=2.d0**n*binom(n,k)*binom_gen(arg,n)
end


!  Y_l^m(theta,phi) = i^(m+|m|) ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2  P_l^|m|(cos(theta))  exp(i m phi)
!  l=0,1,2,....
!  m=0,1,...,l
! Here:
!  n=l (n=0,1,...)  
!  m=0,1,...,n 
!  x=cos(theta) 0 < x < 1
!
!  
!  This routine computes:   PM(m,n) for n=0,...,N (number N in input) and m=0,..,n

!   Exemples (see 'Associated Legendre Polynomilas wikipedia')
!    P_{0}^{0}(x)=1
!    P_{1}^{-1}(x)=-1/2 P_{1}^{1}(x)
!    P_{1}^{0}(x)=x
!    P_{1}^{1}(x)=-(1-x^2)^{1/2}
!    P_{2}^{-2}(x)=1/24 P_{2}^{2}(x)
!    P_{2}^{-1}(x)=-1/6 P_{2}^{1}(x)
!    P_{2}^{0}(x)=1/2 (3x^{2}-1)
!    P_{2}^{1}(x)=-3x(1-x^2)^{1/2}
!    P_{2}^{2}(x)=3(1-x^2)


!c  Computation of real spherical harmonics Ylm(theta,phi)
!c
!c  l=0,1,....
!c  m=-l,l
!c
!c  m>0: Y_lm = sqrt(2) ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2  P_l^|m|(cos(theta)) cos(m phi)
!c  m=0: Y_l0 =         ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2  P_l^0  (cos(theta)) 
!c  m<0: Y_lm = sqrt(2) ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2  P_l^|m|(cos(theta)) sin(|m|phi)

!Examples(wikipedia http://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics)

! l = 0

! Y_00 = \sqrt{1/(4pi)}

! l = 1

! Y_1-1= \sqrt{3/(4pi)} y/r
! Y_10 = \sqrt{3/(4pi)} z/r
! Y_11 = \sqrt{3/(4pi)} x/r
!
! l = 2
!
! Y_2,-2= 1/2 \sqrt{15/pi}  xy/r^2 
! Y_2,-1= 1/2 \sqrt{15/pi}  yz/r^2 
! Y_20  = 1/4 \sqrt{15/pi} (-x^2-y^2 +2z^2)/r^2 
! Y_21  = 1/2 \sqrt{15/pi}  zx/r^2 
! Y_22  = 1/4 \sqrt{15/pi} (x^2-y^2)/r^2 
!
!c
double precision function ylm(l,m,theta,phi)
implicit none
integer l,m
double precision theta,phi,pm,factor,pi,x,fact,sign
DIMENSION PM(0:100,0:100)
pi=dacos(-1.d0)
x=dcos(theta)
sign=(-1.d0)**m
CALL LPMN(100,l,l,X,PM)
factor=dsqrt( (2*l+1)*fact(l-iabs(m)) /(4.d0*pi*fact(l+iabs(m))) )
if(m.gt.0)ylm=sign*dsqrt(2.d0)*factor*pm(m,l)*dcos(dfloat(m)*phi)
if(m.eq.0)ylm=factor*pm(m,l)
if(m.lt.0)ylm=sign*dsqrt(2.d0)*factor*pm(iabs(m),l)*dsin(dfloat(iabs(m))*phi)
end

double precision function binom_gen(alpha,n)
 implicit none
 integer :: n,k
 double precision :: fact,alpha,prod
 prod=1.d0
 do k=0,n-1
   prod=prod*(alpha-k)
   binom_gen = prod/(fact(n))
 enddo
end

!!!   coef_nk(n,k)= 1/2**k 1/[k! [2(n+k)+1)]!!
!!!

      double precision function coef_nk(n,k)
      implicit none
      integer n,k
      double precision gam,dblefact,fact
      gam=dblefact(2*(n+k)+1)
      coef_nk=1.d0/(2.d0**k*fact(k)*gam)
      end

      double precision function doublefact(n)
      implicit double precision(a-h,o-z)
      doublefact=1.d0
      if(n.le.2)return
      d=1.d0
      do i=n,1,-2
       d=d*dfloat(i)
      enddo
      doublefact=d
      end

     double precision function one_electron_Kij(na,ra,gamA,nb,rb,gamB)
     implicit double precision(a-h,o-z)
     dimension na(3),ra(3),nb(3),rb(3),rp(3)
     dimension na2(3)

     sab=one_electron_Sij(na,ra,gamA,nb,rb,gamB)
     term1=gamA*(3.d0+2.d0*(na(1)+na(2)+na(3)))*sab
    
     sab_2=0.d0 
     do l=1,3
      do li=1,3
       if(li.eq.l)then
        na2(li)=na(li)+2
       else  
        na2(li)=na(li)
       endif
      enddo
      sab_2=sab_2+one_electron_Sij(na2,ra,gamA,nb,rb,gamB)
     enddo
     
     term2=-2.d0*gamA**2*sab_2

     sab_12=0.d0 
     do l=1,3
      do li=1,3
       if(li.eq.l.and.na(li).ge.2)then
        na2(li)=na(li)-2
       else  
        na2(li)=na(li)
       endif
      enddo
      sab_12=sab_12+na(l)*(na(l)-1)*one_electron_Sij(na2,ra,gamA,nb,rb,gamB)
     enddo

     term3=-0.5d0*sab_12

     one_electron_Kij=term1+term2+term3

     end


!!cc Calculation of  I = <phiA|phiB> = /int dxdydz  phiA(r) phiB(r)
!!cc *************************************************************
!!cc  r=(x,y,z)
!!cc  A=[ra(1),ra(2),ra(3]
!!cc  B=[rb(1),rb(2),rb(3)]
!!cc phiA(r)=[x-ra(1)]**na(1)*[y-ra(2)]**na(2)*[z-ra(3)]**na(3)*exp[-gamA*(x-ra(1))**2-gamA(y-ra(2))**2-gamA(z-ra(3))**2)]
!!cc phiB(r)=[x-rb(1)]**nb(1)*[y-rb(2)]**nb(2)*[z-rb(3)]**nb(3)*exp[-gamB*(x-rb(1))**2-gamB(y-rb(2))**2-gamB(z-rb(3))**2)]

!!cc  mu=gamA*gamB/(gamA+gamB)
!!cc  rp(l)=(gamA*ra(l)+gamB*rb(l))/(gamA+gamB)

!!c  I = exp[-mu*(ra-rb)**2] * Ix*Iy*Iz
!!cc    where:
!!cc  Ix = /int dx [x-ra(1)]**na(1)*[x-rb(1)]**nb(1)*exp[-(gamA+gamB)*(x-rp(1))**2)]
!!cc  Iy = /int dy [y-ra(2)]**na(2)*[y-rb(2)]**nb(2)*exp[-(gamA+gamB)*(y-rp(2))**2)]
!!cc  Iz = /int dz [z-ra(3)]**na(3)*[z-rb(3)]**nb(3)*exp[-(gamA+gamB)*(z-rp(3))**2)]
!!cc
!!cc  one_electron_sij = I

      double precision function one_electron_Sij(na,ra,gamA,nb,rb,gamB)
      implicit double precision(a-h,o-z)
      dimension na(3),ra(3),nb(3),rb(3),rp(3)
      pi=dacos(-1.d0)

      gamtot=gamA+gamB
      rmu=gamA*gamB/gamtot

      do l=1,3
       rp(l)=(gamA*ra(l)+gamB*rb(l))/gamtot
      enddo
      AmB=(ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2

      one_electron_sij=0.d0
      arg=rmu*AmB
      if(arg.gt.-dlog(1.d-15))return

      one_electron_sij=1.d0
      do l=1,3
       one_electron_sij=one_electron_sij &
     * one_electron_I(na(l),ra(l),nb(l),rb(l),gamtot,rp(l))

      enddo
      one_electron_sij=dexp(-arg)*one_electron_sij
      end

      double precision function check_one_electron_Sij(na,ra,gamA,nb,rb,gamB)
      implicit double precision(a-h,o-z)
      dimension na(3),ra(3),nb(3),rb(3),rp(3)
      pi=dacos(-1.d0)

      bigM=30.d0

      do kk=2,4

       npts=4**kk
       dx=2.d0*bigM/npts
       accu=0.d0

       do ix=1,npts
        x=-bigM+(ix-1)*dx
        do iy=1,npts
         y=-bigM+(iy-1)*dx
         do iz=1,npts
          z=-bigM+(iz-1)*dx

          pA=(x-ra(1))**na(1)*(y-ra(2))**na(2)*(z-ra(3))**na(3)
          pB=(x-rb(1))**nb(1)*(y-rb(2))**nb(2)*(z-rb(3))**nb(3)

          rA2=(x-ra(1))**2+(y-ra(2))**2+(z-ra(3))**2
          rB2=(x-rb(1))**2+(y-rb(2))**2+(z-rb(3))**2
          arg=gamA*rA2+gamB*rB2

          if(arg.lt.-dlog(1.d-15))then
           accu=accu+pA*pB*exp(-arg)
          endif

          enddo
         enddo
        enddo

        rint=(dx)**3*accu
!        print*,'N',npts**3,rint

       enddo

       check_one_electron_Sij=(dx)**3*accu

      end

      double precision function one_electron_Vij(na,ra,gamA,nb,rb,gamB)
      include 'j.inc'
      dimension na(3),ra(3),nb(3),rb(3)

      double precision NAI_pol_mult

      integer dim_int_tab(4)
      integer dim_int

      dim_int_tab(1) = maxval(na)
      dim_int_tab(2) = maxval(nb)
      dim_int_tab(3) = 0
      dim_int_tab(4) = 0
      dim_int = (maxval(dim_int_tab) *24 + 4)*8

      ac=0.d0
      do kn=1,number_nuclei
       ac=ac-charge(kn)*NAI_pol_mult(ra,rb,na,nb,gamA,gamB,centers_nuclei(1,kn),dim_int)
      enddo

      one_electron_Vij=ac

      end

      double precision function check_one_electron_Vij(na,ra,gamA,nb,rb,gamB)
      include 'j.inc'
      dimension na(3),ra(3),nb(3),rb(3),rp(3)

      pi=dacos(-1.d0)

      bigM=30.d0

      do kk=2,4

       npts=4**kk
       dx=2.d0*bigM/npts
       accu=0.d0

       do ix=1,npts
        x=-bigM+(ix-1)*dx+ 0.23*dx
        do iy=1,npts
         y=-bigM+(iy-1)*dx + 0.22*dx
         do iz=1,npts
          z=-bigM+(iz-1)*dx + 0.29*dx

          pA=(x-ra(1))**na(1)*(y-ra(2))**na(2)*(z-ra(3))**na(3)
          pB=(x-rb(1))**nb(1)*(y-rb(2))**nb(2)*(z-rb(3))**nb(3)

          pot=0.d0
          do kn=1,number_nuclei
           rC=dsqrt( &
   (x-centers_nuclei(1,kn))**2 + (y-centers_nuclei(2,kn))**2 + (z-centers_nuclei(3,kn))**2 )
           pot=pot-charge(kn)/rC
          enddo
 
          rA2=(x-ra(1))**2+(y-ra(2))**2+(z-ra(3))**2
          rB2=(x-rb(1))**2+(y-rb(2))**2+(z-rb(3))**2
          arg=gamA*rA2+gamB*rB2


          if(arg.lt.-dlog(1.d-15))then
           accu=accu+pot*pA*pB*dexp(-arg)
          endif

          enddo
         enddo
        enddo

        rint=(dx)**3*accu
!        print*,'N',npts**3,rint

       enddo

       check_one_electron_Vij=(dx)**3*accu

      end

      double precision function check_one_electron_Kij(na,ra,gamA,nb,rb,gamB)
      implicit double precision(a-h,o-z)
      dimension na(3),ra(3),nb(3),rb(3),rp(3)

      pi=dacos(-1.d0)

      bigM=30.d0

      do kk=2,4

       npts=4**kk
       dx=2.d0*bigM/npts
       accu=0.d0

       do ix=1,npts
        x=-bigM+(ix-1)*dx+0.21*dx
        do iy=1,npts
         y=-bigM+(iy-1)*dx+0.11*dx
         do iz=1,npts
          z=-bigM+(iz-1)*dx+0.27*dx

          pA=(x-ra(1))**na(1)*(y-ra(2))**na(2)*(z-ra(3))**na(3)
          pB=(x-rb(1))**nb(1)*(y-rb(2))**nb(2)*(z-rb(3))**nb(3)

          rA2=(x-ra(1))**2+(y-ra(2))**2+(z-ra(3))**2
          rB2=(x-rb(1))**2+(y-rb(2))**2+(z-rb(3))**2
          arg=gamA*rA2+gamB*rB2

          rkin=-0.5d0*( -6.d0*gamA-4.d0*gamA*(na(1)+na(2)+na(3))  &
       + 4.d0*gamA**2*((x-ra(1))**2+ (y-ra(2))**2+(z-ra(3))**2) &
       + na(1)*(na(1)-1.d0)/(x-ra(1))**2 &
       + na(2)*(na(2)-1.d0)/(y-ra(2))**2 &
       + na(3)*(na(3)-1.d0)/(z-ra(3))**2 )


          if(arg.lt.-dlog(1.d-15))then
           accu=accu+rkin*pA*pB*exp(-arg)
          endif

          enddo
         enddo
        enddo

        rint=(dx)**3*accu

       enddo

       check_one_electron_Kij=(dx)**3*accu

      end

!! Here gam=gamA+gamB
!! one_electron_I= Ix
      double precision function one_electron_I(na,ra,nb,rb,gam,rp)
      implicit double precision(a-h,o-z)
      pi=dacos(-1.d0)
      one_electron_I=0.d0

      rint=0.d0
      do kx=0,na
       do lx=0,nb
        iskip=0
        kbig=na+nb-kx-lx
        if(mod(kbig,2).eq.1)iskip=1
        if(iskip.eq.0)then
         kbig=kbig/2
         f3=doublefact(2*kbig-1) 
         f2=(1.d0/(2.d0*gam))**kbig
         rint=rint+f2*f3*binom(na,kx)*binom(nb,lx)*(rp-ra)**kx*(rp-rb)**lx
        endif
       enddo
      enddo

      one_electron_I=dsqrt(pi/gam)*rint
      end

      double precision function check_one_electron_I(na,ra,nb,rb,gam,rp)
      implicit double precision(a-h,o-z)

      pi=dacos(-1.d0)

      bigM=8.d0

      do kk=2,6

       npts=4**kk
       dx=2.d0*BigM/npts
       accu=0.d0

       do ix=1,npts
        x=-bigM+(ix-1)*dx
        pA=(x-ra)**na
        pB=(x-rb)**nb
        arg=gam*(x-rp)**2
        if(arg.lt.-dlog(1.d-15))then
         accu=accu+pA*pB*exp(-arg)
        endif
       enddo

        rint=dx*accu

       enddo

       check_one_electron_I=dx*accu

      end

      subroutine erreurXoverY(x,y,n,rmoy,error)
      implicit double precision (a-h,o-z)
      dimension x(n),y(n)
      ansum=0.d0
      avsum=0.d0
      do i=1,n
        if(i.ne.1) avcu0=avsum/ansum
        ansum=ansum+y(i)
        avsum=avsum+x(i)
        avbl=x(i)/y(i)
        if(ansum.ne.0.0) avcu=avsum/ansum
        if(i.eq.1)then
          avsq=0.d0
        else
          avsq=avsq+(1.d0-y(i)/ansum)*(avbl-avcu0)**2*y(i)
        endif
      enddo
      arg=dabs(avsq/(ansum*(dfloat(n)-1.d0)))
      error=dsqrt(arg)
      rmoy=avcu
      end

!  Y_l^m(theta,phi) = i^(m+|m|) ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2  P_l^|m|(cos(theta))  exp(i m phi)
!  l=0,1,2,....
!  m=0,1,...,l
! Here:
!  n=l (n=0,1,...)  
!  m=0,1,...,n 
!  x=cos(theta) 0 < x < 1
!
!  
!  This routine computes:   PM(m,n) for n=0,...,N (number N in input) and m=0,..,n

!   Exemples (see 'Associated Legendre Polynomilas wikipedia')
!    P_{0}^{0}(x)=1
!    P_{1}^{-1}(x)=-1/2 P_{1}^{1}(x)
!    P_{1}^{0}(x)=x
!    P_{1}^{1}(x)=-(1-x^2)^{1/2}
!    P_{2}^{-2}(x)=1/24 P_{2}^{2}(x)
!    P_{2}^{-1}(x)=-1/6 P_{2}^{1}(x)
!    P_{2}^{0}(x)=1/2 (3x^{2}-1)
!    P_{2}^{1}(x)=-3x(1-x^2)^{1/2}
!    P_{2}^{2}(x)=3(1-x^2)


!!   ntot_fit_slater(ng,n,l,m)
!!   n_fit_slater(ng,n,l,m,1,i) [ i=1,ntot_fit_slater(ng,n,l,m) ]
!!   n_fit_slater(ng,n,l,m,2,i) [ i=1,ntot_fit_slater(ng,n,l,m) ]
!!   n_fit_slater(ng,n,l,m,3,i) [ i=1,ntot_fit_slater(ng,n,l,m) ]
!!   g_fit_slater(ng,n,l,m,i)   [ i=1,ntot_fit_slater(ng,n,l,m) ]
!!   c_fit_slater(ng,n,l,m,i)   [ i=1,ntot_fit_slater(ng,n,l,m) ]
!!
!!   phi_S(n,l,m,g=1,vec_r)= N_n(1) r^{n-1} exp(-r) Y_lm(theta,phi)
!!
!! is written as  sum{i=1}^{ntot} c(i)*phi_g(nx(i),ny(i),nz(i),g(i),vec_r)
!!

double precision function determinant(a)
implicit double precision(a-h,o-z)
dimension a(6,6),indx(6)
      n=6
      np=6
      call ludcmp(a,n,np,indx,d)
      d=1.d0
      do j=1,n
       d=d*a(j,j)
      enddo
      determinant=d
      end

      subroutine ludcmp(a,n,np,indx,d)
      implicit double precision(a-h,o-z)
      parameter (nmax=6,tiny=1.0d-20)
      dimension a(np,np),indx(n),vv(nmax)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0.) stop 'singular matrix.'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.)a(j,j)=tiny
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.)a(n,n)=tiny
      return
      end

!! crochet = int_0^infty x**n exp(-g*x**2)
!! crochet= (n-1)!!/(2g)**(n+1)/2
!! crochet = crochet * sqrt(pi/2) if n is even
!! crochet = crochet if n is odd

double precision function crochet(n,g)
implicit none
integer n
double precision g,pi,dblefact,expo
pi=dacos(-1.d0)
expo=0.5d0*dfloat(n+1)
crochet=dblefact(n-1)/(2.d0*g)**expo
if(mod(n,2).eq.0)crochet=crochet*dsqrt(pi/2.d0)
end

double precision function binom(i,j)
 implicit none
 integer :: i,j,k
 double precision :: fact,prod
 if(i.ge.0)then
  binom=fact(i)/(fact(j)*fact(i-j))
  return
 else
  if(j.eq.0)then
   binom=1.d0
   return
  endif
  prod=1.d0
  do k=0,j-1
   prod=prod*(i-k)
  enddo
  binom=prod/fact(j)
 endif
end

!!    M_n(x) modified spherical bessel function 
!!
      double precision function int_prod_bessel_loc(l,gam,n,a)
      implicit none
      integer n,k,l,kcp
      double precision gam,a,gam_save,a_save
      double precision int,intold,coef_nk,crochet
      logical done

      a_save=a
      gam_save=gam
      a=a/dsqrt(gam)
      gam=1.d0

      k=0
      intold=-1.d0
      int=0.d0
      done=.false.
      kcp=0
      do while (.not.done)
       kcp=kcp+1
       int=int+coef_nk(n,k)*a**(n+2*k)*crochet(2*k+n+l,gam)
       if(dabs(int-intold).lt.1d-15)then
        done=.true.
       else
        k=k+1
        intold=int
       endif
      enddo
      int=int*(1.d0/dsqrt(gam_save))**(l+1)
      int_prod_bessel_loc=int

      a=a_save
      gam=gam_save

      end

double precision function der(g1,g2,g3,g4,c1,c2,c3,c4,n1,n2,n3,n4)
implicit double precision (a-h,o-z)
dimension c1(3),c2(3),c3(3),c4(3)
dimension sig(4)

!do neps=2,6
! eps=10.d0**(-neps)

eps=1.d-04

 ac=0.d0
 do i1=0,n1
  sig(1)=(-1.d0)**i1
  do i2=0,n2
   sig(2)=(-1.d0)**i2
   do i3=0,n3
    sig(3)=(-1.d0)**i3
    do i4=0,n4
     sig(4)=(-1.d0)**i4
     ac=ac+sig(1)*sig(2)*sig(3)*sig(4)*   &
           refi( g1+(1-2*i1)*eps,g2+(1-2*i2)*eps, &
                 g3+(1-2*i3)*eps,g4+(1-2*i4)*eps, &
                 c1,c2,c3,c4 )
    enddo
   enddo
  enddo
 enddo

 ntot=n1+n2+n3+n4
 denom=(2.d0*eps)**ntot
 der=(-1.d0)**ntot*ac/denom

!enddo

end

double precision function refi(g1,g2,g3,g4,c1,c2,c3,c4)
implicit double precision (a-h,o-z)
dimension c1(3),c2(3),c3(3),c4(3)
dimension vector_P(3),vector_Q(3)

p=g1+g2
q=g3+g4
rho=p*q/(p+q)
pmq=0.d0
amb=0.d0
cmd=0.d0
do l=1,3
 vector_P(l)=(g1*c1(l)+g2*c2(l))/p
 vector_Q(l)=(g3*c3(l)+g4*c4(l))/q
 pmq=pmq+(vector_P(l)-vector_Q(l))**2
 amb=amb+(c1(l)-c2(l))**2
 cmd=cmd+(c3(l)-c4(l))**2
enddo
eab=dexp(-g1*g2/p*amb)
ecd=dexp(-g3*g4/q*cmd)

!write(40,*)' p q rh amb cmd ead ecd',p,q,rho,amb,cmd,eab,ecd

pi=dacos(-1.d0)

factor=2.d0*pi**2.5d0
factor=factor/(p*q*dsqrt(p+q))

u=dsqrt(rho*pmq)

if(u.ne.0.d0)then
 refi=factor*eab*ecd*dsqrt(pi/4.d0)* (erf(u)/u)
else
 refi=factor*eab*ecd*dsqrt(pi/4.d0)* (2.d0/dsqrt(pi))
endif

end

recursive subroutine I_x1_new(a,c,B_10,B_01,B_00,res,n_pt) 
  !EGIN_DOC
  !  recursive function involved in the bielectronic integral
  !ND_DOC
  implicit none
  integer, intent(in)            :: a,c,n_pt
  double precision, intent(in)   :: B_10(1000),B_01(1000),B_00(1000)
  double precision, intent(out)  :: res(1000)
  double precision               :: res2(1000)
  integer                        :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter             :: max_dim = 511
  double precision, parameter    :: pi =  dacos(-1.d0)
  double precision, parameter    :: sqpi =  dsqrt(dacos(-1.d0))
  double precision, parameter    :: pi_5_2 =  34.9868366552d0
  double precision, parameter    :: dfour_pi =  4.d0*dacos(-1.d0)
  double precision, parameter    :: dtwo_pi =  2.d0*dacos(-1.d0)
  double precision, parameter    :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: thresh = 1.d-15
  double precision, parameter    :: cx_lda = -0.73855876638202234d0
  double precision, parameter    :: c_2_4_3 = 2.5198420997897464d0
  double precision, parameter    :: cst_lda = -0.93052573634909996d0
  double precision, parameter    :: c_4_3 = 1.3333333333333333d0
  double precision, parameter    :: c_1_3 = 0.3333333333333333d0
  
  if(c<0)then
    do i=1,n_pt
      res(i) = 0.d0
    enddo
  else if (a==0) then
    call I_x2_new(c,B_10,B_01,B_00,res,n_pt)
  else if (a==1) then
    call I_x2_new(c-1,B_10,B_01,B_00,res,n_pt)
    do i=1,n_pt
      res(i) = c * B_00(i) * res(i)
    enddo
  else
    call I_x1_new(a-2,c,B_10,B_01,B_00,res,n_pt)
    call I_x1_new(a-1,c-1,B_10,B_01,B_00,res2,n_pt)
    do i=1,n_pt
      res(i) = (a-1) * B_10(i) * res(i) &
               + c * B_00(i) * res2(i)
    enddo
  endif
end
  
!*****************************************************************************

!******************************************************************************
subroutine CalcBoysF(maxm,t,Fm)

! Comute the generalized Boys function Fm(t) using Slatec library

  implicit none

! Input variables

  double precision,intent(in)   :: t
  integer,intent(in)            :: maxm

! Local variables

  integer                       :: m
  double precision              :: dm
  double precision              :: dgami

! Output variables

  double precision,intent(inout):: Fm(0:maxm)

  if(t == 0d0) then 
    do m=0,maxm
      dm = dble(m)
      Fm(m) = 1d0/(2d0*dm+1d0)
     enddo
   else
     do m=0,maxm
       dm = dble(m)
       Fm(m) = dgami(dm+0.5d0,t)/(2d0*t**(dm+0.5d0))
     enddo
   endif

end subroutine CalcBoysF
!******************************************************************************

!******************************************************************************
subroutine CalcOm(maxm,ExpPQi,NormPQSq,Om)

! Comute the 0^m: (00|00)^m

  implicit none

! Input variables

  integer,intent(in)            :: maxm
  double precision,intent(in)   :: ExpPQi,NormPQSq

! Local variables 

  integer                       :: m
  double precision              :: pi,dm,t
  double precision,allocatable  :: Fm(:)

! Output variables

  double precision,intent(inout):: Om (0:maxm)

  allocate(Fm(0:maxm))

  pi = 4d0*atan(1d0)

! Campute generalized Boys functions

  t = NormPQSq/ExpPQi
  call CalcBoysF(maxm,t,Fm)

! Compute (00|00)^m

  do m=0,maxm
    dm =dble(m)
    Om(m) = (2d0/sqrt(pi))*(-1d0)**dm*(1d0/ExpPQi)**(dm+0.5d0)*Fm(m)
  enddo

  deallocate(Fm)

end subroutine CalcOm

!******************************************************************************
recursive function GVRR(m,AngMomA,AngMomC,TotAngMomA,TotAngMomC,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ) &
                   result(Gac)

! Compute two-electron integrals over Gaussian geminals

  implicit none

! Input variables 

  integer,intent(in)            :: m
  integer,intent(in)            :: AngMomA(3),AngMomC(3)
  integer,intent(in)            :: TotAngMomA,TotAngMomC
  integer,intent(in)            :: maxm
  double precision,intent(in)   :: Om(0:maxm)
  double precision,intent(in)   :: ExpB,ExpD,ExpPi,ExpQi
  double precision,intent(in)   :: CenterAB(3),CenterCD(3),CenterPQ(3)

! Local variables

  logical                       :: NegAngMomA,NegAngMomC
  integer                       :: xyz,am(3),amm(3),cm(3),cmm(3)
  integer                       :: i

! Output variables

  double precision              :: Gac
  
!------------------------------------------------------------------------
! Termination condition
!------------------------------------------------------------------------
  NegAngMomA = AngMomA(1) < 0 .or. AngMomA(2) < 0 .or. AngMomA(3) < 0
  NegAngMomC = AngMomC(1) < 0 .or. AngMomC(2) < 0 .or. AngMomC(3) < 0

  if(NegAngMomA .or. NegAngMomC) then

    Gac = 0d0

  else
!------------------------------------------------------------------------
! Fundamental integral: (00|00)^m
!------------------------------------------------------------------------
    if(TotAngMomA == 0 .and. TotAngMomC == 0) then

      Gac = Om(m)

    else
!------------------------------------------------------------------------
! 1st vertical recurrence relation (4 terms): (a+0|00)^m
!------------------------------------------------------------------------
      if(TotAngMomC == 0) then
        do i=1,3
          am(i)  = AngMomA(i)
          amm(i) = AngMomA(i)
        enddo
! Loop over cartesian directions
        if    (AngMomA(1) > 0) then
          xyz = 1
        elseif(AngMomA(2) > 0) then
          xyz = 2
        elseif(AngMomA(3) > 0) then
          xyz = 3
        endif
! End loop over cartesian directions
        am(xyz)  = am(xyz)  - 1
        amm(xyz) = amm(xyz) - 2
        Gac = - ExpB*ExpPi*CenterAB(xyz)*                                                                                         &
                GVRR(m,am,AngMomC,TotAngMomA-1,TotAngMomC,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)               &
              + ExpPi*CenterPQ(xyz)*                                                                                              &
                GVRR(m+1,am,AngMomC,TotAngMomA-1,TotAngMomC,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)             &
              + dble(AngMomA(xyz)-1)*ExpPi/2d0*(                                                                                  &
                    GVRR(m,amm,AngMomC,TotAngMomA-2,TotAngMomC,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)          &
                  + ExpPi*GVRR(m+1,amm,AngMomC,TotAngMomA-2,TotAngMomC,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ))
!------------------------------------------------------------------------
! 2nd vertical recurrence relation (5 terms): (a0|c+0)^m
!------------------------------------------------------------------------
      else
        do i=1,3
          am(i)  = AngMomA(i)
          cm(i)  = AngMomC(i)
          cmm(i) = AngMomC(i)
        enddo
! Loop over cartesian directions
        if    (AngMomC(1) > 0) then
          xyz = 1
        elseif(AngMomC(2) > 0) then
          xyz = 2
        elseif(AngMomC(3) > 0) then
          xyz = 3
        endif
! End loop over cartesian directions
        am(xyz)  = am(xyz)  - 1
        cm(xyz)  = cm(xyz)  - 1
        cmm(xyz) = cmm(xyz) - 2
        Gac = - ExpD*ExpQi*CenterCD(xyz)*                                                                                         &
                GVRR(m,AngMomA,cm,TotAngMomA,TotAngMomC-1,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)               &
              - ExpQi*CenterPQ(xyz)*                                                                                              &
                GVRR(m+1,AngMomA,cm,TotAngMomA,TotAngMomC-1,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)             &
              + dble(AngMomC(xyz)-1)*ExpQi/2d0*(                                                                                  &
                    GVRR(m,AngMomA,cmm,TotAngMomA,TotAngMomC-2,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)          &
                  + ExpQi*GVRR(m+1,AngMomA,cmm,TotAngMomA,TotAngMomC-2,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)) &
              - dble(AngMomA(xyz))*ExpPi*ExpQi/2d0*                                                                               &
                GVRR(m+1,am,cm,TotAngMomA-1,TotAngMomC-1,maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)
      endif
    endif
  endif
   
end function GVRR
!*****************************************************************************

!*****************************************************************************
recursive function GHRR(m,AngMomA,AngMomB,AngMomC,AngMomD,TotAngMomA,TotAngMomB,TotAngMomC,TotAngMomD, & 
                        maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)                      &
                   result(Gabcd)

! Compute two-electron integrals over Gaussian geminals

  implicit none

! Input variables 

  integer,intent(in)            :: m
  integer,intent(in)            :: AngMomA(3),AngMomB(3),AngMomC(3),AngMomD(3)
  integer,intent(in)            :: TotAngMomA,TotAngMomB,TotAngMomC,TotAngMomD
  integer,intent(in)            :: maxm
  double precision,intent(in)   :: Om(0:maxm)
  double precision,intent(in)   :: ExpB,ExpD,ExpPi,ExpQi
  double precision,intent(in)   :: CenterAB(3),CenterCD(3),CenterPQ(3)

! Local variables

  logical                       :: NegAngMomB,NegAngMomD
  integer                       :: xyz,ap(3),bm(3),cp(3),dm(3)
  integer                       :: i
  double precision              :: GVRR

! Output variables

  double precision              :: Gabcd


!------------------------------------------------------------------------
! Termination condition
!------------------------------------------------------------------------
  NegAngMomB = AngMomB(1) < 0 .or. AngMomB(2) < 0 .or. AngMomB(3) < 0
  NegAngMomD = AngMomD(1) < 0 .or. AngMomD(2) < 0 .or. AngMomD(3) < 0

  if(NegAngMomB .or. NegAngMomD) then

    Gabcd = 0d0

  else
!------------------------------------------------------------------------
! 1st and 2nd vertical recurrence relations: (a0|c0)
!------------------------------------------------------------------------
    if(TotAngMomB == 0 .and. TotAngMomD == 0) then
      Gabcd = GVRR(m,AngMomA,AngMomC,TotAngMomA,TotAngMomC, &
                   maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)
    else
!------------------------------------------------------------------------
! 1st horizontal recurrence relation (2 terms): (ab+|c0)
!------------------------------------------------------------------------
      if(TotAngMomD == 0) then
        do i=1,3
          ap(i) = AngMomA(i)
          bm(i) = AngMomB(i)
        enddo
! Loop over cartesian directions 
        if    (AngMomB(1) > 0) then
          xyz = 1
        elseif(AngMomB(2) > 0) then
          xyz = 2
        elseif(AngMomB(3) > 0) then
          xyz = 3
        endif
! End loop over cartesian directions
        ap(xyz) = ap(xyz) + 1
        bm(xyz) = bm(xyz) - 1
        Gabcd = GHRR(m,ap,bm,AngMomC,AngMomD,TotAngMomA+1,TotAngMomB-1,TotAngMomC,TotAngMomD,      &
                     maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)                     &
              + CenterAB(xyz)*                                                                     &
                  GHRR(m,AngMomA,bm,AngMomC,AngMomD,TotAngMomA,TotAngMomB-1,TotAngMomC,TotAngMomD, &
                       maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ) 
!------------------------------------------------------------------------
! 2nd horizontal recurrence relation (2 terms): (ab|cd+)
!------------------------------------------------------------------------
      else
        do i=1,3
          cp(i) = AngMomC(i)
          dm(i) = AngMomD(i)
        enddo
! Loop over cartesian directions
        if    (AngMomD(1) > 0) then
          xyz = 1
        elseif(AngMomD(2) > 0) then
          xyz = 2
        elseif(AngMomD(3) > 0) then
          xyz = 3
        endif
! End loop over cartesian directions
        cp(xyz) = cp(xyz) + 1
        dm(xyz) = dm(xyz) - 1
        Gabcd = GHRR(m,AngMomA,AngMomB,cp,dm,TotAngMomA,TotAngMomB,TotAngMomC+1,TotAngMomD-1,      &
                     maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ)                     &
              + CenterCD(xyz)*                                                                     & 
                  GHRR(m,AngMomA,AngMomB,AngMomC,dm,TotAngMomA,TotAngMomB,TotAngMomC,TotAngMomD-1, &
                       maxm,Om,ExpB,ExpD,ExpPi,ExpQi,CenterAB,CenterCD,CenterPQ) 
      endif
    endif
  endif
   
end function GHRR
!*****************************************************************************

!     DOUBLE PRECISION FUNCTION DGAMI (A, X)
!
!***BEGIN PROLOGUE  DGAMI
!***PURPOSE  Evaluate the incomplete Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (GAMI-S, DGAMI-D)
!***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate the incomplete gamma function defined by
!
! DGAMI = integral from T = 0 to X of EXP(-T) * T**(A-1.0) .
!
! DGAMI is evaluated for positive values of A and non-negative values
! of X.  A slight deterioration of 2 or 3 digits accuracy will occur
! when DGAMI is very large or very small, because logarithmic variables
! are used.  The function and both arguments are double precision.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DGAMIT, DLNGAM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DGAMI

! $Id$
!DECK D9GMIT
      DOUBLE PRECISION FUNCTION D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
!***BEGIN PROLOGUE  D9GMIT
!***SUBSIDIARY
!***PURPOSE  Compute Tricomi's incomplete Gamma function for small
!            arguments.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9GMIT-S, D9GMIT-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
!             SPECIAL FUNCTIONS, TRICOMI
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute Tricomi's incomplete gamma function for small X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DLNGAM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9GMIT


      DOUBLE PRECISION A, X, ALGAP1, SGNGAM, ALX, AE, AEPS, ALGS, ALG2, &
       BOT, EPS, FK, S, SGNG2, T, TE, D1MACH, DLNGAM
      LOGICAL FIRST
      SAVE EPS, BOT, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9GMIT
      IF (FIRST) THEN
         EPS = 0.5D0*D1MACH(3)
         BOT = LOG (D1MACH(1))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LE. 0.D0) CALL XERMSG ('SLATEC', 'D9GMIT', &
        'X SHOULD BE GT 0', 1, 2)
!
      MA = A + 0.5D0
      IF (A.LT.0.D0) MA = A - 0.5D0
      AEPS = A - MA
!
      AE = A
      IF (A.LT.(-0.5D0)) AE = AEPS
!
      T = 1.D0
      TE = AE
      S = T
      DO 20 K=1,200
        FK = K
        TE = -X*TE/FK
        T = TE/(AE+FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'D9GMIT', &
        'NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES', 2, 2)
!
 30   IF (A.GE.(-0.5D0)) ALGS = -ALGAP1 + LOG(S)
      IF (A.GE.(-0.5D0)) GO TO 60
!
      ALGS = -DLNGAM(1.D0+AEPS) + LOG(S)
      S = 1.0D0
      M = -MA - 1
      IF (M.EQ.0) GO TO 50
      T = 1.0D0
      DO 40 K=1,M
        T = X*T/(AEPS-(M+1-K))
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40   CONTINUE
!
 50   D9GMIT = 0.0D0
      ALGS = -MA*LOG(X) + ALGS
      IF (S.EQ.0.D0 .OR. AEPS.EQ.0.D0) GO TO 60
!
      SGNG2 = SGNGAM * SIGN (1.0D0, S)
      ALG2 = -X - ALGAP1 + LOG(ABS(S))
!
      IF (ALG2.GT.BOT) D9GMIT = SGNG2 * EXP(ALG2)
      IF (ALGS.GT.BOT) D9GMIT = D9GMIT + EXP(ALGS)
      RETURN
!
 60   D9GMIT = EXP (ALGS)
      RETURN
!
      END
!DECK D9LGIC
      DOUBLE PRECISION FUNCTION D9LGIC (A, X, ALX)
!***BEGIN PROLOGUE  D9LGIC
!***SUBSIDIARY
!***PURPOSE  Compute the log complementary incomplete Gamma function
!            for large X and for A .LE. X.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9LGIC-S, D9LGIC-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
!             LOGARITHM, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log complementary incomplete gamma function for large X
! and for A .LE. X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9LGIC
      DOUBLE PRECISION A, X, ALX, EPS, FK, P, R, S, T, XMA, XPA, D1MACH
      SAVE EPS
      DATA EPS / 0.D0 /
!***FIRST EXECUTABLE STATEMENT  D9LGIC
      IF (EPS.EQ.0.D0) EPS = 0.5D0*D1MACH(3)
!
      XPA = X + 1.0D0 - A
      XMA = X - 1.D0 - A
!
      R = 0.D0
      P = 1.D0
      S = P
      DO 10 K=1,300
        FK = K
        T = FK*(A-FK)*(1.D0+R)
        R = -T/((XMA+2.D0*FK)*(XPA+2.D0*FK)+T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'D9LGIC', &
        'NO CONVERGENCE IN 300 TERMS OF CONTINUED FRACTION', 1, 2)
!
 20   D9LGIC = A*ALX - X + LOG(S/XPA)
!
      RETURN
      END
!DECK D9LGIT
      DOUBLE PRECISION FUNCTION D9LGIT (A, X, ALGAP1)
!***BEGIN PROLOGUE  D9LGIT
!***SUBSIDIARY
!***PURPOSE  Compute the logarithm of Tricomi's incomplete Gamma
!            function with Perron's continued fraction for large X and
!            A .GE. X.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9LGIT-S, D9LGIT-D)
!***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
!             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log of Tricomi's incomplete gamma function with Perron's
! continued fraction for large X and for A .GE. X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9LGIT
      DOUBLE PRECISION A, X, ALGAP1, AX, A1X, EPS, FK, HSTAR, P, R, S, &
       SQEPS, T, D1MACH
      LOGICAL FIRST
      SAVE EPS, SQEPS, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9LGIT
      IF (FIRST) THEN
         EPS = 0.5D0*D1MACH(3)
         SQEPS = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LE. 0.D0 .OR. A .LT. X) CALL XERMSG ('SLATEC', 'D9LGIT', &
        'X SHOULD BE GT 0.0 AND LE A', 2, 2)
!
      AX = A + X
      A1X = AX + 1.0D0
      R = 0.D0
      P = 1.D0
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.D0+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'D9LGIT', &
        'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 3, 2)
!
 30   HSTAR = 1.0D0 - X*S/A1X
      IF (HSTAR .LT. SQEPS) CALL XERMSG ('SLATEC', 'D9LGIT', &
        'RESULT LESS THAN HALF PRECISION', 1, 1)
!
      D9LGIT = -X - ALGAP1 - LOG(HSTAR)
      RETURN
!
      END
!DECK D9LGMC
      DOUBLE PRECISION FUNCTION D9LGMC (X)
!***BEGIN PROLOGUE  D9LGMC
!***SUBSIDIARY
!***PURPOSE  Compute the log Gamma correction factor so that
!            LOG(DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X
!            + D9LGMC(X).
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
!             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log gamma correction factor for X .GE. 10. so that
! LOG (DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X)
!
! Series for ALGM       on the interval  0.          to  1.00000E-02
!                                        with weighted error   1.28E-31
!                                         log weighted error  30.89
!                               significant figures required  29.81
!                                    decimal places required  31.48
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9LGMC
      DOUBLE PRECISION X, ALGMCS(15), XBIG, XMAX, DCSEVL, D1MACH
      LOGICAL FIRST
      SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
      DATA ALGMCS(  1) / +.1666389480451863247205729650822D+0      /
      DATA ALGMCS(  2) / -.1384948176067563840732986059135D-4      /
      DATA ALGMCS(  3) / +.9810825646924729426157171547487D-8      /
      DATA ALGMCS(  4) / -.1809129475572494194263306266719D-10     /
      DATA ALGMCS(  5) / +.6221098041892605227126015543416D-13     /
      DATA ALGMCS(  6) / -.3399615005417721944303330599666D-15     /
      DATA ALGMCS(  7) / +.2683181998482698748957538846666D-17     /
      DATA ALGMCS(  8) / -.2868042435334643284144622399999D-19     /
      DATA ALGMCS(  9) / +.3962837061046434803679306666666D-21     /
      DATA ALGMCS( 10) / -.6831888753985766870111999999999D-23     /
      DATA ALGMCS( 11) / +.1429227355942498147573333333333D-24     /
      DATA ALGMCS( 12) / -.3547598158101070547199999999999D-26     /
      DATA ALGMCS( 13) / +.1025680058010470912000000000000D-27     /
      DATA ALGMCS( 14) / -.3401102254316748799999999999999D-29     /
      DATA ALGMCS( 15) / +.1276642195630062933333333333333D-30     /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9LGMC
      IF (FIRST) THEN
         NALGM = INITDS (ALGMCS, 15, REAL(D1MACH(3)) )
         XBIG = 1.0D0/SQRT(D1MACH(3))
         XMAX = EXP (MIN(LOG(D1MACH(2)/12.D0), -LOG(12.D0*D1MACH(1))))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LT. 10.D0) CALL XERMSG ('SLATEC', 'D9LGMC', &
        'X MUST BE GE 10', 1, 2)
      IF (X.GE.XMAX) GO TO 20
!
      D9LGMC = 1.D0/(12.D0*X)
      IF (X.LT.XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS, &
       NALGM) / X
      RETURN
!
 20   D9LGMC = 0.D0
      CALL XERMSG ('SLATEC', 'D9LGMC', 'X SO BIG D9LGMC UNDERFLOWS', 2, &
        1)
      RETURN
!
      END
!DECK DCSEVL
      DOUBLE PRECISION FUNCTION DCSEVL (X, CS, N)
!***BEGIN PROLOGUE  DCSEVL
!***PURPOSE  Evaluate a Chebyshev series.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
!***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!  Evaluate the N-term Chebyshev series CS at X.  Adapted from
!  a method presented in the paper by Broucke referenced below.
!
!       Input Arguments --
!  X    value at which the series is to be evaluated.
!  CS   array of N terms of a Chebyshev series.  In evaluating
!       CS, only half the first coefficient is summed.
!  N    number of terms in array CS.
!
!***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
!                 Chebyshev series, Algorithm 446, Communications of
!                 the A.C.M. 16, (1973) pp. 254-256.
!               L. Fox and I. B. Parker, Chebyshev Polynomials in
!                 Numerical Analysis, Oxford University Press, 1968,
!                 page 56.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900329  Prologued revised extensively and code rewritten to allow
!           X to be slightly outside interval (-1,+1).  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCSEVL
      DOUBLE PRECISION B0, B1, B2, CS(*), ONEPL, TWOX, X, D1MACH
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DCSEVL
      IF (FIRST) ONEPL = 1.0D0 + D1MACH(4)
      FIRST = .FALSE.
      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'DCSEVL', &
        'NUMBER OF TERMS .LE. 0', 2, 2)
      IF (N .GT. 1000) CALL XERMSG ('SLATEC', 'DCSEVL', &
        'NUMBER OF TERMS .GT. 1000', 3, 2)
      IF (ABS(X) .GT. ONEPL) CALL XERMSG ('SLATEC', 'DCSEVL', &
        'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
!
      B1 = 0.0D0
      B0 = 0.0D0
      TWOX = 2.0D0*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
!
      DCSEVL = 0.5D0*(B0-B2)
!
      RETURN
      END
!DECK DGAMI
      DOUBLE PRECISION FUNCTION DGAMI (A, X)
!***BEGIN PROLOGUE  DGAMI
!***PURPOSE  Evaluate the incomplete Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (GAMI-S, DGAMI-D)
!***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate the incomplete gamma function defined by
!
! DGAMI = integral from T = 0 to X of EXP(-T) * T**(A-1.0) .
!
! DGAMI is evaluated for positive values of A and non-negative values
! of X.  A slight deterioration of 2 or 3 digits accuracy will occur
! when DGAMI is very large or very small, because logarithmic variables
! are used.  The function and both arguments are double precision.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DGAMIT, DLNGAM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DGAMI
      DOUBLE PRECISION A, X, FACTOR, DLNGAM, DGAMIT
!***FIRST EXECUTABLE STATEMENT  DGAMI
      IF (A .LE. 0.D0) CALL XERMSG ('SLATEC', 'DGAMI', &
        'A MUST BE GT ZERO', 1, 2)
      IF (X .LT. 0.D0) CALL XERMSG ('SLATEC', 'DGAMI', &
        'X MUST BE GE ZERO', 2, 2)
!
      DGAMI = 0.D0
      IF (X.EQ.0.0D0) RETURN
!
! THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
      FACTOR = EXP (DLNGAM(A) + A*LOG(X))
!
      DGAMI = FACTOR * DGAMIT (A, X)
!
      RETURN
      END
!DECK DGAMIT
      DOUBLE PRECISION FUNCTION DGAMIT (A, X)
!***BEGIN PROLOGUE  DGAMIT
!***PURPOSE  Calculate Tricomi's form of the incomplete Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (GAMIT-S, DGAMIT-D)
!***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
!             SPECIAL FUNCTIONS, TRICOMI
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!   Evaluate Tricomi's incomplete Gamma function defined by
!
!   DGAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) *
!              T**(A-1.)
!
!   for A .GT. 0.0 and by analytic continuation for A .LE. 0.0.
!   GAMMA(X) is the complete gamma function of X.
!
!   DGAMIT is evaluated for arbitrary real values of A and for non-
!   negative values of X (even though DGAMIT is defined for X .LT.
!   0.0), except that for X = 0 and A .LE. 0.0, DGAMIT is infinite,
!   which is a fatal error.
!
!   The function and both arguments are DOUBLE PRECISION.
!
!   A slight deterioration of 2 or 3 digits accuracy will occur when
!   DGAMIT is very large or very small in absolute value, because log-
!   arithmic variables are used.  Also, if the parameter  A  is very
!   close to a negative integer (but not a negative integer), there is
!   a loss of accuracy, which is reported if the result is less than
!   half machine precision.
!
!***REFERENCES  W. Gautschi, A computational procedure for incomplete
!                 gamma functions, ACM Transactions on Mathematical
!                 Software 5, 4 (December 1979), pp. 466-481.
!               W. Gautschi, Incomplete gamma functions, Algorithm 542,
!                 ACM Transactions on Mathematical Software 5, 4
!                 (December 1979), pp. 482-489.
!***ROUTINES CALLED  D1MACH, D9GMIT, D9LGIC, D9LGIT, DGAMR, DLGAMS,
!                    DLNGAM, XERCLR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!***END PROLOGUE  DGAMIT
      DOUBLE PRECISION A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNG, ALX, &
       BOT, H, SGA, SGNGAM, SQEPS, T, D1MACH, DGAMR, D9GMIT, D9LGIT, &
       DLNGAM, D9LGIC
      LOGICAL FIRST
      SAVE ALNEPS, SQEPS, BOT, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DGAMIT
      IF (FIRST) THEN
         ALNEPS = -LOG (D1MACH(3))
         SQEPS = SQRT(D1MACH(4))
         BOT = LOG (D1MACH(1))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LT. 0.D0) CALL XERMSG ('SLATEC', 'DGAMIT', 'X IS NEGATIVE' &
        , 2, 2)
!
      IF (X.NE.0.D0) ALX = LOG (X)
      SGA = 1.0D0
      IF (A.NE.0.D0) SGA = SIGN (1.0D0, A)
      AINTA = AINT (A + 0.5D0*SGA)
      AEPS = A - AINTA
!
      IF (X.GT.0.D0) GO TO 20
      DGAMIT = 0.0D0
      IF (AINTA.GT.0.D0 .OR. AEPS.NE.0.D0) DGAMIT = DGAMR(A+1.0D0)
      RETURN
!
 20   IF (X.GT.1.D0) GO TO 30
      IF (A.GE.(-0.5D0) .OR. AEPS.NE.0.D0) CALL DLGAMS (A+1.0D0, ALGAP1, &
       SGNGAM)
      DGAMIT = D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      RETURN
!
 30   IF (A.LT.X) GO TO 40
      T = D9LGIT (A, X, DLNGAM(A+1.0D0))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = EXP (T)
      RETURN
!
 40   ALNG = D9LGIC (A, X, ALX)
!
! EVALUATE DGAMIT IN TERMS OF LOG (DGAMIC (A, X))
!
      H = 1.0D0
      IF (AEPS.EQ.0.D0 .AND. AINTA.LE.0.D0) GO TO 50
!
      CALL DLGAMS (A+1.0D0, ALGAP1, SGNGAM)
      T = LOG (ABS(A)) + ALNG - ALGAP1
      IF (T.GT.ALNEPS) GO TO 60
!
      IF (T.GT.(-ALNEPS)) H = 1.0D0 - SGA * SGNGAM * EXP(T)
      IF (ABS(H).GT.SQEPS) GO TO 50
!
      CALL XERCLR
      CALL XERMSG ('SLATEC', 'DGAMIT', 'RESULT LT HALF PRECISION', 1, &
        1)
!
 50   T = -A*ALX + LOG(ABS(H))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = SIGN (EXP(T), H)
      RETURN
!
 60   T = T - A*ALX
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = -SGA * SGNGAM * EXP(T)
      RETURN
!
      END
!DECK DGAMLM
      SUBROUTINE DGAMLM (XMIN, XMAX)
!***BEGIN PROLOGUE  DGAMLM
!***PURPOSE  Compute the minimum and maximum bounds for the argument in
!            the Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A, R2
!***TYPE      DOUBLE PRECISION (GAMLIM-S, DGAMLM-D)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Calculate the minimum and maximum legal bounds for X in gamma(X).
! XMIN and XMAX are not the only bounds, but they are the only non-
! trivial ones to calculate.
!
!             Output Arguments --
! XMIN   double precision minimum legal value of X in gamma(X).  Any
!        smaller value of X might result in underflow.
! XMAX   double precision maximum legal value of X in gamma(X).  Any
!        larger value of X might cause overflow.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DGAMLM
      DOUBLE PRECISION XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD, D1MACH
!***FIRST EXECUTABLE STATEMENT  DGAMLM
      ALNSML = LOG(D1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML) &
         / (XMIN*XLN+0.5D0)
        IF (ABS(XMIN-XOLD).LT.0.005D0) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMIN', 1, 2)
!
 20   XMIN = -XMIN + 0.01D0
!
      ALNBIG = LOG (D1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG) &
         / (XMAX*XLN-0.5D0)
        IF (ABS(XMAX-XOLD).LT.0.005D0) GO TO 40
 30   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMAX', 2, 2)
!
 40   XMAX = XMAX - 0.01D0
      XMIN = MAX (XMIN, -XMAX+1.D0)
!
      RETURN
      END
!DECK UTIL_DGAMMA
      DOUBLE PRECISION FUNCTION UTIL_DGAMMA (X)
!***BEGIN PROLOGUE  UTIL_DGAMMA
!***PURPOSE  Compute the complete Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! UTIL_DGAMMA(X) calculates the double precision complete Gamma function
! for double precision argument X.
!
! Series for GAM        on the interval  0.          to  1.00000E+00
!                                        with weighted error   5.79E-32
!                                         log weighted error  31.24
!                               significant figures required  30.00
!                                    decimal places required  32.05
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9LGMC, DCSEVL, DGAMLM, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920618  Removed space from variable name.  (RWC, WRB)
!***END PROLOGUE  UTIL_DGAMMA
      DOUBLE PRECISION X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX, &
       XMIN, Y, D9LGMC, DCSEVL, D1MACH
      LOGICAL FIRST
!
      SAVE GAMCS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL, FIRST
      DATA GAMCS(  1) / +.8571195590989331421920062399942D-2      /
      DATA GAMCS(  2) / +.4415381324841006757191315771652D-2      /
      DATA GAMCS(  3) / +.5685043681599363378632664588789D-1      /
      DATA GAMCS(  4) / -.4219835396418560501012500186624D-2      /
      DATA GAMCS(  5) / +.1326808181212460220584006796352D-2      /
      DATA GAMCS(  6) / -.1893024529798880432523947023886D-3      /
      DATA GAMCS(  7) / +.3606925327441245256578082217225D-4      /
      DATA GAMCS(  8) / -.6056761904460864218485548290365D-5      /
      DATA GAMCS(  9) / +.1055829546302283344731823509093D-5      /
      DATA GAMCS( 10) / -.1811967365542384048291855891166D-6      /
      DATA GAMCS( 11) / +.3117724964715322277790254593169D-7      /
      DATA GAMCS( 12) / -.5354219639019687140874081024347D-8      /
      DATA GAMCS( 13) / +.9193275519859588946887786825940D-9      /
      DATA GAMCS( 14) / -.1577941280288339761767423273953D-9      /
      DATA GAMCS( 15) / +.2707980622934954543266540433089D-10     /
      DATA GAMCS( 16) / -.4646818653825730144081661058933D-11     /
      DATA GAMCS( 17) / +.7973350192007419656460767175359D-12     /
      DATA GAMCS( 18) / -.1368078209830916025799499172309D-12     /
      DATA GAMCS( 19) / +.2347319486563800657233471771688D-13     /
      DATA GAMCS( 20) / -.4027432614949066932766570534699D-14     /
      DATA GAMCS( 21) / +.6910051747372100912138336975257D-15     /
      DATA GAMCS( 22) / -.1185584500221992907052387126192D-15     /
      DATA GAMCS( 23) / +.2034148542496373955201026051932D-16     /
      DATA GAMCS( 24) / -.3490054341717405849274012949108D-17     /
      DATA GAMCS( 25) / +.5987993856485305567135051066026D-18     /
      DATA GAMCS( 26) / -.1027378057872228074490069778431D-18     /
      DATA GAMCS( 27) / +.1762702816060529824942759660748D-19     /
      DATA GAMCS( 28) / -.3024320653735306260958772112042D-20     /
      DATA GAMCS( 29) / +.5188914660218397839717833550506D-21     /
      DATA GAMCS( 30) / -.8902770842456576692449251601066D-22     /
      DATA GAMCS( 31) / +.1527474068493342602274596891306D-22     /
      DATA GAMCS( 32) / -.2620731256187362900257328332799D-23     /
      DATA GAMCS( 33) / +.4496464047830538670331046570666D-24     /
      DATA GAMCS( 34) / -.7714712731336877911703901525333D-25     /
      DATA GAMCS( 35) / +.1323635453126044036486572714666D-25     /
      DATA GAMCS( 36) / -.2270999412942928816702313813333D-26     /
      DATA GAMCS( 37) / +.3896418998003991449320816639999D-27     /
      DATA GAMCS( 38) / -.6685198115125953327792127999999D-28     /
      DATA GAMCS( 39) / +.1146998663140024384347613866666D-28     /
      DATA GAMCS( 40) / -.1967938586345134677295103999999D-29     /
      DATA GAMCS( 41) / +.3376448816585338090334890666666D-30     /
      DATA GAMCS( 42) / -.5793070335782135784625493333333D-31     /
      DATA PI / 3.14159265358979323846264338327950D0 /
      DATA SQ2PIL / 0.91893853320467274178032973640562D0 /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  UTIL_DGAMMA
      IF (FIRST) THEN
         NGAM = INITDS (GAMCS, 42, 0.1*REAL(D1MACH(3)) )
!
         CALL DGAMLM (XMIN, XMAX)
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
!
      Y = ABS(X)
      IF (Y.GT.10.D0) GO TO 50
!
! COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
! GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
!
      N = X
      IF (X.LT.0.D0) N = N - 1
      Y = X - N
      N = N - 1
      UTIL_DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
!
      IF (N.GT.0) GO TO 30
!
! COMPUTE GAMMA(X) FOR X .LT. 1.0
!
      N = -N
      IF (X .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA', 'X IS 0', 4, 2)
      IF (X .LT. 0.0 .AND. X+N-2 .EQ. 0.D0) CALL XERMSG ('SLATEC', &
        'DGAMMA', 'X IS A NEGATIVE INTEGER', 4, 2)
      IF (X .LT. (-0.5D0) .AND. ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) &
        CALL XERMSG ('SLATEC', 'DGAMMA', &
        'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', &
        1, 1)
!
      DO 20 I=1,N
        UTIL_DGAMMA = UTIL_DGAMMA/(X+I-1 )
 20   CONTINUE
      RETURN
!
! GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
!
 30   DO 40 I=1,N
        UTIL_DGAMMA = (Y+I) * UTIL_DGAMMA
 40   CONTINUE
      RETURN
!
! GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
 50   IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DGAMMA', &
        'X SO BIG GAMMA OVERFLOWS', 3, 2)
!
      UTIL_DGAMMA = 0.D0
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DGAMMA', &
        'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
      IF (X.LT.XMIN) RETURN
!
      UTIL_DGAMMA = EXP ((Y-0.5D0)*LOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
      IF (X.GT.0.D0) RETURN
!
      IF (ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) CALL XERMSG ('SLATEC', &
        'DGAMMA', &
        'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
      SINPIY = SIN (PI*Y)
      IF (SINPIY .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA', &
        'X IS A NEGATIVE INTEGER', 4, 2)
!
      UTIL_DGAMMA = -PI/(Y*SINPIY*UTIL_DGAMMA)
!
      RETURN
      END
!DECK DGAMR
      DOUBLE PRECISION FUNCTION DGAMR (X)
!***BEGIN PROLOGUE  DGAMR
!***PURPOSE  Compute the reciprocal of the Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C)
!***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DGAMR(X) calculates the double precision reciprocal of the
! complete Gamma function for double precision argument X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DGAMMA, DLGAMS, XERCLR, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  DGAMR
      DOUBLE PRECISION X, ALNGX, SGNGX, UTIL_DGAMMA
      EXTERNAL UTIL_DGAMMA
!***FIRST EXECUTABLE STATEMENT  DGAMR
      DGAMR = 0.0D0
      IF (X.LE.0.0D0 .AND. AINT(X).EQ.X) RETURN
!
      CALL XGETF (IROLD)
      CALL XSETF (1)
      IF (ABS(X).GT.10.0D0) GO TO 10
      DGAMR = 1.0D0/UTIL_DGAMMA(X)
      CALL XERCLR
      CALL XSETF (IROLD)
      RETURN
!
 10   CALL DLGAMS (X, ALNGX, SGNGX)
      CALL XERCLR
      CALL XSETF (IROLD)
      DGAMR = SGNGX * EXP(-ALNGX)
      RETURN
!
      END
!DECK DLGAMS
      SUBROUTINE DLGAMS (X, DLGAM, SGNGAM)
!***BEGIN PROLOGUE  DLGAMS
!***PURPOSE  Compute the logarithm of the absolute value of the Gamma
!            function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (ALGAMS-S, DLGAMS-D)
!***KEYWORDS  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION,
!             FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DLGAMS(X,DLGAM,SGNGAM) calculates the double precision natural
! logarithm of the absolute value of the Gamma function for
! double precision argument X and stores the result in double
! precision argument DLGAM.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DLNGAM
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DLGAMS
      DOUBLE PRECISION X, DLGAM, SGNGAM, DLNGAM
!***FIRST EXECUTABLE STATEMENT  DLGAMS
      DLGAM = DLNGAM(X)
      SGNGAM = 1.0D0
      IF (X.GT.0.D0) RETURN
!
      INT = MOD (-AINT(X), 2.0D0) + 0.1D0
      IF (INT.EQ.0) SGNGAM = -1.0D0
!
      RETURN
      END
!DECK DLNGAM
      DOUBLE PRECISION FUNCTION DLNGAM (X)
!***BEGIN PROLOGUE  DLNGAM
!***PURPOSE  Compute the logarithm of the absolute value of the Gamma
!            function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C)
!***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DLNGAM(X) calculates the double precision logarithm of the
! absolute value of the Gamma function for double precision
! argument X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9LGMC, DGAMMA, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  DLNGAM
      DOUBLE PRECISION X, DXREL, PI, SINPIY, SQPI2L, SQ2PIL, XMAX, &
       Y, UTIL_DGAMMA, D9LGMC, D1MACH, TEMP
      LOGICAL FIRST
      EXTERNAL UTIL_DGAMMA
      SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL, FIRST
      DATA SQ2PIL / 0.91893853320467274178032973640562D0 /
      DATA SQPI2L / +.225791352644727432363097614947441D+0    /
      DATA PI / 3.14159265358979323846264338327950D0 /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DLNGAM
      IF (FIRST) THEN
         TEMP = 1.D0/LOG(D1MACH(2))
         XMAX = TEMP*D1MACH(2)
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
!
      Y = ABS (X)
      IF (Y.GT.10.D0) GO TO 20
!
! LOG (ABS (DGAMMA(X)) ) FOR ABS(X) .LE. 10.0
!
      DLNGAM = LOG (ABS (UTIL_DGAMMA(X)) )
      RETURN
!
! LOG ( ABS (DGAMMA(X)) ) FOR ABS(X) .GT. 10.0
!
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'DLNGAM', &
        'ABS(X) SO BIG DLNGAM OVERFLOWS', 2, 2)
!
      IF (X.GT.0.D0) DLNGAM = SQ2PIL + (X-0.5D0)*LOG(X) - X + D9LGMC(Y)
      IF (X.GT.0.D0) RETURN
!
      SINPIY = ABS (SIN(PI*Y))
      IF (SINPIY .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DLNGAM', &
        'X IS A NEGATIVE INTEGER', 3, 2)
!
      IF (ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) CALL XERMSG ('SLATEC', &
        'DLNGAM', &
        'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', &
        1, 1)
!
      DLNGAM = SQPI2L + (X-0.5D0)*LOG(Y) - X - LOG(SINPIY) - D9LGMC(Y)
      RETURN
!
      END
!DECK INITDS
      FUNCTION INITDS (OS, NOS, ETA)
!***BEGIN PROLOGUE  INITDS
!***PURPOSE  Determine the number of terms needed in an orthogonal
!            polynomial series so that it meets a specified accuracy.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
!***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
!             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!  Initialize the orthogonal series, represented by the array OS, so
!  that INITDS is the number of terms needed to insure the error is no
!  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
!  machine precision.
!
!             Input Arguments --
!   OS     double precision array of NOS coefficients in an orthogonal
!          series.
!   NOS    number of coefficients in OS.
!   ETA    single precision scalar containing requested accuracy of
!          series.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891115  Modified error message.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  INITDS
      DOUBLE PRECISION OS(*)
!***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS .LT. 1) CALL XERMSG ('SLATEC', 'INITDS', &
        'Number of coefficients is less than 1', 2, 1)
!
      ERR = 0.
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(REAL(OS(I)))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
!
   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'INITDS', &
        'Chebyshev series too short for specified accuracy', 1, 1)
      INITDS = I
!
      RETURN
      END
!DECK XERCLR
      SUBROUTINE XERCLR
!***BEGIN PROLOGUE  XERCLR
!***PURPOSE  Reset current error number to zero.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERCLR-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        This routine simply resets the current error number to zero.
!        This may be necessary in order to determine that a certain
!        error has occurred again since the last time NUMXER was
!        referenced.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERCLR
!***FIRST EXECUTABLE STATEMENT  XERCLR
      JUNK = J4SAVE(1,0,.TRUE.)
      RETURN
      END
!DECK J4SAVE
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
!***BEGIN PROLOGUE  J4SAVE
!***SUBSIDIARY
!***PURPOSE  Save or recall global variables needed by error
!            handling routines.
!***LIBRARY   SLATEC (XERROR)
!***TYPE      INTEGER (J4SAVE-I)
!***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.
!
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                = 6 Refers to the 2nd unit for error messages
!                = 7 Refers to the 3rd unit for error messages
!                = 8 Refers to the 4th unit for error messages
!                = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
!
!***SEE ALSO  XERMSG
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900205  Minor modifications to prologue.  (WRB)
!   900402  Added TYPE section.  (WRB)
!   910411  Added KEYWORDS section.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
      implicit none
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      integer nerr,level
!!!      call errquit(SUBROU//MESSG,NERR,0)
      return
      end
!DECK XSETF
      SUBROUTINE XSETF (KONTRL)
!***BEGIN PROLOGUE  XSETF
!***PURPOSE  Set the error control flag.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3A
!***TYPE      ALL (XSETF-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        XSETF sets the error control flag value to KONTRL.
!        (KONTRL is an input parameter only.)
!        The following table shows how each message is treated,
!        depending on the values of KONTRL and LEVEL.  (See XERMSG
!        for description of LEVEL.)
!
!        If KONTRL is zero or negative, no information other than the
!        message itself (including numeric values, if any) will be
!        printed.  If KONTRL is positive, introductory messages,
!        trace-backs, etc., will be printed in addition to the message.
!
!              ABS(KONTRL)
!        LEVEL        0              1              2
!        value
!          2        fatal          fatal          fatal
!
!          1     not printed      printed         fatal
!
!          0     not printed      printed        printed
!
!         -1     not printed      printed        printed
!                                  only           only
!                                  once           once
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900510  Change call to XERRWV to XERMSG.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XSETF
      CHARACTER *8 XERN1
!***FIRST EXECUTABLE STATEMENT  XSETF
      IF (ABS(KONTRL) .GT. 2) THEN
         write(*,*)'HERE!'
         WRITE (XERN1, '(I8)') KONTRL
         CALL XERMSG ('SLATEC', 'XSETF', &
           'INVALID ARGUMENT = ' // XERN1, 1, 2)
         RETURN
      ENDIF
!
      JUNK = J4SAVE(2,KONTRL,.TRUE.)
      RETURN
      END
!DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH_OLD (I)
!***BEGIN PROLOGUE  D1MACH
!***PURPOSE  Return floating point machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        D = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of D above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH( 3) = B**(-T), the smallest relative spacing.
!   D1MACH( 4) = B**(1-T), the largest relative spacing.
!   D1MACH( 5) = LOG10(B)
!
!   Assume double precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!   EMIN .LE. E .LE. EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(14) = T, the number of base-B digits.
!   I1MACH(15) = EMIN, the smallest exponent E.
!   I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890213  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   900911  Added SUN 386i constants.  (WRB)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!***END PROLOGUE  D1MACH
!
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
!
      DOUBLE PRECISION DMACH(5)
      SAVE DMACH
!
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
!
!***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'D1MACH', &
        'I OUT OF BOUNDS', 1, 2)
!
      D1MACH_OLD = DMACH(I)
      RETURN
!
      END
!DECK XGETF
      SUBROUTINE XGETF (KONTRL)
!***BEGIN PROLOGUE  XGETF
!***PURPOSE  Return the current value of the error control flag.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XGETF-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!   Abstract
!        XGETF returns the current value of the error control flag
!        in KONTRL.  See subroutine XSETF for flag value meanings.
!        (KONTRL is an output parameter only.)
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XGETF
!***FIRST EXECUTABLE STATEMENT  XGETF
      KONTRL = J4SAVE(2,0,.FALSE.)
      RETURN
      END
!
! Adapted from:
! https://github.com/certik/fortran-utils/blob/master/src/legacy/amos/d1mach.f90
!
      DOUBLE PRECISION FUNCTION D1MACH (I)
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION B, X
!***BEGIN PROLOGUE  D1MACH
!***PURPOSE  Return floating point machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (D1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH(3) = B**(-T), the smallest relative spacing.
!   D1MACH(4) = B**(1-T), the largest relative spacing.
!   D1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!   EMIN .LE. E .LE. EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by EHG)      
!***END PROLOGUE  D1MACH
!      
      X = 1.0D0
      B = RADIX(X)
!
      if (I.eq.1) then 
        D1MACH = TINY(X)               ! the smallest positive magnitude.
      else if (I.eq.2) then
        D1MACH = HUGE(X)               ! the largest magnitude.
      else if (I.eq.3) then
        D1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
      else if (I.eq.4) then
        D1MACH = B**(1-DIGITS(X))      ! the largest relative spacing.
      else if (I.eq.5) then
        D1MACH = LOG10(B)
      else
          WRITE (*, FMT = 9000)
 9000     FORMAT ('1ERROR    1 IN D1MACH - I OUT OF BOUNDS')
          STOP
      end if
!
      RETURN
      END
!! We have

!! 1/x = 1/sqrt(pi A) \int dy exp-y**2  exp(-y**2/A *x**2) exp(y**2)
!!  for any A> 0
!!
!! Using Gauss-Hermite we have
!!
!!  1/x = \sum_i  c_i exp(-gam_i x**2)
!!
!!
!!     with   c_i= w_i exp(y_i**2) 1/sqrt(pi A) gam_i=y_i**2/A
!!
double precision function fit_12(r12)
implicit double precision (a-h,o-z)

common/int_grid/npts_gr
common/real_grid/w_grid(100,100),x_grid(100,100),a_grid

pi=dacos(-1.d0)
a=a_grid

if(a.eq.0.d0)then
 fit_12=1.d0/r12
else
sum=0.d0
do i=1,npts_gr
 sum=sum +w_grid(i,npts_gr)*&
 dexp(-x_grid(i,npts_gr)**2/a**2*(r12**2-a**2))/dsqrt(pi*a**2)
enddo
fit_12=sum

if(r12.lt.0.1)fit_12=1.d0/r12

endif
end

      double precision function fit(r)
      implicit double precision(a-h,o-z)
      dimension x(10)
      n=7
      x(1)=1.5374608736368700     
      x(2)= -0.30168769375627269     
      x(3)=   2.5068485742098797E-002
      x(4)=  -9.8060522151395828E-004
      x(5)=   1.9390298869266076E-005
      x(6)=  -1.8777097417190198E-007
      x(7)=   7.0743805081128674E-010

      fit=0.d0
      do i=1,n
       new=2*(i-1)
       fit=fit+x(i)*r**new
      enddo

      end

         double precision function rinteg_1(n)
         implicit double precision (a-h,o-z)
         rinteg_1=1.d0
         if(n.eq.0)return
         rinteg_1=0.d0
         if(mod(n,2).eq.0)rinteg_1=dblefact(n-1)
         end

!! Computation of intd1d2 pi(1)pi(2) f0(1,2) x1^n(1) y1^n(2) z1^n(3) x2^n(4) y2^n(5) z2^n(6)
!!
!!  if f0=1      intd1 pi(1) x1^n(1) y1^n(2) z1^n(3)  * intd2 pi(2) x2^n(4) y2^n(5) z2^n(6)
!!          =   int dx norm exp(-0.5 x**2) x^n(1) int dx norm exp(-0.5 x**2) x^n(2) ....

!         double precision function rinteg(id,n)
!         implicit double precision (a-h,o-z)
!         dimension n(6)
!         rinteg=1.d0
!         do l=1,id
!          rinteg=rinteg*rinteg_1(n(l))
!         enddo
!         end

double precision function rinteg(id,n,ng)
implicit double precision (a-h,o-z)
dimension n(6),n_p(3),n_q(3),n_r(3),n_s(3)
dimension a(3),b(3),c(3),d(3)
double precision integral_pqrs

gamma_p=0.5d0
gamma_q=0.5d0
gamma_r=0.5d0
gamma_s=0.5d0

nA=1
nB=1
nC=1
nD=1

do l=1,3
 a(l)=0.d0
 b(l)=0.d0
 c(l)=0.d0
 d(l)=0.d0
enddo

n_p(1)=0
n_p(2)=0
n_p(3)=0
n_q(1)=n(1)
n_q(2)=n(2)
n_q(3)=n(3)
n_r(1)=0
n_r(2)=0
n_r(3)=0
n_s(1)=n(4)
n_s(2)=n(5)
n_s(3)=n(6)

pi=dacos(-1.d0)

!rnorm=(1.d0/dsqrt(2.d0*pi))**6
!rinteg=integral_pqrs(n_p,gamma_p,n_q,gamma_q,n_r,gamma_r,n_s,gamma_s)*rnorm

!rinteg=bielec_integral(gamma_p,gamma_q,gamma_r,gamma_s,a,b,c,d,n_p,n_q,n_r,n_s)*rnorm

prefactor=1.d0/(4.d0*pi)**2

rinteg=1.d0

end

      double precision function ranf()
      implicit double precision(a-h,o-z)
      logical logic
      common/gen1/fnorm
      common/gen2/n1,n2
      nn1=0
      nn2=n1
      np1=n1*2
      np2=n2*2
      ntemp=np1*8
      if(ntemp.lt.0) then
      logic=.true.
      else
      logic=.false.
      endif
      if(logic) then
      np2=np2+1
      ntemp=np1*8
      if(ntemp.lt.0) then
      np1=np1-268435456
      endif
      endif
      ntemp=np2*8
      if(ntemp.lt.0) then
      np2=np2-268435456
      endif
      ns1=np1
      ns2=nn2+np2
      ntemp=ns2*8
      if(ntemp.lt.0) then
      ns2=ns2-268435456
      endif
      nt1=np1+n1
      nt2=ns2+n2
      ntemp=nt1*8
      if(ntemp.lt.0) then
      logic=.true.
      else 
      logic=.false.
      endif
      if(logic) then
      nt2=nt2+1
      ntemp=nt1*8
      if(ntemp.lt.0) then
      nt1=nt1-268435456
      endif
      endif
      ntemp=nt2*8
      if(ntemp.lt.0) then
      nt2=nt2-268435456
      endif
      n1=nt1
      n2=nt2
      xal1=dfloat(n1)
      xal2=dfloat(n2)*268435456.d0
      xal=xal1+xal2
      ranf=xal*fnorm
      end

         double precision function fref(x1,x2)
         implicit double precision (a-h,o-z)
         dimension x1(3),x2(3)
         fref=1.d0/dsqrt( (x1(1)-x2(1))**2+ (x1(2)-x2(2))**2+(x1(3)-x2(3))**2)
         end


         double precision function rho(r)
         implicit double precision (a-h,o-z)
         common/parameters_rho/gamma_a,gamma_b,center_A(3),center_B(3)
         dimension r(3)
         rA=0.d0
         rB=0.d0
         do l=1,3
          rA=rA+(r(l)-center_A(l))**2
          rB=rB+(r(l)-center_B(l))**2
         enddo
         rho=gamma_a*dsqrt(rA)+gamma_b*dsqrt(rB)
         rho=dexp(-rho)
         end

      double precision function fkt(rci,r)
      implicit double precision(a-h,o-z)
      common/parameters_gauss_int/ng
      if(npow.eq.2)then
       fkt=dexp(-r)
       return
      endif

      argmax=-dlog(1.d-20)
      fkt=0.d0
      if(r.lt.rci)then
       arg=r/dsqrt(1.d0-r/rci)
       if(arg.gt.argmax)return
       fkt=dexp(-arg)
      endif
      end

double precision function NAI_pol_mult(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in)
! function that calculate the folowing integral :
!       int{dr} of (x-A_x)^ax (x-B_X)^bx exp(-alpha (x-A_x)^2 - beta (x-B_x)^2 ) 1/(r-R_c)

implicit none
double precision,intent(in) :: C_center(3),A_center(3),B_center(3),alpha,beta
integer :: power_A(3),power_B(3)
integer :: i,j,k,l,n_pt
double precision :: P_center(3)
double precision :: d(0:n_pt_in),pouet,coeff,rho,dist,const,pouet_2,p,p_inv,factor
double precision :: I_n_special_exact,integrate_bourrin,I_n_bibi
double precision ::  V_e_n,const_factor,dist_integral,two_pi
double precision :: accu,epsilo,rint
integer :: n_pt_in,n_pt_out,lmax

  two_pi = 2.d0 * dacos(-1.d0)
  p = alpha + beta
  p_inv = 1.d0/p
  rho = alpha * beta * p_inv
  dist = 0.d0
  dist_integral = 0.d0
  do i = 1, 3
   P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
   dist = dist+ (A_center(i) - B_center(i))*(A_center(i) - B_center(i))
   dist_integral =  dist_integral +(P_center(i) - C_center(i))*(P_center(i) - C_center(i))
  enddo
  const_factor = dist*rho
  const = p * dist_integral
  factor = dexp(-const_factor)
  coeff = two_pi * factor * p_inv
  lmax = 20

  do i = 0, n_pt_in
    d(i) = 0.d0
  enddo
  n_pt =  2 * ( (power_A(1) + power_B(1)) +(power_A(2) + power_B(2)) +(power_A(3) + power_B(3)) )
  if (n_pt == 0) then
   epsilo = 1.d0
   pouet = rint(0,const)
   NAI_pol_mult = coeff * pouet
   return
  endif

  call give_polynom_mult_center_mono_elec(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out)
  if(n_pt_out<0)then
   NAI_pol_mult = 0.d0
   return
  endif
  accu = 0.d0

! 1/r1 standard attraction integral
  epsilo = 1.d0
! sum of integrals of type : int {t,[0,1]}  exp-(rho.(P-Q)^2 * t^2) * t^i
  do i =0 ,n_pt_out,2
   if(abs(d(i)).lt.0.000000001d0)cycle
   accu = accu + d(i) * rint(i/2,const)
  enddo
  NAI_pol_mult = accu * coeff

end

subroutine give_polynom_mult_center_mono_elec(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out)
!!!! subroutine that returns the explicit polynom in term of the "t" variable of the following polynomw ::
!!!!         I_x1(a_x, d_x,p,q) * I_x1(a_y, d_y,p,q) * I_x1(a_z, d_z,p,q)
!!!! it is for the nuclear electron atraction
implicit none
integer, intent(in) :: n_pt_in
integer,intent(out) :: n_pt_out
double precision, intent(in) :: A_center(3), B_center(3),C_center(3)
double precision, intent(in) :: alpha,beta
integer, intent(in) :: power_A(3), power_B(3)
integer :: a_x,b_x,a_y,b_y,a_z,b_z
double precision :: d(0:n_pt_in)
double precision :: d1(0:n_pt_in)
double precision :: d2(0:n_pt_in)
double precision :: d3(0:n_pt_in)
double precision :: accu,  pq_inv, p10_1, p10_2, p01_1, p01_2
double precision :: p,P_center(3),rho,p_inv,p_inv_2
double precision :: R1x(0:2), B01(0:2), R1xp(0:2),R2x(0:2)
integer :: n_pt1,n_pt2,n_pt3,dim,i
integer :: n_pt_tmp
!print*,'n_pt_in = ',n_pt_in
 accu = 0.d0
!COMPTEUR irp_rdtsc1 = irp_rdtsc()
 p = alpha+beta
 p_inv = 1.d0/p
 p_inv_2 = 0.5d0/p
 do i =1, 3
  P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
 enddo
! print*,'passed the P_center'
 
 R1x(0)  = (P_center(1) - A_center(1))
 R1x(1)  = 0.d0
 R1x(2)  = -(P_center(1) - C_center(1))
 ! R1x = (P_x - A_x) - (P_x - C_x) t^2
 R1xp(0)  = (P_center(1) - B_center(1))
 R1xp(1)  = 0.d0
 R1xp(2)  =-(P_center(1) - C_center(1))
 !R1xp = (P_x - B_x) - (P_x - C_x) t^2
 R2x(0)  =  p_inv_2
 R2x(1)  = 0.d0
 R2x(2)  = -p_inv_2
 !R2x  = 0.5 / p - 0.5/p t^2
 do i = 0,n_pt_in
  d(i) = 0.d0
 enddo
 do i = 0,n_pt_in
  d1(i) = 0.d0
 enddo
 do i = 0,n_pt_in
  d2(i) = 0.d0
 enddo
 do i = 0,n_pt_in
  d3(i) = 0.d0
 enddo
 n_pt1 = n_pt_in
 n_pt2 = n_pt_in              
 n_pt3 = n_pt_in              
 a_x = power_A(1)
 b_x = power_B(1)
  call I_x1_pol_mult_mono_elec(a_x,b_x,R1x,R1xp,R2x,d1,n_pt1,n_pt_in) 
! print*,'passed the first I_x1'
  if(n_pt1<0)then
   n_pt_out = -1
   do i = 0,n_pt_in
    d(i) = 0.d0
   enddo
   return
  endif


 R1x(0)  = (P_center(2) - A_center(2))
 R1x(1)  = 0.d0
 R1x(2)  = -(P_center(2) - C_center(2))
 ! R1x = (P_x - A_x) - (P_x - C_x) t^2
 R1xp(0)  = (P_center(2) - B_center(2))
 R1xp(1)  = 0.d0
 R1xp(2)  =-(P_center(2) - C_center(2))
 !R1xp = (P_x - B_x) - (P_x - C_x) t^2
 a_y = power_A(2)
 b_y = power_B(2)
  call I_x1_pol_mult_mono_elec(a_y,b_y,R1x,R1xp,R2x,d2,n_pt2,n_pt_in) 
! print*,'passed the second I_x1'
  if(n_pt2<0)then
   n_pt_out = -1
   do i = 0,n_pt_in
    d(i) = 0.d0
   enddo
   return
  endif


 R1x(0)  = (P_center(3) - A_center(3))
 R1x(1)  = 0.d0
 R1x(2)  = -(P_center(3) - C_center(3))
 ! R1x = (P_x - A_x) - (P_x - C_x) t^2
 R1xp(0)  = (P_center(3) - B_center(3))
 R1xp(1)  = 0.d0
 R1xp(2)  =-(P_center(3) - C_center(3))
 !R2x  = 0.5 / p - 0.5/p t^2
 a_z = power_A(3)
 b_z = power_B(3)

! print*,'a_z = ',a_z
! print*,'b_z = ',b_z
  call I_x1_pol_mult_mono_elec(a_z,b_z,R1x,R1xp,R2x,d3,n_pt3,n_pt_in) 
! print*,'passed the third I_x1'
  if(n_pt3<0)then
   n_pt_out = -1
   do i = 0,n_pt_in
    d(i) = 0.d0
   enddo
   return
  endif
 n_pt_tmp = 0
 call multiply_poly_2(d1,n_pt1,d2,n_pt2,d,n_pt_tmp,n_pt_in)
 do i = 0,n_pt_tmp
  d1(i) = 0.d0
 enddo
 n_pt_out = 0
 call multiply_poly_2(d ,n_pt_tmp ,d3,n_pt3,d1,n_pt_out,n_pt_in)
 do i = 0, n_pt_out
  d(i) = d1(i)
 enddo

end

recursive subroutine I_x1_pol_mult_mono_elec(a,c,R1x,R1xp,R2x,d,nd,n_pt_in) 
!!!!  recursive function involved in the electron nucleus potential
 implicit none
 integer , intent(in) :: n_pt_in
 double precision,intent(inout) :: d(0:n_pt_in)
 integer,intent(inout) :: nd
 integer, intent(in):: a,c
 double precision, intent(in) :: R1x(0:2),R1xp(0:2),R2x(0:2)
  double precision :: X(0:n_pt_in)
  double precision :: Y(0:n_pt_in)
  integer :: nx, ix,dim,iy,ny
  dim = n_pt_in
! print*,'a,c = ',a,c
! print*,'nd_in = ',nd

  if( (a==0) .and. (c==0))then
!  print*,'coucou !'
   nd = 0
   d(0) = 1.d0
   return
  elseif( (c<0).or.(nd<0) )then
     nd = -1
     return
  else if ((a==0).and.(c.ne.0)) then
     call I_x2_pol_mult_mono_elec(c,R1x,R1xp,R2x,d,nd,n_pt_in)
!    print*,'nd 0,c',nd
  else if (a==1) then
     nx = nd
     do ix=0,n_pt_in
       X(ix) = 0.d0
       Y(ix) = 0.d0
     enddo
     call I_x2_pol_mult_mono_elec(c-1,R1x,R1xp,R2x,X,nx,n_pt_in)
       do ix=0,nx
         X(ix) =X(ix)* c
       enddo
       call multiply_poly_2(X,nx,R2x,2,d,nd,n_pt_in)
     ny=0
     call I_x2_pol_mult_mono_elec(c,R1x,R1xp,R2x,Y,ny,n_pt_in)
     call multiply_poly_2(Y,ny,R1x,2,d,nd,n_pt_in)
  else
     do ix=0,n_pt_in
       X(ix) = 0.d0
       Y(ix) = 0.d0
     enddo
     nx = 0
     call I_x1_pol_mult_mono_elec(a-2,c,R1x,R1xp,R2x,X,nx,n_pt_in)
!    print*,'nx a-2,c= ',nx
       do ix=0,nx
         X(ix) = X(ix)*(a-1)
       enddo
       call multiply_poly_2(X,nx,R2x,2,d,nd,n_pt_in)
!    print*,'nd out = ',nd

     nx = nd
     do ix=0,n_pt_in
       X(ix) = 0.d0
     enddo
     call I_x1_pol_mult_mono_elec(a-1,c-1,R1x,R1xp,R2x,X,nx,n_pt_in)
!      print*,'nx a-1,c-1 = ',nx
       do ix=0,nx
         X(ix) = X(ix)*c
       enddo
       call multiply_poly_2(X,nx,R2x,2,d,nd,n_pt_in)
     ny=0
      call I_x1_pol_mult_mono_elec(a-1,c,R1x,R1xp,R2x,Y,ny,n_pt_in)
      call multiply_poly_2(Y,ny,R1x,2,d,nd,n_pt_in)
  endif
end

recursive subroutine I_x2_pol_mult_mono_elec(c,R1x,R1xp,R2x,d,nd,dim) 
 implicit none
 integer , intent(in) :: dim
 double precision :: d(0:dim)
 integer,intent(inout) :: nd
 integer, intent(in):: c
 double precision, intent(in) :: R1x(0:2),R1xp(0:2),R2x(0:2)
 integer :: i,ix,nx,ny
 double precision :: X(0:dim),Y(0:dim)
!print*,'X2,c = ',c
!print*,'nd_in = ',nd

 if(c==0) then
   nd = 0
   d(0) = 1.d0
!  print*,'nd  IX2 = ',nd
   return
 elseif ((nd<0).or.(c<0))then
   nd = -1
    return
 else
     do ix=0,dim
       X(ix) = 0.d0
       Y(ix) = 0.d0
     enddo
     nx = 0
     call I_x1_pol_mult_mono_elec(0,c-2,R1x,R1xp,R2x,X,nx,dim)
!      print*,'nx 0,c-2 = ',nx
       do ix=0,nx
         X(ix) = X(ix)*(c-1)
       enddo
       call multiply_poly_2(X,nx,R2x,2,d,nd,dim)
!      print*,'nd = ',nd
       ny = 0
       do ix=0,dim
        Y(ix) = 0.d0
       enddo
       
       call I_x1_pol_mult_mono_elec(0,c-1,R1x,R1xp,R2x,Y,ny,dim)
!      print*,'ny = ',ny
!      do ix=0,ny
!        print*,'Y(ix) = ',Y(ix)
!      enddo
       if(ny.ge.0)then
        call multiply_poly_2(Y,ny,R1xp,2,d,nd,dim)
       endif
 endif
end

double precision function int_coul_bare_large_Q(nx,ny,nz,q,x_c)
implicit none
integer nx,ny,nz,ntot,big_m,i,n,m,l,n_triplet(10000),m_triplet(10000),l_triplet(10000),number_triplet
double precision q(3),x_c,accu,alpha,u1,u2,u3,term,sin_k_cos_l_twopi,sin_k_cos_l_pi,h_q
double precision accu_save(3)
logical iyes

ntot=nx+ny+nz

accu=0.d0
do big_m=0,500
 call find_triplet(big_m,number_triplet,n_triplet,m_triplet,l_triplet)
 if(number_triplet.gt.10000)stop 'increase 10000 in int_coul_bare_large_Q'
 do i=1,number_triplet
  n=n_triplet(i)
  m=m_triplet(i)
  l=l_triplet(i)
  alpha=h_q(n,q(1))*h_q(m,q(2))*h_q(l,q(3))
  u1=sin_k_cos_l_twopi(ny+m,nx+n)
  u2=sin_k_cos_l_pi(nx+n+ny+m+1,nz+l)
  u3=x_c**(ntot+n+m+l+2)/dfloat(ntot+n+m+l+2)
  term=u1*u2*u3*alpha
  accu=accu+term
 enddo
if(big_m.le.2)then
accu_save(big_m+1)=accu
else
accu_save(1)=accu_save(2)
accu_save(2)=accu_save(3)
accu_save(3)=accu
call test_convergence(accu_save,iyes)
if(iyes)then
!write(*,*)'big_m=',big_m
goto 100
endif
endif

enddo

stop 'danger!! big_m larger than 500'
100 continue
int_coul_bare_large_Q=accu
end

!
! int_coul_bare_small =  int dr 1/r  x**nx * y**ny * z**nz exp[-gam(r^2 - 2 p.r)]
!
!  
double precision function int_coul_bare_small(nx,ny,nz,gam,pt,r_c)
implicit none
double precision gam,x_c,r_c,pt(3)
double precision sin_k_cos_l_pi,sin_k_cos_l_twopi
double precision r_gauss,accu_save(3),r_gauss_comp
integer :: nx,ny,nz,ntot,n,p,q
double precision factor,accu,ux,uy,uz,coef,accunew,accuold,eps,t1,t2,t3,t4
logical run,iyes

ux=2.d0*gam*pt(1)
uy=2.d0*gam*pt(2)
uz=2.d0*gam*pt(3)

ntot=nx+ny+nz
x_c=dsqrt(gam)*r_c
factor=1.d0/dsqrt(gam)**(ntot+2)
accu=0.d0

do n=0,300
do p=0,n
do q=0,n-p
t1=coef(n,p,q,gam,ux,uy,uz)
t2=sin_k_cos_l_twopi(ny+q,n-p-q+nx)
t3=sin_k_cos_l_pi(nx+ny+1+n-p,nz+p)
!t4=r_gauss(x_c,ntot+n+1)
t4=r_gauss_comp(x_c,ntot+n+1)
accu=accu+t1*t2*t3*t4
enddo
enddo

if(n.le.2)then
accu_save(n+1)=accu
else
accu_save(1)=accu_save(2)
accu_save(2)=accu_save(3)
accu_save(3)=accu
call test_convergence(accu_save,iyes)
if(iyes)then
!write(*,*)'convergence reached in small for n=',n
goto 100
endif
endif

enddo
write(*,*)'danger!! n larger than 300'
100 continue
int_coul_bare_small=factor*accu

end

double precision function h_q(n,q)
implicit none
integer n
double precision q,factmichel,hermite
h_q=(-1.d0)**n/factmichel(n)*hermite(n,-q)
end

!
! sin_k_cos_l= int_0^{pi/2} dphi sin^k(phi) cos^l(phi)
!
double precision function sin_k_cos_l(k,l)
implicit none
integer :: k,l,n,m,k_parity,l_parity
double precision pi
double precision dblefact,factmichel
pi=dacos(-1.d0)
k_parity=mod(k,2)
l_parity=mod(l,2)
if(k_parity.eq.0.and.l_parity.eq.0)then
 n=k/2
 m=l/2
 sin_k_cos_l=pi/2.d0**(n+m+1)*dblefact(k-1)*dblefact(l-1)/factmichel(n+m)
 return
endif
if(k_parity.eq.0.and.l_parity.eq.1)then
 n=k/2
 m=(l-1)/2
 sin_k_cos_l=dblefact(2*m)*dblefact(k-1)/dblefact(k+l)
 return
endif
if(k_parity.eq.1.and.l_parity.eq.0)then
 n=(k-1)/2
 m=l/2
 sin_k_cos_l=dblefact(2*n)*dblefact(l-1)/dblefact(k+l)
 return
endif
if(k_parity.eq.1.and.l_parity.eq.1)then
 n=(k-1)/2
 m=(l-1)/2
 sin_k_cos_l=0.5d0*factmichel(n)*factmichel(m)/factmichel(n+m+1)
 return
endif
end

!
! sin_k_cos_l_pi= int_0^{pi} dphi sin^k(phi) cos^l(phi)
!
double precision function sin_k_cos_l_pi(k,l)
implicit none
integer :: k,l
double precision sin_k_cos_l
double precision ikl,ilk
ikl=sin_k_cos_l(k,l)
ilk=sin_k_cos_l(l,k)
sin_k_cos_l_pi=ikl+(-1.d0)**l*ilk
end

!
! sin_k_cos_l_twopi= int_0^{twopi} dphi sin^k(phi) cos^l(phi)
!
double precision function sin_k_cos_l_twopi(k,l)
implicit none
integer :: k,l
double precision sin_k_cos_l_pi
double precision jkl
jkl=sin_k_cos_l_pi(k,l)
sin_k_cos_l_twopi=jkl*(1.d0+(-1.d0)**(k+l))
end

subroutine find_triplet(big_m,number_triplet,n_triplet,m_triplet,l_triplet)
implicit none
integer big_m,number_triplet,n,m,l,i
integer n_triplet(10000),m_triplet(10000),l_triplet(10000)
i=0
do n=0,big_m
do m=0,big_m
do l=0,big_m
if(n+m+l.eq.big_m)then
i=i+1
if(i.gt.10000)stop 'number of triplets too large'
n_triplet(i)=n
m_triplet(i)=m
l_triplet(i)=l
endif
enddo
enddo
enddo
number_triplet=i
end

double precision function factmichel(n)
  implicit none
  integer :: n
  double precision, save :: memo(1:100)
  integer, save :: memomax = 1
  integer :: i

  if (n<=memomax) then
    if (n<2) then
      factmichel = 1.d0
    else
      factmichel = memo(n)
    endif
    return
  endif

  memo(1) = 1.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)*float(i)
  enddo
  memomax = min(n,100)
  factmichel = memo(memomax)
  do i=101,n
    factmichel = factmichel*float(i)
  enddo
end function

double precision function binommichel(i,j)
 implicit none
 integer :: i,j
 double precision :: factmichel
 binommichel = factmichel(i)/(factmichel(j)*factmichel(i-j))
end

subroutine test_convergence(accu,iyes)
implicit none
dimension accu(3)
double precision accu,d1,d2,d3,eps
logical iyes
iyes=.false.
eps=1.d-12
d1=dabs(accu(2)-accu(1))
d2=dabs(accu(3)-accu(1))
d3=dabs(accu(3)-accu(2))
if(d1.lt.eps.and.d2.lt.eps.and.d3.lt.eps)iyes=.true.
end


!      double precision function hermite(n,x)
!      implicit none
!      integer :: n,k
!      double precision :: h0,x,h1,h2
!      h0=1.d0
!      if(n.eq.0)then
!       hermite=h
!       return
!      endif
!      h1=2.d0*x
!      if(n.eq.1)then
!       hermite=h1
!       return
!      endif
!      do k=1,n-1
!       h2=2.d0*x*h1-2.d0*dfloat(k)*h0
!       h0=h1
!       h1=h2
!      enddo
!      hermite=h2
!      end


!* The routine  rint_c(n,rho) computes:
!**
!**   int_0^1 dt  exp(-rho*t**2) t^n  n=0,1,....
!*!*
!**   rho positif ou nul

      double precision function rint_c(n,rho)
      implicit double precision(a-h,o-z)
      if(mod(n,2).eq.0)then 
        rint_c=rint(n/2,rho)
      else
        rint_c=0.5d0*rintp((n-1)/2,rho)
      endif
      end

!** The routine  rintp(n,rho) computes:
!**
!**   int_0^1 dt  exp(-rho*t) t^n  n=0,1,....
!**
!**   rho positif ou nul
      double precision function rintp(n,rho)
      implicit double precision(a-h,o-z)
      if(rho.eq.0.d0)then
       rintp=1.d0/(n+1)
       return
      endif
      rintp=(1.d0-dexp(-rho))/rho
      if(n.eq.0)return
      do k=1,n
       rintp=(k*rintp-dexp(-rho))/rho
      enddo 
      end


! r_gauss_comp= int_0^{x_c} dx x^n exp(-x**2)
!
double precision function r_gauss_comp(x_c,n)
implicit none
double precision x_c,rint_c
integer :: n
r_gauss_comp=x_c**(n+1)*rint_c(n,x_c**2)
end


double precision function coef(n,p,q,gam,ux,uy,uz)
implicit none
integer ::n,p,q,i
double precision gam,f,ux,uy,uz
coef=1.d0
f=dsqrt(gam)
do i=1,p
coef=coef*uz/dfloat(i)/f
enddo
do i=1,q
coef=coef*uy/dfloat(i)/f
enddo
do i=1,n-p-q
coef=coef*ux/dfloat(i)/f
enddo
end

subroutine multiply_poly_2(b,nb,c,nc,d,nd,dim_int)
 implicit none
 ! D(t) += B(t)*C(t)
 integer, intent(in) :: dim_int
 integer, intent(in) :: nb, nc
 integer, intent(inout) :: nd
 double precision, intent(in) :: b(0:dim_int), c(0:dim_int)
 double precision, intent(inout) :: d(0:dim_int)

 integer :: ndtmp
 integer :: ib, ic, id
 if(nc==-1.or.nb==-1)then
  return
 endif
 ndtmp = nb+nc
 do ib=0,nb
 !print*,ib
  if ( b(ib)==0.d0) then
    cycle
  endif
  do ic = 0,nc
   d(ib+ic) = d(ib+ic) + c(ic) * b(ib)
  enddo
 enddo
 ndtmp = max(ndtmp,nd)
!do while (d(ndtmp) == 0.d0)
!   ndtmp -= 1
!   if(ndtmp.lt.0)then
!    exit
!   endif
!enddo
 do ndtmp = max(ndtmp,nd),0,-1
  if (d(ndtmp) /= 0.d0) then
    exit
  endif
 enddo
 nd = ndtmp

end

integer function ifind_kcp(i,j,k,l)
include 'j.inc'

!!  i(1) j(2) k(1) l(2)

if( (i.ge.k).and.(j.ge.l) )then

!! ijkl

 if(i.ge.j)then
  ni=i
  nj=j
  nk=k
  nl=l
 else
!! permutation of 1 and 2
  ni=j
  nj=i
  nk=l
  nl=k
 endif

endif

if( (i.ge.k).and.(j.lt.l) )then
!! permutation  : ilkj

 if(i.ge.l)then
  ni=i
  nj=l
  nk=k
  nl=j
 else
  ni=l
  nj=i
  nk=j
  nl=k
 endif
endif

if( (i.lt.k).and.(j.ge.l) )then
!! kjil

 if(k.ge.j)then
  ni=k
  nj=j
  nk=i
  nl=l
 else
  ni=j
  nj=k
  nk=l
  nl=i
 endif

endif

if( (i.lt.k).and.(j.lt.l) )then
!! klij

if(k.ge.l)then
  ni=k
  nj=l
  nk=i
  nl=j
 else
  ni=l
  nj=k
  nk=j
  nl=i
 endif
endif

isave=i
jsave=j
ksave=k
lsave=l

i=ni
j=nj
k=nk
l=nl

if(i.eq.j)then

 if(i.ge.k.and.i.ge.l.and.k.ge.l)then

!  i(1)i(2)k(1)l(2)
  ni1=i
  nj1=j
  nk1=k
  nl1=l
 endif

 if(i.ge.k.and.i.ge.l.and.k.lt.l)then

!  i(1)i(2)k(1)l(2)
!  i(2)i(1)k(2)l(1)
!  i(1)i(2)l(1)k(2)
  ni1=i
  nj1=i
  nk1=l
  nl1=k
 endif

 if(i.lt.k.and.i.ge.l)then

!  i(1)i(2)k(1)l(2)
!  k(1)i(2)i(1)l(2)
  ni1=k
  nj1=i
  nk1=i
  nl1=l
 endif

 if(i.lt.k.and.i.lt.l)then

!  i(1)i(2)k(1)l(2)
!  k(1)i(2)i(1)l(2)
!  k(1)l(2)i(1)i(2)

  if(k.ge.l)then
!  k(1)l(2)i(1)i(2)
   ni1=k
   nj1=l
   nk1=i
   nl1=i
  else
!  k(1)l(2)i(1)i(2)
!  k(2)l(1)i(2)i(1)
!  l(1)k(2)i(2)i(1)
   ni1=l
   nj1=k
   nk1=i
   nl1=i
  endif
 endif

 if(i.ge.k.and.i.lt.l)then
!  i(1)i(2)k(1)l(2)
!  i(1)l(2)k(1)i(2)
!  l(1)i(2)i(1)k(2)
  ni1=l
  nj1=i
  nk1=i
  nl1=k
 endif

 ni=ni1
 nj=nj1
 nk=nk1
 nl=nl1

endif

i=isave
j=jsave
k=ksave
l=lsave

ifind_kcp=kcp_ijkl(ni,nj,nk,nl)

if(ifind_kcp.eq.0)then
 write(*,*)' i j k l',i,j,k,l
 write(*,*)'ni nj nk nl',ni,nj,nk,nl
 write(*,*)'ifind_kcp=0',ni,nj,nk,nl,kcp_ijkl(ni,nj,nk,nl)
 stop
endif
end

      subroutine qcksrt(n,arr)
      implicit double precision (a-h,o-z)
      parameter (m=7,nstack=50,fm=7875.,fa=211.,fc=1663. &
         ,fmi=1.2698413e-4)
      dimension arr(n),istack(nstack)
      jstack=0
      l=1
      ir=n
      fx=0.
10    if(ir-l.lt.m)then
        do 13 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)go to 12
            arr(i+1)=arr(i)
11        continue
          i=0
12        arr(i+1)=a
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        i=l
        j=ir
        fx=mod(fx*fa+fc,fm)
        iq=l+(ir-l+1)*(fx*fmi)
        a=arr(iq)
        arr(iq)=arr(l)
20      continue
21        if(j.gt.0)then
            if(a.lt.arr(j))then
              j=j-1
              go to 21
            endif
          endif
          if(j.le.i)then
            arr(i)=a
            go to 30
          endif
          arr(i)=arr(j)
          i=i+1
22        if(i.le.n)then
            if(a.gt.arr(i))then
              i=i+1
              go to 22
            endif
          endif
          if(j.le.i)then
            arr(j)=a
            i=j
            go to 30
          endif
          arr(j)=arr(i)
          j=j-1
        go to 20
30      jstack=jstack+2
        if(jstack.gt.nstack)stop 'nstack must be made larger.'
        if(ir-i.ge.i-l)then
          istack(jstack)=ir
          istack(jstack-1)=i+1
          ir=i-1
        else
          istack(jstack)=i-1
          istack(jstack-1)=l
          l=i+1
        endif
      endif
      go to 10
      end

double precision function bielec_integral(alpha,beta,delta,gama,A_center,B_center, &
                          C_center,D_center,power_A,power_B,power_C,power_D)

  implicit none
  !DOC
  !  integral of the AO basis <ik|jl> =  (ij|kl)
  !                           <AC|BD> =  (AB|CD) 
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  !DOC

  integer, parameter             :: max_dim = 511
  integer,intent(in)             :: power_A(3),power_B(3),power_C(3),power_D(3)
  double precision, intent(in)   :: alpha,beta,delta,gama
  double precision, intent(in)   :: A_center(3), B_center(3), C_center(3), D_center(3)
  integer                        :: i,j,k,l
  integer                        :: p,q,r,s
  integer                        :: num_i,num_j,num_k,num_l,dim1
  double precision               :: integral
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: ao_bielec_integral_schwartz_accel
  double precision               :: p_inv,q_inv
  double precision               :: general_primitive_integral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision, parameter    :: pi =  dacos(-1.d0)
  double precision, parameter    :: sqpi =  dsqrt(dacos(-1.d0))
  double precision, parameter    :: pi_5_2 =  34.9868366552d0
  double precision, parameter    :: dfour_pi =  4.d0*dacos(-1.d0)
  double precision, parameter    :: dtwo_pi =  2.d0*dacos(-1.d0)
  double precision, parameter    :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: thresh = 1.d-15
  double precision, parameter    :: cx_lda = -0.73855876638202234d0
  double precision, parameter    :: c_2_4_3 = 2.5198420997897464d0
  double precision, parameter    :: cst_lda = -0.93052573634909996d0
  double precision, parameter    :: c_4_3 = 1.3333333333333333d0
  double precision, parameter    :: c_1_3 = 0.3333333333333333d0
  
  dim1 = 511 
  
  call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
          alpha,beta,power_A,power_B,A_center,B_center,dim1)
  p_inv = 1.d0/pp

  call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
           delta,gama,power_C,power_D,C_center,D_center,dim1)
  q_inv = 1.d0/qq

  bielec_integral = general_primitive_integral(dim1,              &
          P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
          Q_new,Q_center,fact_q,qq,q_inv,iorder_q)

end

!! p=alpha+beta
!! q=delta+gam
!! p_inv=1/p
!! q_inv=1/q
!!
!! fact_p=E_AB
!! fact_a=E_CD

double precision function general_primitive_integral(dim,            &
      P_new,P_center,fact_p,p,p_inv,iorder_p,                        &
      Q_new,Q_center,fact_q,q,q_inv,iorder_q)
  implicit none
  !EGIN_DOC
  ! Computes the integral <pq|rs> where p,q,r,s are Gaussian primitives
  !ND_DOC
  integer, parameter             :: max_dim = 511
  integer,intent(in)             :: dim
  double precision, intent(in)   :: P_new(0:max_dim,3),P_center(3),fact_p,p,p_inv
  double precision, intent(in)   :: Q_new(0:max_dim,3),Q_center(3),fact_q,q,q_inv
  integer, intent(in)            :: iorder_p(3)
  integer, intent(in)            :: iorder_q(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision, parameter    :: pi =  dacos(-1.d0)
  double precision, parameter    :: sqpi =  dsqrt(dacos(-1.d0))
  double precision, parameter    :: pi_5_2 =  34.9868366552d0
  double precision, parameter    :: dfour_pi =  4.d0*dacos(-1.d0)
  double precision, parameter    :: dtwo_pi =  2.d0*dacos(-1.d0)
  double precision, parameter    :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: thresh = 1.d-15
  double precision, parameter    :: cx_lda = -0.73855876638202234d0
  double precision, parameter    :: c_2_4_3 = 2.5198420997897464d0
  double precision, parameter    :: cst_lda = -0.93052573634909996d0
  double precision, parameter    :: c_4_3 = 1.3333333333333333d0
  double precision, parameter    :: c_1_3 = 0.3333333333333333d0

  
  double precision               :: r_cut,gama_r_cut,rho,dist
  double precision               :: dx(0:max_dim),Ix_pol(0:max_dim),dy(0:max_dim),Iy_pol(0:max_dim),dz(0:max_dim),Iz_pol(0:max_dim)
  integer                        :: n_Ix,n_Iy,n_Iz,nx,ny,nz
  double precision               :: bla
  integer                        :: ix,iy,iz,jx,jy,jz,i
  double precision               :: a,b,c,d,e,f,accu,pq,const
  double precision               :: pq_inv, p10_1, p10_2, p01_1, p01_2,pq_inv_2
  integer                        :: n_pt_tmp,n_pt_out, iorder
  double precision               :: d1(0:max_dim),d_poly(0:max_dim),rint,d1_screened(0:max_dim)
  double precision               :: rint_sum
  
  general_primitive_integral = 0.d0
  
  ! Gaussian Product
  ! ----------------
  
  pq = p_inv*0.5d0*q_inv
  pq_inv = 0.5d0/(p+q)
  p10_1 = q*pq  ! 1/(2p)
  p01_1 = p*pq  ! 1/(2q)
  pq_inv_2 = pq_inv+pq_inv
  p10_2 = pq_inv_2 * p10_1*q !0.5d0*q/(pq + p*p)
  p01_2 = pq_inv_2 * p01_1*p !0.5d0*p/(q*q + pq)
  
!! Ix
!****  
  accu = 0.d0
  iorder = iorder_p(1)+iorder_q(1)+iorder_p(1)+iorder_q(1)
  do ix=0,iorder
    Ix_pol(ix) = 0.d0
  enddo
  n_Ix = 0
  do ix = 0, iorder_p(1)
    if (abs(P_new(ix,1)) < thresh) cycle
    a = P_new(ix,1)
    do jx = 0, iorder_q(1)
      d = a*Q_new(jx,1)
      if (abs(d) < thresh) cycle
      call give_polynom_mult_center_x &
      (P_center(1),Q_center(1),ix,jx,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dx,nx)
      call add_poly_multiply(dx,nx,d,Ix_pol,n_Ix)
    enddo
  enddo
  if (n_Ix == -1) then
    return
  endif

!end Ix
!******

!! Iy
!****
  iorder = iorder_p(2)+iorder_q(2)+iorder_p(2)+iorder_q(2)
  do ix=0, iorder
    Iy_pol(ix) = 0.d0
  enddo
  n_Iy = 0
  do iy = 0, iorder_p(2)
    if (abs(P_new(iy,2)) > thresh) then
      b = P_new(iy,2)
      do jy = 0, iorder_q(2)
        e = b*Q_new(jy,2)
        if (abs(e) < thresh) cycle
        call give_polynom_mult_center_x &
        (P_center(2),Q_center(2),iy,jy,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dy,ny)
        call add_poly_multiply(dy,ny,e,Iy_pol,n_Iy)
      enddo
    endif
  enddo
  if (n_Iy == -1) then
    return
  endif

!end Iy
!******

!! Iz
!****
  iorder = iorder_p(3)+iorder_q(3)+iorder_p(3)+iorder_q(3)
  do ix=0,iorder
    Iz_pol(ix) = 0.d0
  enddo
  n_Iz = 0
  do iz = 0, iorder_p(3)
    if (abs(P_new(iz,3)) > thresh) then
      c = P_new(iz,3)
      do jz = 0, iorder_q(3)
        f = c*Q_new(jz,3)
        if (abs(f) < thresh) cycle

        call give_polynom_mult_center_x &
        (P_center(3),Q_center(3),iz,jz,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dz,nz)

        call add_poly_multiply(dz,nz,f,Iz_pol,n_Iz)
      enddo
    endif
  enddo
  if (n_Iz == -1) then
    return
  endif

!end Iz
!******

  rho = p*q *pq_inv_2
  dist =  (P_center(1) - Q_center(1))*(P_center(1) - Q_center(1)) +  &
          (P_center(2) - Q_center(2))*(P_center(2) - Q_center(2)) +  &
          (P_center(3) - Q_center(3))*(P_center(3) - Q_center(3))
  const = dist*rho
  
  n_pt_tmp = n_Ix + n_Iy
  do i=0,n_pt_tmp
   d_poly(i)=0.d0
  enddo
  call multiply_poly(Ix_pol,n_Ix,Iy_pol,n_Iy,d_poly,n_pt_tmp)
  if (n_pt_tmp == -1) then
   return
  endif

  n_pt_out = n_pt_tmp + n_Iz
  do i=0,n_pt_out
   d1(i)=0.d0
  enddo
  call multiply_poly(d_poly ,n_pt_tmp ,Iz_pol,n_Iz,d1,n_pt_out)

!!  rint_sum = accu= int_0^1 dt exp(-const t^2) Ix Iy Iz
!!
!!  const=rho*(P-Q)**2
!!
!! n_pt_out = n_Ix+n_Iy+n_Iz = nxA+nxB+ nyA+nyB nzA+nzB
!!
!! d1(i) i=0 to n_pt_out   d1(0:max_dim) max_dim=511
!!
!! Ix(t**2) Iy(t**2) Iz(t**2)  = sum_i=0^n_pt_out d(i)  (t^2)^i 

  accu = accu + rint_sum(n_pt_out,const,d1)

!! pi_5_2 = 2 pi^5/2

!! p=alpha+beta
!! q=gama+delta

!!  p_inv=1/p
!!  q_inv=1/q

!! fact_p = E_AB  fact_q= E_CD
!!
!! E_AB= exp(-alpha*beta/(alpha+beta)*(A-B)**2)
!! E_CD= exp(-gama*delta/(gama+delta)*(C-D)**2)

! accu= int_0^1 dt exp(-rho*(P-Q)**2 t^2) Ix Iy Iz

! ERI= 2 pi**(5/2)/( pq*sqrt(p+q) ) E_AB E_CD  int_0^1....

 general_primitive_integral = fact_p * fact_q * accu *pi_5_2*p_inv*q_inv/dsqrt(p+q)

end

subroutine give_polynom_mult_center_x(P_center,Q_center,a_x,d_x,p,q,n_pt_in,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,d,n_pt_out)
  implicit none
!  BEGIN_DOC
  ! subroutine that returns the explicit polynom in term of the "t"
  ! variable of the following polynomw :
  !         I_x1(a_x, d_x,p,q) * I_x1(a_y, d_y,p,q) * I_x1(a_z, d_z,p,q)
!  END_DOC

  integer, intent(in)            :: n_pt_in
  integer,intent(out)            :: n_pt_out
  integer, intent(in)            :: a_x,d_x
  double precision, intent(in)   :: P_center, Q_center
  double precision, intent(in)   :: p,q,pq_inv,p10_1,p01_1,p10_2,p01_2,pq_inv_2

!!  include 'Utils/constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: SIMD_vector = 32
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15


  double precision               :: B10(0:2), B01(0:2), B00(0:2),C00(0:2),D00(0:2)
  integer                        :: n_pt1,dim,i

  double precision,intent(out)   :: d(0:max_dim)
  double precision               :: accu
  accu = 0.d0
  ! pq_inv = 0.5d0/(p+q)
  ! pq_inv_2 = 1.d0/(p+q)
  ! p10_1 = 0.5d0/p
  ! p01_1 = 0.5d0/q
  ! p10_2 = 0.5d0 *  q /(p * q + p * p)
  ! p01_2 = 0.5d0 *  p /(q * q + q * p)
  B10(0)  = p10_1
  B10(1)  = 0.d0
  B10(2)  = - p10_2
  ! B10 = p01_1 - t**2 * p10_2
  B01(0)  = p01_1
  B01(1)  = 0.d0
  B01(2)  = - p01_2
  ! B01 = p01_1- t**2 * pq_inv
  B00(0)  = 0.d0
  B00(1)  = 0.d0
  B00(2)  = pq_inv
  ! B00 = t**2 * pq_inv
  do i = 0,n_pt_in
    d(i) = 0.d0
  enddo
  n_pt1 = n_pt_in
  ! C00 = -q/(p+q)*(Px-Qx) * t^2
  C00(0) = 0.d0
  C00(1) = 0.d0
  C00(2) =  -q*(P_center-Q_center) * pq_inv_2
  ! D00 = -p/(p+q)*(Px-Qx) * t^2
  D00(0) = 0.d0
  D00(1) = 0.d0
  D00(2) =  -p*(Q_center-P_center) * pq_inv_2
  !D00(2) =  -p*(Q_center(1)-P_center(1)) /(p+q)
  call I_x1_pol_mult(a_x,d_x,B10,B01,B00,C00,D00,d,n_pt1,n_pt_in)
  n_pt_out = n_pt1
  if(n_pt1<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif
  
end

subroutine I_x1_pol_mult(a,c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
  implicit none
!  BEGIN_DOC
  ! recursive function involved in the bielectronic integral
!  END_DOC

  integer , intent(in)           :: n_pt_in

!!  include 'Utils/constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: SIMD_vector = 32
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15


  double precision,intent(inout) :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: a,c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  if( (c>=0).and.(nd>=0) )then
    
    if (a==1) then
      call I_x1_pol_mult_a1(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
    else if (a==2) then
      call I_x1_pol_mult_a2(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
    else if (a>2) then
      call I_x1_pol_mult_recurs(a,c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
    else  ! a == 0
      
      if( c==0 )then
        nd = 0
        d(0) = 1.d0
        return
      endif
      
      call I_x2_pol_mult(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
    endif
  else
    nd = -1
  endif
end

recursive subroutine I_x1_pol_mult_recurs(a,c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
  implicit none
!  BEGIN_DOC
  ! recursive function involved in the bielectronic integral
!  END_DOC

  integer , intent(in)           :: n_pt_in

!!  include 'Utils/constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: SIMD_vector = 32
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15


  double precision,intent(inout) :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: a,c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  double precision               :: X(0:max_dim)
  double precision               :: Y(0:max_dim)
  integer                        :: nx, ix,iy,ny
  
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo
  nx = 0
  if (a==3) then
    call I_x1_pol_mult_a1(c,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
  else if (a==4) then
    call I_x1_pol_mult_a2(c,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
  else
    call I_x1_pol_mult_recurs(a-2,c,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
  endif
  
  do ix=0,nx
    X(ix) = X(ix) * dble(a-1)
  enddo
  
  call multiply_poly(X,nx,B_10,2,d,nd)
  
  nx = nd
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo
  
  if (c>0) then
    if (a==3) then
      call I_x1_pol_mult_a2(c-1,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
    else
      call I_x1_pol_mult_recurs(a-1,c-1,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
    endif
    if (c>1) then
      do ix=0,nx
        X(ix) = X(ix)* c
      enddo
    endif
    call multiply_poly(X,nx,B_00,2,d,nd)
  endif
  
  ny=0
  
  do ix=0,n_pt_in
    Y(ix) = 0.d0
  enddo
  if (a==3) then
    call I_x1_pol_mult_a2(c,B_10,B_01,B_00,C_00,D_00,Y,ny,n_pt_in)
  else
    call I_x1_pol_mult_recurs(a-1,c,B_10,B_01,B_00,C_00,D_00,Y,ny,n_pt_in)
  endif
  
  call multiply_poly(Y,ny,C_00,2,d,nd)
  
end

recursive subroutine I_x1_pol_mult_a1(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
  implicit none
!  BEGIN_DOC
  ! recursive function involved in the bielectronic integral
!  END_DOC

  integer , intent(in)           :: n_pt_in

!!  include 'Utils/constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15



  double precision,intent(inout) :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  double precision               :: X(0:max_dim)
  double precision               :: Y(0:max_dim)
  integer                        :: nx, ix,iy,ny
  
  if( (c<0).or.(nd<0) )then
    nd = -1
    return
  endif
  
  nx = nd
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo
  call I_x2_pol_mult(c-1,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
  
  if (c>1) then
    do ix=0,nx
      X(ix) = X(ix) * dble(c)
    enddo
  endif
  
  call multiply_poly(X,nx,B_00,2,d,nd)
  
  ny=0
  
  do ix=0,n_pt_in
    Y(ix) = 0.d0
  enddo
  call I_x2_pol_mult(c,B_10,B_01,B_00,C_00,D_00,Y,ny,n_pt_in)
  
  call multiply_poly(Y,ny,C_00,2,d,nd)
  
end

recursive subroutine I_x1_pol_mult_a2(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
  implicit none
!  BEGIN_DOC
  ! recursive function involved in the bielectronic integral
!  END_DOC

  integer , intent(in)           :: n_pt_in

!!  include 'Utils/constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15


  double precision,intent(inout) :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  double precision               :: X(0:max_dim)
  double precision               :: Y(0:max_dim)
  integer                        :: nx, ix,iy,ny
  
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo
  nx = 0
  call I_x2_pol_mult(c,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
  
  call multiply_poly(X,nx,B_10,2,d,nd)
  
  nx = nd
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo
  
  call I_x1_pol_mult_a1(c-1,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
  
  if (c>1) then
    do ix=0,nx
      X(ix) = X(ix) * dble(c)
    enddo
  endif
  
  call multiply_poly(X,nx,B_00,2,d,nd)
  
  ny=0
  do ix=0,n_pt_in
    Y(ix) = 0.d0
  enddo
  call I_x1_pol_mult_a1(c,B_10,B_01,B_00,C_00,D_00,Y,ny,n_pt_in)
  
  call multiply_poly(Y,ny,C_00,2,d,nd)
  
end

recursive subroutine I_x2_pol_mult(c,B_10,B_01,B_00,C_00,D_00,d,nd,dim)
  implicit none
!  BEGIN_DOC
  ! recursive function involved in the bielectronic integral
!  END_DOC

  integer , intent(in)           :: dim

!!  include 'Utils/constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15


  double precision               :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  integer                        :: nx, ix,ny
  double precision               :: X(0:max_dim),Y(0:max_dim)
  integer                        :: i
  
  select case (c)
    case (0)
      nd = 0
      d(0) = 1.d0
      return
      
    case (:-1)
      nd = -1
      return
      
    case (1)
      nd = 2
      d(0) = D_00(0)
      d(1) = D_00(1)
      d(2) = D_00(2)
      return
      
    case (2)
      nd = 2
      d(0) = B_01(0)
      d(1) = B_01(1)
      d(2) = B_01(2)
      
      ny = 2
      Y(0) = D_00(0)
      Y(1) = D_00(1)
      Y(2) = D_00(2)
      
      call multiply_poly(Y,ny,D_00,2,d,nd)
      return
      
      case default
      
      do ix=0,c+c
        X(ix) = 0.d0
      enddo
      nx = 0
      call I_x2_pol_mult(c-2,B_10,B_01,B_00,C_00,D_00,X,nx,dim)
      
      do ix=0,nx
        X(ix) = X(ix) * dble(c-1)
      enddo
      
      call multiply_poly(X,nx,B_01,2,d,nd)
      
      ny = 0
      do ix=0,c+c
        Y(ix) = 0.d0
      enddo
      call I_x2_pol_mult(c-1,B_10,B_01,B_00,C_00,D_00,Y,ny,dim)
      
      call multiply_poly(Y,ny,D_00,2,d,nd)
      
  end select
end

subroutine give_explicit_poly_and_gaussian_x(P_new,P_center,p,fact_k,iorder,alpha,beta,a,b,A_center,B_center,dim)
!  BEGIN_DOC
! Transform the product of
!  (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) 
!  exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
!
! into
!        fact_k  (x-x_P)^iorder(1)  (y-y_P)^iorder(2)  (z-z_P)^iorder(3) exp(-p(r-P)^2)
!  END_DOC
  implicit none

!!  include 'constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15

  integer, intent(in)            :: dim
  integer, intent(in)            :: a,b               ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)   :: alpha, beta       ! exponents
  double precision, intent(in)   :: A_center          ! A center
  double precision, intent(in)   :: B_center          ! B center
  double precision, intent(out)  :: P_center          ! new center
  double precision, intent(out)  :: p                 ! new exponent
  double precision, intent(out)  :: fact_k            ! constant factor
  double precision, intent(out)  :: P_new(0:max_dim)  ! polynomial
  integer, intent(out)           :: iorder            ! order of the polynomials

  double precision               :: P_a(0:max_dim), P_b(0:max_dim)
  integer                        :: n_new,i,j
  double precision               :: p_inv,ab,d_AB

  ! Do the gaussian product to get the new center and the new exponent
  P_new = 0.d0
  p = alpha+beta
  p_inv = 1.d0/p
  ab = alpha * beta
  d_AB = (A_center - B_center) * (A_center - B_center)
  P_center = (alpha * A_center + beta * B_center) * p_inv
  fact_k = exp(-ab*p_inv * d_AB)

  ! Recenter the polynomials P_a and P_b on x
  call recentered_poly2(P_a(0),A_center,P_center,a,P_b(0),B_center,P_center,b)
  n_new = 0

  call multiply_poly(P_a(0),a,P_b(0),b,P_new(0),n_new)
  iorder = a + b
end


subroutine give_explicit_poly_and_gaussian(P_new,P_center,p,fact_k,iorder,alpha,beta,a,b,A_center,B_center,dim)
!  BEGIN_DOC
! Transforms the product of
! (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) 
! exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
! into
! fact_k *[sum (l_x = 0,i_order(1)) P_new(l_x,1)*(x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
!        *[sum (l_y = 0,i_order(2)) P_new(l_y,2)*(y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
!        *[sum (l_z = 0,i_order(3)) P_new(l_z,3)*(z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
!  END_DOC
  implicit none

!!  include 'constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15

  integer, intent(in)            :: dim
  integer, intent(in)            :: a(3),b(3)         ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)   :: alpha, beta       ! exponents
  double precision, intent(in)   :: A_center(3)       ! A center
  double precision, intent(in)   :: B_center (3)      ! B center
  double precision, intent(out)  :: P_center(3)       ! new center
  double precision, intent(out)  :: p                 ! new exponent
  double precision, intent(out)  :: fact_k            ! constant factor
  double precision, intent(out)  :: P_new(0:max_dim,3)! polynomial
  integer, intent(out)           :: iorder(3)         ! i_order(i) = order of the polynomials

  double precision               :: P_a(0:max_dim,3), P_b(0:max_dim,3)
  integer                        :: n_new,i,j
  
  iorder(1) = 0
  iorder(2) = 0
  iorder(3) = 0
  P_new(0,1) = 0.d0
  P_new(0,2) = 0.d0
  P_new(0,3) = 0.d0
  
  call gaussian_product(alpha,A_center,beta,B_center,fact_k,p,P_center)
  if (fact_k < thresh) then
    fact_k = 0.d0
    return
  endif
  
  call recentered_poly2(P_a(0,1),A_center(1),P_center(1),a(1),P_b(0,1),B_center(1),P_center(1),b(1))
  iorder(1) = a(1) + b(1)
  do i=0,iorder(1)
    P_new(i,1) = 0.d0
  enddo
  n_new=0
  call multiply_poly(P_a(0,1),a(1),P_b(0,1),b(1),P_new(0,1),n_new)
  
  call recentered_poly2(P_a(0,2),A_center(2),P_center(2),a(2),P_b(0,2),B_center(2),P_center(2),b(2))
  iorder(2) = a(2) + b(2)
  do i=0,iorder(2)
    P_new(i,2) = 0.d0
  enddo
  n_new=0
  call multiply_poly(P_a(0,2),a(2),P_b(0,2),b(2),P_new(0,2),n_new)
  
  call recentered_poly2(P_a(0,3),A_center(3),P_center(3),a(3),P_b(0,3),B_center(3),P_center(3),b(3))
  iorder(3) = a(3) + b(3)
  do i=0,iorder(3)
    P_new(i,3) = 0.d0
  enddo
  n_new=0
  call multiply_poly(P_a(0,3),a(3),P_b(0,3),b(3),P_new(0,3),n_new)
  
end


subroutine give_explicit_poly_and_gaussian_double(P_new,P_center,p,fact_k,iorder, &
                            alpha,beta,gama,a,b,A_center,B_center,Nucl_center,dim)
!  BEGIN_DOC
  ! Transforms the product of
  !          (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3) 
  !          exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta) exp(-(r-Nucl_center)^2 gama 
  !
  ! into
  !        fact_k * [ sum (l_x = 0,i_order(1)) P_new(l_x,1) * (x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
  !               * [ sum (l_y = 0,i_order(2)) P_new(l_y,2) * (y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
  !               * [ sum (l_z = 0,i_order(3)) P_new(l_z,3) * (z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
!  END_DOC
  implicit none

!!  include 'constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15

  integer, intent(in)            :: dim
  integer, intent(in)            :: a(3),b(3)         ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)   :: alpha, beta, gama ! exponents
  double precision, intent(in)   :: A_center(3)       ! A center
  double precision, intent(in)   :: B_center (3)      ! B center
  double precision, intent(in)   :: Nucl_center(3)    ! B center
  double precision, intent(out)  :: P_center(3)       ! new center
  double precision, intent(out)  :: p                 ! new exponent
  double precision, intent(out)  :: fact_k            ! constant factor
  double precision, intent(out)  :: P_new(0:max_dim,3)! polynomial
  integer         , intent(out)  :: iorder(3)         ! i_order(i) = order of the polynomials

  double precision  :: P_center_tmp(3)       ! new center
  double precision  :: p_tmp                 ! new exponent
  double precision  :: fact_k_tmp,fact_k_bis ! constant factor
  double precision  :: P_new_tmp(0:max_dim,3)! polynomial
  integer :: i,j
  double precision :: binom_func

  ! First you transform the two primitives into a sum of primitive with the same center P_center_tmp and gaussian exponent p_tmp
  call give_explicit_poly_and_gaussian(P_new_tmp,P_center_tmp,p_tmp,fact_k_tmp,iorder,alpha,beta,a,b,A_center,B_center,dim)
  ! Then you create the new gaussian from the product of the new one per the Nuclei one
  call gaussian_product(p_tmp,P_center_tmp,gama,Nucl_center,fact_k_bis,p,P_center)
  fact_k = fact_k_bis * fact_k_tmp

  ! Then you build the coefficient of the new polynom 
  do i = 0, iorder(1)
   P_new(i,1) = 0.d0
   do j = i,iorder(1)
    P_new(i,1) = P_new(i,1) + P_new_tmp(j,1) * binom_func(j,j-i) * (P_center(1) - P_center_tmp(1))**(j-i)
   enddo
  enddo
  do i = 0, iorder(2)
   P_new(i,2) = 0.d0
   do j = i,iorder(2)
    P_new(i,2) = P_new(i,2) + P_new_tmp(j,2) * binom_func(j,j-i) * (P_center(2) - P_center_tmp(2))**(j-i)
   enddo
  enddo
  do i = 0, iorder(3)
   P_new(i,3) = 0.d0
   do j = i,iorder(3)
    P_new(i,3) = P_new(i,3) + P_new_tmp(j,3) * binom_func(j,j-i) * (P_center(3) - P_center_tmp(3))**(j-i)
   enddo
  enddo
  
end



subroutine gaussian_product(a,xa,b,xb,k,p,xp)
  implicit none
!  BEGIN_DOC
  ! Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
!  END_DOC
  
  double precision, intent(in)   :: a,b         ! Exponents
  double precision, intent(in)   :: xa(3),xb(3) ! Centers
  double precision, intent(out)  :: p           ! New exponent
  double precision, intent(out)  :: xp(3)       ! New center
  double precision, intent(out)  :: k           ! Constant
  
  double precision               :: p_inv
  double precision               :: xab(3), ab
  
  p = a+b
  p_inv = 1.d0/(a+b)
  ab = a*b
  xab(1) = xa(1)-xb(1)
  xab(2) = xa(2)-xb(2)
  xab(3) = xa(3)-xb(3)
  ab = ab*p_inv
  k = ab*(xab(1)*xab(1)+xab(2)*xab(2)+xab(3)*xab(3))
  if (k > 40.d0) then
    k=0.d0
    return
  endif
  k = dexp(-k)
  xp(1) = (a*xa(1)+b*xb(1))*p_inv
  xp(2) = (a*xa(2)+b*xb(2))*p_inv
  xp(3) = (a*xa(3)+b*xb(3))*p_inv
end subroutine




subroutine gaussian_product_x(a,xa,b,xb,k,p,xp)
  implicit none
!  BEGIN_DOC
  ! Gaussian product in 1D.
  ! e^{-a (x-x_A)^2} e^{-b (x-x_B)^2} = K_{ab}^x e^{-p (x-x_P)^2}
!  END_DOC
  
  double precision  , intent(in) :: a,b      ! Exponents
  double precision  , intent(in) :: xa,xb    ! Centers
  double precision  , intent(out) :: p       ! New exponent
  double precision  , intent(out) :: xp      ! New center
  double precision  , intent(out) :: k       ! Constant
  
  double precision               :: p_inv
  double precision               :: xab, ab
  
  p = a+b
  p_inv = 1.d0/(a+b)
  ab = a*b
  xab = xa-xb
  ab = ab*p_inv
  k = ab*xab*xab
  if (k > 40.d0) then
    k=0.d0
    return
  endif
  k = exp(-k)
  xp = (a*xa+b*xb)*p_inv
end subroutine

subroutine multiply_poly(b,nb,c,nc,d,nd)
  implicit none
!  BEGIN_DOC
  ! Multiply two polynomials
  ! D(t) += B(t)*C(t)
!  END_DOC
  
  integer, intent(in)            :: nb, nc
  integer, intent(out)           :: nd
  double precision, intent(in)   :: b(0:nb), c(0:nc)
  double precision, intent(inout) :: d(0:nb+nc)
  
  integer                        :: ndtmp
  integer                        :: ib, ic, id, k
  if(ior(nc,nb) >= 0) then ! True if nc>=0 and nb>=0
    continue
  else
    return
  endif
  ndtmp = nb+nc
  
  do ic = 0,nc
    d(ic) = d(ic) + c(ic) * b(0)
  enddo
  
  do ib=1,nb
    d(ib) = d(ib) + c(0) * b(ib)
    do ic = 1,nc
      d(ib+ic) = d(ib+ic) + c(ic) * b(ib)
    enddo
  enddo
  
  do nd = ndtmp,0,-1
    if (d(nd) == 0.d0) then
      cycle
    endif
    exit
  enddo
  
end

subroutine add_poly(b,nb,c,nc,d,nd)
  implicit none
!  BEGIN_DOC
  ! Add two polynomials
  ! D(t) += B(t)+C(t)
!  END_DOC

  integer                        :: ib, ic, id

  integer, intent(inout)         :: nb, nc
  integer, intent(out)           :: nd
  double precision, intent(in)   :: b(0:nb), c(0:nc)
  double precision, intent(out)  :: d(0:nb+nc)
  
  nd = nb+nc
  do ib=0,max(nb,nc)
    d(ib) = d(ib) + c(ib) + b(ib)
  enddo
  do while ( (d(nd) == 0.d0).and.(nd>=0) )
    nd = nd -1
    if (nd < 0) then
      exit
    endif
  enddo
  
end

subroutine add_poly_multiply(b,nb,cst,d,nd)
  implicit none
!  BEGIN_DOC
  ! Add a polynomial multiplied by a constant
  ! D(t) += cst * B(t)
!  END_DOC

    integer                        :: ib, ic, id

  integer, intent(in)            :: nb
  integer, intent(inout)         :: nd
  double precision, intent(in)   :: b(0:nb),cst
  double precision, intent(inout) :: d(0:max(nb,nd))
  
  nd = max(nd,nb)
  if (nd /= -1) then
    do ib=0,nb
      d(ib) = d(ib) + cst*b(ib)
    enddo
    do while ( d(nd) == 0.d0 )
      nd = nd - 1
      if (nd < 0) then
        exit
      endif
    enddo
  endif
  
end



subroutine recentered_poly2(P_new,x_A,x_P,a,P_new2,x_B,x_Q,b)
  implicit none
!  BEGIN_DOC
  ! Recenter two polynomials
!  END_DOC
  integer, intent(in)            :: a,b
  double precision, intent(in)   :: x_A,x_P,x_B,x_Q
  double precision, intent(out)  :: P_new(0:a),P_new2(0:b)
  double precision               :: pows_a(-2:a+b+4), pows_b(-2:a+b+4)
  double precision               :: binom_func
  double precision               :: binom_transp
  integer                        :: i,j,k,l, minab, maxab
  if ((a<0).or.(b<0) ) return
  maxab = max(a,b)
  minab = max(min(a,b),0)
  pows_a(0) = 1.d0
  pows_a(1) = (x_P - x_A)
  pows_b(0) = 1.d0
  pows_b(1) = (x_Q - x_B)
  do i =  2,maxab
    pows_a(i) = pows_a(i-1)*pows_a(1)
    pows_b(i) = pows_b(i-1)*pows_b(1)
  enddo
  P_new (0) =  pows_a(a)
  P_new2(0) =  pows_b(b)
  do i =  1,min(minab,20)
    P_new (i) =  binom_transp(a-i,a) * pows_a(a-i)
    P_new2(i) =  binom_transp(b-i,b) * pows_b(b-i)
  enddo
  do i =  minab+1,min(a,20)
    P_new (i) =  binom_transp(a-i,a) * pows_a(a-i)
  enddo
  do i =  minab+1,min(b,20)
    P_new2(i) =  binom_transp(b-i,b) * pows_b(b-i)
  enddo
  do i =  101,a
    P_new(i) =  binom_func(a,a-i) * pows_a(a-i)
  enddo
  do i =  101,b
    P_new2(i) =  binom_func(b,b-i) * pows_b(b-i)
  enddo
end




double precision function F_integral(n,p)
!  BEGIN_DOC
  ! function that calculates the following integral
  ! \int_{\-infty}^{+\infty} x^n \exp(-p x^2) dx
!  END_DOC
  implicit none
  integer                        :: n
  double precision               :: p
  integer                        :: i,j
  double precision               :: accu,sqrt_p,fact_ratio,tmp,fact

!!  include 'constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15


  if(n < 0)then
    F_integral = 0.d0
  endif
  if(iand(n,1).ne.0)then
    F_integral = 0.d0
    return
  endif
  sqrt_p = 1.d0/dsqrt(p)
  if(n==0)then
    F_integral = sqpi * sqrt_p
    return
  endif
  F_integral = sqpi * 0.5d0**n * sqrt_p**(n+1) * fact(n)/fact(ishft(n,-1))
end

double precision function rint(n,rho)
  implicit none
!  BEGIN_DOC
!.. math::
!
!  \int_0^1 dx \exp(-p x^2) x^n
!
!  END_DOC

!!  include 'constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15



  double precision               :: rho,u,rint1,v,val0,rint_large_n,u_inv
  integer                        :: n,k
  double precision               :: two_rho_inv

  if(rho.lt.-100.)then
   write(*,*)'rho=',rho
   stop
  endif
  
  if(n.eq.0)then
    if(rho == 0.d0)then
      rint=1.d0
    else
      u_inv=1.d0/dsqrt(rho)
      u=rho*u_inv
      rint=0.5d0*u_inv*sqpi*erf(u)
    endif
    return
  endif
  if(rho.lt.1.d0)then
    rint=rint1(n,rho)
  else
    if(n.le.20)then
      u_inv=1.d0/dsqrt(rho)
      if(rho.gt.80.d0)then
       v=0.d0
      else
       v=dexp(-rho)
      endif
      u=rho*u_inv
      two_rho_inv = 0.5d0*u_inv*u_inv
      val0=0.5d0*u_inv*sqpi*erf(u)
      rint=(val0-v)*two_rho_inv
      do k=2,n
        rint=(rint*dfloat(k+k-1)-v)*two_rho_inv
      enddo
    else
      rint=rint_large_n(n,rho)
    endif
  endif
end



double precision function rint_sum(n_pt_out,rho,d1)
  implicit none
!  BEGIN_DOC
  ! Needed for the calculation of two-electron integrals.
!  END_DOC

!!  include 'constants.include.F'
integer, parameter :: max_dim = 511
integer, parameter :: N_int_max = 32

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15


  integer, intent(in)            :: n_pt_out
  double precision, intent(in)   :: rho,d1(0:n_pt_out)
  double precision               :: u,rint1,v,val0,rint_large_n,u_inv
  integer                        :: n,k,i
  double precision               :: two_rho_inv, rint_tmp, di
  
  
  if(rho < 1.d0)then
    
    if(rho == 0.d0)then
      rint_sum=d1(0)
    else
      u_inv=1.d0/dsqrt(rho)
      u=rho*u_inv
      rint_sum=0.5d0*u_inv*sqpi*erf(u) *d1(0)
    endif
    
    do i=2,n_pt_out,2
      n = shiftr(i,1)
      rint_sum = rint_sum + d1(i)*rint1(n,rho)
    enddo
    
  else
    
    if(rho.gt.80.d0)then
     v=0.d0
    else
     v=dexp(-rho)
    endif
     
    u_inv=1.d0/dsqrt(rho)
    u=rho*u_inv
    two_rho_inv = 0.5d0*u_inv*u_inv
    val0=0.5d0*u_inv*sqpi*erf(u)
    rint_sum=val0*d1(0)
    rint_tmp=(val0-v)*two_rho_inv
    di = 3.d0
    do i=2,min(n_pt_out,40),2
      rint_sum = rint_sum + d1(i)*rint_tmp
      rint_tmp = (rint_tmp*di-v)*two_rho_inv
      di = di+2.d0
    enddo
    do i=42,n_pt_out,2
      n = shiftr(i,1)
      rint_sum = rint_sum + d1(i)*rint_large_n(n,rho)
    enddo
    
  endif
end

double precision function hermite(n,x)
  implicit none
!  BEGIN_DOC
! Hermite polynomial
!  END_DOC
  integer                        :: n,k
  double precision               :: h0,x,h1,h2
  h0=1.d0
  if(n.eq.0)then
    hermite=h0
    return
  endif
  h1=x+x
  if(n.eq.1)then
    hermite=h1
    return
  endif
  do k=1,n-1
    h2=(x+x)*h1-dfloat(k+k)*h0
    h0=h1
    h1=h2
  enddo
  hermite=h2
end

double precision function rint_large_n(n,rho)
  implicit none
!  BEGIN_DOC
! Version of rint for large values of n
!  END_DOC
  integer                        :: n,k,l
  double precision               :: rho,u,accu,eps,t1,t2,fact,alpha_k,rajout,hermite
  u=dsqrt(rho)
  accu=0.d0
  k=0
  eps=1.d0
  do while (eps.gt.1.d-15)
    t1=1.d0
    do l=0,k
      t1=t1*(n+n+l+1.d0)
    enddo
    t2=0.d0
    do l=0,k
      t2=t2+(-1.d0)**l/(fact(l+1)*fact(k-l))
    enddo
    alpha_k=t2*fact(k+1)*fact(k)*(-1.d0)**k
    alpha_k= alpha_k/t1
    rajout=(-1.d0)**k*u**k*hermite(k,u)*alpha_k/fact(k)
    accu=accu+rajout
    eps=dabs(rajout)/accu
    k=k+1
  enddo
  rint_large_n=dexp(-rho)*accu
end


double precision function rint1(n,rho)
  implicit none
!  BEGIN_DOC
! Standard version of rint
!  END_DOC
  integer, intent(in)            :: n
  double precision, intent(in)   :: rho
!!!!!!  double precision, parameter    :: eps=1.d-13
  double precision, parameter    :: eps=1.d-15
  double precision               :: rho_tmp, diff
  double precision               :: fact_inv
  integer                        :: k
  double precision, parameter, dimension(100) :: inv_int = (/ 1.d0, 1.d0/2.d0, 1.d0/3.d0, &
      1.d0/4.d0, 1.d0/5.d0, 1.d0/6.d0, 1.d0/7.d0, 1.d0/8.d0, 1.d0/9.d0, 1.d0/10.d0, &
      1.d0/11.d0, 1.d0/12.d0, 1.d0/13.d0, 1.d0/14.d0, 1.d0/15.d0, 1.d0/16.d0, 1.d0/17.d0, &
      1.d0/18.d0, 1.d0/19.d0, 1.d0/20.d0, 1.d0/21.d0, 1.d0/22.d0, 1.d0/23.d0, 1.d0/24.d0, &
      1.d0/25.d0, 1.d0/26.d0, 1.d0/27.d0, 1.d0/28.d0, 1.d0/29.d0, 1.d0/30.d0, 1.d0/31.d0, &
      1.d0/32.d0, 1.d0/33.d0, 1.d0/34.d0, 1.d0/35.d0, 1.d0/36.d0, 1.d0/37.d0, 1.d0/38.d0, &
      1.d0/39.d0, 1.d0/40.d0, 1.d0/41.d0, 1.d0/42.d0, 1.d0/43.d0, 1.d0/44.d0, 1.d0/45.d0, &
      1.d0/46.d0, 1.d0/47.d0, 1.d0/48.d0, 1.d0/49.d0, 1.d0/50.d0, 1.d0/51.d0, 1.d0/52.d0, &
      1.d0/53.d0, 1.d0/54.d0, 1.d0/55.d0, 1.d0/56.d0, 1.d0/57.d0, 1.d0/58.d0, 1.d0/59.d0, &
      1.d0/60.d0, 1.d0/61.d0, 1.d0/62.d0, 1.d0/63.d0, 1.d0/64.d0, 1.d0/65.d0, 1.d0/66.d0, &
      1.d0/67.d0, 1.d0/68.d0, 1.d0/69.d0, 1.d0/70.d0, 1.d0/71.d0, 1.d0/72.d0, 1.d0/73.d0, &
      1.d0/74.d0, 1.d0/75.d0, 1.d0/76.d0, 1.d0/77.d0, 1.d0/78.d0, 1.d0/79.d0, 1.d0/80.d0, &
      1.d0/81.d0, 1.d0/82.d0, 1.d0/83.d0, 1.d0/84.d0, 1.d0/85.d0, 1.d0/86.d0, 1.d0/87.d0, &
      1.d0/88.d0, 1.d0/89.d0, 1.d0/90.d0, 1.d0/91.d0, 1.d0/92.d0, 1.d0/93.d0, 1.d0/94.d0, &
      1.d0/95.d0, 1.d0/96.d0, 1.d0/97.d0, 1.d0/98.d0, 1.d0/99.d0, 1.d0/100.d0 /)

  if (n<50) then
    rint1=inv_int(n+n+1)
    rho_tmp = 1.d0
    do k=1,20
      rho_tmp = -rho_tmp*rho
      diff=rho_tmp*fact_inv(k)*inv_int(ishft(k+n,1)+1)

      rint1=rint1+diff
      if (dabs(diff) > eps) then
        cycle
      endif
      return
    enddo
  else
    rint1=1.d0/dfloat(n+n+1)
    rho_tmp = 1.d0
    do k=1,20
      rho_tmp = -rho_tmp*rho
      diff=rho_tmp*fact_inv(k)/dfloat(shiftl(k+n,1)+1)
      rint1=rint1+diff
      if (dabs(diff) > eps) then
        cycle
      endif
      return
    enddo
  endif

   write(*,*)'n rhod=',n,rho
   write(*,*)'diff=',diff,' pb in rint1 k too large!'
  stop 1
end

double precision function gauleg_t2(m,n)
   implicit none
   integer, intent(in) :: m,n
   !EGIN_DOC
   ! t_w(i,1,k) = w(i)
   ! t_w(i,2,k) = t(i)
   !ND_DOC
   integer                        :: i,j,l
   double precision :: tmp_1(1000,1000),tmp_2(1000,1000)
   integer :: imax
   l=0
   imax = max(m,n)
   imax = imax * 2
   do i = 2,imax,2
     l = l+1
!    call gauleg(0.d0,1.d0,gauleg_t2(1,l),gauleg_w(1,l),i)
     call gauleg(0.d0,1.d0,tmp_1(1,l),tmp_2(1,l),i)
     do j=1,i
       tmp_1(j,l) = tmp_1(j,l) * tmp_1(j,l)
     enddo
   enddo
   gauleg_t2 = tmp_1(m,n)

end


double precision function gauleg_w(m,n)
   implicit none
   integer, intent(in) :: m,n
   !EGIN_DOC
   ! t_w(i,1,k) = w(i)
   ! t_w(i,2,k) = t(i)
   !ND_DOC
   integer                        :: i,j,l
   double precision :: tmp_1(1000,1000),tmp_2(1000,1000)
   integer :: imax
   l=0
   imax = max(m,n)
   imax = imax * 2
   do i = 2,imax,2
     l = l+1
!    call gauleg(0.d0,1.d0,gauleg_t2(1,l),gauleg_w(1,l),i)
     call gauleg(0.d0,1.d0,tmp_1(1,l),tmp_2(1,l),i)
     do j=1,i
       tmp_1(j,l) = tmp_1(j,l) * tmp_1(j,l)
     enddo
   enddo
   gauleg_w= tmp_2(m,n)

end


subroutine gauleg(x1,x2,x,w,n)
   implicit none
   !EGIN_DOC
   ! Gauss-Legendre
   !ND_DOC
   integer, intent(in)            :: n
   double precision, intent(in)   :: x1, x2
   double precision, intent (out) :: x(n),w(n)
   double precision, parameter    :: eps=3.d-14

   integer                        :: m,i,j
   double precision               :: xm, xl, z, z1, p1, p2, p3, pp, dn
   m=(n+1)/2
   xm=0.5d0*(x2+x1)
   xl=0.5d0*(x2-x1)
   dn = dble(n)
   do i=1,m
     z=dcos(3.141592654d0*(dble(i)-.25d0)/(dble(n)+.5d0))
     z1 = z+1.d0
     do while (dabs(z-z1) > eps)
       p1=1.d0
       p2=0.d0
       do j=1,n
         p3=p2
         p2=p1
         p1=(dble(j+j-1)*z*p2-dble(j-1)*p3)/j
       enddo
       pp=dn*(z*p1-p2)/(z*z-1.d0)
       z1=z
       z=z1-p1/pp
     end do
     x(i)=xm-xl*z
     x(n+1-i)=xm+xl*z
     w(i)=(xl+xl)/((1.d0-z*z)*pp*pp)
     w(n+1-i)=w(i)
   enddo
end

!***** various simple routines

double precision function binom_func(i,j)
  implicit none
  !EGIN_DOC
  !.. math                       :: 
  !
  !  \frac{i!}{j!(i-j)!}
  !
  !ND_DOC
  integer,intent(in)             :: i,j
  double precision               :: logfact
  integer, save                  :: ifirst
  double precision, save         :: memo(0:15,0:15)
  integer                        :: k,l
  if (ifirst == 0) then
    ifirst = 1
    do k=0,15
      do l=0,15
        memo(k,l) = dexp( logfact(k)-logfact(l)-logfact(k-l) )
      enddo
    enddo
  endif
  if ( (i<=15).and.(j<=15) ) then
    binom_func = memo(i,j)
  else
    binom_func = dexp( logfact(i)-logfact(j)-logfact(i-j) )
  endif
end

double precision function binom_transp(i,j)
implicit none
integer,intent(in)             :: i,j
double precision               :: binom_func
binom_transp=binom_func(j,i)
end

double precision function logfact(n)
  implicit none
  !EGIN_DOC
  ! n!
  !ND_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1
  integer                        :: i
  
  if (n<=memomax) then
    if (n<2) then
      logfact = 0.d0
    else
      logfact = memo(n)
    endif
    return
  endif
  
  memo(1) = 0.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)+dlog(dble(i))
  enddo
  memomax = min(n,100)
  logfact = memo(memomax)
  do i=101,n
    logfact = logfact + dlog(dble(i))
  enddo
end function

double precision function fact_inv(i)
 implicit none
 integer :: i
 double precision :: fact
 fact_inv = 1.d0/fact(i)
end

double precision function dble_fact(n)
  implicit none
  integer :: n
  double precision :: dble_fact_even, dble_fact_odd

  dble_fact = 1.d0

  if(n.lt.0) return

  if(iand(n,1).eq.0)then
    dble_fact = dble_fact_even(n)
  else
    dble_fact= dble_fact_odd(n)
  endif

end function

double precision function dble_fact_even(n) result(fact2)
  implicit none
  !EGIN_DOC
  ! n!!
  !ND_DOC
  integer                        :: n,k
  double precision, save         :: memo(0:100)
  integer, save                  :: memomax = 0
  double precision               :: prod
  integer                        :: i
  double precision :: dble_logfact

  if (n <= memomax) then
    if (n < 2) then
      fact2 = 1.d0
    else
      fact2 = memo(n)
    endif
    return
  endif

  memo(0)=1.d0
  memo(1)=1.d0
  do i=memomax+2,min(n,100),2
    memo(i) = memo(i-2)* dble(i)
  enddo
  memomax = min(n,100)
  fact2 = memo(memomax)
  
  if (n > 100) then
    fact2 = dexp(dble_logfact(n))
  endif

end function

double precision function dble_fact_odd(n) result(fact2)
  implicit none
  !EGIN_DOC
  ! n!!
  !ND_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1
  integer                        :: i
  double precision :: dble_logfact
  
  if (n<=memomax) then
    if (n<3) then
      fact2 = 1.d0
    else
      fact2 = memo(n)
    endif
    return
  endif
  
  memo(1) = 1.d0
  do i=memomax+2,min(n,99),2
    memo(i) = memo(i-2)* dble(i)
  enddo
  memomax = min(n,99)
  fact2 = memo(memomax)
  
  if (n > 99) then
    fact2 = dexp(dble_logfact(n))
  endif

end function

double precision function dble_logfact(n) result(logfact2)
  implicit none
  !EGIN_DOC
  ! n!!
  !ND_DOC
  integer                        :: n
  integer :: k
  double precision :: prod
  prod=0.d0
  do k=2,n,2
   prod=prod+dlog(dfloat(k))
  enddo
  logfact2=prod
  return
  
end function


double precision function fact(n)
  implicit none
  !EGIN_DOC
  ! n!
  !ND_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1
  integer                        :: i
  double precision :: logfact
  
  if (n<=memomax) then
    if (n<2) then
      fact = 1.d0
    else
      fact = memo(n)
    endif
    return
  endif
  
  memo(1) = 1.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)*dble(i)
  enddo
  memomax = min(n,100)
  fact = dexp(logfact(n))
end function

!******* 
! what follows erf function needed in rint

      double precision FUNCTION ERF(X)
      implicit double precision(a-h,o-z)
      IF(X.LT.0.d0)THEN
        ERF=-GAMMP(.5d0,X**2)
      ELSE
        ERF=GAMMP(.5d0,X**2)
      ENDIF
      RETURN
      END

      double precision FUNCTION GAMMLN(XX)
      implicit double precision(a-h,o-z)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0, &
          -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
      FUNCTION GAMMP(A,X)
      implicit double precision(a-h,o-z)
      IF(X.LT.0.d0.OR.A.LE.0.d0)STOP
      IF(X.LT.A+1.d0)THEN
        CALL GSER(GAMMP,A,X,GLN)
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.d0-GAMMCF
      ENDIF
      RETURN
      END
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      implicit double precision(a-h,o-z)
      PARAMETER (ITMAX=100,EPS=3.D-7)
      GLN=GAMMLN(A)
      GOLD=0.d0
      A0=1.d0
      A1=X
      B0=0.d0
      B1=1.d0
      FAC=1.d0
      DO 11 N=1,ITMAX
        AN=DFLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.d0)THEN
          FAC=1.d0/A1
          G=B1*FAC
          IF(DABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      STOP 'A TOO LARGE, ITMAX TOO SMALL'
1     GAMMCF=DEXP(-X+A*LOG(X)-GLN)*G
      RETURN
      END
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      implicit double precision(a-h,o-z)
      PARAMETER (ITMAX=100,EPS=3.D-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.d0)THEN
        IF(X.LT.0.d0)STOP ' X LT 0'
        GAMSER=0.d0
        RETURN
      ENDIF
      AP=A
      SUM=1.d0/A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.d0
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(DABS(DEL).LT.DABS(SUM)*EPS)GO TO 1
11    CONTINUE
      STOP 'A TOO LARGE, ITMAX TOO SMALL'
1     GAMSER=SUM*DEXP(-X+A*LOG(X)-GLN)
      RETURN
      END

!******************
!*****************
! Monte Carlo part for checking


subroutine integral_pqrs_via_MC  &
(nmc,a,nA,gA,b,nB,gB,c,nC,gC,d,nD,gD,rmoy,error)

implicit none
double precision gA,gB,gC,gD,gApB,gCpD,sig1,sig2,gauss,rmoy,error
integer :: nmc,nblocks,i,l,kk,n1,n2
double precision xpi,dpi,fact,fnorm,pnorm,r12
double precision x1,y1,z1,x2,y2,z2
double precision x10,y10,z10,x20,y20,z20
double precision a(3),b(3),c(3),d(3),GAB(3),GCD(3),DG(3)
double precision accu(50),factor
double precision Poly1,Poly2,eab,ecd,amb2,cmd2
integer nA(3),nB(3),nC(3),nD(3)

!******************************************************************
! Data for gaussian generator:
      xpi=dacos(-1.d0)
      dpi=2.d0*xpi
      fact=1.d0/dsqrt(2.d0)
!*******************************************************************

gApB=gA+gB
gCpD=gC+gD
sig1=1.d0/dsqrt(2.d0*gApB)
sig2=1.d0/dsqrt(2.d0*gCpD)
factor=sig1**3*sig2**3*dpi**3

do l=1,3
 GAB(l)=(gA*a(l)+gB*b(l))/gApB
 GCD(l)=(gC*c(l)+gD*d(l))/gCpD
 DG(l)=GAB(l)-GCD(l)
enddo

amb2=(a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2
cmd2=(c(1)-d(1))**2+(c(2)-d(2))**2+(c(3)-d(3))**2
eab=dexp(-gA*gB/gApB *amb2)
ecd=dexp(-gC*gD/gCpD *cmd2)
factor=factor*eab*ecd

nblocks=50

do kk=1,nblocks

accu(kk)=0.d0

do i=1,nmc

x1=gauss(sig1)
y1=gauss(sig1)
z1=gauss(sig1)
x2=gauss(sig2)
y2=gauss(sig2)
z2=gauss(sig2)
Poly1=( x1-(a(1)-GAB(1)) )**nA(1) *( x1-(b(1)-GAB(1)))**nB(1)* &
      ( y1-(a(2)-GAB(2)) )**nA(2) *( y1-(b(2)-GAB(2)))**nB(2)* &
      ( z1-(a(3)-GAB(3)) )**nA(3) *( z1-(b(3)-GAB(3)))**nB(3)

Poly2=( x2-(c(1)-GCD(1)) )**nC(1) *( x2-(d(1)-GCD(1)))**nD(1)* &
      ( y2-(c(2)-GCD(2)) )**nC(2) *( y2-(d(2)-GCD(2)))**nD(2)* &
      ( z2-(c(3)-GCD(3)) )**nC(3) *( z2-(d(3)-GCD(3)))**nD(3)

r12=dsqrt((x1-x2 +DG(1) )**2+(y1-y2 +DG(2))**2+(z1-z2 +DG(3))**2)

accu(kk)=accu(kk)+Poly1*Poly2/r12

enddo
accu(kk)=accu(kk)/nmc*factor
enddo

call erreur(accu,nblocks,rmoy,error)

end

      double precision function gauss(sigma)
      implicit double precision(a-h,o-z)
      dpi=2.d0*dacos(-1.d0)
      call random_number(r1)
      call random_number(r2)
      gauss=sigma*dsqrt(-2.d0*dlog(r1)*dcos(dpi*r2))
      end

!c  calcul de la moyenne / erreur
!c ------------------------------
      subroutine erreur(x,n,rmoy,error)
      implicit double precision(a-h,o-z)
      dimension x(n)
! calcul de la moyenne
      rmoy=0.d0
      do i=1,n
        rmoy=rmoy+x(i)
      enddo
      rmoy=rmoy/dfloat(n)
! calcul de l'erreur
      error=0.d0
      do i=1,n
       error=error+(x(i)-rmoy)**2
      enddo
      if(n.gt.1)then
        rn=dfloat(n)
        rn1=dfloat(n-1)
        error=dsqrt(error)/dsqrt(rn*rn1)
      else
        write(2,*)'Seulement un block Erreur nondefinie'
        error=0.d0
      endif
      end

subroutine read_geometry(MOLECULE)
include 'j.inc'
character*80 MOLECULE
character*80 charabid,filename

conversion_ang_to_ua=1.d0/0.52917721092d0 ! drawn from toto
filename=trim(MOLECULE)//'.xyz'

open(2,file=filename)
read(2,*)number_nuclei
if(number_nuclei.gt.number_nuclei_max)stop 'increase number_nuclei_max'
read(2,'(a80)')charabid

do i=1,number_nuclei
 read(2,*)ATOM(i),centers_nuclei(1,i),centers_nuclei(2,i),centers_nuclei(3,i)
 charge(i)=number_atom(ATOM(i))
 do l=1,3
  centers_nuclei(l,i)=conversion_ang_to_ua*centers_nuclei(l,i)
 enddo
enddo

!! Calculation of the nuclear energy'
!*******************************
enucl=0.d0
do i=1,number_nuclei
 do j=i+1,number_nuclei

  rij=dsqrt( (centers_nuclei(1,i)-centers_nuclei(1,j))**2    &
            +(centers_nuclei(2,i)-centers_nuclei(2,j))**2   &
            +(centers_nuclei(3,i)-centers_nuclei(3,j))**2 )
  enucl=enucl+charge(i)*charge(j)/rij
 enddo
enddo
end

!! iread_c0=0: initial vectors are built from H0
!! iread_c0=1: initial vectors are read in file 'c0'
!! At the end of SCF final vectors are written in 'c0'

!! bielec used to make SCF:
!!
!! ifind_kcp_red(i,j,k,l) gives the kcp in ijkl_red() that gives the bielec for i,j,k,l
!! biel_c=ijkl_red(ifind_kcp_red(i,k,j,l))
!! biel_e=ijkl_red(ifind_kcp_red(i,j,k,l))


subroutine SCF(nocc,S_ij,V_ij, K_ij,ijkl,iread_c0,ehf,niter_SCF)
include 'j.inc'

double precision ijkl(nint_max)
double precision S_ij(nbasis_max,nbasis_max),V_ij(nbasis_max,nbasis_max), &
                 K_ij(nbasis_max,nbasis_max)

double precision h(nbasis_max,nbasis_max),s(nbasis_max,nbasis_max)
double precision h0(nbasis_max,nbasis_max),s0(nbasis_max,nbasis_max)
double precision pot(nbasis_max,nbasis_max)
double precision hs(nbasis_max,nbasis_max)
double precision cold(nbasis_max,nbasis_max)
double precision cnew(nbasis_max,nbasis_max)
double precision w(nbasis_max)

dimension rho(nbasis_max,nbasis_max)

logical ZV
character*80 charabid

call cpu_time(t0)

do i=1,nbasis
 do j=1,nbasis
  h0(i,j)=K_ij(i,j)+V_ij(i,j)
  s0(i,j)=S_ij(i,j)
 enddo
enddo

if(iread_c0.eq.0)then
 print*,'BUILD INITIAL VECTORS'
 call diag_nonortho_basis(s0,h0,w,cold)
else
 print*,'READ INITIAL VECTORS in file c0'
 open(55,file='c0')
 rewind 55
 do i=1,nocc
  read(55,'(a80)')charabid
  do k=1,nbasis
   read(55,*)cold(k,i)
  enddo
 enddo
 close(55)
endif

do k=1,nbasis
 do l=1,nbasis
  rho(k,l)=0.d0
  do i=1,nocc
   rho(k,l)=rho(k,l)+cold(k,i)*cold(l,i)
  enddo
 enddo
enddo

do iter=1,niter_SCF

 do i=1,nbasis
  do j=1,nbasis

   coulomb=0.d0
   exchang=0.d0

   do k=1,nbasis
    do l=1,nbasis

     if(ifind_kcp(i,k,j,l).eq.0)then
      biel_c=0.d0
     else
      biel_c=ijkl(ifind_kcp(i,k,j,l))
     endif
     if(ifind_kcp(i,j,k,l).eq.0)then
      biel_e=0.d0
     else
      biel_e=ijkl(ifind_kcp(i,j,k,l))
     endif

     coulomb=coulomb+rho(k,l)*biel_c
     exchang=exchang+rho(k,l)*biel_e

    enddo !l
   enddo !k

   pot(i,j)=2.d0*coulomb-exchang
   h(i,j)=h0(i,j)+pot(i,j)

  enddo!j
 enddo !i

 ehf=e_hartree_fock(h0,ijkl,rho)
 write(*,*)'ehf=',ehf

 do i=1,nbasis
  do j=1,nbasis
   hs(i,j)=0.5d0*(h(i,j)+h(j,i))
  enddo
 enddo

 call diag_nonortho_basis(s0,hs,w,cnew)

 do k=1,nbasis
  do l=1,nbasis
   rho(k,l)=0.d0
   do i=1,nocc
    rho(k,l)=rho(k,l)+0.8d0*cnew(k,i)*cnew(l,i)+0.2d0*cold(k,i)*cold(l,i)
   enddo
  enddo
 enddo

 do k=1,nbasis
  do i=1,nocc
   cold(k,i)=cnew(k,i)
  enddo
 enddo

enddo !iter

open(55,file='c0')
rewind 55
print*,'WRITE FINAL VECTORS in file c0'
do i=1,nocc
 write(55,*)'occupied orbital',' #',i, 'epsi=',w(i)
 do k=1,nbasis
  write(55,*)cnew(k,i)
 enddo
enddo
close(55)

do k=1,nbasis
 do l=1,nbasis
  rho(k,l)=0.d0
  do i=1,nocc
   rho(k,l)=rho(k,l)+cnew(k,i)*cnew(l,i)
  enddo
 enddo
enddo

end

double precision function e_hartree_fock(h0,ijkl,rho)
include 'j.inc'

double precision h0(nbasis_max,nbasis_max)
dimension rho(nbasis_max,nbasis_max)
double precision ijkl(nint_max)

e0=0.d0
do i=1,nbasis
 do j=1,nbasis
  e0=e0+2.d0*rho(i,j)*h0(i,j)
 enddo
enddo

e1=0.d0
do i=1,nbasis
 do j=1,nbasis
  do k=1,nbasis
   do l=1,nbasis
    e1=e1+(2.d0*rho(i,j)*rho(k,l)-rho(i,k)*rho(j,l))*ijkl(ifind_kcp(i,k,j,l))
   enddo
  enddo
 enddo
enddo

e_hartree_fock=e0+e1+enucl

end

      subroutine diag_nonortho_basis(s,h,d,z)
      include 'j.inc'
      dimension s(nbasis_max,nbasis_max)
      dimension prod(nbasis_max,nbasis_max)
      dimension htilde(nbasis_max,nbasis_max)
      dimension smud(nbasis_max,nbasis_max)
      dimension d(nbasis_max),v(nbasis_max,nbasis_max)
      dimension h(nbasis_max,nbasis_max),z(nbasis_max,nbasis_max)
      dimension h0(nbasis_max,nbasis_max),s0(nbasis_max,nbasis_max)

      if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'

      do i=1,nbasis
       do j=1,nbasis
        s0(i,j)=s(i,j)
        h0(i,j)=h(i,j)
       enddo
      enddo

!! Diagonalisation of overlap S
!!
      call jacobi(s,nbasis,nbasis_max,d,v,nrot)
      call eigsrt(d,v,nbasis,nbasis_max)

!      print*,'vp S min=',d(1)
!      print*,'vp S max=',d(nbasis)

!! Construction of matrix S-1/2
!!
      do i=1,nbasis
       do j=1,nbasis
        smud(i,j)=0.d0
         do k=1,nbasis
          smud(i,j)=smud(i,j)+1.d0/dsqrt(dabs(d(k)))*v(i,k)*v(j,k)
         enddo
       enddo
      enddo

!! Calculation of   H_tilde= S-1/2 H S-1/2
!!
      call multiply(nbasis,nbasis_max,smud,h,prod)
      call multiply(nbasis,nbasis_max,prod,smud,htilde)

      call jacobi(htilde,nbasis,nbasis_max,d,v,nrot)

      call multiply(nbasis,nbasis_max,smud,v,z)
      call eigsrt(d,z,nbasis,nbasis_max)

     do i=1,nbasis
       do j=1,nbasis
        s(i,j)=s0(i,j)
        h(i,j)=h0(i,j)
       enddo
      enddo

      end

      subroutine multiply(n,np,a,b,c)
      implicit double precision(a-h,o-z)
      dimension a(np,np),b(np,np),c(np,np)
      do i=1,n
       do j=1,n
        c(i,j)=0.d0
        do k=1,n
          c(i,j)=c(i,j)+a(i,k)*b(k,j)
        enddo
       enddo
      enddo
      end

      subroutine jacobi(a,n,np,d,v,nrot)
      implicit double precision (a-h,o-z)
      parameter (nmax=1000)
      dimension a(np,np),d(np),v(np,np),b(nmax),z(nmax)
      if(nmax.lt.np)stop 'increase nmax'

      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+dabs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*dabs(a(ip,iq))
            if((i.gt.4).and.(dabs(d(ip))+g.eq.dabs(d(ip))) &
              .and.(dabs(d(iq))+g.eq.dabs(d(iq))))then
              a(ip,iq)=0.
            else if(dabs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(dabs(h)+g.eq.dabs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(dabs(theta)+dsqrt(1.d0+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./dsqrt(1.d0+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      stop '50 iterations should never happen'
      return
      end

      subroutine eigsrt(d,v,n,np)
      implicit double precision (a-h,o-z)
      dimension d(np),v(np,np)
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).lt.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      end

      subroutine indexx(n,arrin,indx)
      implicit double precision (a-h,o-z)
      dimension arrin(n),indx(n)
      do 11 j=1,n
        indx(j)=j
11    continue
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
        else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1)then
            indx(1)=indxt
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
          endif
          if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        indx(i)=indxt
      go to 10
      end

     SUBROUTINE qsort(INDEX,N,DATA)
     implicit none
!===================================================================
!
!     SORTRX -- SORT, Real input, indeX output
!
!
!     Input:  N     INTEGER
!             DATA  double precision
!
!     Output: INDEX INTEGER (DIMENSION N)
!
! This routine performs an in-memory sort of the first N elements of
! array DATA, returning into array INDEX the indices of elements of
! DATA arranged in ascending order.  Thus,
!
!    DATA(INDEX(1)) will be the smallest number in array DATA;
!    DATA(INDEX(N)) will be the largest number in DATA.
!
! The original data is not physically rearranged.  The original order
! of equal input values is not necessarily preserved.
!
!===================================================================
!
! SORTRX uses a hybrid QuickSort algorithm, based on several
! suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
! "pivot key" [my term] for dividing each subsequence is chosen to be
! the median of the first, last, and middle values of the subsequence;
! and the QuickSort is cut off when a subsequence has 9 or fewer
! elements, and a straight insertion sort of the entire array is done
! at the end.  The result is comparable to a pure insertion sort for
! very short arrays, and very fast for very large arrays (of order 12
! micro-sec/element on the 3081K for arrays of 10K elements).  It is
! also not subject to the poor performance of the pure QuickSort on
! partially ordered data.
!
! Created:  15 Jul 1986  Len Moss
!
!===================================================================
 
!      INTEGER   N,INDEX(:) ! (N)
!      double precision  DATA(:) ! (N)

      INTEGER   N,INDEX(N)
      double precision  DATA(N)

 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      double precision       DATAP
 
!     QuickSort Cutoff
!
!     Quit QuickSort-ing when a subsequence contains M or fewer
!     elements and finish off at end with straight insertion sort.
!     According to Knuth, V.3, the optimum value of M is around 9.
 
      INTEGER   M
      PARAMETER (M=9)
 
!===================================================================
!
!     Make initial guess for INDEX
 
      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE
 
!     If array is short, skip QuickSort and go directly to
!     the straight insertion sort.
 
      IF (N.LE.M) GOTO 900
 
!===================================================================
!
!     QuickSort
!
!     The "Qn:"s correspond roughly to steps in Algorithm Q,
!     Knuth, V.3, PP.116-117, modified to select the median
!     of the first, last, and middle elements as the "pivot
!     key" (in Knuth's notation, "K").  Also modified to leave
!     data in place and produce an INDEX array.  To simplify
!     comments, let DATA[I]=DATA(INDEX(I)).
 
! Q1: Initialize
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
! Q2: Sort the subsequence DATA[L]..DATA[R].
!
!     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
!     r > R, and L <= m <= R.  (First time through, there is no
!     DATA for l < L or r > R.)
 
      I=L
      J=R
 
! Q2.5: Select pivot key
!
!     Let the pivot, P, be the midpoint of this subsequence,
!     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
!     so the corresponding DATA values are in increasing order.
!     The pivot key, DATAP, is then DATA[P].
 
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
 
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
!     Now we swap values between the right and left sides and/or
!     move DATAP until all smaller values are on the left and all
!     larger values are on the right.  Neither the left or right
!     side will be internally ordered yet; however, DATAP will be
!     in its final position.
 
  300 CONTINUE
 
! Q3: Search for datum on left >= DATAP
!
!     At this point, DATA[L] <= DATAP.  We can therefore start scanning
!     up from L, looking for a value >= DATAP (this scan is guaranteed
!     to terminate since we initially placed DATAP near the middle of
!     the subsequence).
 
         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
! Q4: Search for datum on right <= DATAP
!
!     At this point, DATA[R] >= DATAP.  We can therefore start scanning
!     down from R, looking for a value <= DATAP (this scan is guaranteed
!     to terminate since we initially placed DATAP near the middle of
!     the subsequence).
 
         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
 
! Q5: Have the two scans collided?
 
      IF (I.LT.J) THEN
 
! Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
 
! Q7: Yes, select next subsequence to sort
!
!     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
!     for all L <= l < I and J < r <= R.  If both subsequences are
!     more than M elements long, push the longer one on the stack and
!     go back to QuickSort the shorter; if only one is more than M
!     elements long, go back and QuickSort it; otherwise, pop a
!     subsequence off the stack and QuickSort it.
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
! Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
!===================================================================
!
! Q9: Straight Insertion sort
 
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
!===================================================================
!
!     All done
 
      END SUBROUTINE qsort

! Computation of I= (nlm_1 nlm_2 | nlm_3 nlm_4)= /int d vector_1 d vector_2  (n1 l1 m1)*(r1)(n2 l2 m2)(r1) 1/r12 (n3 l3 m3)*(r2) (n4 l4 m4)(r2) 
!!
!! where (nlm)*(vector r) = N_n r^{n-1} exp(-g r) Y_l^m(theta,phi)   N_n=sqrt((2g)**(2n+1)/(2n)!) 
!! with Y_l^m(theta,phi) = i^(m+|m|) ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2 P_l^|m|(cos(theta))  exp(i m phi)
!!
!! l=0,1,2,...,n-1
!! m=-l...,l
!!
!! Here P_l^|m|(cos(theta)) = "Associated Legendre Polynomials wikipedia" computed with routine:  SUBROUTINE LPMN(MM,M,N,X,PM)
!!
!! In input of LPMN:  n=l (n=0,1,...)  m=0,1,...,n x=cos(theta) 0 (< or =)  x (< ! or =) 1
!!
!! the routine computes:   PM(m,n) for n=0,...,N (number N in input) and  m=0,..,n
!!
!! Exemples (see 'Associated Legendre Polynomilas wikipedia')
!!   P_{0}^{0}(x)=1
!!   P_{1}^{-1}(x)=-1/2 P_{1}^{1}(x)
!!   P_{1}^{0}(x)=x
!!   P_{1}^{1}(x)=-(1-x^2)^{1/2}
!!   P_{2}^{-2}(x)=1/24 P_{2}^{2}(x)
!!   P_{2}^{-1}(x)=-1/6 P_{2}^{1}(x)
!!   P_{2}^{0}(x)=1/2 (3x^{2}-1)
!!   P_{2}^{1}(x)=-3x(1-x^2)^{1/2}
!!   P_{2}^{2}(x)=3(1-x^2)
!!
!!
!!*************
!!  Conditions: l1+l2+l3+l4 even and  m2-m1=m3-m4  otherwise I=0
!!*************
!!
!! I= N_n1*N_n2*N_n3*N_n4  \sum_{lmin}^{lmax},2  4*pi/(2l+1)  [G_{l1m1 lm}^{l2m2} G_{l1m1 lm}^{l2m2}] * bigF(n1+n2,g1+g2,n3+n4,g3+g4,l)
!!
!!  Gaunt coefficient:  G_{l1m1 l2m2}^lm =  int_dOmega  Ylm^*  Yl1m1  Yl2m2       m=m1+m2 (0 otherwise)
!!
!!  lmin=max(l12_min, l34_min)  
!!  lmax=min(l12_max, l34_max)  
!!
!!  lij_max=li+lj   
!!  lij_min=  max(|li-lj|,|mi-mj|)     if [lij_max + max(|li-lj|,|mi-mj|)] = even
!!  lij_min=  max(|li-lj|,|mi-mj|)+1   if [lij_max + max(|li-lj|,|mi-mj|)] = odd
!!  
!!
!! bifF defined as:
!!
!! bigF(k1,g1,k2,g2,l)= int_0^inf r1^k1 exp(-g1 r1) [1/r1**(l+1) int_0^r1 r2**(k2+l) exp(-g2 r2) + r1**l int_r1^inf r2**(k2-l-1) exp(-g2 r2)]
!!
!! k1 integer larger or equal to 0 , k2 integer larger or equal to 0, k1 and k2 larger and equal to (l+1)
!!
!! We can show that
!!
!! bigF=  [(k1-l-1)!(k2+l)!/g1**(k1-l)/g2**(k2+l+1)]
!!      + \sum_{i=0}^{k2-l-1} [(k2-l-1)!(k1+k2-1-i)!/(k2-l-1-i)!/g2**(i+1)/(g1+g2)**(k1+k2-i)]
!!      - \sum_{i=0}^{k2+l  } [(k2+l)!(k1+k2-1-i)!/(k2+l-i)!/g2**(i+1)/(g1+g2)**(k1+k2-i)]
!!
!! This formula can been checked with bigF_num
!!
!!
double precision function ijkl_slater(nslat,lslat,mslat,gamma)
implicit none
integer nslat(4),lslat(4),mslat(4),i,l,lmin,lmax
integer l1,l2,l3,l4,m1,m2,m3,m4,n12,n34,m
double precision gamma(4),pi,int,norm,norm_sto,gaunt,gamm12,gamm34,bigF

ijkl_slater=0.d0

l1=lslat(1)
l2=lslat(2)
l3=lslat(3)
l4=lslat(4)
if(mod(l1+l2+l3+l4,2).eq.1)return

m1=mslat(1)
m2=mslat(2)
m3=mslat(3)
m4=mslat(4)
if((m2-m1).ne.(m3-m4))then
 return
else
 m=m2-m1
endif

pi=dacos(-1.d0)

norm=1.d0
do i=1,4
 norm=norm*norm_sto(gamma(i),nslat(i))
enddo

call compute_lmin_lmax(l1,l2,l3,l4,m1,m2,m3,m4,lmin,lmax)

n12=nslat(1)+nslat(2)
n34=nslat(3)+nslat(4)
gamm12=gamma(1)+gamma(2)
gamm34=gamma(3)+gamma(4)

int=0.d0
do l=lmin,lmax,2
 int=int+gaunt(l2,l,l1,m2,m,m1)*gaunt(l3,l,l4,m3,m,m4)*bigF(n12,gamm12,n34,gamm34,l)/(2.d0*l+1.d0)
enddo
ijkl_slater=4.d0*pi*int*norm

end

!! Computation of overlap (n1l1m1 | n2l2m2)=  N_n1*N_n2 !(n1+n2)!(g1+g2)**(n1+n2+1)   l1=l1 m1=m2
!!
double precision function overlap_slater(nslat,lslat,mslat,gamma)
implicit none
integer nslat(4),lslat(4),mslat(4),i
integer l1,l2,m1,m2,n12
double precision gamma(4),int,norm,norm_sto,g12,fact


overlap_slater=0.d0

l1=lslat(1)
l2=lslat(2)
if(l1.ne.l2)return
m1=mslat(1)
m2=mslat(2)
if(m1.ne.m2)return

norm=1.d0
do i=1,2
 norm=norm*norm_sto(gamma(i),nslat(i))
enddo

n12=nslat(1)+nslat(2)
g12=gamma(1)+gamma(2)

overlap_slater=norm*fact(n12)/g12**(n12+1)

end

!! Computation of nuclear integral I= (n1l1m1 |-1/r| n2l2m2)= - N_n1*N_n2 !(n1+n2-1)!(g1+g2)**(n1+n2)   l1=l1 m1=m2
!!
double precision function nuclear_slater(nslat,lslat,mslat,gamma)
implicit none
integer nslat(4),lslat(4),mslat(4),i
integer l1,l2,m1,m2,n12
double precision gamma(4),int,norm,norm_sto,g12,fact

nuclear_slater=0.d0

l1=lslat(1)
l2=lslat(2)
if(l1.ne.l2)return
m1=mslat(1)
m2=mslat(2)
if(m1.ne.m2)return

norm=1.d0
do i=1,2
 norm=norm*norm_sto(gamma(i),nslat(i))
enddo

n12=nslat(1)+nslat(2)
g12=gamma(1)+gamma(2)

nuclear_slater=-norm*fact(n12-1)/g12**n12

end

!! Computation of kinetic integral I=-1/2 (n1 l1 m1 | Laplacian | n2 l2 m2)
!!
!! I= -1/2  N_n1*N_n2 [ g2**2 *(n1 l1 m1 | Lapl| n2 l2 m2)/N_n1*N_n2 
!!                     -2*n2*g2*(n1 l1 m1|Lap| n2-1 l2 m2)/N_n1*N_(n2-1)
!!                     +2*(n2+l2)*(n2-l2-1)*(n1 l1 m1|Lap| n2-2 l2 m2)/N_n1*N_(n2-2)
!!
!! Be careful, error in Berlu's formula
!!
double precision function kinetic_slaterII(nslat,lslat,mslat,gamma)
implicit none
integer nslat(4),lslat(4),mslat(4),i,nslatemp(4)
integer l1,l2,m1,m2,n2,ntot
double precision gamma(4),int,t0,t1,t2,g2,overlap_slater,g1,g12
double precision norm,norm_sto,fact,normtemp

kinetic_slaterII=0.d0

l1=lslat(1)
l2=lslat(2)
if(l1.ne.l2)return
m1=mslat(1)
m2=mslat(2)
if(m1.ne.m2)return
g1=gamma(2)
g2=gamma(2)
g12=g1+g2
n2=nslat(2)

norm=1.d0
do i=1,2
 norm=norm*norm_sto(gamma(i),nslat(i))
enddo

normtemp=norm_sto(gamma(1),nslat(1))*norm_sto(gamma(2),nslat(2))
t0=overlap_slater(nslat,lslat,mslat,gamma)/normtemp

nslatemp(1)=nslat(1)
nslatemp(2)=nslat(2)-1
t1=0.d0
if(nslatemp(2).ge.1)then
normtemp=norm_sto(gamma(2),nslat(2))*norm_sto(gamma(2),nslatemp(2))
t1=overlap_slater(nslatemp,lslat,mslat,gamma)/normtemp
else
ntot=nslatemp(1)+nslatemp(2)
t1=fact(ntot)/g12**(ntot+1)
endif

nslatemp(1)=nslat(1)
nslatemp(2)=nslat(2)-2
t2=0.d0
if(nslatemp(2).ge.1)then
normtemp=norm_sto(gamma(2),nslat(2))*norm_sto(gamma(2),nslatemp(2))
t2=overlap_slater(nslatemp,lslat,mslat,gamma)/normtemp
else
ntot=nslatemp(1)+nslatemp(2)-1
t2=fact(ntot)/g12**(ntot+1)
endif

kinetic_slaterII=-0.5d0*norm*(g2**2*t0- 2.d0*n2*g2*t1+(n2+l2)*(n2-l2-1)*t2)

end

double precision function ijkl_slater_cartesian(nslat,lslat,mslat,gamma)
implicit none
integer nslat(4),lslat(4),mslat(4)
double precision gamma(4),ijkl_slater,rac2
integer sig(4),epsk1,epsk2,epsk3,epsk4,eps(4),i,mslat_k(4)
double complex int,coef(-1:1,-1:1,-1:1),imag
parameter(rac2=1.41421356237310d0)

do i=1,4
 sig(i)=(-1)**mslat(i)
 eps(i)=mslat(i)*sig(i)
enddo

coef(-1,-1,-1)=(0.d0, 1.d0)
coef(-1,-1,+1)=(0.d0, 1.d0)

coef(-1, 0,-1)=(0.d0, 0.d0)
coef(-1, 0,+1)=(0.d0, 0.d0)

coef(-1,+1,-1)=(0.d0, 1.d0)
coef(-1,+1,+1)=(0.d0,-1.d0)

coef( 0,-1,-1)=(0.d0, 0.d0)
coef( 0,-1,+1)=(0.d0, 0.d0)

coef( 0, 0,-1)=(rac2, 0.d0)
coef( 0, 0,+1)=(rac2, 0.d0)

coef( 0,+1,-1)=(0.d0, 0.d0)
coef( 0,+1,+1)=(0.d0, 0.d0)

coef(+1,-1,-1)=(1.d0, 0.d0)
coef(+1,-1,+1)=(1.d0, 0.d0)

coef(+1, 0,-1)=(0.d0, 0.d0)
coef(+1, 0,+1)=(0.d0, 0.d0)

coef(+1,+1,-1)=(0.d0,-1.d0)
coef(+1,+1,+1)=(0.d0, 1.d0)

int=(0.d0,0.d0)

do epsk1=-1,1
 do epsk2=-1,1
  do epsk3=-1,1
   do epsk4=-1,1

    mslat_k(1)=epsk1*iabs(mslat(1))
    mslat_k(2)=epsk2*iabs(mslat(2))
    mslat_k(3)=epsk3*iabs(mslat(3))
    mslat_k(4)=epsk4*iabs(mslat(4))

    int=int +conjg(coef(eps(1),sig(1),epsk1))*coef(eps(2),sig(2),epsk2) &
    *conjg(coef(eps(3),sig(3),epsk3))*coef(eps(4),sig(4),epsk4)* &
    ijkl_slater(nslat,lslat,mslat_k,gamma)
   enddo
  enddo
 enddo
enddo

print*,'int=',int

ijkl_slater_cartesian=int
end


double precision function kinetic_slater(nslat,lslat,mslat,gamma)
implicit none
integer nslat(4),lslat(4),mslat(4),i,nslatemp(4)
integer l1,l2,m1,m2,n2,n1,n12
double precision gamma(4),int,norm,norm_sto,g2,overlap_slater,g1,g12,fact
double precision t1,t2,t3,t4

kinetic_slater=0.d0

l1=lslat(1)
l2=lslat(2)
if(l1.ne.l2)return
m1=mslat(1)
m2=mslat(2)
if(m1.ne.m2)return
g1=gamma(1)
g2=gamma(2)
g12=g1+g2
n1=nslat(1)
n2=nslat(2)
n12=n1+n2

norm=1.d0
do i=1,2
 norm=norm*norm_sto(gamma(i),nslat(i))
enddo

t1=n2*(n2-1)*fact(n12-2)/g12**(n12-1)
t2=2.d0*n2*g2*fact(n12-1)/g12**n12
t3=g2**2*fact(n12)/g12**(n12+1)
t4=l2*(l2+1)*fact(n12-2)/g12**(n12-1)

kinetic_slater=-0.5d0*norm*(t1-t2+t3-t4)

end

double precision function norm_sto(gam,n)
implicit none
double precision gam,g,fact
integer n,p
p=2*n
g=2.d0*gam
norm_sto=dsqrt( g**(p+1)/fact(p))
end

subroutine compute_lmin_lmax(l1,l2,l3,l4,m1,m2,m3,m4,lmin,lmax)
implicit none
integer l1,l2,l3,l4,m1,m2,m3,m4,lmin,lmax
integer mx12,l12min,mx34,l34min

lmax=min(l1+l2,l3+l4)

mx12=max(iabs(l1-l2),iabs(m1-m2))
if(mod(lmax+mx12,2).eq.0)then
l12min=mx12
else
l12min=mx12+1
endif

mx34=max(iabs(l3-l4),iabs(m3-m4))
if(mod(lmax+mx34,2).eq.0)then
l34min=mx34
else
l34min=mx34+1
endif

lmin=max(l12min,l34min)

end

! F is defined as:

! F(k1,g1,k2,g2,l)= int_0^inf r1^k1 exp(-g1 r1) 
!                                             [ 1/r1**(l+1) int_0^r1 r2**(k2+l) exp(-g2 r2)
!                                             +   r1**l     int_r1^inf  r2**(k2-l-1) exp(-g2 r2) ]
!
! k1 integer ge 0 
! k2 integer ge 0 
! k1 and k2 ge (l+1)
!
! We can show:
!
!  F= (k1-l-1)! (k2+l)!/g1**(k1-l)/g2**(k2+l+1)  + \sum_{i=0}^{k2-l-1} (k2-l-1)!(k1+k2-1-i)!/(k2-l-1-i)!/g2**(i+1)/(g1+g2)**(k1+k2-i)
!                                                 - sum_{i=0}^{k2+l  } (k2+l  )!(k1+k2-1-i)!/(k2+l-i  )!/g2**(i+1)/(g1+g2)**(k1+k2-i)
!
! This formula can been checked with bigF_num
!
! Here F must be computed with  k1=n12  k2=n34 g1 -->g12 et g2 --->g13
!
!
double precision function bigF(n12,g12,n34,g34,l)
implicit none
integer n12,n34,l,l12,l34
double precision g12,g34,g,first_term,smallg,fact

l12=-l
l34=l+1
g=g12+g34

first_term= fact(n34+l34-1)/g34**(n34+l34)*fact(n12+l12-1)/g12**(n12+l12)

bigF=first_term+ smallg(n34,n12,g34,g,l12-1)-smallg(n34,n12,g34,g,l34-1) 
end

double precision function bigF_num(np,rmax,k1,g1,k2,g2,l)
implicit none
integer k1,k2,l,np,i1,i2
double precision g1,g2,rmax,dr1,dr2,r1min,accu1,accu2,accu,r1,r2,ijkl_only_s_num

dr1=rmax/np
r1min=-dr1/2.d0

write(*,*)'l=',l
accu=0.d0
do i1=1,np
 r1=r1min+(i1-1)*dr1
 accu1=0.d0
 dr2=r1/np
 do i2=1,np
  r2=(i2-1)*dr2
  accu1=accu1+r2**(k2+l)*dexp(-g2*r2)/r1**(l+1)
 enddo
 accu1=accu1*dr2
 dr2=(rmax-r1)/np
 accu2=0.d0
 do i2=1,np
  r2=r1+(i2-1)*dr2
  accu2=accu2+r2**(k2-l-1)*dexp(-g2*r2)*r1**l
 enddo
 accu2=accu2*dr2
 accu=accu+r1**k1*dexp(-g1*r1)*(accu1+accu2)
enddo

bigF_num=accu*dr1
end


double precision function smallg(n,m,g,gp,l)
implicit none
integer n,m,l,i
double precision g,gp,fact
smallg=0.d0
do i=0,n+l
 smallg=smallg+fact(n+l)*fact(n+m-1-i)/fact(n+l-i)/g**(i+1)/gp**(n+m-i)
enddo
end


!  Y_l^m(theta,phi) = i^(m+|m|) ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2  P_l^|m|(cos(theta))  exp(i m phi)
!  l=0,1,2,....
!  m=0,1,...,l
! Here:
!  n=l (n=0,1,...)  
!  m=0,1,...,n 
!  x=cos(theta) 0 < x < 1
!
!  
!  This routine computes:   PM(m,n) for n=0,...,N (number N in input) and m=0,..,n

!   Exemples (see 'Associated Legendre Polynomilas wikipedia')
!    P_{0}^{0}(x)=1
!    P_{1}^{-1}(x)=-1/2 P_{1}^{1}(x)
!    P_{1}^{0}(x)=x
!    P_{1}^{1}(x)=-(1-x^2)^{1/2}
!    P_{2}^{-2}(x)=1/24 P_{2}^{2}(x)
!    P_{2}^{-1}(x)=-1/6 P_{2}^{1}(x)
!    P_{2}^{0}(x)=1/2 (3x^{2}-1)
!    P_{2}^{1}(x)=-3x(1-x^2)^{1/2}
!    P_{2}^{2}(x)=3(1-x^2)


!  Y_l^m(theta,phi) = i^(m+|m|) ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2
!  P_l^|m|(cos(theta))  exp(i m phi)


        double precision function ylm_no_phi(l,m,x)
        implicit double precision (a-h,o-z)
        DIMENSION PM(0:100,0:100)
        MM=100
        pi=dacos(-1.d0)
        iabs_m=iabs(m)
!        if(iabs_m.gt.l)stop 'm must be between -l and l'
        if(m.ge.0)then
         phase=1.d0
        else
         phase=(-1.d0)**m
        endif
        factor= dsqrt( ((2*l+1)*fact(l-iabs_m))/(4.d0*pi*fact(l+iabs_m)) )

        if(dabs(x).gt.1.d0)then
         print*,'pb. in ylm_no'
         print*,'x=',x
         stop 
        endif

        call LPMN(MM,l,l,X,PM)
        plm=PM(iabs_m,l)
        ylm_no_phi=phase*factor*plm
        end


! Definition of Gaunt coefficients:
!
!  G(l1m1,l2l2,lm) =  int_dOmega  Ylm^*  Yl1m1  Yl2m2
!
!                  =  int_{-1}^{+1} dx Ylm_no_phi(x) Yl1m1_no_phi(x) Yl2m2_no_phi(x)   with m=m1+m2, if not G=0
!    x=cos(theta)

        double precision function gaunt_num(l1,m1,l2,m2,l,m)
        implicit double precision (a-h,o-z)
        pi=dacos(-1.d0)
        do knp=5,5
         np=10**knp
         dx=2.d0/np 
         accu=0.d0
         x=-1.d0
         do i=1,np
          x=x+dx
          accu=accu+ ylm_no_phi(l,m,x)* ylm_no_phi(l1,m1,x)* ylm_no_phi(l2,m2,x)
         enddo
         gaunt_num=2.d0*pi*accu*dx
        enddo
        end

!      This routine computes    I= int_{0}^x du u**n  exp(-alpha*u)
! 
!       I= n!/alpha**(n+1) - x**n/alpha exp(-alpha*x) \sum_{i=0}^n  n!/(n-i)!  (1/(alpha*x)^i

        double precision function expo_uncomplete(x,n,alpha)
        implicit none
        double precision x,alpha,accu,fact
        integer n,i
        accu=0.d0
        do i=0,n
          accu=accu+fact(n)/fact(n-i)/(alpha*x)**i
        enddo
        expo_uncomplete= fact(n)/alpha**(n+1)-x**n/alpha*dexp(-alpha*x)*accu
        end

        double precision function expo_uncomplete_num(x,n,alpha)
        implicit none
        double precision x,alpha,accu,fact,du,u
        integer n,i,np,k
        np=10**6
        du=x/np
        u=0.d0
        accu=0.d0
        do i=1,np
          u=u+du
          accu=accu+u**n*dexp(-alpha*u)
        enddo
        expo_uncomplete_num= accu*du
        end
!  expo_complete= int_0^inf  r**n exp(-gamma*r) = fact(n)/gamma**(n+1)
!
        double precision function expo_complete(n,gamma)
        implicit none
        double precision gamma,fact
        integer n
        expo_complete= fact(n)/gamma**(n+1)
        end
!  ijkl_only_s = (n1n2|n3n4)  where all STO are s type
!
        double precision function ijkl_only_s(nslat,gslat)
        implicit none
        integer nslat(4)
        double precision gslat(4)
        double precision g1,g2,g3,g4,g12,g34,norm,t1,t2,t3,norm_sto,fact,expo_complete
        integer n1,n2,n3,n4,n12,n34,i

        n1=nslat(1)
        n2=nslat(2)
        n3=nslat(3)
        n4=nslat(4)
        g1=gslat(1)
        g2=gslat(2)
        g3=gslat(3)
        g4=gslat(4)

        n34=n3+n4
        n12=n1+n2
        g12=g1+g2
        g34=g3+g4

        norm=norm_sto(g1,n1)*norm_sto(g2,n2)*norm_sto(g3,n3)*norm_sto(g4,n4)

        t1= fact(n34)/g34**(n34+1)*expo_complete(n12-1,g12)
        t2=0.d0
        do i=0,n34
         t2=t2+ 1.d0/fact(n34-i)/g34**(i+1)*expo_complete(n12+n34-1-i,g12+g34)
        enddo
        t2=t2*fact(n34)
        t3=0.d0
        do i=0,n34-1
         t3=t3+ 1.d0/fact(n34-i-1)/g34**(i+1)*expo_complete(n12+n34-1-i,g12+g34)
        enddo
        t3=t3*fact(n34-1)

         ijkl_only_s=norm*(t1-t2+t3)

!        print*,'norm=',norm,t1,t2,t3
!        print*,ijkl_only_s
        end

      double precision function repart(x,g)
      implicit none
      double precision x,g,y
      y=g*x
      repart=1.d0+y+0.5d0*y**2
      repart=repart*dexp(-y)
      repart=1.d0-repart
      end

      double precision function x_from_repart(y,g)
      implicit none
      double precision x0,x1,y,g,diff,xmid,repartnew,repartold,u
      double precision repart
      integer iter
      x0=0.d0
      x1=0.d0
      do while (repart(x1,g).lt.y) 
       x1=x1+1.d0
      enddo
      diff=1.d0
      repartold=repart(x1,g)
      iter=0
      do while (diff.gt.1.d-12) 
       xmid=(x0+x1)/2.d0
       u=repart(xmid,g)
       if(u.lt.y)then
        x0=xmid
       else
        x1=xmid
       endif
       diff=dabs(repartold-u)
       repartold=u
       iter=iter+1
      enddo
      x_from_repart=x1

!      print*,'x0=',x0
!      print*,'x1=',x1
!      print*,'y=',y
!      print*,'repart0=',repart(x0,g)
!      print*,'repart1=',repart(x1,g)
!      print*,'iter=',iter
      end 

      double precision function draw_r(g)
      implicit none
      double precision g,y,ranf,x_from_repart
      y=ranf()
      draw_r=x_from_repart(y,g)
      end

!    IMPORTANT:    call as follows:  gaunt(l,l1,l2,m,m1,m2) = /int d_Omega Y*_l^m  Y_l1^m1  Y_l2^m2
!
!      Y_l^m = (-1)**[(m+|m|)/2] C_lm  P_l^|m| ( cos[theta] ) exp(i m phi)
!
!    C_lm =sqrt(  (2*l+1)/(4pi) (l-|m|)!/(l+|m|)!  )
!
!     m=m1+m2
!
!     Be careful in the routine  l1 --> m  l2---> l1  l3 --> l2
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: gaunt
! !INTERFACE:
Real (8) Function gaunt (l1, l2, l3, m1, m2, m3)
! !INPUT/OUTPUT PARAMETERS:
! l1, l2, l3 : angular momentum quantum numbers (in,integer)
! m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
! Returns the Gaunt coefficient given by
! $$ \langle Y^{l_1}_{m_1}|Y^{l_2}_{m_2}|Y^{l_3}_{m_3} \rangle
! = (-1)^{m_1}\left[\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi} \right]
! ^{\frac{1}{2}}
! \begin{pmatrix} l_1 & l_2 & l_3 \\ 0 & 0 & 0 \end{pmatrix}
! \begin{pmatrix} l_1 & l_2 & l_3 \\ -m_1 & m_2 & m_3 \end{pmatrix}. $$
! Suitable for $l_i$ less than 50.
!
! !REVISION HISTORY:
! Created November 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer :: l1
      Integer :: l2
      Integer :: l3
      Integer :: m1
      Integer :: m2
      Integer :: m3
! local variables
      Integer :: j, j1, j2, j3, jh
      Real (8) :: t1
! real constant 1/sqrt(4*pi)
      Real (8), Parameter :: c1 = 0.28209479177387814347d0
! external functions
      Real (8) :: wigner3j, factr, factnm
      External wigner3j, factr, factnm
      If ((l1 .Lt. 0) .Or. (l2 .Lt. 0) .Or. (l3 .Lt. 0) .Or. (Abs(m1) &
     & .Gt. l1) .Or. (Abs(m2) .Gt. l2) .Or. (Abs(m3) .Gt. l3)) Then
Write (*,*)
         Write (*, '("Error(gaunt): non-physical arguments :")')
         Write (*, '("l1 = ", I8, " l2 = ", I8, " l3 = ", I8)') l1, l2, &
        & l3
         Write (*, '("m1 = ", I8, " m2 = ", I8, " m3 = ", I8)') m1, m2, &
        & m3
         Write (*,*)
         Stop
End If
If ((l1 .Gt. 50) .Or. (l2 .Gt. 50) .Or. (l3 .Gt. 50)) Then
Write (*,*)
         Write (*, '("Error(gaunt): angular momenta out of range : ", 3&
&I8)') l1, l2, l3
         Write (*,*)
         Stop
End If
If (m1-m2-m3 .Ne. 0) Then
gaunt = 0.d0
         Return
End If
j1 = l2 - l1 + l3
      j2 = l1 - l2 + l3
      j3 = l1 + l2 - l3
      If ((j1 .Lt. 0) .Or. (j2 .Lt. 0) .Or. (j3 .Lt. 0)) Then
gaunt = 0.d0
         Return
End If
j = l1 + l2 + l3
      If (Mod(j, 2) .Ne. 0) Then
gaunt = 0.d0
         Return
End If
jh = j / 2
      t1 = Sqrt (dble((2*l1+1)*(2*l2+1)*(2*l3+1))*factr(j1, &
     & j+1)*factnm(j2, 1)*factnm(j3, 1))
      t1 = t1 * factr (jh, jh-l1) / (factnm(jh-l2, 1)*factnm(jh-l3, 1))
      gaunt = t1 * c1 * wigner3j (l1, l2, l3,-m1, m2, m3)
      If (Mod(m1+jh, 2) .Ne. 0) gaunt = - gaunt
      Return
End Function
!EOC

!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: factr
! !INTERFACE:
Real (8) Function factr (n, d)
! !INPUT/OUTPUT PARAMETERS:
! n : numerator (in,integer)
! d : denominator (in,integer)
! !DESCRIPTION:
! Returns the ratio $n!/d!$ for $n,d\ge 0$. Performs no under- or overflow
! checking.
!
! !REVISION HISTORY:
! Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Integer, Intent (In) :: d
! local variables
      Integer :: i
      If (n .Lt. 0) Then
Write (*,*)
         Write (*, '("Error(factr): n < 0 : ", I8)') n
         Write (*,*)
         Stop
End If
If (d .Lt. 0) Then
Write (*,*)
         Write (*, '("Error(factr): d < 0 : ", I8)') d
         Write (*,*)
         Stop
End If
If (n .Lt. d) Then
factr = 1.d0 / dble (n+1)
         Do i = n + 2, d
            factr = factr / dble (i)
         End Do
Else If (n .Eq. d) Then
factr = 1.d0
      Else
factr = dble (d+1)
         Do i = d + 2, n
            factr = factr * dble (i)
         End Do
End If
Return
End Function
!EOC

!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: factnm
! !INTERFACE:
Real (8) Function factnm (n, m)
! !INPUT/OUTPUT PARAMETERS:
! n : input (in,integer)
! m : order of multifactorial (in,integer)
! !DESCRIPTION:
! Returns the multifactorial
! $$ n\underbrace{!!\,...\,!}_{m\,{\rm times}}=\prod_{i\ge 0,\,n-im>0}n-im $$
! for $n,\,m \ge 0$. $n$ should be less than 150.
!
! !REVISION HISTORY:
! Created January 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Integer, Intent (In) :: m
! local variables
      Integer :: i, j
      Real (8) :: f1 (24), f2 (38)
      Data f1 / 1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0, 5040.d0, &
     & 40320.d0, 362880.d0, 3628800.d0, 39916800.d0, 479001600.d0, &
     & 6227020800.d0, 87178291200.d0, 1307674368000.d0, &
     & 20922789888000.d0, 355687428096000.d0, 6402373705728000.d0, &
     & 121645100408832000.d0, 2432902008176640000.d0, &
     & 51090942171709440000.d0, 1124000727777607680000.d0, &
     & 25852016738884976640000.d0, 620448401733239439360000.d0 /
      Data f2 / 1.d0, 2.d0, 3.d0, 8.d0, 15.d0, 48.d0, 105.d0, 384.d0, &
     & 945.d0, 3840.d0, 10395.d0, 46080.d0, 135135.d0, 645120.d0, &
     & 2027025.d0, 10321920.d0, 34459425.d0, 185794560.d0, &
     & 654729075.d0, 3715891200.d0, 13749310575.d0, 81749606400.d0, &
     & 316234143225.d0, 1961990553600.d0, 7905853580625.d0, &
     & 51011754393600.d0, 213458046676875.d0, 1428329123020800.d0, &
     & 6190283353629375.d0, 42849873690624000.d0, &
     & 191898783962510625.d0, 1371195958099968000.d0, &
     & 6332659870762850625.d0, 46620662575398912000.d0, &
     & 221643095476699771875.d0, 1678343852714360832000.d0, &
     & 8200794532637891559375.d0, 63777066403145711616000.d0 /
! fast return if possible
      If (n .Eq. 0) Then
factnm = 1.d0
         Return
End If
If (m .Eq. 1) Then
If ((n .Ge. 1) .And. (n .Le. 24)) Then
factnm = f1 (n)
            Return
End If
End If
If (m .Eq. 2) Then
If ((n .Ge. 1) .And. (n .Le. 38)) Then
factnm = f2 (n)
            Return
End If
End If
If (n .Lt. 0) Then
Write (*,*)
         Write (*, '("Error(factnm): n < 0 : ", I8)') n
         Write (*,*)
         Stop
End If
If (m .Le. 0) Then
Write (*,*)
         Write (*, '("Error(factnm): m <= 0 : ", I8)') m
         Write (*,*)
         Stop
End If
If (n .Gt. 150) Then
Write (*,*)
         Write (*, '("Error(factnm): n out of range : ", I8)') n
         Write (*,*)
         Stop
End If
If (m .Eq. 1) Then
factnm = f1 (24)
         Do i = 25, n
            factnm = factnm * dble (i)
         End Do
Else
j = n / m
         If (Mod(n, m) .Eq. 0) j = j - 1
         factnm = dble (n)
         Do i = 1, j
            factnm = factnm * dble (n-i*m)
         End Do
End If
Return
End Function
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: wigner3j
! !INTERFACE:
Real (8) Function wigner3j (j1, j2, j3, m1, m2, m3)
! !INPUT/OUTPUT PARAMETERS:
! j1, j2, j3 : angular momentum quantum numbers (in,integer)
! m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
! Returns the Wigner $3j$-symbol. There are many equivalent definitions for
! the $3j$-symbols, the following provides high accuracy for $j\le 50$
! \begin{align*}
! &\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix}=(-1)^
! {j1+j2+m3}\\
! &\times\sqrt{\frac{(j_1+m_1)!(j_2+m_2)!(j_3+m_3)!(j_3-m_3)!(j_1-m_1)!
! (j_2-m_2)!}{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!(1+j_1+j_2+j_3)!}}
! \times\sum_{\max(0,j_2-j_3-m_1,j_1-j_3+m_2)}^
! {\min(j_1+j_2-j_3,j_1-m_1,j_2+m_2)}\\
! &(-1)^k\frac{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!}
! {(j_3-j_1-m_2+k)!(j_3-j_2+m_1+k)!(j_1+j_2-j_3-k)!k!(j_1-m_1-k)!
! (j_2+m_2-k)}.
! \end{align*}
!
! !REVISION HISTORY:
! Created November 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: j1
      Integer, Intent (In) :: j2
      Integer, Intent (In) :: j3
      Integer, Intent (In) :: m1
      Integer, Intent (In) :: m2
      Integer, Intent (In) :: m3
! local variables
      Integer :: k, k1, k2, l1, l2, l3, n1, n2
      Real (8) :: sgn, sum, t1
! external functions
      Real (8) :: factnm, factr
      External factnm, factr
! check input variables
      If ((j1 .Lt. 0) .Or. (j2 .Lt. 0) .Or. (j3 .Lt. 0) .Or. (Abs(m1) &
     & .Gt. j1) .Or. (Abs(m2) .Gt. j2) .Or. (Abs(m3) .Gt. j3)) Then
Write (*,*)
         Write (*, '("Error(wigner3j): invalid arguments :")')
         Write (*, '("j1 = ", I8, " j2 = ", I8, " j3 = ", I8)') j1, j2, &
        & j3
         Write (*, '("m1 = ", I8, " m2 = ", I8, " m3 = ", I8)') m1, m2, &
        & m3
         Write (*,*)
         Stop
End If
If ((j1 .Eq. 0) .And. (j2 .Eq. 0) .And. (j3 .Eq. 0)) Then
wigner3j = 1.d0
         Return
End If
If ((j1 .Gt. 50) .Or. (j2 .Gt. 50) .Or. (j3 .Gt. 50)) Then
Write (*,*)
         Write (*, '("Error(wigner3j): angular momenta out of range : "&
&, 3I8)') j1, j2, j3
         Write (*,*)
         Stop
End If
l1 = j2 - j1 + j3
      l2 = j1 - j2 + j3
      l3 = j1 + j2 - j3
      If ((m1+m2+m3 .Ne. 0) .Or. (l1 .Lt. 0) .Or. (l2 .Lt. 0) .Or. (l3 &
     & .Lt. 0)) Then
wigner3j = 0.d0
         Return
End If
n1 = j1 - m1
      n2 = j2 + m2
      k1 = Max (0, j2-j3-m1, j1-j3+m2)
      k2 = Min (l3, n1, n2)
      If (Mod(k1+j1+j2+m3, 2) .Ne. 0) Then
sgn = - 1.d0
      Else
sgn = 1.d0
      End If
sum = 0.d0
      Do k = k1, k2
         t1 = sgn * factr (l1, l1-n2+k) * factr (l2, l2-n1+k) * factr &
        & (l3, l3-k)
         sum = sum + t1 / (factnm(k, 1)*factnm(n1-k, 1)*factnm(n2-k, &
        & 1))
         sgn = - sgn
      End Do
t1 = factr (j1+m1, l1) * factr (j2+m2, l2) * factr (j3+m3, l3)
      t1 = t1 * factr (j3-m3, 1+j1+j2+j3) * factnm (j1-m1, 1) * factnm &
     & (j2-m2, 1)
      wigner3j = sum * Sqrt (t1)
      Return
End Function
!EOC

double precision function ijkl_slater_xnynzn_num(npts,nix,niy,niz,power,gamma)
implicit none

integer nix(4),niy(4),niz(4),npts
integer power(4)

double precision gamma(4),pi,dphi1,dphi2,du1,du2,phi1,phi2,u1,u2,r1,r2,dr
double precision rmax,v1,v2,x1,y1,z1,x2,y2,z2,r12,f,rint,gtot
integer i1phi1,i2phi2,iu1,ir1,iu2,ir2

pi=dacos(-1.d0)

gtot=(gamma(1)+gamma(2)+gamma(3)+gamma(4))/2.d0

rmax=-1.d0/gtot*dlog(1.d-10)
print*,'rmax?'
read(*,*)rmax

dphi1=2.d0*pi/npts
du1=2.d0/npts
dphi2=2.d0*pi/npts
du2=2.d0/npts
dr=rmax/npts

do i1phi1=1,npts
 phi1=(i1phi1-1)*dphi1+ 0.33*dphi1
 do iu1=1,npts
  u1=-1.d0+(iu1-1)*du1+ 0.33*du1
  do ir1=1,npts
   r1=(ir1-1)*dr+ 0.33*dr

   v1=dsqrt(1.d0-u1**2)
   x1=r1*dcos(phi1)*v1
   y1=r1*dsin(phi1)*v1
   z1=r1*u1

   do i2phi2=1,npts
    phi2=(i2phi2-1)*dphi2+ 0.77*dphi2
    do iu2=1,npts
     u2=-1.d0+(iu2-1)*du2+ 0.77*du2
     do ir2=1,npts
      r2=(ir2-1)*dr+ 0.77*dr

      v2=dsqrt(1.d0-u2**2)
      x2=r2*dcos(phi2)*v2
      y2=r2*dsin(phi2)*v2
      z2=r2*u2

      r12=dsqrt( (x1-x2)**2+(y1-y2)**2+(z1-z2)**2 )

      f= x1**(nix(1)+nix(2))*y1**(niy(1)+niy(2))*z1**(niz(1)+niz(2)) * &
         x2**(nix(3)+nix(4))*y2**(niy(3)+niy(4))*z2**(niz(3)+niz(4))
      f=f*r1**(power(1)+power(2))*r2**(power(3)+power(4))
      f=f*r1**2*r2**2*dexp(-(gamma(1)+gamma(2))*r1-(gamma(3)+gamma(4))*r2)
      f=f/r12
      rint=rint+f
     enddo
    enddo
   enddo
  enddo
 enddo
enddo

ijkl_slater_xnynzn_num=rint*dr**2*dphi1*dphi2*du1*du2

end

double precision function ijkl_slater_xnynzn(nix,niy,niz,power,gamma)
implicit none

integer nix(4),niy(4),niz(4)
integer power(4)
double precision gamma(4)

integer max_ni_max
parameter(max_ni_max=2)

integer nslat(4),lslat(4),mslat(4)
integer lslat_max(4)
double precision sij,tij,vij,vijkl
double precision overlap_slater,kinetic_slater,nuclear_slater,ijkl_slater

integer nx,ny,nz,l,m,max_ni,k,ntot
integer m1,m2,m3,m4,l1,l2,l3,l4
double complex imag
double complex c(0:max_ni_max,0:max_ni_max,0:max_ni_max,0:3*max_ni_max,-3*max_ni_max:3*max_ni_max)
double complex rint,coef
double precision norm,norm_sto
double precision pi,rac4pi,rac4pi_3,racpi_15,c0

!sij=overlap_slater(nslat,lslat,mslat,gamma)
!tij=kinetic_slater(nslat,lslat,mslat,gamma)
!vij=nuclear_slater(nslat,lslat,mslat,gamma)
!write(*,*)sij,tij,vij
!vijkl=ijkl_slater(nslat,lslat,mslat,gamma)/sij**2
!write(*,*)vijkl

max_ni=2

do k=1,4
 ntot=nix(k)+niy(k)+niz(k)
 if(ntot.gt.max_ni)stop ' x^nx x^nx x^nx Slater no yet implemented'
enddo

pi=dacos(-1.d0)
rac4pi=dsqrt(4.d0*pi)
rac4pi_3=dsqrt(4.d0*pi/3.d0)
racpi_15=dsqrt(pi/15.d0)

imag=(0.d0,1.d0)

do nx=0,max_ni
 do ny=0,max_ni
  do nz=0,max_ni
   do l=0,nx+ny+nz
    do m=-l,l
     c(nx,ny,nz,l,m)=0.d0
    enddo
   enddo
  enddo
 enddo
enddo

!! Decomposition of x^nx y^ny z^nz over spherical harmonics (complex)
!!
!!******************************************************************
!! x^nx y^ny z^nz = \sum_{l=0,lmax} {m=-l,l}  c(nx,ny,nz,l,m) Y(l,m)
!! lmax=nx+ny+nz
!!*****************************************************************

! 1
c(0,0,0,0,0) = rac4pi

!x
c(1,0,0,1,-1)=  rac4pi_3*1.d0/dsqrt(2.d0)
c(1,0,0,1, 1)= -rac4pi_3*1.d0/dsqrt(2.d0)

!y
c(0,1,0,1,-1)=  imag*rac4pi_3*1.d0/dsqrt(2.d0)
c(0,1,0,1, 1)=  imag*rac4pi_3*1.d0/dsqrt(2.d0)

!z
c(0,0,1,1, 0)=  rac4pi_3

!xy
c(1,1,0,2,-2)=  imag*racpi_15*dsqrt(2.d0)
c(1,1,0,2, 2)= -imag*racpi_15*dsqrt(2.d0)

!yz
c(0,1,1,2,-1)=  imag*racpi_15*dsqrt(2.d0) 
c(0,1,1,2, 1)=  imag*racpi_15*dsqrt(2.d0)

!xz
c(1,0,1,2,-1)=  racpi_15*dsqrt(2.d0)
c(1,0,1,2, 1)= -racpi_15*dsqrt(2.d0)

!!x2
c(2,0,0,0, 0)=  rac4pi/3.d0
c(2,0,0,2,-2)=  racpi_15*dsqrt(2.d0)
c(2,0,0,2, 2)=  racpi_15*dsqrt(2.d0)
c(2,0,0,2, 0)= -racpi_15*2.d0/dsqrt(3.d0)

!!y2
c(0,2,0,0, 0)=  rac4pi/3.d0
c(0,2,0,2,-2)= -racpi_15*dsqrt(2.d0)
c(0,2,0,2, 2)= -racpi_15*dsqrt(2.d0)
c(0,2,0,2, 0)= -racpi_15*2.d0/dsqrt(3.d0)

!!z2
c(0,0,2,0, 0)=  rac4pi/3.d0
c(0,0,2,2, 0)=  racpi_15*4.d0/dsqrt(3.d0)

do k=1,4
 nslat(k)=1+nix(k)+niy(k)+niz(k)+power(k)
 lslat_max(k)=nix(k)+niy(k)+niz(k)
enddo

rint=(0.d0,0.d0)

do l1=0,lslat_max(1)
do m1=-l1,l1
 do l2=0,lslat_max(2)
 do m2=-l2,l2
  do l3=0,lslat_max(3)
  do m3=-l3,l3
   do l4=0,lslat_max(4)
   do m4=-l4,l4

    lslat(1)=l1
    lslat(2)=l2
    lslat(3)=l3
    lslat(4)=l4

    mslat(1)=-m1
    mslat(2)=m2
    mslat(3)=-m3
    mslat(4)=m4

    coef=c(nix(1),niy(1),niz(1),l1,m1)*c(nix(2),niy(2),niz(2),l2,m2) * &
    c(nix(3),niy(3),niz(3),l3,m3)*c(nix(4),niy(4),niz(4),l4,m4)

    if(abs(coef).gt.1.d-20)then
     rint=rint + coef*(-1.d0)**(m1+m3)*ijkl_slater(nslat,lslat,mslat,gamma)
    endif

   enddo
   enddo
  enddo
  enddo
 enddo
 enddo
enddo
enddo

norm=1.d0
do k=1,4
 norm=norm*norm_sto(gamma(k),nslat(k))
enddo

rint=rint/norm

ijkl_slater_xnynzn=dreal(rint)

end

subroutine ijkl_slater_mc(np,nslat,lslat,mslat,gslat,rmoy,error)
implicit none
integer n1,n2,n3,n4,l1,l2,l3,l4,m1,m2,m3,m4
integer i1,i2,i3,i4,i5,i6,kk
integer nslat(4),lslat(4),mslat(4),np
double precision gslat(4),g1,g2,g3,g4,arg,arg1,arg2,g12,g34
double precision r1,r2,phi1,phi2,x1,x2,norm,norm_sto,ylm_no_phi
double precision r12,xp1,xp2,yp1,yp2,zp1,zp2,accu,term1,term2
double precision draw_r,ranf,rmoy,error,res(20),pi

pi=dacos(-1.d0)

rmoy=0.d0
error=0.d0

n1=nslat(1)
n2=nslat(2)
n3=nslat(3)
n4=nslat(4)

l1=lslat(1)
l2=lslat(2)
l3=lslat(3)
l4=lslat(4)
if(mod(l1+l2+l3+l4,2).eq.1)return

m1=mslat(1)
m2=mslat(2)
m3=mslat(3)
m4=mslat(4)
if((m2-m1).ne.(m3-m4))then
 return
endif

g1=gslat(1)
g2=gslat(2)
g12=g1+g2
g3=gslat(3)
g4=gslat(4)
g34=g3+g4

norm=norm_sto(g1,n1)*norm_sto(g2,n2)*norm_sto(g3,n3)*norm_sto(g4,n4)

do kk=1,20
accu=0.d0

do i1=1,np

r1=draw_r(g12)
r2=draw_r(g34)
phi1=2.d0*pi*ranf()
phi2=2.d0*pi*ranf()
x1= 2.d0*(ranf()-0.5d0)
x2= 2.d0*(ranf()-0.5d0)

term1=r1**(n1+n2-2)* ylm_no_phi(l1,m1,x1)*ylm_no_phi(l2,m2,x1)
term2=r2**(n3+n4-2)* ylm_no_phi(l3,m3,x2)*ylm_no_phi(l4,m4,x2)

arg1=1.d0-x1**2
if(arg1.lt.0.d0)then
write(*,*)'arg1=',arg1
stop
endif

xp1=r1*dsqrt(dabs(1.d0-x1**2))*dcos(phi1)
yp1=r1*dsqrt(dabs(1.d0-x1**2))*dsin(phi1)
zp1=r1*x1

arg2=1.d0-x2**2
if(arg2.lt.0.d0)then
write(*,*)'arg2=',arg2
stop
endif

xp2=r2*dsqrt(dabs(1.d0-x2**2))*dcos(phi2)
yp2=r2*dsqrt(dabs(1.d0-x2**2))*dsin(phi2)
zp2=r2*x2

arg= (xp1-xp2)**2 + (yp1-yp2)**2+ (zp1-zp2)**2

if(arg.lt.0.d0)then
print*,'arg=',arg
print*,xp1,yp1,zp1,xp2,yp2,zp2
stop
endif

r12=dsqrt(arg)

if(r12.eq.0.d0)then
print*,'r12=',r12
stop
endif

accu=accu+term1*term2/r12
enddo

res(kk)=16.d0*pi**2*norm*accu*4.d0/g12**3/g34**3/np

enddo

call erreur(res,20,rmoy,error)

end

subroutine overlap_slater_mc(np,nslat,lslat,mslat,gslat,rmoy,error)
implicit none
integer n1,n2,l1,l2,m1,m2
integer i1,kk
integer nslat(4),lslat(4),mslat(4),np
double precision gslat(4),g1,g2,arg,g12
double precision r,phi,x,ylm_no_phi
double precision accu,term,norm,norm_sto
double precision draw_r,ranf,rmoy,error,res(20),pi

pi=dacos(-1.d0)

l1=lslat(1)
l2=lslat(2)
if(l1.ne.l2)then
print*,'l1 ne l2'
return
endif
m1=mslat(1)
m2=mslat(2)
if(m1.ne.m2)then
print*,'m1 ne m2'
return
endif

pi=dacos(-1.d0)

rmoy=0.d0
error=0.d0

n1=nslat(1)
n2=nslat(2)
g1=gslat(1)
g2=gslat(2)
g12=g1+g2

norm=norm_sto(g1,n1)*norm_sto(g2,n2)

do kk=1,20
accu=0.d0

do i1=1,np

r=draw_r(g12)
x= 2.d0*(ranf()-0.5d0)

term=r**(n1+n2-2)* ylm_no_phi(l1,m1,x)*ylm_no_phi(l2,m2,x)

arg=1.d0-x**2
if(arg.lt.0.d0)then
write(*,*)'arg=',arg
stop
endif

accu=accu+term

enddo

res(kk)=norm*4.d0*pi*accu*2.d0/g12**3/np

enddo

call erreur(res,20,rmoy,error)

end

subroutine nuclear_slater_mc(np,nslat,lslat,mslat,gslat,rmoy,error)
implicit none
integer n1,n2,l1,l2,m1,m2
integer i1,kk
integer nslat(4),lslat(4),mslat(4),np
double precision gslat(4),g1,g2,arg,g12
double precision r,phi,x,ylm_no_phi
double precision accu,term,norm,norm_sto
double precision draw_r,ranf,rmoy,error,res(20),pi

pi=dacos(-1.d0)

l1=lslat(1)
l2=lslat(2)
if(l1.ne.l2)return
m1=mslat(1)
m2=mslat(2)
if(m1.ne.m2)return

pi=dacos(-1.d0)

rmoy=0.d0
error=0.d0

n1=nslat(1)
n2=nslat(2)
g1=gslat(1)
g2=gslat(2)
g12=g1+g2

norm=norm_sto(g1,n1)*norm_sto(g2,n2)

do kk=1,20
accu=0.d0

do i1=1,np

r=draw_r(g12)
x= 2.d0*(ranf()-0.5d0)

term=r**(n1+n2-3)* ylm_no_phi(l1,m1,x)*ylm_no_phi(l2,m2,x)

arg=1.d0-x**2
if(arg.lt.0.d0)then
write(*,*)'arg=',arg
stop
endif

accu=accu+term

enddo

res(kk)=-norm*4.d0*pi*accu*2.d0/g12**3/np

enddo

call erreur(res,20,rmoy,error)

end

double precision function ijkl_only_s_num(np,rmax,nslat,gslat)
implicit none
integer n1,n2,n3,n4
integer i1,i2
integer nslat(4),np,pi
double precision gslat(4),g1,g2,g3,g4
double precision r1,r2,dr1,dr2,norm_sto
double precision rmax,accu,term1,term2,r12max
double precision r1min,r2min,norm

pi=dacos(-1.d0)

n1=nslat(1)
n2=nslat(2)
n3=nslat(3)
n4=nslat(4)

g1=gslat(1)
g2=gslat(2)
g3=gslat(3)
g4=gslat(4)

dr1=rmax/np
dr2=rmax/np

r1min=-dr2/3.d0
r2min=-dr2/2.d0

norm=norm_sto(g1,n1)*norm_sto(g2,n2)*norm_sto(g3,n3)*norm_sto(g4,n4)

accu=0.d0

do i1=1,np
do i2=1,np
r1=r1min+(i1-1)*dr1
r2=r2min+(i2-1)*dr2

term1=r1**(n1+n2)* dexp(-(g1+g2)*r1)
term2=r2**(n3+n4)* dexp(-(g3+g4)*r2)

r12max=max(r1,r2)

accu=accu+term1*term2/r12max

enddo
enddo

ijkl_only_s_num=norm*accu*dr1*dr2

end

subroutine read_basis(filename_basis)
include 'j.inc'
character*(128) :: filename_basis
character*80 ATOM_READ
character*80 orb

!! read basis set for the 36 first atoms
!! n_b(i) number of basis (S,P,D,F, or G) centered on nucleus i=1,36
!! orb_b(k,i) name of the basis function, that is S,P,D,F or G for k=1,n_b(i)
!! n_cont_b(k,i) number of contracted primitives for the radial part
!! for the moment STO have only one primtive
!! u_STO(r)=r**n_sto* exp(-gi r)
!! u_GTO(r) = sum_i c_i exp(-gi r**2)
!! g_i and c_i  gamma_b(k,m,i),coef_b(k,m,i) m=1,n_cont_b(k,i)


if(basis_type.eq.'STO')then
 open(1,file=filename_basis)
 do i=1,n_atoms_max
  n_b(i)=0
 enddo
 i=0
 do
  read(1,*,end=1000)ATOM_READ,orb,gam
  i = i + 1
  i_atom=number_atom(ATOM_READ)
  n_b(i_atom)=n_b(i_atom)+1
  if(n_b(i_atom).gt.nbasis_max)stop 'in read_basis: increase nbasis_max'
  gamma_b(n_b(i_atom),1,i_atom)=gam
  orb_b(n_b(i_atom),i_atom)=orb
  n_cont_b(n_b(i_atom),i_atom)=1
  coef_b(n_b(i_atom),1,i_atom)=1.d0
 enddo
 1000  continue
endif

if(basis_type.eq.'GTO')then
 open(1,file=filename_basis)

 do i=1,36
  n_b(i)=0
  i_error=0
  read(1,'(a80)',iostat=i_error)ATOM_READ
  do while (i_error.eq.0)
   read(1,'(a1,i5)',iostat=i_error)ORB,n_contr
   if(n_contr.eq.0)then
    i_error=-1
   endif
   if(i_error.eq.0)then
     n_b(i)=n_b(i)+1
     if(n_b(i).gt.nbasis_max)stop 'in read_basis: increase nbasis_max'
     orb_b(n_b(i),i)=orb
     n_cont_b(n_b(i),i)=n_contr
     do m=1,n_contr
      read(1,*,iostat=i_error)ibid,gamma_b(n_b(i),m,i),coef_b(n_b(i),m,i)
     enddo
   endif
  enddo
 enddo !i=1,36
endif

end
!! Build the gaussian expansion of each atomic orbitals 
!!*****************************************************
!!
!! For i=1,nbasis
!!
!!  phi_i(r) centered at [ center(1,i),center(2,i),center(3,i) ]
!!  phi_i(r) = x^npower(1,i) y^npower(2,i) z^npower(3,i)  u_i(r)
!!
!!  u_i(r)= sum_{m=1,n_gauss(i)} [ c_gauss(1,m,i) + c_gauss(2,m,i) r^2 ] exp[-g_gauss(m,i) r^2]
!!
!!  i_type(i)=1  :     c_gauss(2,...)=0 and only first part is computed
!!  i_type(i)=2  :     c_gauss(1...)=0 and only second part is computed
!!  i_type(i)=3  :     c_gauss(1,...) and c_gauss(2,...) computed
!!
!!  For STO g_min(i)=g_slater 
!!  For GTO = smallest exponent of the expansion (to be used for computing G_center in Monte Carlo)

subroutine build_gaussian_expansion_of_orbitals(mpi_rank)
include 'j.inc'
character*80,orb

open(21,file='info_basis')
open(22,file='orbital_coefficients_gaussian_expansion')
open(23,file='info_simplified')

if(basis_type.eq.'STO')then

 if(finite_range)then
  call read_fit_exp_rc
 else
  call read_fit_SMILES
 endif

 nbasis=0

 do i=1,number_nuclei

  i_atom=number_atom(ATOM(i))

  do k=1,n_b(i_atom)

   orb=orb_b(k,i_atom)

   if(orb.eq.'1S')then
    nbasis=nbasis+1
    if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
    orb_name(nbasis)='1S'
    orb_name_full(nbasis)='1S'
    n_STO(nbasis)=0 
    i_type(nbasis)=1 
    do l=1,3
     npower(l,nbasis)=0
     center(l,nbasis)=centers_nuclei(l,i)
    enddo
    nucleus_number(nbasis)=i

    ng(1,nbasis)=ng0
    if(gamma_b(k,1,i_atom).gt.g_thr)then
     ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
    endif
    n_contract(nbasis)=n_cont_b(k,i_atom)
    do mm=1,n_contract(nbasis)
     g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
     c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
    enddo

    if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
    g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
    g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

    if(finite_range)then
     call compute_rc_of_orbital(nbasis)
    endif
    call build_c_g_gauss_STO(nbasis)
    if(mpi_rank.eq.0)call write_STO_in_file_info_basis(nbasis)

   endif  !1S

   if(orb.eq.'2S')then

    nbasis=nbasis+1
    if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
    orb_name(nbasis)='2S'
    orb_name_full(nbasis)='2S'
    n_STO(nbasis)=1
    i_type(nbasis)=2

    do l=1,3
     npower(l,nbasis)=0
     center(l,nbasis)=centers_nuclei(l,i)
    enddo
    nucleus_number(nbasis)=i

    ng(1,nbasis)=ng0
    if(gamma_b(k,1,i_atom).gt.g_thr)then
     ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
    endif
    n_contract(nbasis)=n_cont_b(k,i_atom)
    do mm=1,n_contract(nbasis)
     g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
     c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
    enddo

    if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
    g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
    g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

    if(finite_range)then
     call compute_rc_of_orbital(nbasis)
    endif
    call build_c_g_gauss_STO(nbasis)
    if(mpi_rank.eq.0)call write_STO_in_file_info_basis(nbasis)

   endif !2S

   if(orb.eq.'3S')then

    nbasis=nbasis+1
    if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
    orb_name(nbasis)='3S'
    orb_name_full(nbasis)='3S'
    n_STO(nbasis)=2
    i_type(nbasis)=2

    do l=1,3
     npower(l,nbasis)=0
     center(l,nbasis)=centers_nuclei(l,i)
    enddo
    nucleus_number(nbasis)=i

    ng(1,nbasis)=ng0
    if(gamma_b(k,1,i_atom).gt.g_thr)then
     ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
    endif
    n_contract(nbasis)=n_cont_b(k,i_atom)
    do mm=1,n_contract(nbasis)
     g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
     c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
    enddo

    if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
    g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
    g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

    if(finite_range)then
     call compute_rc_of_orbital(nbasis)
    endif
    call build_c_g_gauss_STO(nbasis)
    if(mpi_rank.eq.0)call write_STO_in_file_info_basis(nbasis)

   endif !3S

   if(orb.eq.'2P')then

    do m=1,3
     nbasis=nbasis+1
     if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
     n_STO(nbasis)=0
     i_type(nbasis)=1
     orb_name(nbasis)='2P'
     do l=1,3
      center(l,nbasis)=centers_nuclei(l,i)
      if(l.eq.m)then
       npower(l,nbasis)=1
      else
       npower(l,nbasis)=0
      endif 
      nucleus_number(nbasis)=i
     enddo
     if(m.eq.1)orb_name_full(nbasis)='2P_X'
     if(m.eq.2)orb_name_full(nbasis)='2P_Y'
     if(m.eq.3)orb_name_full(nbasis)='2P_Z'

     ng(1,nbasis)=ng0+i_add_2p
     if(gamma_b(k,1,i_atom).gt.g_thr)then
      ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
     endif
     n_contract(nbasis)=n_cont_b(k,i_atom)
     do mm=1,n_contract(nbasis)
      g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
      c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
     enddo

     if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
     g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
     g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

     if(finite_range)then
      call compute_rc_of_orbital(nbasis)
     endif
     call build_c_g_gauss_STO(nbasis)
     if(mpi_rank.eq.0)call write_STO_in_file_info_basis(nbasis)

    enddo ! m

   endif !2P

   if(orb.eq.'3P')then

    do m=1,3
     nbasis=nbasis+1
     if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
     if(m.eq.1)orb_name_full(nbasis)='3P_X'
     if(m.eq.2)orb_name_full(nbasis)='3P_Y'
     if(m.eq.3)orb_name_full(nbasis)='3P_Z'
     n_STO(nbasis)=1
     i_type(nbasis)=2
     orb_name(nbasis)='3P'
     do l=1,3
      center(l,nbasis)=centers_nuclei(l,i)
      if(l.eq.m)then
       npower(l,nbasis)=1
      else
       npower(l,nbasis)=0
      endif
     enddo
     nucleus_number(nbasis)=i

     ng(1,nbasis)=ng0+i_add_2p
     if(gamma_b(k,1,i_atom).gt.g_thr)then
      ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
     endif
     n_contract(nbasis)=n_cont_b(k,i_atom)
     do mm=1,n_contract(nbasis)
      g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
      c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
     enddo

     if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
     g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
     g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

     if(finite_range)then
      call compute_rc_of_orbital(nbasis)
     endif
     call build_c_g_gauss_STO(nbasis)
     if(mpi_rank.eq.0)call write_STO_in_file_info_basis(nbasis)

    enddo ! m

   endif !3P

   if(orb.eq.'3D')then

    do m=1,6
     nbasis=nbasis+1
     if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
     n_STO(nbasis)=0
     i_type(nbasis)=1
     orb_name(nbasis)='3D'
     do l=1,3
      center(l,nbasis)=centers_nuclei(l,i)
     enddo
     if(m.eq.1)then
      npower(1,nbasis)=2
      npower(2,nbasis)=0
      npower(3,nbasis)=0
      orb_name_full(nbasis)='3D_XX'
     endif
     if(m.eq.2)then
      npower(1,nbasis)=1
      npower(2,nbasis)=1
      npower(3,nbasis)=0
      orb_name_full(nbasis)='3D_XY'
     endif
     if(m.eq.3)then
      npower(1,nbasis)=1
      npower(2,nbasis)=0
      npower(3,nbasis)=1
      orb_name_full(nbasis)='3D_XZ'
     endif
     if(m.eq.4)then
      npower(1,nbasis)=0
      npower(2,nbasis)=2
      npower(3,nbasis)=0
      orb_name_full(nbasis)='3D_YY'
     endif
     if(m.eq.5)then
      npower(1,nbasis)=0
      npower(2,nbasis)=1
      npower(3,nbasis)=1
      orb_name_full(nbasis)='3D_YZ'
     endif
     if(m.eq.6)then
      npower(1,nbasis)=0
      npower(2,nbasis)=0
      npower(3,nbasis)=2
      orb_name_full(nbasis)='3D_ZZ'
     endif
     nucleus_number(nbasis)=i

     ng(1,nbasis)=ng0+i_add_3d
     if(gamma_b(k,1,i_atom).gt.g_thr)then
      ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
     endif
     n_contract(nbasis)=n_cont_b(k,i_atom)
     do mm=1,n_contract(nbasis)
      g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
      c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
     enddo

     if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
     g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
     g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

     if(finite_range)then
      call compute_rc_of_orbital(nbasis)
     endif
     call build_c_g_gauss_STO(nbasis)
     if(mpi_rank.eq.0)call write_STO_in_file_info_basis(nbasis)

    enddo ! m 

   endif  !3D

   if(orb.eq.'4F')then

    do m=1,10
     nbasis=nbasis+1
     if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
     n_STO(nbasis)=0
     i_type(nbasis)=1
     orb_name(nbasis)='4F'
     do l=1,3
      center(l,nbasis)=centers_nuclei(l,i)
     enddo
     if(m.eq.1)then
      npower(1,nbasis)=3
      npower(2,nbasis)=0
      npower(3,nbasis)=0
      orb_name_full(nbasis)='4F_XXX'
     endif
     if(m.eq.2)then
      npower(1,nbasis)=2
      npower(2,nbasis)=1
      npower(3,nbasis)=0
      orb_name_full(nbasis)='4F_XXY'
     endif
     if(m.eq.3)then
      npower(1,nbasis)=2
      npower(2,nbasis)=0
      npower(3,nbasis)=1
      orb_name_full(nbasis)='4F_XXZ'
     endif
     if(m.eq.4)then
      npower(1,nbasis)=1
      npower(2,nbasis)=2
      npower(3,nbasis)=0
      orb_name_full(nbasis)='4F_XYY'
     endif
     if(m.eq.5)then
      npower(1,nbasis)=1
      npower(2,nbasis)=1
      npower(3,nbasis)=1
      orb_name_full(nbasis)='4F_XYZ'
     endif
     if(m.eq.6)then
      npower(1,nbasis)=1
      npower(2,nbasis)=0
      npower(3,nbasis)=2
      orb_name_full(nbasis)='4F_XZZ'
     endif
     if(m.eq.7)then
      npower(1,nbasis)=0
      npower(2,nbasis)=3
      npower(3,nbasis)=0
      orb_name_full(nbasis)='4F_YYY'
     endif
     if(m.eq.8)then
      npower(1,nbasis)=0
      npower(2,nbasis)=2
      npower(3,nbasis)=1
      orb_name_full(nbasis)='4F_YYZ'
     endif
     if(m.eq.9)then
      npower(1,nbasis)=0
      npower(2,nbasis)=1
      npower(3,nbasis)=2
      orb_name_full(nbasis)='4F_YZZ'
     endif
     if(m.eq.10)then
      npower(1,nbasis)=0
      npower(2,nbasis)=0
      npower(3,nbasis)=3
      orb_name_full(nbasis)='4F_ZZZ'
     endif
     nucleus_number(nbasis)=i

     ng(1,nbasis)=ng0+i_add_3d
     if(gamma_b(k,1,i_atom).gt.g_thr)then
      ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
     endif
     n_contract(nbasis)=n_cont_b(k,i_atom)
     do mm=1,n_contract(nbasis)
      g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
      c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
     enddo

     if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
     g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
     g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

     if(finite_range)then
      call compute_rc_of_orbital(nbasis)
     endif
     call build_c_g_gauss_STO(nbasis)
     if(mpi_rank.eq.0)call write_STO_in_file_info_basis(nbasis)

    enddo ! m

   endif  !4F

   if(orb.eq.'5G')then

    do m=1,15
     nbasis=nbasis+1
     if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
     n_STO(nbasis)=0
     i_type(nbasis)=1
     orb_name(nbasis)='5G'
     do l=1,3
      center(l,nbasis)=centers_nuclei(l,i)
     enddo
     if(m.eq.1)then
      npower(1,nbasis)=4
      npower(2,nbasis)=0
      npower(3,nbasis)=0
      orb_name_full(nbasis)='5G_XXXX'
     endif
     if(m.eq.2)then
      npower(1,nbasis)=3
      npower(2,nbasis)=1
      npower(3,nbasis)=0
      orb_name_full(nbasis)='5G_XXXY'
     endif
     if(m.eq.3)then
      npower(1,nbasis)=2
      npower(2,nbasis)=2
      npower(3,nbasis)=0
      orb_name_full(nbasis)='5G_XXYY'
     endif
     if(m.eq.4)then
      npower(1,nbasis)=1
      npower(2,nbasis)=3
      npower(3,nbasis)=0
      orb_name_full(nbasis)='5G_XYYY'
     endif
     if(m.eq.5)then
      npower(1,nbasis)=0
      npower(2,nbasis)=4
      npower(3,nbasis)=0
      orb_name_full(nbasis)='5G_YYYY'
     endif
     if(m.eq.6)then
      npower(1,nbasis)=3
      npower(2,nbasis)=0
      npower(3,nbasis)=1
      orb_name_full(nbasis)='5G_XXXZ'
     endif
     if(m.eq.7)then
      npower(1,nbasis)=2
      npower(2,nbasis)=1
      npower(3,nbasis)=1
      orb_name_full(nbasis)='5G_XXYZ'
     endif
     if(m.eq.8)then
      npower(1,nbasis)=1
      npower(2,nbasis)=2
      npower(3,nbasis)=1
      orb_name_full(nbasis)='5G_XYYZ'
     endif
     if(m.eq.9)then
      npower(1,nbasis)=0
      npower(2,nbasis)=3
      npower(3,nbasis)=1
      orb_name_full(nbasis)='5G_YYYZ'
     endif
     if(m.eq.10)then
      npower(1,nbasis)=2
      npower(2,nbasis)=0
      npower(3,nbasis)=2
      orb_name_full(nbasis)='5G_XXZZ'
     endif
     if(m.eq.11)then
      npower(1,nbasis)=1
      npower(2,nbasis)=1
      npower(3,nbasis)=2
      orb_name_full(nbasis)='5G_XYZZ'
     endif
     if(m.eq.12)then
      npower(1,nbasis)=0
      npower(2,nbasis)=2
      npower(3,nbasis)=2
      orb_name_full(nbasis)='5G_YYZZ'
     endif
     if(m.eq.13)then
      npower(1,nbasis)=1
      npower(2,nbasis)=0
      npower(3,nbasis)=3
      orb_name_full(nbasis)='5G_XZZZ'
     endif
     if(m.eq.14)then
      npower(1,nbasis)=0
      npower(2,nbasis)=1
      npower(3,nbasis)=3
      orb_name_full(nbasis)='5G_YZZZ'
     endif
     if(m.eq.15)then
      npower(1,nbasis)=0
      npower(2,nbasis)=0
      npower(3,nbasis)=4
      orb_name_full(nbasis)='5G_ZZZZ'
     endif
     nucleus_number(nbasis)=i

     ng(1,nbasis)=ng0+i_add_3d
     if(gamma_b(k,1,i_atom).gt.g_thr)then
      ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
     endif
     n_contract(nbasis)=n_cont_b(k,i_atom)
     do mm=1,n_contract(nbasis)
      g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
      c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
     enddo

     if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
     g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
     g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

     if(finite_range)then
      call compute_rc_of_orbital(nbasis)
     endif
     call build_c_g_gauss_STO(nbasis)
     if(mpi_rank.eq.0)call write_STO_in_file_info_basis(nbasis)

    enddo !m

   endif  !5G

  enddo ! k=1,n_b(i_atom)
 enddo ! i=1,number_nuclei

 return

endif

if(basis_type.eq.'GTO')then

 if(finite_range)call read_fit_gaussian_rc

 nbasis=0

 do i=1,number_nuclei
  i_atom=number_atom(ATOM(i))
  do k=1,n_b(i_atom)

   orb=orb_b(k,i_atom)

   if(orb.ne.'S'.and.orb.ne.'P'.and.orb.ne.'D'.and.orb.ne.'F'.and.orb.ne.'G')then
    write(*,*)'WARNING: orb=',orb,' not yet coded!'
   endif

   if(orb.eq.'S')then
    nbasis=nbasis+1
    if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
    orb_name(nbasis)='S'
    orb_name_full(nbasis)='S'
    i_type(nbasis)=1
    do l=1,3
     npower(l,nbasis)=0
     center(l,nbasis)=centers_nuclei(l,i)
    enddo
    nucleus_number(nbasis)=i

    n_contract(nbasis)=n_cont_b(k,i_atom)
    nx=npower(1,nbasis)
    ny=npower(2,nbasis)
    nz=npower(3,nbasis)
    do mm=1,n_contract(nbasis)
     g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
     c_contract(mm,nbasis)=coef_b(k,mm,i_atom)/rnorm_prim(nx,ny,nz,g_contract(mm,nbasis))
    enddo
    g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)

    if(finite_range)then
     call compute_rc_of_orbital(nbasis)
     call compute_n_star_of_primitive(nbasis)
    endif
    call build_c_g_gauss_GTO(nbasis)
    if(mpi_rank.eq.0)call write_GTO_in_file_info_basis(nbasis)
   endif  !S

   if(orb_b(k,i_atom).eq.'P')then

    do m=1,3
     nbasis=nbasis+1
     if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
     orb_name(nbasis)='P'
     if(m.eq.1)orb_name_full(nbasis)='X'
     if(m.eq.2)orb_name_full(nbasis)='Y'
     if(m.eq.3)orb_name_full(nbasis)='Z'
     i_type(nbasis)=1
     do l=1,3
      center(l,nbasis)=centers_nuclei(l,i)
      if(l.eq.m)then
       npower(l,nbasis)=1
      else
       npower(l,nbasis)=0
      endif
     enddo
     nucleus_number(nbasis)=i

     n_contract(nbasis)=n_cont_b(k,i_atom)
     nx=npower(1,nbasis)
     ny=npower(2,nbasis)
     nz=npower(3,nbasis)
     do mm=1,n_contract(nbasis)
      g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
      c_contract(mm,nbasis)=coef_b(k,mm,i_atom)/rnorm_prim(nx,ny,nz,g_contract(mm,nbasis))
     enddo
     g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)

     if(finite_range)then
      call compute_rc_of_orbital(nbasis)
      call compute_n_star_of_primitive(nbasis)
     endif
     call build_c_g_gauss_GTO(nbasis)
     if(mpi_rank.eq.0)call write_GTO_in_file_info_basis(nbasis)

    enddo ! m
   endif

   if(orb.eq.'D')then
    do m=1,6
     nbasis=nbasis+1
     if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
     do l=1,3
      center(l,nbasis)=centers_nuclei(l,i)
     enddo
     orb_name(nbasis)='D'
     if(m.eq.1)then
      npower(1,nbasis)=2
      npower(2,nbasis)=0
      npower(3,nbasis)=0
      orb_name_full(nbasis)='XX'
     endif
     if(m.eq.2)then
      npower(1,nbasis)=1
      npower(2,nbasis)=1
      npower(3,nbasis)=0
      orb_name_full(nbasis)='XY'
     endif
     if(m.eq.3)then
      npower(1,nbasis)=1
      npower(2,nbasis)=0
      npower(3,nbasis)=1
      orb_name_full(nbasis)='XZ'
     endif
     if(m.eq.4)then
      npower(1,nbasis)=0
      npower(2,nbasis)=2
      npower(3,nbasis)=0
      orb_name_full(nbasis)='YY'
     endif
     if(m.eq.5)then
      npower(1,nbasis)=0
      npower(2,nbasis)=1
      npower(3,nbasis)=1
      orb_name_full(nbasis)='YZ'
     endif
     if(m.eq.6)then
      npower(1,nbasis)=0
      npower(2,nbasis)=0
      npower(3,nbasis)=2
      orb_name_full(nbasis)='ZZ'
     endif
     i_type(nbasis)=1
     nucleus_number(nbasis)=i

     n_contract(nbasis)=n_cont_b(k,i_atom)
     nx=npower(1,nbasis)
     ny=npower(2,nbasis)
     nz=npower(3,nbasis)
     do mm=1,n_contract(nbasis)
      g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
      c_contract(mm,nbasis)=coef_b(k,mm,i_atom)/rnorm_prim(nx,ny,nz,g_contract(mm,nbasis))
     enddo
     g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)

     if(finite_range)then
      call compute_rc_of_orbital(nbasis)
      call compute_n_star_of_primitive(nbasis)
     endif
     call build_c_g_gauss_GTO(nbasis)
     if(mpi_rank.eq.0)call write_GTO_in_file_info_basis(nbasis)

    enddo ! m
   endif  !D

   if(orb.eq.'F')then
    do m=1,10
     nbasis=nbasis+1
     if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
     do l=1,3
      center(l,nbasis)=centers_nuclei(l,i)
     enddo
     orb_name(nbasis)='F'
     if(m.eq.1)then
      npower(1,nbasis)=3
      npower(2,nbasis)=0
      npower(3,nbasis)=0
      orb_name_full(nbasis)='XXX'
     endif
     if(m.eq.2)then
      npower(1,nbasis)=2
      npower(2,nbasis)=1
      npower(3,nbasis)=0
      orb_name_full(nbasis)='XXY'
     endif
     if(m.eq.3)then
      npower(1,nbasis)=2
      npower(2,nbasis)=0
      npower(3,nbasis)=1
      orb_name_full(nbasis)='XXZ'
     endif
     if(m.eq.4)then
      npower(1,nbasis)=1
      npower(2,nbasis)=2
      npower(3,nbasis)=0
      orb_name_full(nbasis)='XYY'
     endif
     if(m.eq.5)then
      npower(1,nbasis)=1
      npower(2,nbasis)=1
      npower(3,nbasis)=1
      orb_name_full(nbasis)='XYZ'
     endif
     if(m.eq.6)then
      npower(1,nbasis)=1
      npower(2,nbasis)=0
      npower(3,nbasis)=2
      orb_name_full(nbasis)='XZZ'
     endif
     if(m.eq.7)then
      npower(1,nbasis)=0
      npower(2,nbasis)=3
      npower(3,nbasis)=0
      orb_name_full(nbasis)='YYY'
     endif
     if(m.eq.8)then
      npower(1,nbasis)=0
      npower(2,nbasis)=2
      npower(3,nbasis)=1
      orb_name_full(nbasis)='YYZ'
     endif
     if(m.eq.9)then
      npower(1,nbasis)=0
      npower(2,nbasis)=1
      npower(3,nbasis)=2
      orb_name_full(nbasis)='YZZ'
     endif
     if(m.eq.10)then
      npower(1,nbasis)=0
      npower(2,nbasis)=0
      npower(3,nbasis)=3
      orb_name_full(nbasis)='ZZZ'
     endif
     i_type(nbasis)=1
     nucleus_number(nbasis)=i

     n_contract(nbasis)=n_cont_b(k,i_atom)
     nx=npower(1,nbasis)
     ny=npower(2,nbasis)
     nz=npower(3,nbasis)
     do mm=1,n_contract(nbasis)
      g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
      c_contract(mm,nbasis)=coef_b(k,mm,i_atom)/rnorm_prim(nx,ny,nz,g_contract(mm,nbasis))
     enddo
     g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)

     if(finite_range)then
      call compute_rc_of_orbital(nbasis)
      call compute_n_star_of_primitive(nbasis)
     endif
     call build_c_g_gauss_GTO(nbasis)
     if(mpi_rank.eq.0)call write_GTO_in_file_info_basis(nbasis)

    enddo ! m

   endif  !F

   if(orb.eq.'G')then
    do m=1,15
     nbasis=nbasis+1
     if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
     do l=1,3
      center(l,nbasis)=centers_nuclei(l,i)
     enddo
     orb_name(nbasis)='G'
     if(m.eq.1)then
      npower(1,nbasis)=4
      npower(2,nbasis)=0
      npower(3,nbasis)=0
      orb_name_full(nbasis)='XXXX'
     endif
     if(m.eq.2)then
      npower(1,nbasis)=3
      npower(2,nbasis)=1
      npower(3,nbasis)=0
      orb_name_full(nbasis)='XXXY'
     endif
     if(m.eq.3)then
      npower(1,nbasis)=2
      npower(2,nbasis)=2
      npower(3,nbasis)=0
      orb_name_full(nbasis)='XXYY'
     endif
     if(m.eq.4)then
      npower(1,nbasis)=1
      npower(2,nbasis)=3
      npower(3,nbasis)=0
      orb_name_full(nbasis)='XYYY'
     endif
     if(m.eq.5)then
      npower(1,nbasis)=0
      npower(2,nbasis)=4
      npower(3,nbasis)=0
      orb_name_full(nbasis)='YYYY'
     endif
     if(m.eq.6)then
      npower(1,nbasis)=3
      npower(2,nbasis)=0
      npower(3,nbasis)=1
      orb_name_full(nbasis)='XXXZ'
     endif
     if(m.eq.7)then
      npower(1,nbasis)=2
      npower(2,nbasis)=1
      npower(3,nbasis)=1
      orb_name_full(nbasis)='XXYZ'
     endif
     if(m.eq.8)then
      npower(1,nbasis)=1
      npower(2,nbasis)=2
      npower(3,nbasis)=1
      orb_name_full(nbasis)='XYYZ'
     endif
     if(m.eq.9)then
      npower(1,nbasis)=0
      npower(2,nbasis)=3
      npower(3,nbasis)=1
      orb_name_full(nbasis)='YYYZ'
     endif
     if(m.eq.10)then
      npower(1,nbasis)=2
      npower(2,nbasis)=0
      npower(3,nbasis)=2
      orb_name_full(nbasis)='XXZZ'
     endif
     if(m.eq.11)then
      npower(1,nbasis)=1
      npower(2,nbasis)=1
      npower(3,nbasis)=2
      orb_name_full(nbasis)='XYZZ'
     endif
     if(m.eq.12)then
      npower(1,nbasis)=0
      npower(2,nbasis)=2
      npower(3,nbasis)=2
      orb_name_full(nbasis)='YYZZ'
     endif
     if(m.eq.13)then
      npower(1,nbasis)=1
      npower(2,nbasis)=0
      npower(3,nbasis)=3
      orb_name_full(nbasis)='XZZZ'
     endif
     if(m.eq.14)then
      npower(1,nbasis)=0
      npower(2,nbasis)=1
      npower(3,nbasis)=3
      orb_name_full(nbasis)='YZZZ'
     endif
     if(m.eq.15)then
      npower(1,nbasis)=0
      npower(2,nbasis)=0
      npower(3,nbasis)=4
      orb_name_full(nbasis)='ZZZZ'
     endif
     i_type(nbasis)=1
     nucleus_number(nbasis)=i

     n_contract(nbasis)=n_cont_b(k,i_atom)
     nx=npower(1,nbasis)
     ny=npower(2,nbasis)
     nz=npower(3,nbasis)
     do mm=1,n_contract(nbasis)
      g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
      c_contract(mm,nbasis)=coef_b(k,mm,i_atom)/rnorm_prim(nx,ny,nz,g_contract(mm,nbasis))
     enddo
     g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)

     if(finite_range)then
      call compute_rc_of_orbital(nbasis)
      call compute_n_star_of_primitive(nbasis)
     endif
     call build_c_g_gauss_GTO(nbasis)
     if(mpi_rank.eq.0)call write_GTO_in_file_info_basis(nbasis)

    enddo !m
   endif  !G

  enddo !k loop over all basis of atom
 enddo !number of nuclei

 if(nbasis.gt.nbasis_max)stop 'nbasis_max too small'

endif  ! basis_type.eq.'GTO'

do i=1,nbasis
 do j=1,nbasis
  dist_ij(i,j)=  &
  dsqrt((center(1,i)-center(1,j))**2+(center(2,i)-center(2,j))**2+(center(3,i)-center(3,j))**2)
 enddo
enddo

end

!! ROUTINE GAUSS_IJKL
!********************

! A_x=center(1,i)
! A_y=center(2,i)
! A_z=center(3,i)

! nA_x=npower(1,i)
! nA_y=npower(2,i)
! nA_z=npower(3,i)

! idem for B,C,D for orbitals i,k,j,l
!
! The following quantity is computed:

! gauss_ijkl_new = int d1 d2 (x1-A_x)**nA_x * (y1-A_y)**nA_y * (z1-A_z)**nA_z *  u_i(r1)
!                            (x1-B_x)**nB_x * (y1-B_y)**nB_y * (z1-B_z)**nB_z *  u_k(r1)     
! 1/r_12
!                            (x2-C_x)**nC_x * (y2-C_y)**nC_y * (z2-C_z)**nC_z * u_j(r2)
!                            (x2-D_x)**nD_x * (y2-D_y)**nD_y * (z2-D_z)**nD_z * u_l(r2)
!
! with the radial parts given by
!
!  u_i(r1) = sum_{i1=1,ngauss(1)}  [ c0(i1,1)+ c2(i1,1)*(r1-A)**2] exp[-gg(i1,1)*(r1-A)**2]
!  u_k(r1) = sum_{i2=1,ngauss(2)}  [ c0(i2,1)+ c2(i2,1)*(r1-A)**2] exp[-gg(i2,1)*(r1-B)**2]
!
!  u_j(r2) = sum_{i3=1,ngauss(3)}  [ c0(i3,1)+ c2(i3,1)*(r2-C)**2] exp[-gg(i3,1)*(r2-C)**2]
!  u_l(r2) = sum_{i4=1,ngauss(4)}  [ c0(i4,1)+ c2(i4,1)*(r2-D)**2] exp[-gg(i4,1)*(r2-D)**2]

! To save computational cost:
!
! it=1 the r**2-component is not calculated (c2=0)
! it=2 only the r**2-component is calculated (c1=0)
! it=3 both r**0 and r**2-component are calculated 

double precision function gauss_ijkl(i,k,j,l)
include 'j.inc'

!! local arrays:
dimension n_orb(4),nc(4),d(n_gauss_max,4,4),n_c(3,4,4)

n_orb(1)=i
n_orb(2)=k
n_orb(3)=j
n_orb(4)=l

do kk=1,4

 i_o=n_orb(kk)

 if(i_type(i_o).eq.1)then
  nc(kk)=1
  do m=1,n_gauss(i_o)
   do ii=1,nc(kk)
    d(m,ii,kk)=c_gauss(1,m,i_o)
   enddo
  enddo
  do ll=1,3
   do ii=1,nc(kk)
    n_c(ll,ii,kk)=npower(ll,i_o)
   enddo
  enddo
 endif

 if(i_type(i_o).eq.2)then
  nc(kk)=3
  do m=1,n_gauss(i_o)
   do ii=1,nc(kk)
   d(m,ii,kk)=c_gauss(2,m,i_o)
   enddo
  enddo
  do ll=1,3
   do ii=1,nc(kk)
    if(ll.eq.ii)then
     n_c(ll,ii,kk)=npower(ll,i_o)+2
    else
     n_c(ll,ii,kk)=npower(ll,i_o)
    endif
   enddo
  enddo
 endif

 if(i_type(i_o).eq.3)then
  nc(kk)=4
  do m=1,n_gauss(i_o)
   do ii=1,nc(kk)
    if(ii.eq.1)d(m,ii,kk)=c_gauss(1,m,i_o)
    if(ii.gt.1)d(m,ii,kk)=c_gauss(2,m,i_o)
   enddo
  enddo
  do ll=1,3
   do ii=1,nc(kk)
    if(ll.eq.(ii-1))then
     n_c(ll,ii,kk)=npower(ll,i_o)+2
    else
     n_c(ll,ii,kk)=npower(ll,i_o)
    endif
   enddo
  enddo
 endif

enddo ! i_o=1,4

rint=0.d0

do i1=1,n_gauss(n_orb(1))
 do i2=1,n_gauss(n_orb(2))
  do i3=1,n_gauss(n_orb(3))
   do i4=1,n_gauss(n_orb(4))

    do ii1=1,nc(1)
     do ii2=1,nc(2)
      do ii3=1,nc(3)
       do ii4=1,nc(4)

        rint=rint+   &
        d(i1,ii1,1)*d(i2,ii2,2)*d(i3,ii3,3)*d(i4,ii4,4)* &
        bielec_integral(g_gauss(i1,i),g_gauss(i2,k),g_gauss(i3,j),g_gauss(i4,l), &
                        center(1,i),center(1,k),center(1,j),center(1,l), &
                        n_c(1,ii1,1),n_c(1,ii2,2),n_c(1,ii3,3),n_c(1,ii4,4))

       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
enddo

gauss_ijkl=rint

end

subroutine one_elect(i,k,Sij,Vij,Kij)
include 'j.inc'

double precision Kij

!! local arrays:
dimension n_orb(4),nc(4),d(n_gauss_max,4,4),n_c(3,4,4)

n_orb(1)=i
n_orb(2)=k

do kk=1,2

 i_o=n_orb(kk)

 if(i_type(i_o).eq.1)then
  
  nc(kk)=1
  do m=1,n_gauss(i_o)
   do ii=1,nc(kk)
    d(m,ii,kk)=c_gauss(1,m,i_o)
   enddo
  enddo
  do ll=1,3
   do ii=1,nc(kk)
    n_c(ll,ii,kk)=npower(ll,i_o)
   enddo
  enddo
 endif

 if(i_type(i_o).eq.2)then
  nc(kk)=3
  do m=1,n_gauss(i_o)
   do ii=1,nc(kk)
   d(m,ii,kk)=c_gauss(2,m,i_o)
   enddo
  enddo
  do ll=1,3
   do ii=1,nc(kk)
    if(ll.eq.ii)then
     n_c(ll,ii,kk)=npower(ll,i_o)+2
    else
     n_c(ll,ii,kk)=npower(ll,i_o)
    endif
   enddo
  enddo
 endif

 if(i_type(i_o).eq.3)then
  nc(kk)=4
  do m=1,n_gauss(i_o)
   do ii=1,nc(kk)
    if(ii.eq.1)d(m,ii,kk)=c_gauss(1,m,i_o)
    if(ii.gt.1)d(m,ii,kk)=c_gauss(2,m,i_o)
   enddo
  enddo
  do ll=1,3
   do ii=1,nc(kk)
    if(ll.eq.(ii-1))then
     n_c(ll,ii,kk)=npower(ll,i_o)+2
    else
     n_c(ll,ii,kk)=npower(ll,i_o)
    endif
   enddo
  enddo
 endif

enddo ! i_o=1,2

Sij=0.d0
Vij=0.d0
Kij=0.d0

do i1=1,n_gauss(n_orb(1))
 do i2=1,n_gauss(n_orb(2))
  do ii1=1,nc(1)
   do ii2=1,nc(2)

!    write(37,*)'in one-lect=',d(i1,ii1,1)*d(i2,ii2,2)

    Sij=Sij+   &
    d(i1,ii1,1)*d(i2,ii2,2)* &
    one_electron_Sij(n_c(1,ii1,1),center(1,i),g_gauss(i1,i),n_c(1,ii2,2),center(1,k),g_gauss(i2,k))

    Vij=Vij+   &
    d(i1,ii1,1)*d(i2,ii2,2)* &
    one_electron_Vij(n_c(1,ii1,1),center(1,i),g_gauss(i1,i),n_c(1,ii2,2),center(1,k),g_gauss(i2,k))

    Kij=Kij+   &
    d(i1,ii1,1)*d(i2,ii2,2)* &
    one_electron_Kij(n_c(1,ii1,1),center(1,i),g_gauss(i1,i),n_c(1,ii2,2),center(1,k),g_gauss(i2,k))

   enddo
  enddo
 enddo
enddo

end

integer function number_atom(ATOM)
character*80 ATOM
character*80 HYDROGEN
character*80 HELIUM
character*80 LITHIUM
character*80 BERYLLIUM
character*80 BORON
character*80 CARBON
character*80 NITROGEN
character*80 OXYGEN
character*80 FLUORINE
character*80 NEON
character*80 SODIUM
character*80 MAGNESIUM
character*80 ALUMINUM
character*80 SILICON
character*80 PHOSPHORUS
character*80 SULFUR
character*80 CHLORINE
character*80 ARGON

if(ATOM.eq.'H')then
 number_atom=1
 return
endif
if(ATOM.eq.'He')then
 number_atom=2
 return
endif
if(ATOM.eq.'Li')then
 number_atom=3
 return
endif
if(ATOM.eq.'Be')then
 number_atom=4
 return
endif
if(ATOM.eq.'B')then
 number_atom=5
 return
endif
if(ATOM.eq.'C')then
 number_atom=6
 return
endif
if(ATOM.eq.'N')then
 number_atom=7
 return
endif
if(ATOM.eq.'O')then
 number_atom=8
 return
endif
if(ATOM.eq.'F')then
 number_atom=9
 return
endif
if(ATOM.eq.'Ne')then
 number_atom=10
 return
endif
if(ATOM.eq.'Na')then
 number_atom=11
 return
endif
if(ATOM.eq.'Mg')then
 number_atom=12
 return
endif
if(ATOM.eq.'Al')then
 number_atom=13
 return
endif
if(ATOM.eq.'Si')then
 number_atom=14
 return
endif
if(ATOM.eq.'P')then
 number_atom=15
 return
endif
if(ATOM.eq.'S')then
 number_atom=16
 return
endif
if(ATOM.eq.'Cl')then
 number_atom=17
 return
endif
if(ATOM.eq.'Ar')then
 i_atom=18
 return
endif
stop 'ATOM not defined'
end

subroutine read_fit_SMILES
include 'j.inc'

c_fit_exp(1, 1)= 0.770925400000000D+00
g_fit_exp(1, 1)= 0.607017700000000D+00

c_fit_exp(1, 2)= 0.481012932300000D+00
g_fit_exp(1, 2)= 0.851818663537752D+00

c_fit_exp(2, 2)= 0.208059273700000D+00
g_fit_exp(2, 2)= 0.151623292722915D+00

c_fit_exp(1, 3)= 0.355424585100000D+00
g_fit_exp(1, 3)= 0.222766058294164D+01

c_fit_exp(2, 3)= 0.343751184800000D+00
g_fit_exp(2, 3)= 0.405771156247865D+00

c_fit_exp(3, 3)= 0.107132229600000D+00
g_fit_exp(3, 3)= 0.109817510390447D+00

c_fit_exp(1, 4)= 0.247466170504506D+00
g_fit_exp(1, 4)= 0.521684453382758D+01

c_fit_exp(2, 4)= 0.317363686360385D+00
g_fit_exp(2, 4)= 0.954618276011928D+00

c_fit_exp(3, 4)= 0.248749179481746D+00
g_fit_exp(3, 4)= 0.265203410216609D+00

c_fit_exp(4, 4)= 0.595294511501881D-01
g_fit_exp(4, 4)= 0.880186277439521D-01

c_fit_exp(1, 5)= 0.172442089510056D+00
g_fit_exp(1, 5)= 0.113056369555928D+02

c_fit_exp(2, 5)= 0.247677417058381D+00
g_fit_exp(2, 5)= 0.207172817827251D+01

c_fit_exp(3, 5)= 0.278094350350624D+00
g_fit_exp(3, 5)= 0.578648483322853D+00

c_fit_exp(4, 5)= 0.180650313016640D+00
g_fit_exp(4, 5)= 0.197572457261057D+00

c_fit_exp(5, 5)= 0.348527400685348D-01
g_fit_exp(5, 5)= 0.744527174597236D-01

c_fit_exp(1, 6)= 0.121983877649254D+00
g_fit_exp(1, 6)= 0.231030314879614D+02

c_fit_exp(2, 6)= 0.184112601688539D+00
g_fit_exp(2, 6)= 0.423591553412138D+01

c_fit_exp(3, 6)= 0.241817473442379D+00
g_fit_exp(3, 6)= 0.118505651867278D+01

c_fit_exp(4, 6)= 0.238573127640704D+00
g_fit_exp(4, 6)= 0.407098898180031D+00

c_fit_exp(5, 6)= 0.131906430347041D+00
g_fit_exp(5, 6)= 0.158088415115159D+00

c_fit_exp(6, 6)= 0.212215075497035D-01
g_fit_exp(6, 6)= 0.651095395441374D-01

c_fit_exp(1, 7)= 0.878033213887358D-01
g_fit_exp(1, 7)= 0.450587089159263D+02

c_fit_exp(2, 7)= 0.135628748317659D+00
g_fit_exp(2, 7)= 0.826349795861767D+01

c_fit_exp(3, 7)= 0.191668575025147D+00
g_fit_exp(3, 7)= 0.231329257358592D+01

c_fit_exp(4, 7)= 0.229568572069155D+00
g_fit_exp(4, 7)= 0.796146258719258D+00

c_fit_exp(5, 7)= 0.201809727849932D+00
g_fit_exp(5, 7)= 0.311276302660432D+00

c_fit_exp(6, 7)= 0.969446140814731D-01
g_fit_exp(6, 7)= 0.132280698472294D+00

c_fit_exp(7, 7)= 0.133287719635348D-01
g_fit_exp(7, 7)= 0.582363437713665D-01

c_fit_exp(1, 8)= 0.642516313276343D-01
g_fit_exp(1, 8)= 0.845781563304623D+02

c_fit_exp(2, 8)= 0.100385624210043D+00
g_fit_exp(2, 8)= 0.155129625353873D+02

c_fit_exp(3, 8)= 0.147063430327926D+00
g_fit_exp(3, 8)= 0.434394734098412D+01

c_fit_exp(4, 8)= 0.193752560197877D+00
g_fit_exp(4, 8)= 0.149607112612222D+01

c_fit_exp(5, 8)= 0.213469026225290D+00
g_fit_exp(5, 8)= 0.586147458716441D+00

c_fit_exp(6, 8)= 0.169178581226784D+00
g_fit_exp(6, 8)= 0.250952839465497D+00

c_fit_exp(7, 8)= 0.717426601419014D-01
g_fit_exp(7, 8)= 0.114110939145779D+00

c_fit_exp(8, 8)= 0.858707209406920D-02
g_fit_exp(8, 8)= 0.529406321961345D-01

c_fit_exp(1, 9)= 0.477206905701700D-01
g_fit_exp(1, 9)= 0.153728624497012D+03

c_fit_exp(2, 9)= 0.749947953184247D-01
g_fit_exp(2, 9)= 0.281979359285725D+02

c_fit_exp(3, 9)= 0.111913605320956D+00
g_fit_exp(3, 9)= 0.789709680353634D+01

c_fit_exp(4, 9)= 0.154827763712346D+00
g_fit_exp(4, 9)= 0.272064070449719D+01

c_fit_exp(5, 9)= 0.191491258069609D+00
g_fit_exp(5, 9)= 0.106675720199037D+01

c_fit_exp(6, 9)= 0.195508945402579D+00
g_fit_exp(6, 9)= 0.457779546063757D+00

c_fit_exp(7, 9)= 0.141014529442202D+00
g_fit_exp(7, 9)= 0.209837187888388D+00

c_fit_exp(8, 9)= 0.534591749350126D-01
g_fit_exp(8, 9)= 0.100627425452267D+00

c_fit_exp(9, 9)= 0.565173408345841D-02
g_fit_exp(9, 9)= 0.487180150656087D-01

c_fit_exp(1,10)= 0.359129962819811D-01
g_fit_exp(1,10)= 0.271813821974153D+03

c_fit_exp(2,10)= 0.566140388656530D-01
g_fit_exp(2,10)= 0.498594956327534D+02

c_fit_exp(3,10)= 0.853172206215175D-01
g_fit_exp(3,10)= 0.139646272657160D+02

c_fit_exp(4,10)= 0.121137987302091D+00
g_fit_exp(4,10)= 0.481171730584777D+01

c_fit_exp(5,10)= 0.159315391058145D+00
g_fit_exp(5,10)= 0.188731442181522D+01

c_fit_exp(6,10)= 0.186027678035103D+00
g_fit_exp(6,10)= 0.810611495451433D+00

c_fit_exp(7,10)= 0.177048733885195D+00
g_fit_exp(7,10)= 0.372519490811607D+00

c_fit_exp(8,10)= 0.117123692038386D+00
g_fit_exp(8,10)= 0.180179052359733D+00

c_fit_exp(9,10)= 0.401022877801723D-01
g_fit_exp(9,10)= 0.902214251178216D-01

c_fit_exp(10,10)= 0.378856726819921D-02
g_fit_exp(10,10)= 0.452609263385258D-01

c_fit_exp(1,11)= 0.273444184322038D-01
g_fit_exp(1,11)= 0.469212624533611D+03

c_fit_exp(2,11)= 0.431796908811632D-01
g_fit_exp(2,11)= 0.860703463180347D+02

c_fit_exp(3,11)= 0.654222807364076D-01
g_fit_exp(3,11)= 0.241074885797090D+02

c_fit_exp(4,11)= 0.942206675295390D-01
g_fit_exp(4,11)= 0.830726346088361D+01

c_fit_exp(5,11)= 0.128154921716897D+00
g_fit_exp(5,11)= 0.325893549836361D+01

c_fit_exp(6,11)= 0.161055324631454D+00
g_fit_exp(6,11)= 0.140026148154793D+01

c_fit_exp(7,11)= 0.178316132712128D+00
g_fit_exp(7,11)= 0.644115393750965D+00

c_fit_exp(8,11)= 0.158966346575069D+00
g_fit_exp(8,11)= 0.312414941089569D+00

c_fit_exp(9,11)= 0.970797311049486D-01
g_fit_exp(9,11)= 0.157857443058638D+00

c_fit_exp(10,11)= 0.302755090529202D-01
g_fit_exp(10,11)= 0.819432989927229D-01

c_fit_exp(11,11)= 0.258045035778310D-02
g_fit_exp(11,11)= 0.423706164631500D-01

c_fit_exp(1,12)= 0.210377625472215D-01
g_fit_exp(1,12)= 0.793039446170553D+03

c_fit_exp(2,12)= 0.332525942917884D-01
g_fit_exp(2,12)= 0.145473159894310D+03

c_fit_exp(3,12)= 0.505340129967287D-01
g_fit_exp(3,12)= 0.407465424609267D+02

c_fit_exp(4,12)= 0.733630479964681D-01
g_fit_exp(4,12)= 0.140415863585192D+02

c_fit_exp(5,12)= 0.101691032042973D+00
g_fit_exp(5,12)= 0.550900166024451D+01

c_fit_exp(6,12)= 0.133189763990981D+00
g_fit_exp(6,12)= 0.236748289304553D+01

c_fit_exp(7,12)= 0.160547093545091D+00
g_fit_exp(7,12)= 0.108948879478615D+01

c_fit_exp(8,12)= 0.169116873969973D+00
g_fit_exp(8,12)= 0.528994892587186D+00

c_fit_exp(9,12)= 0.141796964904962D+00
g_fit_exp(9,12)= 0.268099992187749D+00

c_fit_exp(10,12)= 0.803817812532298D-01
g_fit_exp(10,12)= 0.140493996469797D+00

c_fit_exp(11,12)= 0.229957930872552D-01
g_fit_exp(11,12)= 0.751970913050165D-01

c_fit_exp(12,12)= 0.178247441957718D-02
g_fit_exp(12,12)= 0.399126360758679D-01

c_fit_exp(1,13)= 0.163369362051693D-01
g_fit_exp(1,13)= 0.131540702517948D+04

c_fit_exp(2,13)= 0.258366067534699D-01
g_fit_exp(2,13)= 0.241296302024407D+03

c_fit_exp(3,13)= 0.393322430996508D-01
g_fit_exp(3,13)= 0.675871271935986D+02

c_fit_exp(4,13)= 0.573641705576576D-01
g_fit_exp(4,13)= 0.232916593038547D+02

c_fit_exp(5,13)= 0.803847670294060D-01
g_fit_exp(5,13)= 0.913858194234606D+01

c_fit_exp(6,13)= 0.107817759792989D+00
g_fit_exp(6,13)= 0.392767361983254D+01

c_fit_exp(7,13)= 0.136478173904019D+00
g_fit_exp(7,13)= 0.180783719785555D+01

c_fit_exp(8,13)= 0.158234458584262D+00
g_fit_exp(8,13)= 0.878187280586540D+00

c_fit_exp(9,13)= 0.159021520155440D+00
g_fit_exp(9,13)= 0.445590868716678D+00

c_fit_exp(10,13)= 0.125839792603740D+00
g_fit_exp(10,13)= 0.234263318149547D+00

c_fit_exp(11,13)= 0.665338613575184D-01
g_fit_exp(11,13)= 0.126626724804897D+00

c_fit_exp(12,13)= 0.175669900438330D-01
g_fit_exp(12,13)= 0.695902243009532D-01

c_fit_exp(13,13)= 0.124678621204025D-02
g_fit_exp(13,13)= 0.377925974151797D-01

c_fit_exp(1,14)= 0.127932368770716D-01
g_fit_exp(1,14)= 0.214537707477934D+04

c_fit_exp(2,14)= 0.202388392678375D-01
g_fit_exp(2,14)= 0.393546060451553D+03

c_fit_exp(3,14)= 0.308419806231812D-01
g_fit_exp(3,14)= 0.110233113337968D+03

c_fit_exp(4,14)= 0.451032650177116D-01
g_fit_exp(4,14)= 0.379887452963752D+02

c_fit_exp(5,14)= 0.636088831927261D-01
g_fit_exp(5,14)= 0.149054744139986D+02

c_fit_exp(6,14)= 0.865168505250275D-01
g_fit_exp(6,14)= 0.640657955459644D+01

c_fit_exp(7,14)= 0.112706811884299D+00
g_fit_exp(7,14)= 0.294914814559617D+01

c_fit_exp(8,14)= 0.138244899127121D+00
g_fit_exp(8,14)= 0.143291492815084D+01

c_fit_exp(9,14)= 0.154502937420846D+00
g_fit_exp(9,14)= 0.727423635405543D+00

c_fit_exp(10,14)= 0.148481978950622D+00
g_fit_exp(10,14)= 0.382915338818619D+00

c_fit_exp(11,14)= 0.111235294061802D+00
g_fit_exp(11,14)= 0.207692693953825D+00

c_fit_exp(12,14)= 0.550814204181380D-01
g_fit_exp(12,14)= 0.115310760470963D+00

c_fit_exp(13,14)= 0.134927368310292D-01
g_fit_exp(13,14)= 0.648537048783197D-01

c_fit_exp(14,14)= 0.881963752875984D-03
g_fit_exp(14,14)= 0.359421744229747D-01

c_fit_exp(1,15)= 0.100945583030874D-01
g_fit_exp(1,15)= 0.344609816497214D+04

c_fit_exp(2,15)= 0.159726291595280D-01
g_fit_exp(2,15)= 0.632150419226417D+03

c_fit_exp(3,15)= 0.243555714913702D-01
g_fit_exp(3,15)= 0.177067492270140D+03

c_fit_exp(4,15)= 0.356751800975631D-01
g_fit_exp(4,15)= 0.610218749214765D+02

c_fit_exp(5,15)= 0.505053110761937D-01
g_fit_exp(5,15)= 0.239432874126044D+02

c_fit_exp(6,15)= 0.692706168013029D-01
g_fit_exp(6,15)= 0.102914852004936D+02

c_fit_exp(7,15)= 0.918042218103435D-01
g_fit_exp(7,15)= 0.473777505440657D+01

c_fit_exp(8,15)= 0.116465047371226D+00
g_fit_exp(8,15)= 0.230222985655960D+01

c_fit_exp(9,15)= 0.138698887213232D+00
g_fit_exp(9,15)= 0.116901532295553D+01

c_fit_exp(10,15)= 0.149682964158876D+00
g_fit_exp(10,15)= 0.615705406584236D+00

c_fit_exp(11,15)= 0.137836969363124D+00
g_fit_exp(11,15)= 0.334412216082250D+00

c_fit_exp(12,15)= 0.980197748885386D-01
g_fit_exp(12,15)= 0.186343241425286D+00

c_fit_exp(13,15)= 0.456250975215059D-01
g_fit_exp(13,15)= 0.105910013657109D+00

c_fit_exp(14,15)= 0.104165872776284D-01
g_fit_exp(14,15)= 0.607971065105486D-01

c_fit_exp(15,15)= 0.630286148882694D-03
g_fit_exp(15,15)= 0.343106159988495D-01

c_fit_exp(1,16)= 0.802044725896653D-02
g_fit_exp(1,16)= 0.545918839901298D+04

c_fit_exp(2,16)= 0.126922530182093D-01
g_fit_exp(2,16)= 0.100143181846528D+04

c_fit_exp(3,16)= 0.193607407305483D-01
g_fit_exp(3,16)= 0.280505220653120D+03

c_fit_exp(4,16)= 0.283868275774804D-01
g_fit_exp(4,16)= 0.966696144131766D+02

c_fit_exp(5,16)= 0.402807944911994D-01
g_fit_exp(5,16)= 0.379308587975430D+02

c_fit_exp(6,16)= 0.555286741739241D-01
g_fit_exp(6,16)= 0.163040401241013D+02

c_fit_exp(7,16)= 0.743653919186775D-01
g_fit_exp(7,16)= 0.750597003924092D+01

c_fit_exp(8,16)= 0.962953110388442D-01
g_fit_exp(8,16)= 0.364761762693163D+01

c_fit_exp(9,16)= 0.119196557437036D+00
g_fit_exp(9,16)= 0.185240295110498D+01

c_fit_exp(10,16)= 0.138032058078675D+00
g_fit_exp(10,16)= 0.975892064442538D+00

c_fit_exp(11,16)= 0.144054670127166D+00
g_fit_exp(11,16)= 0.530357245983755D+00

c_fit_exp(12,16)= 0.127334835960666D+00
g_fit_exp(12,16)= 0.295960378921370D+00

c_fit_exp(13,16)= 0.861633938491029D-01
g_fit_exp(13,16)= 0.168857591434998D+00

c_fit_exp(14,16)= 0.378227842619708D-01
g_fit_exp(14,16)= 0.979814983667451D-01

c_fit_exp(15,16)= 0.808073645625427D-02
g_fit_exp(15,16)= 0.572818618534217D-01

c_fit_exp(16,16)= 0.454633135926318D-03
g_fit_exp(16,16)= 0.328593868092531D-01

c_fit_exp(1,17)= 0.641300873689383D-02
g_fit_exp(1,17)= 0.853917230601150D+04

c_fit_exp(2,17)= 0.101492390854659D-01
g_fit_exp(2,17)= 0.156642426011542D+04

c_fit_exp(3,17)= 0.154851830899902D-01
g_fit_exp(3,17)= 0.438762676175976D+03

c_fit_exp(4,17)= 0.227182517738975D-01
g_fit_exp(4,17)= 0.151209873247835D+03

c_fit_exp(5,17)= 0.322834450575313D-01
g_fit_exp(5,17)= 0.593315344960561D+02

c_fit_exp(6,17)= 0.446438892357398D-01
g_fit_exp(6,17)= 0.255031135240893D+02

c_fit_exp(7,17)= 0.601753781851046D-01
g_fit_exp(7,17)= 0.117412399145964D+02

c_fit_exp(8,17)= 0.789131676045026D-01
g_fit_exp(8,17)= 0.570601613552084D+01

c_fit_exp(9,17)= 0.100039249761419D+00
g_fit_exp(9,17)= 0.289794190780702D+01

c_fit_exp(10,17)= 0.121001515120965D+00
g_fit_exp(10,17)= 0.152691578391376D+01

c_fit_exp(11,17)= 0.136419142564113D+00
g_fit_exp(11,17)= 0.830051162259839D+00

c_fit_exp(12,17)= 0.137853115085293D+00
g_fit_exp(12,17)= 0.463498194588053D+00

c_fit_exp(13,17)= 0.117152620891234D+00
g_fit_exp(13,17)= 0.264856531148730D+00

c_fit_exp(14,17)= 0.755963354502840D-01
g_fit_exp(14,17)= 0.154302832385906D+00

c_fit_exp(15,17)= 0.313860793335296D-01
g_fit_exp(15,17)= 0.912077609490478D-01

c_fit_exp(16,17)= 0.629741865016941D-02
g_fit_exp(16,17)= 0.542047458651128D-01

c_fit_exp(17,17)= 0.330738280330244D-03
g_fit_exp(17,17)= 0.315586674456672D-01

c_fit_exp(1,18)= 0.515772064542778D-02
g_fit_exp(1,18)= 0.132017859483798D+05

c_fit_exp(2,18)= 0.816298875098312D-02
g_fit_exp(2,18)= 0.242173439991410D+04

c_fit_exp(3,18)= 0.124564536829482D-01
g_fit_exp(3,18)= 0.678339697629839D+03

c_fit_exp(4,18)= 0.182817615829007D-01
g_fit_exp(4,18)= 0.233775289561660D+03

c_fit_exp(5,18)= 0.260024041416456D-01
g_fit_exp(5,18)= 0.917288074500271D+02

c_fit_exp(6,18)= 0.360287185971759D-01
g_fit_exp(6,18)= 0.394290688817472D+02

c_fit_exp(7,18)= 0.487596388886758D-01
g_fit_exp(7,18)= 0.181527712383434D+02

c_fit_exp(8,18)= 0.644513474059740D-01
g_fit_exp(8,18)= 0.882210053608991D+01

c_fit_exp(9,18)= 0.829354670073262D-01
g_fit_exp(9,18)= 0.448070749665194D+01

c_fit_exp(10,18)= 0.103085066271499D+00
g_fit_exp(10,18)= 0.236104847448299D+01

c_fit_exp(11,18)= 0.121975623649862D+00
g_fit_exp(11,18)= 0.128368498863673D+01

c_fit_exp(12,18)= 0.134017978331351D+00
g_fit_exp(12,18)= 0.717025507322824D+00

c_fit_exp(13,18)= 0.131273523185779D+00
g_fit_exp(13,18)= 0.410011636829473D+00

c_fit_exp(14,18)= 0.107411730797403D+00
g_fit_exp(14,18)= 0.239263318473613D+00

c_fit_exp(15,18)= 0.662265932293839D-01
g_fit_exp(15,18)= 0.142018643684173D+00

c_fit_exp(16,18)= 0.260743890920158D-01
g_fit_exp(16,18)= 0.853556764114360D-01

c_fit_exp(17,18)= 0.492895597179743D-02
g_fit_exp(17,18)= 0.514872880447035D-01

c_fit_exp(18,18)= 0.242502872530627D-03
g_fit_exp(18,18)= 0.303849994918093D-01

c_fit_exp(1,19)= 0.417056542116485D-02
g_fit_exp(1,19)= 0.201912902379125D+05

c_fit_exp(2,19)= 0.660083426493615D-02
g_fit_exp(2,19)= 0.370388950634204D+04

c_fit_exp(3,19)= 0.100735740023136D-01
g_fit_exp(3,19)= 0.103747831539269D+04

c_fit_exp(4,19)= 0.147880819060227D-01
g_fit_exp(4,19)= 0.357545212313871D+03

c_fit_exp(5,19)= 0.210452880323935D-01
g_fit_exp(5,19)= 0.140294031366189D+03

c_fit_exp(6,19)= 0.291965183933903D-01
g_fit_exp(6,19)= 0.603048107461182D+02

c_fit_exp(7,19)= 0.396145989757254D-01
g_fit_exp(7,19)= 0.277639929615564D+02

c_fit_exp(8,19)= 0.526269544047791D-01
g_fit_exp(8,19)= 0.134932706020755D+02

c_fit_exp(9,19)= 0.683643209823631D-01
g_fit_exp(9,19)= 0.685334696760377D+01

c_fit_exp(10,19)= 0.864547092236644D-01
g_fit_exp(10,19)= 0.361143936196586D+01

c_fit_exp(11,19)= 0.105481361236077D+00
g_fit_exp(11,19)= 0.196367349206757D+01

c_fit_exp(12,19)= 0.122209754755901D+00
g_fit_exp(12,19)= 0.109701898208572D+01

c_fit_exp(13,19)= 0.130970082701948D+00
g_fit_exp(13,19)= 0.627507515122241D+00

c_fit_exp(14,19)= 0.124476281320858D+00
g_fit_exp(14,19)= 0.366453407241060D+00

c_fit_exp(15,19)= 0.981906062360255D-01
g_fit_exp(15,19)= 0.217893741841034D+00

c_fit_exp(16,19)= 0.579518177444959D-01
g_fit_exp(16,19)= 0.131525673296283D+00

c_fit_exp(17,19)= 0.216883635722516D-01
g_fit_exp(17,19)= 0.802503618512026D-01

c_fit_exp(18,19)= 0.387373815692823D-02
g_fit_exp(18,19)= 0.490687595996008D-01

c_fit_exp(19,19)= 0.179103380996411D-03
g_fit_exp(19,19)= 0.293196533363169D-01

c_fit_exp(1,20)= 0.338925524518914D-02
g_fit_exp(1,20)= 0.305737707306057D+05

c_fit_exp(2,20)= 0.536433869133321D-02
g_fit_exp(2,20)= 0.560845238081001D+04

c_fit_exp(3,20)= 0.818702846043344D-02
g_fit_exp(3,20)= 0.157095673359472D+04

c_fit_exp(4,20)= 0.120204765499523D-01
g_fit_exp(4,20)= 0.541397851053993D+03

c_fit_exp(5,20)= 0.171128956812078D-01
g_fit_exp(5,20)= 0.212434696299232D+03

c_fit_exp(6,20)= 0.237600102172834D-01
g_fit_exp(6,20)= 0.913144457376170D+02

c_fit_exp(7,20)= 0.322912173598550D-01
g_fit_exp(7,20)= 0.420408724571937D+02

c_fit_exp(8,20)= 0.430364681786010D-01
g_fit_exp(8,20)= 0.204320044312754D+02

c_fit_exp(9,20)= 0.562465757844937D-01
g_fit_exp(9,20)= 0.103777516144611D+02

c_fit_exp(10,20)= 0.719231157149912D-01
g_fit_exp(10,20)= 0.546880754521819D+01

c_fit_exp(11,20)= 0.894938900094020D-01
g_fit_exp(11,20)= 0.297373529207063D+01

c_fit_exp(12,20)= 0.107275992448720D+00
g_fit_exp(12,20)= 0.166144190158044D+01

c_fit_exp(13,20)= 0.121789617521563D+00
g_fit_exp(13,20)= 0.950525608188701D+00

c_fit_exp(14,20)= 0.127401418663432D+00
g_fit_exp(14,20)= 0.555286839668447D+00

c_fit_exp(15,20)= 0.117591681569608D+00
g_fit_exp(15,20)= 0.330433600152179D+00

c_fit_exp(16,20)= 0.895350439379999D-01
g_fit_exp(16,20)= 0.199823032282570D+00

c_fit_exp(17,20)= 0.506672131666828D-01
g_fit_exp(17,20)= 0.122468407586625D+00

c_fit_exp(18,20)= 0.180636386890479D-01
g_fit_exp(18,20)= 0.757582532150385D-01

c_fit_exp(19,20)= 0.305632563014643D-02
g_fit_exp(19,20)= 0.469014624281249D-01

c_fit_exp(20,20)= 0.133175125910862D-03
g_fit_exp(20,20)= 0.283474986144323D-01

c_fit_exp(1,21)= 0.276717596651701D-02
g_fit_exp(1,21)= 0.458655224600659D+05

c_fit_exp(2,21)= 0.437979610115597D-02
g_fit_exp(2,21)= 0.841357227966295D+04

c_fit_exp(3,21)= 0.668467725241485D-02
g_fit_exp(3,21)= 0.235668607156192D+04

c_fit_exp(4,21)= 0.981565641944583D-02
g_fit_exp(4,21)= 0.812183710722911D+03

c_fit_exp(5,21)= 0.139773269214216D-01
g_fit_exp(5,21)= 0.318686486188461D+03

c_fit_exp(6,21)= 0.194165554067809D-01
g_fit_exp(6,21)= 0.136986733337473D+03

c_fit_exp(7,21)= 0.264163413612259D-01
g_fit_exp(7,21)= 0.630684570899178D+02

c_fit_exp(8,21)= 0.352802879381937D-01
g_fit_exp(8,21)= 0.306516601441670D+02

c_fit_exp(9,21)= 0.462915992868701D-01
g_fit_exp(9,21)= 0.155686387605128D+02

c_fit_exp(10,21)= 0.596204916131332D-01
g_fit_exp(10,21)= 0.820441065293146D+01

c_fit_exp(11,21)= 0.751374789964347D-01
g_fit_exp(11,21)= 0.446138492939706D+01

c_fit_exp(12,21)= 0.920765642387414D-01
g_fit_exp(12,21)= 0.249272713781956D+01

c_fit_exp(13,21)= 0.108515961357975D+00
g_fit_exp(13,21)= 0.142624496542728D+01

c_fit_exp(14,21)= 0.120795635966311D+00
g_fit_exp(14,21)= 0.833349048962535D+00

c_fit_exp(15,21)= 0.123423293314672D+00
g_fit_exp(15,21)= 0.496087481824107D+00

c_fit_exp(16,21)= 0.110724061838948D+00
g_fit_exp(16,21)= 0.300248834374300D+00

c_fit_exp(17,21)= 0.814660775158014D-01
g_fit_exp(17,21)= 0.184371188353393D+00

c_fit_exp(18,21)= 0.442702216188023D-01
g_fit_exp(18,21)= 0.114577834094565D+00

c_fit_exp(19,21)= 0.150649743817410D-01
g_fit_exp(19,21)= 0.717755301098176D-01

c_fit_exp(20,21)= 0.242035302389690D-02
g_fit_exp(20,21)= 0.449473427872493D-01

c_fit_exp(21,21)= 0.996496112254142D-04
g_fit_exp(21,21)= 0.274561391108200D-01

c_fit_exp(1,22)= 0.226913031761205D-02
g_fit_exp(1,22)= 0.682091371189616D+05

c_fit_exp(2,22)= 0.359153382540182D-02
g_fit_exp(2,22)= 0.125122863577308D+05

c_fit_exp(3,22)= 0.548172544664007D-02
g_fit_exp(3,22)= 0.350475817369770D+04

c_fit_exp(4,22)= 0.804979252715866D-02
g_fit_exp(4,22)= 0.120784373775441D+04

c_fit_exp(5,22)= 0.114645468455698D-01
g_fit_exp(5,22)= 0.473936786485527D+03

c_fit_exp(6,22)= 0.159313390723039D-01
g_fit_exp(6,22)= 0.203721009984637D+03

c_fit_exp(7,22)= 0.216898123883615D-01
g_fit_exp(7,22)= 0.937930102724964D+02

c_fit_exp(8,22)= 0.290074427842301D-01
g_fit_exp(8,22)= 0.455841568024849D+02

c_fit_exp(9,22)= 0.381592527860377D-01
g_fit_exp(9,22)= 0.231533256474799D+02

c_fit_exp(10,22)= 0.493784355047882D-01
g_fit_exp(10,22)= 0.122015459440783D+02

c_fit_exp(11,22)= 0.627515458055122D-01
g_fit_exp(11,22)= 0.663506399890864D+01

c_fit_exp(12,22)= 0.780178359696196D-01
g_fit_exp(12,22)= 0.370735091051748D+01

c_fit_exp(13,22)= 0.942265798031771D-01
g_fit_exp(13,22)= 0.212132300512291D+01

c_fit_exp(14,22)= 0.109247049688351D+00
g_fit_exp(14,22)= 0.123960463639945D+01

c_fit_exp(15,22)= 0.119302737852567D+00
g_fit_exp(15,22)= 0.738073056186813D+00

c_fit_exp(16,22)= 0.119133363977060D+00
g_fit_exp(16,22)= 0.446887144736002D+00

c_fit_exp(17,22)= 0.103955796201880D+00
g_fit_exp(17,22)= 0.274657853812604D+00

c_fit_exp(18,22)= 0.739866428669969D-01
g_fit_exp(18,22)= 0.171028757590205D+00

c_fit_exp(19,22)= 0.386637438851140D-01
g_fit_exp(19,22)= 0.107647183832562D+00

c_fit_exp(20,22)= 0.125814001630370D-01
g_fit_exp(20,22)= 0.682204469331389D-01

c_fit_exp(21,22)= 0.192349510437292D-02
g_fit_exp(21,22)= 0.431757309717551D-01

c_fit_exp(22,22)= 0.750042272901905D-04
g_fit_exp(22,22)= 0.266353445584859D-01

c_fit_exp(1,23)= 0.186832903343599D-02
g_fit_exp(1,23)= 0.100613375749796D+06

c_fit_exp(2,23)= 0.295716895175938D-02
g_fit_exp(2,23)= 0.184565219298485D+05

c_fit_exp(3,23)= 0.451357408924910D-02
g_fit_exp(3,23)= 0.516977088227441D+04

c_fit_exp(4,23)= 0.662837127150985D-02
g_fit_exp(4,23)= 0.178165697824450D+04

c_fit_exp(5,23)= 0.944112205744548D-02
g_fit_exp(5,23)= 0.699091393433676D+03

c_fit_exp(6,23)= 0.131224944853285D-01
g_fit_exp(6,23)= 0.300503623458888D+03

c_fit_exp(7,23)= 0.178739382793442D-01
g_fit_exp(7,23)= 0.138351858962568D+03

c_fit_exp(8,23)= 0.239258148140032D-01
g_fit_exp(8,23)= 0.672402843460585D+02

c_fit_exp(9,23)= 0.315281032926786D-01
g_fit_exp(9,23)= 0.341531561167682D+02

c_fit_exp(10,23)= 0.409248141726735D-01
g_fit_exp(10,23)= 0.179984625127351D+02

c_fit_exp(11,23)= 0.522963197509083D-01
g_fit_exp(11,23)= 0.978747764172284D+01

c_fit_exp(12,23)= 0.656433364815351D-01
g_fit_exp(12,23)= 0.546887361936766D+01

c_fit_exp(13,23)= 0.805752610913946D-01
g_fit_exp(13,23)= 0.312935939592224D+01

c_fit_exp(14,23)= 0.959680141921650D-01
g_fit_exp(14,23)= 0.182876183146045D+01

c_fit_exp(15,23)= 0.109513669104953D+00
g_fit_exp(15,23)= 0.108897986556470D+01

c_fit_exp(16,23)= 0.117380358637583D+00
g_fit_exp(16,23)= 0.659491111496733D+00

c_fit_exp(17,23)= 0.114616686652238D+00
g_fit_exp(17,23)= 0.405498857758017D+00

c_fit_exp(18,23)= 0.973506421104138D-01
g_fit_exp(18,23)= 0.252737499166398D+00

c_fit_exp(19,23)= 0.670864729114690D-01
g_fit_exp(19,23)= 0.159407184945529D+00

c_fit_exp(20,23)= 0.337576501862480D-01
g_fit_exp(20,23)= 0.101515011171007D+00

c_fit_exp(21,23)= 0.105218893066912D-01
g_fit_exp(21,23)= 0.650277363382350D-01

c_fit_exp(22,23)= 0.153378855302899D-02
g_fit_exp(22,23)= 0.415615994482735D-01

c_fit_exp(23,23)= 0.567669834633619D-04
g_fit_exp(23,23)= 0.258765772414374D-01

c_fit_exp(1,24)= 0.154423170516553D-02
g_fit_exp(1,24)= 0.147278106599185D+06

c_fit_exp(2,24)= 0.244420001363844D-02
g_fit_exp(2,24)= 0.270167030914145D+05

c_fit_exp(3,24)= 0.373066221089795D-02
g_fit_exp(3,24)= 0.756752417616781D+04

c_fit_exp(4,24)= 0.547879115887009D-02
g_fit_exp(4,24)= 0.260799456756108D+04

c_fit_exp(5,24)= 0.780425164378925D-02
g_fit_exp(5,24)= 0.102333228248235D+04

c_fit_exp(6,24)= 0.108489825460834D-01
g_fit_exp(6,24)= 0.439878433963367D+03

c_fit_exp(7,24)= 0.147817713765211D-01
g_fit_exp(7,24)= 0.202520214000420D+03

c_fit_exp(8,24)= 0.197986218895006D-01
g_fit_exp(8,24)= 0.984268638893637D+02

c_fit_exp(9,24)= 0.261192118418234D-01
g_fit_exp(9,24)= 0.499938035341068D+02

c_fit_exp(10,24)= 0.339742795967016D-01
g_fit_exp(10,24)= 0.263464913927587D+02

c_fit_exp(11,24)= 0.435746153869190D-01
g_fit_exp(11,24)= 0.143272066815672D+02

c_fit_exp(12,24)= 0.550453989141138D-01
g_fit_exp(12,24)= 0.800560485482314D+01

c_fit_exp(13,24)= 0.683002066065283D-01
g_fit_exp(13,24)= 0.458100574848274D+01

c_fit_exp(14,24)= 0.828214815101769D-01
g_fit_exp(14,24)= 0.267718244935672D+01

c_fit_exp(15,24)= 0.973251029350125D-01
g_fit_exp(15,24)= 0.159429135276150D+01

c_fit_exp(16,24)= 0.109358693164808D+00
g_fit_exp(16,24)= 0.965620511663150D+00

c_fit_exp(17,24)= 0.115092438442037D+00
g_fit_exp(17,24)= 0.593859553173269D+00

c_fit_exp(18,24)= 0.109946675783232D+00
g_fit_exp(18,24)= 0.370306721617428D+00

c_fit_exp(19,24)= 0.909565836042106D-01
g_fit_exp(19,24)= 0.233788459523829D+00

c_fit_exp(20,24)= 0.607458424255352D-01
g_fit_exp(20,24)= 0.149205251688718D+00

c_fit_exp(21,24)= 0.294694320665184D-01
g_fit_exp(21,24)= 0.960534048651219D-01

c_fit_exp(22,24)= 0.881175485936950D-02
g_fit_exp(22,24)= 0.621445939417817D-01

c_fit_exp(23,24)= 0.122696559474345D-02
g_fit_exp(23,24)= 0.400842751926843D-01

c_fit_exp(24,24)= 0.431876572948832D-04
g_fit_exp(24,24)= 0.251726299749625D-01

c_fit_exp(1,25)= 0.128096675493160D-02
g_fit_exp(1,25)= 0.214036509259696D+06

c_fit_exp(2,25)= 0.202751065960036D-02
g_fit_exp(2,25)= 0.392628686732832D+05

c_fit_exp(3,25)= 0.309467840958100D-02
g_fit_exp(3,25)= 0.109977417502483D+05

c_fit_exp(4,25)= 0.454488391382052D-02
g_fit_exp(4,25)= 0.379015026594806D+04

c_fit_exp(5,25)= 0.647424765338986D-02
g_fit_exp(5,25)= 0.148719017672689D+04

c_fit_exp(6,25)= 0.900099953388812D-02
g_fit_exp(6,25)= 0.639267555014913D+03

c_fit_exp(7,25)= 0.122664266496777D-01
g_fit_exp(7,25)= 0.294319239586281D+03

c_fit_exp(8,25)= 0.164362500971814D-01
g_fit_exp(8,25)= 0.143042275189838D+03

c_fit_exp(9,25)= 0.217000164360698D-01
g_fit_exp(9,25)= 0.726553770292016D+02

c_fit_exp(10,25)= 0.282654817957851D-01
g_fit_exp(10,25)= 0.382891501939069D+02

c_fit_exp(11,25)= 0.363425143090276D-01
g_fit_exp(11,25)= 0.208217240880472D+02

c_fit_exp(12,25)= 0.461066868881516D-01
g_fit_exp(12,25)= 0.116346392134755D+02

c_fit_exp(13,25)= 0.576260706132241D-01
g_fit_exp(13,25)= 0.665771926992637D+01

c_fit_exp(14,25)= 0.707266745720698D-01
g_fit_exp(14,25)= 0.389091979605828D+01

c_fit_exp(15,25)= 0.847682863619789D-01
g_fit_exp(15,25)= 0.231717225379771D+01

c_fit_exp(16,25)= 0.983216633596232D-01
g_fit_exp(16,25)= 0.140354306097018D+01

c_fit_exp(17,25)= 0.108823056456259D+00
g_fit_exp(17,25)= 0.863288101133850D+00

c_fit_exp(18,25)= 0.112497638077511D+00
g_fit_exp(18,25)= 0.538438419708802D+00

c_fit_exp(19,25)= 0.105186728416715D+00
g_fit_exp(19,25)= 0.340099759940840D+00

c_fit_exp(20,25)= 0.848093853886320D-01
g_fit_exp(20,25)= 0.217275270033416D+00

c_fit_exp(21,25)= 0.549395821606321D-01
g_fit_exp(21,25)= 0.140187727291110D+00

c_fit_exp(22,25)= 0.257250787594092D-01
g_fit_exp(22,25)= 0.911606444038214D-01

c_fit_exp(23,25)= 0.738996947047390D-02
g_fit_exp(23,25)= 0.595283122066580D-01

c_fit_exp(24,25)= 0.984557822351332D-03
g_fit_exp(24,25)= 0.387267714168891D-01

c_fit_exp(25,25)= 0.330188877862306D-04
g_fit_exp(25,25)= 0.245174914715798D-01

c_fit_exp(1,26)= 0.106621425297179D-02
g_fit_exp(1,26)= 0.308940551230892D+06

c_fit_exp(2,26)= 0.168760365012342D-02
g_fit_exp(2,26)= 0.566720721003190D+05

c_fit_exp(3,26)= 0.257587627318352D-02
g_fit_exp(3,26)= 0.158741543015879D+05

c_fit_exp(4,26)= 0.378301441645845D-02
g_fit_exp(4,26)= 0.547070806199865D+04

c_fit_exp(5,26)= 0.538912255476034D-02
g_fit_exp(5,26)= 0.214661262157819D+04

c_fit_exp(6,26)= 0.749288550954292D-02
g_fit_exp(6,26)= 0.922720016369172D+03

c_fit_exp(7,26)= 0.102126234518251D-01
g_fit_exp(7,26)= 0.424821143011071D+03

c_fit_exp(8,26)= 0.136880442555268D-01
g_fit_exp(8,26)= 0.206467742789176D+03

c_fit_exp(9,26)= 0.180810590655764D-01
g_fit_exp(9,26)= 0.104871171106271D+03

c_fit_exp(10,26)= 0.235738750315597D-01
g_fit_exp(10,26)= 0.552668872399399D+02

c_fit_exp(11,26)= 0.303611536979302D-01
g_fit_exp(11,26)= 0.300543586233681D+02

c_fit_exp(12,26)= 0.386303606887721D-01
g_fit_exp(12,26)= 0.167936894560282D+02

c_fit_exp(13,26)= 0.485200732775694D-01
g_fit_exp(13,26)= 0.960998187253436D+01

c_fit_exp(14,26)= 0.600397965571665D-01
g_fit_exp(14,26)= 0.561636927858683D+01

c_fit_exp(15,26)= 0.729283747605690D-01
g_fit_exp(15,26)= 0.334481424732257D+01

c_fit_exp(16,26)= 0.864284663158991D-01
g_fit_exp(16,26)= 0.202608099691825D+01

c_fit_exp(17,26)= 0.989818327463860D-01
g_fit_exp(17,26)= 0.124628481113614D+01

c_fit_exp(18,26)= 0.107946088866039D+00
g_fit_exp(18,26)= 0.777416451767433D+00

c_fit_exp(19,26)= 0.109649261130234D+00
g_fit_exp(19,26)= 0.491170994393143D+00

c_fit_exp(20,26)= 0.100390133295301D+00
g_fit_exp(20,26)= 0.313947375818574D+00

c_fit_exp(21,26)= 0.789329940938365D-01
g_fit_exp(21,26)= 0.202777225531258D+00

c_fit_exp(22,26)= 0.496373755169979D-01
g_fit_exp(22,26)= 0.132165826214873D+00

c_fit_exp(23,26)= 0.224574827045157D-01
g_fit_exp(23,26)= 0.867533761341587D-01

c_fit_exp(24,26)= 0.620615772120996D-02
g_fit_exp(24,26)= 0.571431520190398D-01

c_fit_exp(25,26)= 0.792363130231384D-03
g_fit_exp(25,26)= 0.374745546077727D-01

c_fit_exp(26,26)= 0.253614816758353D-04
g_fit_exp(26,26)= 0.239058868197446D-01

c_fit_exp(1,27)= 0.890332012605812D-03
g_fit_exp(1,27)= 0.443057524795007D+06

c_fit_exp(2,27)= 0.140921880121019D-02
g_fit_exp(2,27)= 0.812744981645833D+05

c_fit_exp(3,27)= 0.215097047284865D-02
g_fit_exp(3,27)= 0.227654276304345D+05

c_fit_exp(4,27)= 0.315901258843909D-02
g_fit_exp(4,27)= 0.784564721974483D+04

c_fit_exp(5,27)= 0.450029176612941D-02
g_fit_exp(5,27)= 0.307849856411580D+04

c_fit_exp(6,27)= 0.625737120455464D-02
g_fit_exp(6,27)= 0.132329079791175D+04

c_fit_exp(7,27)= 0.852946350261140D-02
g_fit_exp(7,27)= 0.609244482239657D+03

c_fit_exp(8,27)= 0.114342402917285D-01
g_fit_exp(8,27)= 0.296099665456092D+03

c_fit_exp(9,27)= 0.151092700187941D-01
g_fit_exp(9,27)= 0.150398049542812D+03

c_fit_exp(10,27)= 0.197120155931700D-01
g_fit_exp(10,27)= 0.792595710588365D+02

c_fit_exp(11,27)= 0.254165043722815D-01
g_fit_exp(11,27)= 0.431017790173602D+02

c_fit_exp(12,27)= 0.324031754654358D-01
g_fit_exp(12,27)= 0.240843808672152D+02

c_fit_exp(13,27)= 0.408356627511651D-01
g_fit_exp(13,27)= 0.137820730641332D+02

c_fit_exp(14,27)= 0.508140053808003D-01
g_fit_exp(14,27)= 0.805474501526848D+01

c_fit_exp(15,27)= 0.622881563242921D-01
g_fit_exp(15,27)= 0.479705646874746D+01

c_fit_exp(16,27)= 0.749110034392810D-01
g_fit_exp(16,27)= 0.290583240010356D+01

c_fit_exp(17,27)= 0.878147083812176D-01
g_fit_exp(17,27)= 0.178751413237868D+01

c_fit_exp(18,27)= 0.993290646333178D-01
g_fit_exp(18,27)= 0.111511139275615D+01

c_fit_exp(19,27)= 0.106764941414220D+00
g_fit_exp(19,27)= 0.704622950573761D+00

c_fit_exp(20,27)= 0.106595683207391D+00
g_fit_exp(20,27)= 0.450501140973978D+00

c_fit_exp(21,27)= 0.956021800824124D-01
g_fit_exp(21,27)= 0.291131788361478D+00

c_fit_exp(22,27)= 0.733433500688662D-01
g_fit_exp(22,27)= 0.189964062077806D+00

c_fit_exp(23,27)= 0.448075037411844D-01
g_fit_exp(23,27)= 0.124988923960910D+00

c_fit_exp(24,27)= 0.196074337260051D-01
g_fit_exp(24,27)= 0.827640904077737D-01

c_fit_exp(25,27)= 0.521915823447276D-02
g_fit_exp(25,27)= 0.549597616161479D-01

c_fit_exp(26,27)= 0.639487200141191D-03
g_fit_exp(26,27)= 0.363155331907311D-01

c_fit_exp(27,27)= 0.195656662489068D-04
g_fit_exp(27,27)= 0.233333745922059D-01

c_fit_exp(1,28)= 0.745740775205452D-03
g_fit_exp(1,28)= 0.631521964366233D+06

c_fit_exp(2,28)= 0.118036045417660D-02
g_fit_exp(2,28)= 0.115846426810800D+06

c_fit_exp(3,28)= 0.180165522837284D-02
g_fit_exp(3,28)= 0.324492128468479D+05

c_fit_exp(4,28)= 0.264600878786599D-02
g_fit_exp(4,28)= 0.111829695996673D+05

c_fit_exp(5,28)= 0.376952784935628D-02
g_fit_exp(5,28)= 0.438800740005945D+04

c_fit_exp(6,28)= 0.524145785829416D-02
g_fit_exp(6,28)= 0.188618262194621D+04

c_fit_exp(7,28)= 0.714513556435150D-02
g_fit_exp(7,28)= 0.868400653016163D+03

c_fit_exp(8,28)= 0.957970445543414D-02
g_fit_exp(8,28)= 0.422052627683795D+03

c_fit_exp(9,28)= 0.126617559809898D-01
g_fit_exp(9,28)= 0.214373528768774D+03

c_fit_exp(10,28)= 0.165262483220066D-01
g_fit_exp(10,28)= 0.112974675131382D+03

c_fit_exp(11,28)= 0.213255899147088D-01
g_fit_exp(11,28)= 0.614363319953987D+02

c_fit_exp(12,28)= 0.272248105923058D-01
g_fit_exp(12,28)= 0.343294376349108D+02

c_fit_exp(13,28)= 0.343890453571282D-01
g_fit_exp(13,28)= 0.196447962040350D+02

c_fit_exp(14,28)= 0.429567808472476D-01
g_fit_exp(14,28)= 0.114812079623139D+02

c_fit_exp(15,28)= 0.529882234493818D-01
g_fit_exp(15,28)= 0.683777898207794D+01

c_fit_exp(16,28)= 0.643732566130836D-01
g_fit_exp(16,28)= 0.414207434944819D+01

c_fit_exp(17,28)= 0.766807762448570D-01
g_fit_exp(17,28)= 0.254805266496236D+01

c_fit_exp(18,28)= 0.889400303842376D-01
g_fit_exp(18,28)= 0.158963235301454D+01

c_fit_exp(19,28)= 0.993864490113207D-01
g_fit_exp(19,28)= 0.100454526896880D+01

c_fit_exp(20,28)= 0.105314743353556D+00
g_fit_exp(20,28)= 0.642349943155279D+00

c_fit_exp(21,28)= 0.103380455401738D+00
g_fit_exp(21,28)= 0.415228119666079D+00

c_fit_exp(22,28)= 0.908605099027584D-01
g_fit_exp(22,28)= 0.271088576847539D+00

c_fit_exp(23,28)= 0.680490933018638D-01
g_fit_exp(23,28)= 0.178571419715022D+00

c_fit_exp(24,28)= 0.404171221088204D-01
g_fit_exp(24,28)= 0.118534395190990D+00

c_fit_exp(25,28)= 0.171223737357714D-01
g_fit_exp(25,28)= 0.791369209441557D-01

c_fit_exp(26,28)= 0.439509766671479D-02
g_fit_exp(26,28)= 0.529534569036018D-01

c_fit_exp(27,28)= 0.517504170275080D-03
g_fit_exp(27,28)= 0.352393705142203D-01

c_fit_exp(28,28)= 0.151574101444251D-04
g_fit_exp(28,28)= 0.227960903907883D-01

c_fit_exp(1,29)= 0.626448179908404D-03
g_fit_exp(1,29)= 0.894940240138575D+06

c_fit_exp(2,29)= 0.991544329215264D-03
g_fit_exp(2,29)= 0.164167891882899D+06

c_fit_exp(3,29)= 0.151345625743476D-02
g_fit_exp(3,29)= 0.459843179860603D+05

c_fit_exp(4,29)= 0.222275363546510D-02
g_fit_exp(4,29)= 0.158475721577823D+05

c_fit_exp(5,29)= 0.316658684864648D-02
g_fit_exp(5,29)= 0.621831852060851D+04

c_fit_exp(6,29)= 0.440317716688698D-02
g_fit_exp(6,29)= 0.267294108924146D+04

c_fit_exp(7,29)= 0.600266781667194D-02
g_fit_exp(7,29)= 0.123062533413754D+04

c_fit_exp(8,29)= 0.804868118739758D-02
g_fit_exp(8,29)= 0.598098103744377D+03

c_fit_exp(9,29)= 0.106399519832370D-01
g_fit_exp(9,29)= 0.303792574052465D+03

c_fit_exp(10,29)= 0.138916362631564D-01
g_fit_exp(10,29)= 0.160098548507995D+03

c_fit_exp(11,29)= 0.179356381166494D-01
g_fit_exp(11,29)= 0.870626841971220D+02

c_fit_exp(12,29)= 0.229187188257039D-01
g_fit_exp(12,29)= 0.486490351309506D+02

c_fit_exp(13,29)= 0.289961135016513D-01
g_fit_exp(13,29)= 0.278391700666100D+02

c_fit_exp(14,29)= 0.363166243245800D-01
g_fit_exp(14,29)= 0.162704005859194D+02

c_fit_exp(15,29)= 0.449924009406517D-01
g_fit_exp(15,29)= 0.969010998246118D+01

c_fit_exp(16,29)= 0.550427812249893D-01
g_fit_exp(16,29)= 0.586997477179005D+01

c_fit_exp(17,29)= 0.662975259026313D-01
g_fit_exp(17,29)= 0.361105607813852D+01

c_fit_exp(18,29)= 0.782442084028098D-01
g_fit_exp(18,29)= 0.225286279503855D+01

c_fit_exp(19,29)= 0.898175492336992D-01
g_fit_exp(19,29)= 0.142373175203151D+01

c_fit_exp(20,29)= 0.991765050033453D-01
g_fit_exp(20,29)= 0.910471970566278D+00

c_fit_exp(21,29)= 0.103628545054684D+00
g_fit_exp(21,29)= 0.588637940552827D+00

c_fit_exp(22,29)= 0.100042669431459D+00
g_fit_exp(22,29)= 0.384415078334652D+00

c_fit_exp(23,29)= 0.861962251492727D-01
g_fit_exp(23,29)= 0.253369828084075D+00

c_fit_exp(24,29)= 0.630533037142353D-01
g_fit_exp(24,29)= 0.168386145603840D+00

c_fit_exp(25,29)= 0.364335336344937D-01
g_fit_exp(25,29)= 0.112701870387142D+00

c_fit_exp(26,29)= 0.149560174999581D-01
g_fit_exp(26,29)= 0.758254792576768D-01

c_fit_exp(27,29)= 0.370613059906273D-02
g_fit_exp(27,29)= 0.511034537949892D-01

c_fit_exp(28,29)= 0.419880428946407D-03
g_fit_exp(28,29)= 0.342372421423812D-01

c_fit_exp(29,29)= 0.117890562639788D-04
g_fit_exp(29,29)= 0.222906879748054D-01

c_fit_exp(1,30)= 0.527695217548469D-03
g_fit_exp(1,30)= 0.126124103688047D+07

c_fit_exp(2,30)= 0.835238099206235D-03
g_fit_exp(2,30)= 0.231362133071302D+06

c_fit_exp(3,30)= 0.127487769791547D-02
g_fit_exp(3,30)= 0.648057903459222D+05

c_fit_exp(4,30)= 0.187236827915367D-02
g_fit_exp(4,30)= 0.223340151677907D+05

c_fit_exp(5,30)= 0.266743848855275D-02
g_fit_exp(5,30)= 0.876348900987621D+04

c_fit_exp(6,30)= 0.370916253000094D-02
g_fit_exp(6,30)= 0.376698155651932D+04

c_fit_exp(7,30)= 0.505670692247451D-02
g_fit_exp(7,30)= 0.173432307120325D+04

c_fit_exp(8,30)= 0.678071041503699D-02
g_fit_exp(8,30)= 0.842901157277878D+03

c_fit_exp(9,30)= 0.896480913395097D-02
g_fit_exp(9,30)= 0.428135760032203D+03

c_fit_exp(10,30)= 0.117070551932366D-01
g_fit_exp(10,30)= 0.225627458486109D+03

c_fit_exp(11,30)= 0.151208382598782D-01
g_fit_exp(11,30)= 0.122697847692630D+03

c_fit_exp(12,30)= 0.193345749227226D-01
g_fit_exp(12,30)= 0.685613969889904D+02

c_fit_exp(13,30)= 0.244887962784807D-01
g_fit_exp(13,30)= 0.392339969841834D+02

c_fit_exp(14,30)= 0.307281687085362D-01
g_fit_exp(14,30)= 0.229300939435597D+02

c_fit_exp(15,30)= 0.381841927058923D-01
g_fit_exp(15,30)= 0.136564665106698D+02

c_fit_exp(16,30)= 0.469416252611480D-01
g_fit_exp(16,30)= 0.827273396999604D+01

c_fit_exp(17,30)= 0.569781621637300D-01
g_fit_exp(17,30)= 0.508922960104393D+01

c_fit_exp(18,30)= 0.680638359016814D-01
g_fit_exp(18,30)= 0.317512247604159D+01

c_fit_exp(19,30)= 0.796082074416677D-01
g_fit_exp(19,30)= 0.200662782668857D+01

c_fit_exp(20,30)= 0.904605049209562D-01
g_fit_exp(20,30)= 0.128329657523236D+01

c_fit_exp(21,30)= 0.987211190114169D-01
g_fit_exp(21,30)= 0.829749530674091D+00

c_fit_exp(22,30)= 0.101737234831662D+00
g_fit_exp(22,30)= 0.541962791962346D+00

c_fit_exp(23,30)= 0.966170558467152D-01
g_fit_exp(23,30)= 0.357320680535023D+00

c_fit_exp(24,30)= 0.816344253087800D-01
g_fit_exp(24,30)= 0.237615336068643D+00

c_fit_exp(25,30)= 0.583543858572401D-01
g_fit_exp(25,30)= 0.159234179035390D+00

c_fit_exp(26,30)= 0.328246192394495D-01
g_fit_exp(26,30)= 0.107408175273181D+00

c_fit_exp(27,30)= 0.130676101960358D-01
g_fit_exp(27,30)= 0.727907716315045D-01

c_fit_exp(28,30)= 0.312930130529902D-02
g_fit_exp(28,30)= 0.493920325346856D-01

c_fit_exp(29,30)= 0.341524937013461D-03
g_fit_exp(29,30)= 0.333015150233115D-01

c_fit_exp(30,30)= 0.920393285138635D-05
g_fit_exp(30,30)= 0.218142269248978D-01

end

!! exp(-r**2) is replaced by KT(n,alpha(n),beta(n),r,rc) 
!!                          =exp(-alpha(n) r**2 / ( 1-(r/rc(n))**2)**beta(n))
!!
!! where n is such that exp(-rc(n)**2)=10^(-n)
!!
!! KT is expanded as   KT = \sum_{i=1}^ng  c_fit_GTO(i,ng,n)  exp(-g_fit_GTO(i,ng,n) r**2)
!!                                         

subroutine read_fit_gaussian_rc
include 'j.inc'

alpha( 1)= 0.850321097331778E+00
 beta( 1)= 0.200000000000000E+00
c_fit_GTO( 1, 5, 1)= 0.178910187786054E+02
g_fit_GTO( 1, 5, 1)= 0.225528468879624E+01
c_fit_GTO( 2, 5, 1)= 0.552216725541588E+02
g_fit_GTO( 2, 5, 1)= 0.149475199084317E+01
c_fit_GTO( 3, 5, 1)=-0.547308571782568E+02
g_fit_GTO( 3, 5, 1)= 0.184445259542737E+01
c_fit_GTO( 4, 5, 1)= 0.171954682957566E+01
g_fit_GTO( 4, 5, 1)= 0.786739724810603E+00
c_fit_GTO( 5, 5, 1)=-0.190792327226296E+02
g_fit_GTO( 5, 5, 1)= 0.115652226923382E+01
! chi2=  0.421951997213346E-04

alpha( 2)= 0.960148861176769E+00
 beta( 2)= 0.100000000000000E+00
c_fit_GTO( 1, 3, 2)=-0.573272779817103E+00
g_fit_GTO( 1, 3, 2)= 0.115021106170821E+01
c_fit_GTO( 2, 3, 2)= 0.242421640433683E+01
g_fit_GTO( 2, 3, 2)= 0.922943787034219E+00
c_fit_GTO( 3, 3, 2)=-0.852524729102909E+00
g_fit_GTO( 3, 3, 2)= 0.737596385594286E+00
! chi2=  0.788513570912861E-06

alpha( 2)= 0.960148861176769E+00
 beta( 2)= 0.100000000000000E+00
c_fit_GTO( 1, 4, 2)= 0.106824599412030E+01
g_fit_GTO( 1, 4, 2)= 0.162949611413616E+01
c_fit_GTO( 2, 4, 2)= 0.485498464406728E+01
g_fit_GTO( 2, 4, 2)= 0.106883706783677E+01
c_fit_GTO( 3, 4, 2)=-0.381125918326741E+01
g_fit_GTO( 3, 4, 2)= 0.131951143178821E+01
c_fit_GTO( 4, 4, 2)=-0.111038116258298E+01
g_fit_GTO( 4, 4, 2)= 0.831062631687085E+00
! chi2=  0.564615264834507E-06

alpha( 3)= 0.987315089553559E+00
 beta( 3)= 0.500000000000000E-01
c_fit_GTO( 1, 3, 3)=-0.222851069889805E-01
g_fit_GTO( 1, 3, 3)= 0.184129865854256E+01
c_fit_GTO( 2, 3, 3)= 0.105128960551162E+01
g_fit_GTO( 2, 3, 3)= 0.989533144926947E+00
c_fit_GTO( 3, 3, 3)=-0.293199832312457E-01
g_fit_GTO( 3, 3, 3)= 0.508705723257478E+00
! chi2=  0.169375084483946E-06

alpha( 3)= 0.987315089553559E+00
 beta( 3)= 0.500000000000000E-01
c_fit_GTO( 1, 4, 3)= 0.200704923041787E-01
g_fit_GTO( 1, 4, 3)= 0.248190920286532E+01
c_fit_GTO( 2, 4, 3)= 0.108811444316457E+01
g_fit_GTO( 2, 4, 3)= 0.100705724032504E+01
c_fit_GTO( 3, 4, 3)=-0.855945770166696E-01
g_fit_GTO( 3, 4, 3)= 0.169225759971967E+01
c_fit_GTO( 4, 4, 3)=-0.224217765326068E-01
g_fit_GTO( 4, 4, 3)= 0.496529104909361E+00
! chi2=  0.254393451550153E-06

alpha( 4)= 0.994505374495814E+00
 beta( 4)= 0.300000000000000E-01
c_fit_GTO( 1, 2, 4)= 0.107666273895471E+01
g_fit_GTO( 1, 2, 4)= 0.975787574579645E+00
c_fit_GTO( 2, 2, 4)=-0.765310768413567E-01
g_fit_GTO( 2, 2, 4)= 0.721981009969554E+00
! chi2=  0.268413345158232E-08

alpha( 4)= 0.994505374495814E+00
 beta( 4)= 0.300000000000000E-01
c_fit_GTO( 1, 3, 4)= 0.481100088685904E+00
g_fit_GTO( 1, 3, 4)= 0.990360244312386E+00
c_fit_GTO( 2, 3, 4)= 0.962995298541514E+00
g_fit_GTO( 2, 3, 4)= 0.925030648844638E+00
c_fit_GTO( 3, 3, 4)=-0.444065866002273E+00
g_fit_GTO( 3, 3, 4)= 0.838937032972230E+00
! chi2=  0.169775589889501E-09

alpha( 5)= 0.997139483808004E+00
 beta( 5)= 0.200000000000000E-01
c_fit_GTO( 1, 1, 5)= 0.100100367310234E+01
g_fit_GTO( 1, 1, 5)= 0.100100824325331E+01
! chi2=  0.489739443296088E-07

alpha( 5)= 0.997139483808004E+00
 beta( 5)= 0.200000000000000E-01
c_fit_GTO( 1, 2, 5)= 0.113938760661220E+01
g_fit_GTO( 1, 2, 5)= 0.977363198336558E+00
c_fit_GTO( 2, 2, 5)=-0.139366598470653E+00
g_fit_GTO( 2, 2, 5)= 0.834592057922488E+00
! chi2=  0.685996805041962E-10

alpha( 6)= 0.998826308719630E+00
 beta( 6)= 0.100000000000000E-01
c_fit_GTO( 1, 1, 6)= 0.100040764779999E+01
g_fit_GTO( 1, 1, 6)= 0.100040835298401E+01
! chi2=  0.781940613030215E-08

alpha( 6)= 0.998826308719630E+00
 beta( 6)= 0.100000000000000E-01
c_fit_GTO( 1, 2, 6)= 0.107598170469864E+01
g_fit_GTO( 1, 2, 6)= 0.989017141214555E+00
c_fit_GTO( 2, 2, 6)=-0.759753951293319E-01
g_fit_GTO( 2, 2, 6)= 0.859417248852281E+00
! chi2=  0.762912039488209E-11

alpha( 7)= 0.999502681805158E+00
 beta( 7)= 0.500000000000000E-02
c_fit_GTO( 1, 1, 7)= 0.100017157416398E+01
g_fit_GTO( 1, 1, 7)= 0.100017169373977E+01
! chi2=  0.155471995257752E-08

alpha( 8)= 0.999827424967123E+00
 beta( 8)= 0.200000000000000E-02
c_fit_GTO( 1, 1, 8)= 0.100005925346480E+01
g_fit_GTO( 1, 1, 8)= 0.100005926745424E+01
! chi2=  0.159221715995554E-09

alpha( 9)= 0.999923816628816E+00
 beta( 9)= 0.100000000000000E-02
c_fit_GTO( 1, 1, 9)= 0.100002606389866E+01
g_fit_GTO( 1, 1, 9)= 0.100002606670055E+01
! chi2=  0.304008180211781E-10

do nng=10,16
 alpha(nng)= 0.999993179967811E+00
  beta(nng)= 0.100000000000000E-03
 c_fit_GTO( 1, 1,nng)= 0.100000232683378E+01
 g_fit_GTO( 1, 1,nng)= 0.100000232697940E+01
! chi2 for nng=10:   0.440003463760519E-12
enddo

end

!! exp(-r) is replaced by exp(-alpha r / (1-r/rc)**beta)
!!exp(-alpha r / (1-r/rc)**beta) is expanded as sum_{i=1} c(i) exp(-g(i)*r2)
!! eps  determines rc  as  exp(-rc)=eps

!! exp(-r) is replaced by KT(n,alpha(n),beta(n),r,rc) 
!!                          =exp(-alpha(n) r / (1-r/rc(n))**beta(n))
!!
!! where n is such that exp(-rc(n))=10^(-n)
!!
!! KT is expanded as   KT = \sum_{i=1}^ng  c_fit_STO(i,ng,n)  exp(-g_fit_STO(i,ng,n) r**2)
!!                                   

subroutine read_fit_exp_rc
include 'j.inc'

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1, 1, 2)= 0.626707281001719E+00
g_fit_STO( 1, 1, 2)= 0.392124851518778E+00
! chi2=  0.148098034377895E-02

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1, 1, 3)= 0.613519433389094E+00
g_fit_STO( 1, 1, 3)= 0.381687266255794E+00
! chi2=  0.117140302049119E-02

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1, 1, 4)= 0.608527275008729E+00
g_fit_STO( 1, 1, 4)= 0.378272576247815E+00
! chi2=  0.926878581791671E-03

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1, 1, 5)= 0.605878863353211E+00
g_fit_STO( 1, 1, 5)= 0.376491485914202E+00
! chi2=  0.760989690463565E-03

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1, 1, 6)= 0.604234238278698E+00
g_fit_STO( 1, 1, 6)= 0.375386789288498E+00
! chi2=  0.643942820997324E-03

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1, 1, 7)= 0.603112840240676E+00
g_fit_STO( 1, 1, 7)= 0.374632699348452E+00
! chi2=  0.557557341650868E-03

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1, 1, 8)= 0.602298987157683E+00
g_fit_STO( 1, 1, 8)= 0.374084642196233E+00
! chi2=  0.491373195423264E-03

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1, 1, 9)= 0.601681321801279E+00
g_fit_STO( 1, 1, 9)= 0.373668150321058E+00
! chi2=  0.439119752201238E-03

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1, 1,10)= 0.601196495178396E+00
g_fit_STO( 1, 1,10)= 0.373340854983893E+00
! chi2=  0.396850038160585E-03

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1, 1,11)= 0.600805798247970E+00
g_fit_STO( 1, 1,11)= 0.373076843531410E+00
! chi2=  0.361968167005653E-03

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1, 1,12)= 0.600484236627521E+00
g_fit_STO( 1, 1,12)= 0.372859367043491E+00
! chi2=  0.332701356359145E-03

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1, 1,13)= 0.600214951172029E+00
g_fit_STO( 1, 1,13)= 0.372677111549264E+00
! chi2=  0.307799533547637E-03

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1, 1,14)= 0.599986152982410E+00
g_fit_STO( 1, 1,14)= 0.372522161037776E+00
! chi2=  0.286356761736763E-03

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1, 1,15)= 0.599789352386243E+00
g_fit_STO( 1, 1,15)= 0.372388807282695E+00
! chi2=  0.267700860731357E-03

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1, 1,16)= 0.599618278826312E+00
g_fit_STO( 1, 1,16)= 0.372272830575009E+00
! chi2=  0.251322798448900E-03

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1, 2, 2)= 0.443383935177421E+00
g_fit_STO( 1, 2, 2)= 0.183015746867704E+01
c_fit_STO( 2, 2, 2)= 0.384965802273218E+00
g_fit_STO( 2, 2, 2)= 0.250358802425312E+00
! chi2=  0.756874701312236E-04

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1, 2, 3)= 0.473324708288800E+00
g_fit_STO( 1, 2, 3)= 0.158229672358198E+01
c_fit_STO( 2, 2, 3)= 0.337845841307969E+00
g_fit_STO( 2, 2, 3)= 0.219600849351058E+00
! chi2=  0.795463253038728E-04

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1, 2, 4)= 0.481806722676892E+00
g_fit_STO( 1, 2, 4)= 0.152529651382518E+01
c_fit_STO( 2, 2, 4)= 0.323974204795239E+00
g_fit_STO( 2, 2, 4)= 0.210472574712969E+00
! chi2=  0.678791540236921E-04

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1, 2, 5)= 0.485661667089062E+00
g_fit_STO( 1, 2, 5)= 0.150104331971161E+01
c_fit_STO( 2, 2, 5)= 0.317493589835319E+00
g_fit_STO( 2, 2, 5)= 0.206137534797537E+00
! chi2=  0.574768167725585E-04

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1, 2, 6)= 0.487844402567810E+00
g_fit_STO( 1, 2, 6)= 0.148773036650148E+01
c_fit_STO( 2, 2, 6)= 0.313752450143319E+00
g_fit_STO( 2, 2, 6)= 0.203605998828459E+00
! chi2=  0.494449405157718E-04

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1, 2, 7)= 0.489243534664385E+00
g_fit_STO( 1, 2, 7)= 0.147935007493163E+01
c_fit_STO( 2, 2, 7)= 0.311320181341807E+00
g_fit_STO( 2, 2, 7)= 0.201946796630074E+00
! chi2=  0.432513306188401E-04

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1, 2, 8)= 0.490214896949148E+00
g_fit_STO( 1, 2, 8)= 0.147360083284974E+01
c_fit_STO( 2, 2, 8)= 0.309613279089538E+00
g_fit_STO( 2, 2, 8)= 0.200775512227757E+00
! chi2=  0.383824957409414E-04

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1, 2, 9)= 0.490927860678808E+00
g_fit_STO( 1, 2, 9)= 0.146941633622522E+01
c_fit_STO( 2, 2, 9)= 0.308349805085172E+00
g_fit_STO( 2, 2, 9)= 0.199904601874832E+00
! chi2=  0.344732148782867E-04

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1, 2,10)= 0.491473078169944E+00
g_fit_STO( 1, 2,10)= 0.146623650694576E+01
c_fit_STO( 2, 2,10)= 0.307377030628833E+00
g_fit_STO( 2, 2,10)= 0.199231694344579E+00
! chi2=  0.312731554345796E-04

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1, 2,11)= 0.491903338879166E+00
g_fit_STO( 1, 2,11)= 0.146373944846975E+01
c_fit_STO( 2, 2,11)= 0.306605086475296E+00
g_fit_STO( 2, 2,11)= 0.198696183128133E+00
! chi2=  0.286090987926745E-04

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1, 2,12)= 0.492251441102310E+00
g_fit_STO( 1, 2,12)= 0.146172721510596E+01
c_fit_STO( 2, 2,12)= 0.305977653830637E+00
g_fit_STO( 2, 2,12)= 0.198259897021729E+00
! chi2=  0.263587096211075E-04

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1, 2,13)= 0.492538813132767E+00
g_fit_STO( 1, 2,13)= 0.146007150419224E+01
c_fit_STO( 2, 2,13)= 0.305457664556513E+00
g_fit_STO( 2, 2,13)= 0.197897607718327E+00
! chi2=  0.244336455470883E-04

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1, 2,14)= 0.492780044942384E+00
g_fit_STO( 1, 2,14)= 0.145868559024997E+01
c_fit_STO( 2, 2,14)= 0.305019724167033E+00
g_fit_STO( 2, 2,14)= 0.197591969194289E+00
! chi2=  0.227687396714978E-04

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1, 2,15)= 0.492985405377211E+00
g_fit_STO( 1, 2,15)= 0.145750866874169E+01
c_fit_STO( 2, 2,15)= 0.304645849186284E+00
g_fit_STO( 2, 2,15)= 0.197330662443252E+00
! chi2=  0.213149735846012E-04

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1, 2,16)= 0.493162333070465E+00
g_fit_STO( 1, 2,16)= 0.145649694415707E+01
c_fit_STO( 2, 2,16)= 0.304322954345904E+00
g_fit_STO( 2, 2,16)= 0.197104699272989E+00
! chi2=  0.200348202168244E-04

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1, 3, 2)= 0.238801755987333E+00
g_fit_STO( 1, 3, 2)= 0.744122816315128E+01
c_fit_STO( 2, 3, 2)= 0.366340181767888E+00
g_fit_STO( 2, 3, 2)= 0.100422012473975E+01
c_fit_STO( 3, 3, 2)= 0.309762312038216E+00
g_fit_STO( 3, 3, 2)= 0.219761510287146E+00
! chi2=  0.603668685001517E-05

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1, 3, 3)= 0.283070255046330E+00
g_fit_STO( 1, 3, 3)= 0.537987365002368E+01
c_fit_STO( 2, 3, 3)= 0.390840529449356E+00
g_fit_STO( 2, 3, 3)= 0.733904397675014E+00
c_fit_STO( 3, 3, 3)= 0.223740403319561E+00
g_fit_STO( 3, 3, 3)= 0.171951322590657E+00
! chi2=  0.764459164647062E-05

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1, 3, 4)= 0.295910824362574E+00
g_fit_STO( 1, 3, 4)= 0.497140430022392E+01
c_fit_STO( 2, 3, 4)= 0.395939344704345E+00
g_fit_STO( 2, 3, 4)= 0.677841235346170E+00
c_fit_STO( 3, 3, 4)= 0.200618031114837E+00
g_fit_STO( 3, 3, 4)= 0.158089652125428E+00
! chi2=  0.724266873005843E-05

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1, 3, 5)= 0.301524853363032E+00
g_fit_STO( 1, 3, 5)= 0.481687093225282E+01
c_fit_STO( 2, 3, 5)= 0.397522480022960E+00
g_fit_STO( 2, 3, 5)= 0.656343645099150E+00
c_fit_STO( 3, 3, 5)= 0.191106484421993E+00
g_fit_STO( 3, 3, 5)= 0.152022943847156E+00
! chi2=  0.636558176744508E-05

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1, 3, 6)= 0.304666416292376E+00
g_fit_STO( 1, 3, 6)= 0.473679104590938E+01
c_fit_STO( 2, 3, 6)= 0.398169045896779E+00
g_fit_STO( 2, 3, 6)= 0.645122859460358E+00
c_fit_STO( 3, 3, 6)= 0.186005671169438E+00
g_fit_STO( 3, 3, 6)= 0.148647705445970E+00
! chi2=  0.557589248346339E-05

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1, 3, 7)= 0.306674429215314E+00
g_fit_STO( 1, 3, 7)= 0.468804784139509E+01
c_fit_STO( 2, 3, 7)= 0.398477832946384E+00
g_fit_STO( 2, 3, 7)= 0.638257522473799E+00
c_fit_STO( 3, 3, 7)= 0.182841975190729E+00
g_fit_STO( 3, 3, 7)= 0.146503302579725E+00
! chi2=  0.492906550443179E-05

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1, 3, 8)= 0.308069014513749E+00
g_fit_STO( 1, 3, 8)= 0.465533285390012E+01
c_fit_STO( 2, 3, 8)= 0.398640111566363E+00
g_fit_STO( 2, 3, 8)= 0.633630956620546E+00
c_fit_STO( 3, 3, 8)= 0.180692920530453E+00
g_fit_STO( 3, 3, 8)= 0.145021739165712E+00
! chi2=  0.440431879684204E-05

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1, 3, 9)= 0.309094224896215E+00
g_fit_STO( 1, 3, 9)= 0.463188576399222E+01
c_fit_STO( 2, 3, 9)= 0.398730615776559E+00
g_fit_STO( 2, 3, 9)= 0.630304021993779E+00
c_fit_STO( 3, 3, 9)= 0.179139624403351E+00
g_fit_STO( 3, 3, 9)= 0.143937325531597E+00
! chi2=  0.397484295644819E-05

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1, 3,10)= 0.309879757829769E+00
g_fit_STO( 1, 3,10)= 0.461427031787022E+01
c_fit_STO( 2, 3,10)= 0.398782828297581E+00
g_fit_STO( 2, 3,10)= 0.627797551140031E+00
c_fit_STO( 3, 3,10)= 0.177965265218645E+00
g_fit_STO( 3, 3,10)= 0.143109449133913E+00
! chi2=  0.361875146023471E-05

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1, 3,11)= 0.310500924797851E+00
g_fit_STO( 1, 3,11)= 0.460055859801687E+01
c_fit_STO( 2, 3,11)= 0.398813306579097E+00
g_fit_STO( 2, 3,11)= 0.625841839233613E+00
c_fit_STO( 3, 3,11)= 0.177046613761960E+00
g_fit_STO( 3, 3,11)= 0.142456802150364E+00
! chi2=  0.331958274571580E-05

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1, 3,12)= 0.311004461337388E+00
g_fit_STO( 1, 3,12)= 0.458958660769053E+01
c_fit_STO( 2, 3,12)= 0.398830869388092E+00
g_fit_STO( 2, 3,12)= 0.624273593671839E+00
c_fit_STO( 3, 3,12)= 0.176308544292532E+00
g_fit_STO( 3, 3,12)= 0.141929125733873E+00
! chi2=  0.306513422987578E-05

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1, 3,13)= 0.311420914198161E+00
g_fit_STO( 1, 3,13)= 0.458061060067582E+01
c_fit_STO( 2, 3,13)= 0.398840491899702E+00
g_fit_STO( 2, 3,13)= 0.622988208967199E+00
c_fit_STO( 3, 3,13)= 0.175702685919615E+00
g_fit_STO( 3, 3,13)= 0.141493690989748E+00
! chi2=  0.284631268756308E-05

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1, 3,14)= 0.311771081795480E+00
g_fit_STO( 1, 3,14)= 0.457313336478921E+01
c_fit_STO( 2, 3,14)= 0.398845113964235E+00
g_fit_STO( 2, 3,14)= 0.621915611892196E+00
c_fit_STO( 3, 3,14)= 0.175196505149857E+00
g_fit_STO( 3, 3,14)= 0.141128274685919E+00
! chi2=  0.265626146917102E-05

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1, 3,15)= 0.312069637855997E+00
g_fit_STO( 1, 3,15)= 0.456680994294180E+01
c_fit_STO( 2, 3,15)= 0.398846532792762E+00
g_fit_STO( 2, 3,15)= 0.621007082683421E+00
c_fit_STO( 3, 3,15)= 0.174767318041432E+00
g_fit_STO( 3, 3,15)= 0.140817258027741E+00
! chi2=  0.248973898903506E-05

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1, 3,16)= 0.312327211565973E+00
g_fit_STO( 1, 3,16)= 0.456139375597114E+01
c_fit_STO( 2, 3,16)= 0.398845887249489E+00
g_fit_STO( 2, 3,16)= 0.620227728613800E+00
c_fit_STO( 3, 3,16)= 0.174398833327414E+00
g_fit_STO( 3, 3,16)= 0.140549344318740E+00
! chi2=  0.234268250639213E-05

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1, 4, 2)= 0.124821976621491E+00
g_fit_STO( 1, 4, 2)= 0.283176128713753E+02
c_fit_STO( 2, 4, 2)= 0.217225058192056E+00
g_fit_STO( 2, 4, 2)= 0.382227805770788E+01
c_fit_STO( 3, 4, 2)= 0.325301267172421E+00
g_fit_STO( 3, 4, 2)= 0.833945639728744E+00
c_fit_STO( 4, 4, 2)= 0.289030251856363E+00
g_fit_STO( 4, 4, 2)= 0.212651268829803E+00
! chi2=  0.162056564273535E-05

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1, 4, 3)= 0.166955123687561E+00
g_fit_STO( 1, 4, 3)= 0.163215386910652E+02
c_fit_STO( 2, 4, 3)= 0.275302277201921E+00
g_fit_STO( 2, 4, 3)= 0.222010285633075E+01
c_fit_STO( 3, 4, 3)= 0.328593584127601E+00
g_fit_STO( 3, 4, 3)= 0.509346397191463E+00
c_fit_STO( 4, 4, 3)= 0.170392376737730E+00
g_fit_STO( 4, 4, 3)= 0.152284992903861E+00
! chi2=  0.899174885594907E-06

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1, 4, 4)= 0.180294218105403E+00
g_fit_STO( 1, 4, 4)= 0.142026313956870E+02
c_fit_STO( 2, 4, 4)= 0.292388512870731E+00
g_fit_STO( 2, 4, 4)= 0.193200264896463E+01
c_fit_STO( 3, 4, 4)= 0.326545736613179E+00
g_fit_STO( 3, 4, 4)= 0.443583444747303E+00
c_fit_STO( 4, 4, 4)= 0.137155654505549E+00
g_fit_STO( 4, 4, 4)= 0.133319231654188E+00
! chi2=  0.972962916910637E-06

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1, 4, 5)= 0.185750075874510E+00
g_fit_STO( 1, 4, 5)= 0.134962752547686E+02
c_fit_STO( 2, 4, 5)= 0.298974500171981E+00
g_fit_STO( 2, 4, 5)= 0.183572392322593E+01
c_fit_STO( 3, 4, 5)= 0.324649985195423E+00
g_fit_STO( 3, 4, 5)= 0.420994374529283E+00
c_fit_STO( 4, 4, 5)= 0.125005133725648E+00
g_fit_STO( 4, 4, 5)= 0.125519947210107E+00
! chi2=  0.896192564470062E-06

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1, 4, 6)= 0.188684972162474E+00
g_fit_STO( 1, 4, 6)= 0.131531204800517E+02
c_fit_STO( 2, 4, 6)= 0.302380126919540E+00
g_fit_STO( 2, 4, 6)= 0.178890438065433E+01
c_fit_STO( 3, 4, 6)= 0.323247557781601E+00
g_fit_STO( 3, 4, 6)= 0.409885800931790E+00
c_fit_STO( 4, 4, 6)= 0.118983256801475E+00
g_fit_STO( 4, 4, 6)= 0.121376739950117E+00
! chi2=  0.801458712161773E-06

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1, 4, 7)= 0.190519774143868E+00
g_fit_STO( 1, 4, 7)= 0.129516579791826E+02
c_fit_STO( 2, 4, 7)= 0.304451183510694E+00
g_fit_STO( 2, 4, 7)= 0.176139673169500E+01
c_fit_STO( 3, 4, 7)= 0.322213963656911E+00
g_fit_STO( 3, 4, 7)= 0.403312421285979E+00
c_fit_STO( 4, 4, 7)= 0.115431171854294E+00
g_fit_STO( 4, 4, 7)= 0.118822398640250E+00
! chi2=  0.716581776558691E-06

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1, 4, 8)= 0.191776708803555E+00
g_fit_STO( 1, 4, 8)= 0.128194839386571E+02
c_fit_STO( 2, 4, 8)= 0.305841333890574E+00
g_fit_STO( 2, 4, 8)= 0.174333790861836E+01
c_fit_STO( 3, 4, 8)= 0.321431841498837E+00
g_fit_STO( 3, 4, 8)= 0.398973858934155E+00
c_fit_STO( 4, 4, 8)= 0.113099425932816E+00
g_fit_STO( 4, 4, 8)= 0.117093675839843E+00
! chi2=  0.644864084824217E-06

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1, 4, 9)= 0.192692305537813E+00
g_fit_STO( 1, 4, 9)= 0.127262167856831E+02
c_fit_STO( 2, 4, 9)= 0.306838233547577E+00
g_fit_STO( 2, 4, 9)= 0.173058707791354E+01
c_fit_STO( 3, 4, 9)= 0.320823064709236E+00
g_fit_STO( 3, 4, 9)= 0.395897461351293E+00
c_fit_STO( 4, 4, 9)= 0.111455154681435E+00
g_fit_STO( 4, 4, 9)= 0.115847119107583E+00
! chi2=  0.584813571904827E-06

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1, 4,10)= 0.193389313788522E+00
g_fit_STO( 1, 4,10)= 0.126569357200659E+02
c_fit_STO( 2, 4,10)= 0.307587784958877E+00
g_fit_STO( 2, 4,10)= 0.172110973395680E+01
c_fit_STO( 3, 4,10)= 0.320337229595704E+00
g_fit_STO( 3, 4,10)= 0.393602821025969E+00
c_fit_STO( 4, 4,10)= 0.110234920824851E+00
g_fit_STO( 4, 4,10)= 0.114906126498137E+00
! chi2=  0.534299626844052E-06

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1, 4,11)= 0.193937836004044E+00
g_fit_STO( 1, 4,11)= 0.126034716670167E+02
c_fit_STO( 2, 4,11)= 0.308171753445773E+00
g_fit_STO( 2, 4,11)= 0.171379161766525E+01
c_fit_STO( 3, 4,11)= 0.319941188335492E+00
g_fit_STO( 3, 4,11)= 0.391825711149986E+00
c_fit_STO( 4, 4,11)= 0.109294133861597E+00
g_fit_STO( 4, 4,11)= 0.114170803273581E+00
! chi2=  0.491437691144335E-06

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1, 4,12)= 0.194380844868826E+00
g_fit_STO( 1, 4,12)= 0.125609823231919E+02
c_fit_STO( 2, 4,12)= 0.308639486334772E+00
g_fit_STO( 2, 4,12)= 0.170797199357400E+01
c_fit_STO( 3, 4,12)= 0.319612501421180E+00
g_fit_STO( 3, 4,12)= 0.390408874553547E+00
c_fit_STO( 4, 4,12)= 0.108547021708312E+00
g_fit_STO( 4, 4,12)= 0.113580458804870E+00
! chi2=  0.454719218118184E-06

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1, 4,13)= 0.194746161602492E+00
g_fit_STO( 1, 4,13)= 0.125264175670549E+02
c_fit_STO( 2, 4,13)= 0.309022510792885E+00
g_fit_STO( 2, 4,13)= 0.170323454244523E+01
c_fit_STO( 3, 4,13)= 0.319335535008801E+00
g_fit_STO( 3, 4,13)= 0.389252907018447E+00
c_fit_STO( 4, 4,13)= 0.107939564578150E+00
g_fit_STO( 4, 4,13)= 0.113096118237071E+00
! chi2=  0.422969019898140E-06

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1, 4,14)= 0.195052605044258E+00
g_fit_STO( 1, 4,14)= 0.124977606096679E+02
c_fit_STO( 2, 4,14)= 0.309341910434070E+00
g_fit_STO( 2, 4,14)= 0.169930391958081E+01
c_fit_STO( 3, 4,14)= 0.319099093265678E+00
g_fit_STO( 3, 4,14)= 0.388291884521225E+00
c_fit_STO( 4, 4,14)= 0.107436073025065E+00
g_fit_STO( 4, 4,14)= 0.112691614347247E+00
! chi2=  0.395275140478609E-06

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1, 4,15)= 0.195313360548483E+00
g_fit_STO( 1, 4,15)= 0.124736256404293E+02
c_fit_STO( 2, 4,15)= 0.309612313300908E+00
g_fit_STO( 2, 4,15)= 0.169599086534319E+01
c_fit_STO( 3, 4,15)= 0.318894971425306E+00
g_fit_STO( 3, 4,15)= 0.387480377618139E+00
c_fit_STO( 4, 4,15)= 0.107012040681432E+00
g_fit_STO( 4, 4,15)= 0.112348732175302E+00
! chi2=  0.370926349057960E-06

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1, 4,16)= 0.195537945014063E+00
g_fit_STO( 1, 4,16)= 0.124530290426175E+02
c_fit_STO( 2, 4,16)= 0.309844182899421E+00
g_fit_STO( 2, 4,16)= 0.169316102023811E+01
c_fit_STO( 3, 4,16)= 0.318717026301581E+00
g_fit_STO( 3, 4,16)= 0.386786064690199E+00
c_fit_STO( 4, 4,16)= 0.106650091956141E+00
g_fit_STO( 4, 4,16)= 0.112054404587003E+00
! chi2=  0.349363221174002E-06

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1, 5, 2)= 0.132347409623488E+00
g_fit_STO( 1, 5, 2)= 0.251190624427665E+02
c_fit_STO( 2, 5, 2)= 0.226691922016015E+00
g_fit_STO( 2, 5, 2)= 0.341004378340174E+01
c_fit_STO( 3, 5, 2)= 0.321136035017096E+00
g_fit_STO( 3, 5, 2)= 0.766302063011331E+00
c_fit_STO( 4, 5, 2)= 0.273501541403149E+00
g_fit_STO( 4, 5, 2)= 0.206215935333223E+00
c_fit_STO( 5, 5, 2)=-0.396089542697085E-06
g_fit_STO( 5, 5, 2)=-0.426362997722474E+00
! chi2=  0.504528498500195E-06

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1, 5, 3)= 0.999460946953543E-01
g_fit_STO( 1, 5, 3)= 0.463276937833796E+02
c_fit_STO( 2, 5, 3)= 0.174165165552970E+00
g_fit_STO( 2, 5, 3)= 0.629866047863685E+01
c_fit_STO( 3, 5, 3)= 0.260034039554611E+00
g_fit_STO( 3, 5, 3)= 0.143896396210423E+01
c_fit_STO( 4, 5, 3)= 0.284938583818494E+00
g_fit_STO( 4, 5, 3)= 0.421365676744027E+00
c_fit_STO( 5, 5, 3)= 0.146041390307174E+00
g_fit_STO( 5, 5, 3)= 0.143946207446028E+00
! chi2=  0.128856197698321E-06

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1, 5, 4)= 0.112600267499341E+00
g_fit_STO( 1, 5, 4)= 0.371286801020721E+02
c_fit_STO( 2, 5, 4)= 0.194447107408543E+00
g_fit_STO( 2, 5, 4)= 0.504849714385687E+01
c_fit_STO( 3, 5, 4)= 0.279789944871692E+00
g_fit_STO( 3, 5, 4)= 0.115465512121670E+01
c_fit_STO( 4, 5, 4)= 0.271738938081820E+00
g_fit_STO( 4, 5, 4)= 0.340545212416462E+00
c_fit_STO( 5, 5, 4)= 0.102078829632201E+00
g_fit_STO( 5, 5, 4)= 0.119767859516558E+00
! chi2=  0.152566244956621E-06

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1, 5, 5)= 0.117522098901339E+00
g_fit_STO( 1, 5, 5)= 0.344191321128711E+02
c_fit_STO( 2, 5, 5)= 0.202143394161674E+00
g_fit_STO( 2, 5, 5)= 0.467997045304604E+01
c_fit_STO( 3, 5, 5)= 0.286458662096449E+00
g_fit_STO( 3, 5, 5)= 0.107015398894723E+01
c_fit_STO( 4, 5, 5)= 0.265468397287967E+00
g_fit_STO( 4, 5, 5)= 0.315176054709193E+00
c_fit_STO( 5, 5, 5)= 0.873172461946644E-01
g_fit_STO( 5, 5, 5)= 0.110069883398980E+00
! chi2=  0.149176488439571E-06

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1, 5, 6)= 0.120047486186406E+00
g_fit_STO( 1, 5, 6)= 0.331953332495658E+02
c_fit_STO( 2, 5, 6)= 0.206031833229635E+00
g_fit_STO( 2, 5, 6)= 0.451347969023005E+01
c_fit_STO( 3, 5, 6)= 0.289542107637711E+00
g_fit_STO( 3, 5, 6)= 0.103189685217372E+01
c_fit_STO( 4, 5, 6)= 0.261806135709943E+00
g_fit_STO( 4, 5, 6)= 0.303490215064788E+00
c_fit_STO( 5, 5, 6)= 0.805848956900302E-01
g_fit_STO( 5, 5, 6)= 0.105130836644248E+00
! chi2=  0.136711193035428E-06

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1, 5, 7)= 0.121583209270562E+00
g_fit_STO( 1, 5, 7)= 0.325050288433055E+02
c_fit_STO( 2, 5, 7)= 0.208372018499221E+00
g_fit_STO( 2, 5, 7)= 0.441954457489273E+01
c_fit_STO( 3, 5, 7)= 0.291281843346844E+00
g_fit_STO( 3, 5, 7)= 0.101028421292929E+01
c_fit_STO( 4, 5, 7)= 0.259405549081031E+00
g_fit_STO( 4, 5, 7)= 0.296826369313149E+00
c_fit_STO( 5, 5, 7)= 0.768235539161970E-01
g_fit_STO( 5, 5, 7)= 0.102174357223892E+00
! chi2=  0.123782350723932E-06

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1, 5, 8)= 0.122617246794876E+00
g_fit_STO( 1, 5, 8)= 0.320631432277314E+02
c_fit_STO( 2, 5, 8)= 0.209935955502383E+00
g_fit_STO( 2, 5, 8)= 0.435939547724988E+01
c_fit_STO( 3, 5, 8)= 0.292389242952094E+00
g_fit_STO( 3, 5, 8)= 0.996430784365133E+00
c_fit_STO( 4, 5, 8)= 0.257712438663902E+00
g_fit_STO( 4, 5, 8)= 0.292526855173178E+00
c_fit_STO( 5, 5, 8)= 0.744432954058538E-01
g_fit_STO( 5, 5, 8)= 0.100213923162234E+00
! chi2=  0.112239413070087E-06

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1, 5, 9)= 0.123361776773075E+00
g_fit_STO( 1, 5, 9)= 0.317564532314122E+02
c_fit_STO( 2, 5, 9)= 0.211055661532014E+00
g_fit_STO( 2, 5, 9)= 0.431763300761505E+01
c_fit_STO( 3, 5, 9)= 0.293152517931761E+00
g_fit_STO( 3, 5, 9)= 0.986803390089920E+00
c_fit_STO( 4, 5, 9)= 0.256455215067153E+00
g_fit_STO( 4, 5, 9)= 0.289523785589060E+00
c_fit_STO( 5, 5, 9)= 0.728080588819447E-01
g_fit_STO( 5, 5, 9)= 0.988209364488528E-01
! chi2=  0.102299579006731E-06

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1, 5,10)= 0.123923888120342E+00
g_fit_STO( 1, 5,10)= 0.315313361766915E+02
c_fit_STO( 2, 5,10)= 0.211897281831277E+00
g_fit_STO( 2, 5,10)= 0.428696373351694E+01
c_fit_STO( 3, 5,10)= 0.293709000355430E+00
g_fit_STO( 3, 5,10)= 0.979727436091827E+00
c_fit_STO( 4, 5,10)= 0.255485234057325E+00
g_fit_STO( 4, 5,10)= 0.287307547532314E+00
c_fit_STO( 5, 5,10)= 0.716178960356681E-01
g_fit_STO( 5, 5,10)= 0.977809454361573E-01
! chi2=  0.937972961772285E-07

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1, 5,11)= 0.124363510335233E+00
g_fit_STO( 1, 5,11)= 0.313591912302193E+02
c_fit_STO( 2, 5,11)= 0.212553173909835E+00
g_fit_STO( 2, 5,11)= 0.426349713173508E+01
c_fit_STO( 3, 5,11)= 0.294131992177262E+00
g_fit_STO( 3, 5,11)= 0.974308978081455E+00
c_fit_STO( 4, 5,11)= 0.254714434432686E+00
g_fit_STO( 4, 5,11)= 0.285604635162692E+00
c_fit_STO( 5, 5,11)= 0.707139780197426E-01
g_fit_STO( 5, 5,11)= 0.969751858844505E-01
! chi2=  0.865030573575618E-07

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1, 5,12)= 0.124716865989714E+00
g_fit_STO( 1, 5,12)= 0.312233734080808E+02
c_fit_STO( 2, 5,12)= 0.213078820302608E+00
g_fit_STO( 2, 5,12)= 0.424496912821708E+01
c_fit_STO( 3, 5,12)= 0.294464008897023E+00
g_fit_STO( 3, 5,12)= 0.970027601074418E+00
c_fit_STO( 4, 5,12)= 0.254087344150899E+00
g_fit_STO( 4, 5,12)= 0.284255156749901E+00
c_fit_STO( 5, 5,12)= 0.700046747600446E-01
g_fit_STO( 5, 5,12)= 0.963326902676106E-01
! chi2=  0.802056016045019E-07

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1, 5,13)= 0.125007128117468E+00
g_fit_STO( 1, 5,13)= 0.311135528081486E+02
c_fit_STO( 2, 5,13)= 0.213509571090334E+00
g_fit_STO( 2, 5,13)= 0.422997447993210E+01
c_fit_STO( 3, 5,13)= 0.294731357171740E+00
g_fit_STO( 3, 5,13)= 0.966560058630192E+00
c_fit_STO( 4, 5,13)= 0.253567321236885E+00
g_fit_STO( 4, 5,13)= 0.283159407422639E+00
c_fit_STO( 5, 5,13)= 0.694335423409852E-01
g_fit_STO( 5, 5,13)= 0.958084843070333E-01
! chi2=  0.747288425738737E-07

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1, 5,14)= 0.125249837453403E+00
g_fit_STO( 1, 5,14)= 0.310229789792823E+02
c_fit_STO( 2, 5,14)= 0.213869025089058E+00
g_fit_STO( 2, 5,14)= 0.421759481838875E+01
c_fit_STO( 3, 5,14)= 0.294951143035504E+00
g_fit_STO( 3, 5,14)= 0.963695013939300E+00
c_fit_STO( 4, 5,14)= 0.253129194977592E+00
g_fit_STO( 4, 5,14)= 0.282251995395717E+00
c_fit_STO( 5, 5,14)= 0.689639731121195E-01
g_fit_STO( 5, 5,14)= 0.953727056764064E-01
! chi2=  0.699306341636406E-07

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1, 5,15)= 0.125455805216971E+00
g_fit_STO( 1, 5,15)= 0.309470548579346E+02
c_fit_STO( 2, 5,15)= 0.214173546614723E+00
g_fit_STO( 2, 5,15)= 0.420720468851041E+01
c_fit_STO( 3, 5,15)= 0.295134963743015E+00
g_fit_STO( 3, 5,15)= 0.961288458855089E+00
c_fit_STO( 4, 5,15)= 0.252755102745641E+00
g_fit_STO( 4, 5,15)= 0.281488233252035E+00
c_fit_STO( 5, 5,15)= 0.685712030376254E-01
g_fit_STO( 5, 5,15)= 0.950047527683797E-01
! chi2=  0.656972324526213E-07

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1, 5,16)= 0.125632790193197E+00
g_fit_STO( 1, 5,16)= 0.308825416756871E+02
c_fit_STO( 2, 5,16)= 0.214434832469111E+00
g_fit_STO( 2, 5,16)= 0.419836359385679E+01
c_fit_STO( 3, 5,16)= 0.295290941076583E+00
g_fit_STO( 3, 5,16)= 0.959238970173204E+00
c_fit_STO( 4, 5,16)= 0.252432032921804E+00
g_fit_STO( 4, 5,16)= 0.280836564232603E+00
c_fit_STO( 5, 5,16)= 0.682379012495396E-01
g_fit_STO( 5, 5,16)= 0.946899610496582E-01
! chi2=  0.619375499680025E-07

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1, 6, 2)= 0.763999618395128E-01
g_fit_STO( 1, 6, 2)= 0.761625752142724E+02
c_fit_STO( 2, 6, 2)= 0.134988985266797E+00
g_fit_STO( 2, 6, 2)= 0.103510814680310E+02
c_fit_STO( 3, 6, 2)= 0.214176206722078E+00
g_fit_STO( 3, 6, 2)= 0.235408744510055E+01
c_fit_STO( 4, 6, 2)= 0.292338168089450E+00
g_fit_STO( 4, 6, 2)= 0.667259525255384E+00
c_fit_STO( 5, 6, 2)= 0.255497762434741E+00
g_fit_STO( 5, 6, 2)= 0.199769048789301E+00
c_fit_STO( 6, 6, 2)=-0.464522067199094E-05
g_fit_STO( 6, 6, 2)=-0.311019563196220E+00
! chi2=  0.621619900147082E-07

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1, 6, 3)= 0.110721642005023E+00
g_fit_STO( 1, 6, 3)= 0.367530903267014E+02
c_fit_STO( 2, 6, 3)= 0.182469643077358E+00
g_fit_STO( 2, 6, 3)= 0.524152816668969E+01
c_fit_STO( 3, 6, 3)= 0.262864884839872E+00
g_fit_STO( 3, 6, 3)= 0.127142908386752E+01
c_fit_STO( 4, 6, 3)= 0.809791559652845E+01
g_fit_STO( 4, 6, 3)= 0.457601086364510E+00
c_fit_STO( 5, 6, 3)= 0.143750147796058E+00
g_fit_STO( 5, 6, 3)= 0.143118706404996E+00
c_fit_STO( 6, 6, 3)=-0.783715642172542E+01
g_fit_STO( 6, 6, 3)= 0.459822409755285E+00
! chi2=  0.147560278901019E-06

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1, 6, 4)= 0.720962530538626E-01
g_fit_STO( 1, 6, 4)= 0.912159598330879E+02
c_fit_STO( 2, 6, 4)= 0.127380473796650E+00
g_fit_STO( 2, 6, 4)= 0.124017113313697E+02
c_fit_STO( 3, 6, 4)= 0.200762429552592E+00
g_fit_STO( 3, 6, 4)= 0.283398706021713E+01
c_fit_STO( 4, 6, 4)= 0.261831078648656E+00
g_fit_STO( 4, 6, 4)= 0.831417365344359E+00
c_fit_STO( 5, 6, 4)= 0.230639451547105E+00
g_fit_STO( 5, 6, 4)= 0.286342266423976E+00
c_fit_STO( 6, 6, 4)= 0.821887305667896E-01
g_fit_STO( 6, 6, 4)= 0.112040913389373E+00
! chi2=  0.267061394430903E-07

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1, 6, 5)= 0.764847945372521E-01
g_fit_STO( 1, 6, 5)= 0.818930729794197E+02
c_fit_STO( 2, 6, 5)= 0.134839489988886E+00
g_fit_STO( 2, 6, 5)= 0.111340444633476E+02
c_fit_STO( 3, 6, 5)= 0.210716936140758E+00
g_fit_STO( 3, 6, 5)= 0.254422975792375E+01
c_fit_STO( 4, 6, 5)= 0.267600312089050E+00
g_fit_STO( 4, 6, 5)= 0.746275380196253E+00
c_fit_STO( 5, 6, 5)= 0.219104959580674E+00
g_fit_STO( 5, 6, 5)= 0.256805582394793E+00
c_fit_STO( 6, 6, 5)= 0.646152575814164E-01
g_fit_STO( 6, 6, 5)= 0.100257482245589E+00
! chi2=  0.281055659107884E-07

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1, 6, 6)= 0.786274220259672E-01
g_fit_STO( 1, 6, 6)= 0.780087382929410E+02
c_fit_STO( 2, 6, 6)= 0.138455057050563E+00
g_fit_STO( 2, 6, 6)= 0.106057776009500E+02
c_fit_STO( 3, 6, 6)= 0.215394253208816E+00
g_fit_STO( 3, 6, 6)= 0.242341070668619E+01
c_fit_STO( 4, 6, 6)= 0.269813150943781E+00
g_fit_STO( 4, 6, 6)= 0.710616067671822E+00
c_fit_STO( 5, 6, 6)= 0.213083127832209E+00
g_fit_STO( 5, 6, 6)= 0.244087827292023E+00
c_fit_STO( 6, 6, 6)= 0.572380343540142E-01
g_fit_STO( 6, 6, 6)= 0.944555983058103E-01
! chi2=  0.265266603953636E-07

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1, 6, 7)= 0.798906008098516E-01
g_fit_STO( 1, 6, 7)= 0.759130208681763E+02
c_fit_STO( 2, 6, 7)= 0.140576645485091E+00
g_fit_STO( 2, 6, 7)= 0.103206816470674E+02
c_fit_STO( 3, 6, 7)= 0.218081993470394E+00
g_fit_STO( 3, 6, 7)= 0.235818499544008E+01
c_fit_STO( 4, 6, 7)= 0.270887052829049E+00
g_fit_STO( 4, 6, 7)= 0.691325901071497E+00
c_fit_STO( 5, 6, 7)= 0.209381323213018E+00
g_fit_STO( 5, 6, 7)= 0.237122461195813E+00
c_fit_STO( 6, 6, 7)= 0.533509133992430E-01
g_fit_STO( 6, 6, 7)= 0.910799975688582E-01
! chi2=  0.243625551832052E-07

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1, 6, 8)= 0.807248585721399E-01
g_fit_STO( 1, 6, 8)= 0.746069946717378E+02
c_fit_STO( 2, 6, 8)= 0.141973215249485E+00
g_fit_STO( 2, 6, 8)= 0.101429378532731E+02
c_fit_STO( 3, 6, 8)= 0.219824908908628E+00
g_fit_STO( 3, 6, 8)= 0.231750604128397E+01
c_fit_STO( 4, 6, 8)= 0.271492093743700E+00
g_fit_STO( 4, 6, 8)= 0.679277393149327E+00
c_fit_STO( 5, 6, 8)= 0.206874872523457E+00
g_fit_STO( 5, 6, 8)= 0.232737549470056E+00
c_fit_STO( 6, 6, 8)= 0.509865476370802E-01
g_fit_STO( 6, 6, 8)= 0.888866339923007E-01
! chi2=  0.222724080567079E-07

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1, 6, 9)= 0.813178620297484E-01
g_fit_STO( 1, 6, 9)= 0.737165431801422E+02
c_fit_STO( 2, 6, 9)= 0.142963484196151E+00
g_fit_STO( 2, 6, 9)= 0.100216794107542E+02
c_fit_STO( 3, 6, 9)= 0.221046906895077E+00
g_fit_STO( 3, 6, 9)= 0.228974389737265E+01
c_fit_STO( 4, 6, 9)= 0.271868580970644E+00
g_fit_STO( 4, 6, 9)= 0.671044283492284E+00
c_fit_STO( 5, 6, 9)= 0.205065436476186E+00
g_fit_STO( 5, 6, 9)= 0.229723720534680E+00
c_fit_STO( 6, 6, 9)= 0.494069336689969E-01
g_fit_STO( 6, 6, 9)= 0.873508529400038E-01
! chi2=  0.204074844053383E-07

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1, 6,10)= 0.817615000636250E-01
g_fit_STO( 1, 6,10)= 0.730712568355430E+02
c_fit_STO( 2, 6,10)= 0.143702931872240E+00
g_fit_STO( 2, 6,10)= 0.993373532871114E+01
c_fit_STO( 3, 6,10)= 0.221951416290502E+00
g_fit_STO( 3, 6,10)= 0.226960012972878E+01
c_fit_STO( 4, 6,10)= 0.272119966876227E+00
g_fit_STO( 4, 6,10)= 0.665063742955952E+00
c_fit_STO( 5, 6,10)= 0.203697965457864E+00
g_fit_STO( 5, 6,10)= 0.227524364767539E+00
c_fit_STO( 6, 6,10)= 0.482806649369607E-01
g_fit_STO( 6, 6,10)= 0.862167686140238E-01
! chi2=  0.187804084064214E-07

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1, 6,11)= 0.821061174246275E-01
g_fit_STO( 1, 6,11)= 0.725826161440889E+02
c_fit_STO( 2, 6,11)= 0.144276473921535E+00
g_fit_STO( 2, 6,11)= 0.986707195581064E+01
c_fit_STO( 3, 6,11)= 0.222648101503509E+00
g_fit_STO( 3, 6,11)= 0.225432286830450E+01
c_fit_STO( 4, 6,11)= 0.272296906496774E+00
g_fit_STO( 4, 6,11)= 0.660523295491861E+00
c_fit_STO( 5, 6,11)= 0.202628327032600E+00
g_fit_STO( 5, 6,11)= 0.225848259419457E+00
c_fit_STO( 6, 6,11)= 0.474386351153201E-01
g_fit_STO( 6, 6,11)= 0.853455196612103E-01
! chi2=  0.173670242229943E-07

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1, 6,12)= 0.823816373414918E-01
g_fit_STO( 1, 6,12)= 0.722001705596527E+02
c_fit_STO( 2, 6,12)= 0.144734487101466E+00
g_fit_STO( 2, 6,12)= 0.981482709436697E+01
c_fit_STO( 3, 6,12)= 0.223201285338650E+00
g_fit_STO( 3, 6,12)= 0.224234269557282E+01
c_fit_STO( 4, 6,12)= 0.272426656460443E+00
g_fit_STO( 4, 6,12)= 0.656959232429613E+00
c_fit_STO( 5, 6,12)= 0.201768944697846E+00
g_fit_STO( 5, 6,12)= 0.224528329023910E+00
c_fit_STO( 6, 6,12)= 0.467860446822506E-01
g_fit_STO( 6, 6,12)= 0.846554507869780E-01
! chi2=  0.161363731630005E-07

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1, 6,13)= 0.826069986231974E-01
g_fit_STO( 1, 6,13)= 0.718930474300161E+02
c_fit_STO( 2, 6,13)= 0.145108750219382E+00
g_fit_STO( 2, 6,13)= 0.977280447558234E+01
c_fit_STO( 3, 6,13)= 0.223651213000484E+00
g_fit_STO( 3, 6,13)= 0.223269967753217E+01
c_fit_STO( 4, 6,13)= 0.272524972146706E+00
g_fit_STO( 4, 6,13)= 0.654087688546906E+00
c_fit_STO( 5, 6,13)= 0.201063505196432E+00
g_fit_STO( 5, 6,13)= 0.223461889951547E+00
c_fit_STO( 6, 6,13)= 0.462658508171507E-01
g_fit_STO( 6, 6,13)= 0.840955045163437E-01
! chi2=  0.150595175583497E-07

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1, 6,14)= 0.827947511157433E-01
g_fit_STO( 1, 6,14)= 0.716413217856378E+02
c_fit_STO( 2, 6,14)= 0.145420330733573E+00
g_fit_STO( 2, 6,14)= 0.973829515692097E+01
c_fit_STO( 3, 6,14)= 0.224024368748319E+00
g_fit_STO( 3, 6,14)= 0.222477391294886E+01
c_fit_STO( 4, 6,14)= 0.272601539589306E+00
g_fit_STO( 4, 6,14)= 0.651725130516348E+00
c_fit_STO( 5, 6,14)= 0.200474157163108E+00
g_fit_STO( 5, 6,14)= 0.222582281284789E+00
c_fit_STO( 6, 6,14)= 0.458417003809316E-01
g_fit_STO( 6, 6,14)= 0.836321130036139E-01
! chi2=  0.141117127919271E-07

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1, 6,15)= 0.829535783929852E-01
g_fit_STO( 1, 6,15)= 0.714315477825328E+02
c_fit_STO( 2, 6,15)= 0.145683772893193E+00
g_fit_STO( 2, 6,15)= 0.970946871721093E+01
c_fit_STO( 3, 6,15)= 0.224338831549632E+00
g_fit_STO( 3, 6,15)= 0.221814702132134E+01
c_fit_STO( 4, 6,15)= 0.272662518356981E+00
g_fit_STO( 4, 6,15)= 0.649747811857104E+00
c_fit_STO( 5, 6,15)= 0.199974564397384E+00
g_fit_STO( 5, 6,15)= 0.221844425722587E+00
c_fit_STO( 6, 6,15)= 0.454894136910675E-01
g_fit_STO( 6, 6,15)= 0.832423435617031E-01
! chi2=  0.132724674138559E-07

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1, 6,16)= 0.830896974792764E-01
g_fit_STO( 1, 6,16)= 0.712542884075179E+02
c_fit_STO( 2, 6,16)= 0.145909419036777E+00
g_fit_STO( 2, 6,16)= 0.968504477693609E+01
c_fit_STO( 3, 6,16)= 0.224607454403081E+00
g_fit_STO( 3, 6,16)= 0.221252613915335E+01
c_fit_STO( 4, 6,16)= 0.272712030196756E+00
g_fit_STO( 4, 6,16)= 0.648068967637498E+00
c_fit_STO( 5, 6,16)= 0.199545757631710E+00
g_fit_STO( 5, 6,16)= 0.221216643594474E+00
c_fit_STO( 6, 6,16)= 0.451922553093808E-01
g_fit_STO( 6, 6,16)= 0.829099754084022E-01
! chi2=  0.125249964026241E-07

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1, 7, 2)= 0.458885424167736E-01
g_fit_STO( 1, 7, 2)= 0.211809589701046E+03
c_fit_STO( 2, 7, 2)= 0.818886717970142E-01
g_fit_STO( 2, 7, 2)= 0.287937623690715E+02
c_fit_STO( 3, 7, 2)= 0.202602238834114E+00
g_fit_STO( 3, 7, 2)= 0.190379794221872E+01
c_fit_STO( 4, 7, 2)= 0.134521201234128E+00
g_fit_STO( 4, 7, 2)= 0.656896929121867E+01
c_fit_STO( 5, 7, 2)= 0.273979402132560E+00
g_fit_STO( 5, 7, 2)= 0.615908754884964E+00
c_fit_STO( 6, 7, 2)= 0.245178984521786E+00
g_fit_STO( 6, 7, 2)= 0.196088630657315E+00
c_fit_STO( 7, 7, 2)=-0.118038827274348E-04
g_fit_STO( 7, 7, 2)=-0.267903406182100E+00
! chi2=  0.108033177555295E-07

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1, 7, 3)= 0.622053875694260E-01
g_fit_STO( 1, 7, 3)= 0.120279213440970E+03
c_fit_STO( 2, 7, 3)= 0.110377247805451E+00
g_fit_STO( 2, 7, 3)= 0.163528555903519E+02
c_fit_STO( 3, 7, 3)= 0.176986252436944E+00
g_fit_STO( 3, 7, 3)= 0.373581410516049E+01
c_fit_STO( 4, 7, 3)= 0.244334990544326E+00
g_fit_STO( 4, 7, 3)= 0.109351357086875E+01
c_fit_STO( 5, 7, 3)= 0.253506569373002E+00
g_fit_STO( 5, 7, 3)= 0.371893531738258E+00
c_fit_STO( 6, 7, 3)= 0.130946445021309E+00
g_fit_STO( 6, 7, 3)= 0.138829754653350E+00
c_fit_STO( 7, 7, 3)=-0.576454537223164E-06
g_fit_STO( 7, 7, 3)=-0.124211990001521E+00
! chi2=  0.192776976320044E-07

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1, 7, 4)= 0.470755595529027E-01
g_fit_STO( 1, 7, 4)= 0.214559532858958E+03
c_fit_STO( 2, 7, 4)= 0.839238349743176E-01
g_fit_STO( 2, 7, 4)= 0.291700806924872E+02
c_fit_STO( 3, 7, 4)= 0.137047097230960E+00
g_fit_STO( 3, 7, 4)= 0.666445337523973E+01
c_fit_STO( 4, 7, 4)= 0.200205499805010E+00
g_fit_STO( 4, 7, 4)= 0.195255640578433E+01
c_fit_STO( 5, 7, 4)= 0.242902300583743E+00
g_fit_STO( 5, 7, 4)= 0.668177405150732E+00
c_fit_STO( 6, 7, 4)= 0.201371110668777E+00
g_fit_STO( 6, 7, 4)= 0.255949203053193E+00
c_fit_STO( 7, 7, 4)= 0.711083706306515E-01
g_fit_STO( 7, 7, 4)= 0.107746797153786E+00
! chi2=  0.513179049907966E-08

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1, 7, 5)= 0.509991189158313E-01
g_fit_STO( 1, 7, 5)= 0.184777402603910E+03
c_fit_STO( 2, 7, 5)= 0.908048790570700E-01
g_fit_STO( 2, 7, 5)= 0.251206507624501E+02
c_fit_STO( 3, 7, 5)= 0.147557798730488E+00
g_fit_STO( 3, 7, 5)= 0.573928276399635E+01
c_fit_STO( 4, 7, 5)= 0.212247414711706E+00
g_fit_STO( 4, 7, 5)= 0.168158962245635E+01
c_fit_STO( 5, 7, 5)= 0.246698413468600E+00
g_fit_STO( 5, 7, 5)= 0.575662975753481E+00
c_fit_STO( 6, 7, 5)= 0.183499409223293E+00
g_fit_STO( 6, 7, 5)= 0.220935105761592E+00
c_fit_STO( 7, 7, 5)= 0.504597878001231E-01
g_fit_STO( 7, 7, 5)= 0.937750996436013E-01
! chi2=  0.582053806616444E-08

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1, 7, 6)= 0.528261613139224E-01
g_fit_STO( 1, 7, 6)= 0.173400829549175E+03
c_fit_STO( 2, 7, 6)= 0.939973434777755E-01
g_fit_STO( 2, 7, 6)= 0.235733868700883E+02
c_fit_STO( 3, 7, 6)= 0.152362771930565E+00
g_fit_STO( 3, 7, 6)= 0.538567963618199E+01
c_fit_STO( 4, 7, 6)= 0.217463187727225E+00
g_fit_STO( 4, 7, 6)= 0.157786100960521E+01
c_fit_STO( 5, 7, 6)= 0.247561705842765E+00
g_fit_STO( 5, 7, 6)= 0.539913430096705E+00
c_fit_STO( 6, 7, 6)= 0.175009097090308E+00
g_fit_STO( 6, 7, 6)= 0.206773846372478E+00
c_fit_STO( 7, 7, 6)= 0.424102384514667E-01
g_fit_STO( 7, 7, 6)= 0.870270380621985E-01
! chi2=  0.569414627294061E-08

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1, 7, 7)= 0.538675310274197E-01
g_fit_STO( 1, 7, 7)= 0.167559638749923E+03
c_fit_STO( 2, 7, 7)= 0.958128350237810E-01
g_fit_STO( 2, 7, 7)= 0.227785780762005E+02
c_fit_STO( 3, 7, 7)= 0.155069401305638E+00
g_fit_STO( 3, 7, 7)= 0.520399628880096E+01
c_fit_STO( 4, 7, 7)= 0.220292896676975E+00
g_fit_STO( 4, 7, 7)= 0.152453544907468E+01
c_fit_STO( 5, 7, 7)= 0.247714723338341E+00
g_fit_STO( 5, 7, 7)= 0.521477554724784E+00
c_fit_STO( 6, 7, 7)= 0.170085962354483E+00
g_fit_STO( 6, 7, 7)= 0.199346817725303E+00
c_fit_STO( 7, 7, 7)= 0.384249178949077E-01
g_fit_STO( 7, 7, 7)= 0.832026207283928E-01
! chi2=  0.531595412812216E-08

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1, 7, 8)= 0.545408541235768E-01
g_fit_STO( 1, 7, 8)= 0.164025732581544E+03
c_fit_STO( 2, 7, 8)= 0.969848316371578E-01
g_fit_STO( 2, 7, 8)= 0.222973753513922E+02
c_fit_STO( 3, 7, 8)= 0.156805235240649E+00
g_fit_STO( 3, 7, 8)= 0.509396442208794E+01
c_fit_STO( 4, 7, 8)= 0.222059588774538E+00
g_fit_STO( 4, 7, 8)= 0.149222436445052E+01
c_fit_STO( 5, 7, 8)= 0.247667679433417E+00
g_fit_STO( 5, 7, 8)= 0.510283780719632E+00
c_fit_STO( 6, 7, 8)= 0.166873879238786E+00
g_fit_STO( 6, 7, 8)= 0.194793980610871E+00
c_fit_STO( 7, 7, 8)= 0.361024956974110E-01
g_fit_STO( 7, 7, 8)= 0.807669799658148E-01
! chi2=  0.490386623100079E-08

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1, 7, 9)= 0.550128178993745E-01
g_fit_STO( 1, 7, 9)= 0.161662870285406E+03
c_fit_STO( 2, 7, 9)= 0.978054676892177E-01
g_fit_STO( 2, 7, 9)= 0.219752685934464E+02
c_fit_STO( 3, 7, 9)= 0.158014881661515E+00
g_fit_STO( 3, 7, 9)= 0.502027541203286E+01
c_fit_STO( 4, 7, 9)= 0.223265987251499E+00
g_fit_STO( 4, 7, 9)= 0.147057436075736E+01
c_fit_STO( 5, 7, 9)= 0.247562225978737E+00
g_fit_STO( 5, 7, 9)= 0.502771282574874E+00
c_fit_STO( 6, 7, 9)= 0.164613088923120E+00
g_fit_STO( 6, 7, 9)= 0.191718108491208E+00
c_fit_STO( 7, 7, 9)= 0.345968165020813E-01
g_fit_STO( 7, 7, 9)= 0.790864710084219E-01
! chi2=  0.451867922610205E-08

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1, 7,10)= 0.553624721494418E-01
g_fit_STO( 1, 7,10)= 0.159974175813968E+03
c_fit_STO( 2, 7,10)= 0.984129088588709E-01
g_fit_STO( 2, 7,10)= 0.217447410045189E+02
c_fit_STO( 3, 7,10)= 0.158906949495593E+00
g_fit_STO( 3, 7,10)= 0.496750847363692E+01
c_fit_STO( 4, 7,10)= 0.224141736629270E+00
g_fit_STO( 4, 7,10)= 0.145506352373055E+01
c_fit_STO( 5, 7,10)= 0.247444439609349E+00
g_fit_STO( 5, 7,10)= 0.497381299416966E+00
c_fit_STO( 6, 7,10)= 0.162935802733304E+00
g_fit_STO( 6, 7,10)= 0.189499913907209E+00
c_fit_STO( 7, 7,10)= 0.335465359775729E-01
g_fit_STO( 7, 7,10)= 0.778590886253315E-01
! chi2=  0.417449424925956E-08

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1, 7,11)= 0.556320826417495E-01
g_fit_STO( 1, 7,11)= 0.158709469679962E+03
c_fit_STO( 2, 7,11)= 0.988810385269084E-01
g_fit_STO( 2, 7,11)= 0.215717580266297E+02
c_fit_STO( 3, 7,11)= 0.159592509253925E+00
g_fit_STO( 3, 7,11)= 0.492788227752767E+01
c_fit_STO( 4, 7,11)= 0.224806378134338E+00
g_fit_STO( 4, 7,11)= 0.144340792116416E+01
c_fit_STO( 5, 7,11)= 0.247330138095938E+00
g_fit_STO( 5, 7,11)= 0.493325579722833E+00
c_fit_STO( 6, 7,11)= 0.161642115222221E+00
g_fit_STO( 6, 7,11)= 0.187823877944256E+00
c_fit_STO( 7, 7,11)= 0.327742441216164E-01
g_fit_STO( 7, 7,11)= 0.769241047715778E-01
! chi2=  0.387121029528846E-08

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1, 7,12)= 0.558463850433443E-01
g_fit_STO( 1, 7,12)= 0.157728753332027E+03
c_fit_STO( 2, 7,12)= 0.992529862089223E-01
g_fit_STO( 2, 7,12)= 0.214372928916774E+02
c_fit_STO( 3, 7,12)= 0.160136022071467E+00
g_fit_STO( 3, 7,12)= 0.489705023573379E+01
c_fit_STO( 4, 7,12)= 0.225327942478770E+00
g_fit_STO( 4, 7,12)= 0.143433270996441E+01
c_fit_STO( 5, 7,12)= 0.247224574649504E+00
g_fit_STO( 5, 7,12)= 0.490163763324529E+00
c_fit_STO( 6, 7,12)= 0.160614260060245E+00
g_fit_STO( 6, 7,12)= 0.186512628113206E+00
c_fit_STO( 7, 7,12)= 0.321834669335214E-01
g_fit_STO( 7, 7,12)= 0.761885107696956E-01
! chi2=  0.360464226614255E-08

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1, 7,13)= 0.560208433615631E-01
g_fit_STO( 1, 7,13)= 0.156947698803189E+03
c_fit_STO( 2, 7,13)= 0.995556978317208E-01
g_fit_STO( 2, 7,13)= 0.213298790042029E+02
c_fit_STO( 3, 7,13)= 0.160577612035353E+00
g_fit_STO( 3, 7,13)= 0.487239178466311E+01
c_fit_STO( 4, 7,13)= 0.225748156276447E+00
g_fit_STO( 4, 7,13)= 0.142706878996278E+01
c_fit_STO( 5, 7,13)= 0.247128988734662E+00
g_fit_STO( 5, 7,13)= 0.487629907462833E+00
c_fit_STO( 6, 7,13)= 0.159778103666836E+00
g_fit_STO( 6, 7,13)= 0.185458599154354E+00
c_fit_STO( 7, 7,13)= 0.317174484182918E-01
g_fit_STO( 7, 7,13)= 0.755948474600926E-01
! chi2=  0.336984129186804E-08

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1, 7,14)= 0.561656023892144E-01
g_fit_STO( 1, 7,14)= 0.156312697400986E+03
c_fit_STO( 2, 7,14)= 0.998068321210908E-01
g_fit_STO( 2, 7,14)= 0.212422200831667E+02
c_fit_STO( 3, 7,14)= 0.160943527327455E+00
g_fit_STO( 3, 7,14)= 0.485223841854270E+01
c_fit_STO( 4, 7,14)= 0.226093959865862E+00
g_fit_STO( 4, 7,14)= 0.142112615100072E+01
c_fit_STO( 5, 7,14)= 0.247043136292494E+00
g_fit_STO( 5, 7,14)= 0.485554357694885E+00
c_fit_STO( 6, 7,14)= 0.159084803681324E+00
g_fit_STO( 6, 7,14)= 0.184592862806456E+00
c_fit_STO( 7, 7,14)= 0.313407575271321E-01
g_fit_STO( 7, 7,14)= 0.751057664671619E-01
! chi2=  0.316217089608262E-08

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1, 7,15)= 0.562876465297610E-01
g_fit_STO( 1, 7,15)= 0.155787528386156E+03
c_fit_STO( 2, 7,15)= 0.100018550135463E+00
g_fit_STO( 2, 7,15)= 0.211693976952944E+02
c_fit_STO( 3, 7,15)= 0.161251700828152E+00
g_fit_STO( 3, 7,15)= 0.483546894377758E+01
c_fit_STO( 4, 7,15)= 0.226383566073716E+00
g_fit_STO( 4, 7,15)= 0.141617600123439E+01
c_fit_STO( 5, 7,15)= 0.246966199287241E+00
g_fit_STO( 5, 7,15)= 0.483823247781695E+00
c_fit_STO( 6, 7,15)= 0.158500709261028E+00
g_fit_STO( 6, 7,15)= 0.183869032637727E+00
c_fit_STO( 7, 7,15)= 0.310301084924666E-01
g_fit_STO( 7, 7,15)= 0.746959134633559E-01
! chi2=  0.297760257281562E-08

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1, 7,16)= 0.563918250607441E-01
g_fit_STO( 1, 7,16)= 0.155347655130569E+03
c_fit_STO( 2, 7,16)= 0.100199350593988E+00
g_fit_STO( 2, 7,16)= 0.211080783814840E+02
c_fit_STO( 3, 7,16)= 0.161514706577087E+00
g_fit_STO( 3, 7,16)= 0.482131758485964E+01
c_fit_STO( 4, 7,16)= 0.226629475635975E+00
g_fit_STO( 4, 7,16)= 0.141199361310138E+01
c_fit_STO( 5, 7,16)= 0.246897246857529E+00
g_fit_STO( 5, 7,16)= 0.482358770656922E+00
c_fit_STO( 6, 7,16)= 0.158002361548902E+00
g_fit_STO( 6, 7,16)= 0.183255200517401E+00
c_fit_STO( 7, 7,16)= 0.307697754217347E-01
g_fit_STO( 7, 7,16)= 0.743476023490348E-01
! chi2=  0.281273814487675E-08

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1, 8, 2)= 0.660876090839233E-01
g_fit_STO( 1, 8, 2)= 0.846456457803666E+02
c_fit_STO( 2, 8, 2)= 0.899990469074802E-01
g_fit_STO( 2, 8, 2)= 0.153577048549252E+02
c_fit_STO( 3, 8, 2)= 0.130130373183963E+00
g_fit_STO( 3, 8, 2)= 0.467105167124630E+01
c_fit_STO( 4, 8, 2)= 0.195073188256511E+00
g_fit_STO( 4, 8, 2)= 0.160228298404629E+01
c_fit_STO( 5, 8, 2)= 0.261130645877827E+00
g_fit_STO( 5, 8, 2)= 0.565097414562401E+00
c_fit_STO( 6, 8, 2)= 0.231032985054857E+00
g_fit_STO( 6, 8, 2)= 0.190278756579774E+00
c_fit_STO( 7, 8, 2)=-0.839471756387235E-04
g_fit_STO( 7, 8, 2)=-0.169095729577064E+00
c_fit_STO( 8, 8, 2)=-0.601702511903412E-11
g_fit_STO( 8, 8, 2)=-0.882689347151314E+00
! chi2=  0.280367902719398E-07

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1, 8, 3)= 0.674913992216507E-01
g_fit_STO( 1, 8, 3)= 0.847165161739710E+02
c_fit_STO( 2, 8, 3)= 0.928985624760478E-01
g_fit_STO( 2, 8, 3)= 0.152799344088636E+02
c_fit_STO( 3, 8, 3)= 0.127134607777069E+00
g_fit_STO( 3, 8, 3)= 0.468954325557077E+01
c_fit_STO( 4, 8, 3)= 0.171048388624719E+00
g_fit_STO( 4, 8, 3)= 0.174248014236175E+01
c_fit_STO( 5, 8, 3)= 0.208348927248174E+00
g_fit_STO( 5, 8, 3)= 0.706485266412505E+00
c_fit_STO( 6, 8, 3)= 0.200663587255745E+00
g_fit_STO( 6, 8, 3)= 0.300425767796817E+00
c_fit_STO( 7, 8, 3)= 0.105255926601423E+00
g_fit_STO( 7, 8, 3)= 0.129406864942115E+00
c_fit_STO( 8, 8, 3)=-0.768242689253120E-04
g_fit_STO( 8, 8, 3)=-0.228194036028074E-01
! chi2=  0.187335246415179E-07

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1, 8, 4)= 0.470754973083160E-01
g_fit_STO( 1, 8, 4)= 0.214560084782754E+03
c_fit_STO( 2, 8, 4)= 0.839237227071787E-01
g_fit_STO( 2, 8, 4)= 0.291701583747425E+02
c_fit_STO( 3, 8, 4)= 0.137046933483880E+00
g_fit_STO( 3, 8, 4)= 0.666447108408064E+01
c_fit_STO( 4, 8, 4)= 0.200205357168614E+00
g_fit_STO( 4, 8, 4)= 0.195256128260548E+01
c_fit_STO( 5, 8, 4)= 0.242902368249795E+00
g_fit_STO( 5, 8, 4)= 0.668178777306934E+00
c_fit_STO( 6, 8, 4)= 0.201775631118528E+00
g_fit_STO( 6, 8, 4)= 0.255949541982473E+00
c_fit_STO( 7, 8, 4)= 0.711085048142058E-01
g_fit_STO( 7, 8, 4)= 0.107746850697248E+00
c_fit_STO( 8, 8, 4)=-0.404220746721309E-03
g_fit_STO( 8, 8, 4)= 0.255951668417463E+00
! chi2=  0.513179049960053E-08

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1, 8, 5)= 0.686735719277879E-01
g_fit_STO( 1, 8, 5)= 0.847642650656455E+02
c_fit_STO( 2, 8, 5)= 0.942019513333408E-01
g_fit_STO( 2, 8, 5)= 0.152341445530511E+02
c_fit_STO( 3, 8, 5)= 0.124500086154504E+00
g_fit_STO( 3, 8, 5)= 0.474911011728634E+01
c_fit_STO( 4, 8, 5)= 0.161312455988325E+00
g_fit_STO( 4, 8, 5)= 0.183211853991350E+01
c_fit_STO( 5, 8, 5)= 0.191908337903700E+00
g_fit_STO( 5, 8, 5)= 0.778308900822393E+00
c_fit_STO( 6, 8, 5)= 0.184373726631394E+00
g_fit_STO( 6, 8, 5)= 0.351894461378370E+00
c_fit_STO( 7, 8, 5)= 0.116374234122322E+00
g_fit_STO( 7, 8, 5)= 0.167858742587196E+00
c_fit_STO( 8, 8, 5)= 0.310403980205564E-01
g_fit_STO( 8, 8, 5)= 0.842254827332432E-01
! chi2=  0.115664992018419E-07

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1, 8, 6)= 0.687802803767679E-01
g_fit_STO( 1, 8, 6)= 0.847039077383590E+02
c_fit_STO( 2, 8, 6)= 0.945415956779800E-01
g_fit_STO( 2, 8, 6)= 0.152752286754573E+02
c_fit_STO( 3, 8, 6)= 0.127808344137907E+00
g_fit_STO( 3, 8, 6)= 0.470658086450329E+01
c_fit_STO( 4, 8, 6)= 0.168172779238176E+00
g_fit_STO( 4, 8, 6)= 0.177300028056909E+01
c_fit_STO( 5, 8, 6)= 0.199656149608191E+00
g_fit_STO( 5, 8, 6)= 0.734592260856431E+00
c_fit_STO( 6, 8, 6)= 0.185756758715205E+00
g_fit_STO( 6, 8, 6)= 0.324196584885744E+00
c_fit_STO( 7, 8, 6)= 0.105606702802912E+00
g_fit_STO( 7, 8, 6)= 0.151388573527676E+00
c_fit_STO( 8, 8, 6)= 0.219202106774891E-01
g_fit_STO( 8, 8, 6)= 0.750719504754877E-01
! chi2=  0.973457020112165E-08

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1, 8, 7)= 0.688093033173755E-01
g_fit_STO( 1, 8, 7)= 0.955316558629852E+02
c_fit_STO( 2, 8, 7)= 0.108724434075507E+00
g_fit_STO( 2, 8, 7)= 0.146000757333873E+02
c_fit_STO( 3, 8, 7)= 0.152495124199401E+00
g_fit_STO( 3, 8, 7)= 0.392431691827094E+01
c_fit_STO( 4, 8, 7)= 0.190423720442235E+00
g_fit_STO( 4, 8, 7)= 0.138141069417779E+01
c_fit_STO( 5, 8, 7)= 0.202226469439303E+00
g_fit_STO( 5, 8, 7)= 0.568629464377631E+00
c_fit_STO( 6, 8, 7)= 0.162218399121042E+00
g_fit_STO( 6, 8, 7)= 0.258862264022886E+00
c_fit_STO( 7, 8, 7)= 0.769423097422468E-01
g_fit_STO( 7, 8, 7)= 0.126589319904577E+00
c_fit_STO( 8, 8, 7)= 0.128060614695337E-01
g_fit_STO( 8, 8, 7)= 0.657387155227691E-01
! chi2=  0.724240779846748E-08

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1, 8, 8)= 0.690091459086973E-01
g_fit_STO( 1, 8, 8)= 0.847047963005952E+02
c_fit_STO( 2, 8, 8)= 0.949252188148874E-01
g_fit_STO( 2, 8, 8)= 0.152954710613617E+02
c_fit_STO( 3, 8, 8)= 0.130095544838047E+00
g_fit_STO( 3, 8, 8)= 0.468164256354882E+01
c_fit_STO( 4, 8, 8)= 0.173529390752830E+00
g_fit_STO( 4, 8, 8)= 0.173359362461665E+01
c_fit_STO( 5, 8, 8)= 0.205320325629391E+00
g_fit_STO( 5, 8, 8)= 0.704123150042392E+00
c_fit_STO( 6, 8, 8)= 0.185461511455705E+00
g_fit_STO( 6, 8, 8)= 0.305150224960061E+00
c_fit_STO( 7, 8, 8)= 0.973400424491330E-01
g_fit_STO( 7, 8, 8)= 0.139848561544986E+00
c_fit_STO( 8, 8, 8)= 0.164329528720414E-01
g_fit_STO( 8, 8, 8)= 0.673934280323745E-01
! chi2=  0.738367635749682E-08

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1, 8, 9)= 0.690899212013969E-01
g_fit_STO( 1, 8, 9)= 0.847020722031265E+02
c_fit_STO( 2, 8, 9)= 0.950587127522805E-01
g_fit_STO( 2, 8, 9)= 0.152992911077353E+02
c_fit_STO( 3, 8, 9)= 0.130760873371718E+00
g_fit_STO( 3, 8, 9)= 0.467398083892173E+01
c_fit_STO( 4, 8, 9)= 0.174942348872045E+00
g_fit_STO( 4, 8, 9)= 0.172290426694854E+01
c_fit_STO( 5, 8, 9)= 0.206592603411335E+00
g_fit_STO( 5, 8, 9)= 0.696448145403587E+00
c_fit_STO( 6, 8, 9)= 0.185060217933165E+00
g_fit_STO( 6, 8, 9)= 0.300596187832830E+00
c_fit_STO( 7, 8, 9)= 0.952323978598818E-01
g_fit_STO( 7, 8, 9)= 0.137138077288422E+00
c_fit_STO( 8, 8, 9)= 0.153361824524420E-01
g_fit_STO( 8, 8, 9)= 0.654832633275660E-01
! chi2=  0.658818184535291E-08

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1, 8,10)= 0.691529618219387E-01
g_fit_STO( 1, 8,10)= 0.846990508860697E+02
c_fit_STO( 2, 8,10)= 0.951585621143249E-01
g_fit_STO( 2, 8,10)= 0.153028938978494E+02
c_fit_STO( 3, 8,10)= 0.131245946745225E+00
g_fit_STO( 3, 8,10)= 0.466872051309260E+01
c_fit_STO( 4, 8,10)= 0.175937778225017E+00
g_fit_STO( 4, 8,10)= 0.171545712013983E+01
c_fit_STO( 5, 8,10)= 0.207457727375568E+00
g_fit_STO( 5, 8,10)= 0.691138009221625E+00
c_fit_STO( 6, 8,10)= 0.184718129158481E+00
g_fit_STO( 6, 8,10)= 0.297450204011489E+00
c_fit_STO( 7, 8,10)= 0.937451747396456E-01
g_fit_STO( 7, 8,10)= 0.135256430528484E+00
c_fit_STO( 8, 8,10)= 0.146243690276654E-01
g_fit_STO( 8, 8,10)= 0.641345689197744E-01
! chi2=  0.594759136691706E-08

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1, 8,11)= 0.692039641651626E-01
g_fit_STO( 1, 8,11)= 0.846962585870636E+02
c_fit_STO( 2, 8,11)= 0.952365597000290E-01
g_fit_STO( 2, 8,11)= 0.153060473110195E+02
c_fit_STO( 3, 8,11)= 0.131618061757272E+00
g_fit_STO( 3, 8,11)= 0.466486473983092E+01
c_fit_STO( 4, 8,11)= 0.176683392004202E+00
g_fit_STO( 4, 8,11)= 0.170992425563287E+01
c_fit_STO( 5, 8,11)= 0.208084453695333E+00
g_fit_STO( 5, 8,11)= 0.687216473911676E+00
c_fit_STO( 6, 8,11)= 0.184425816533592E+00
g_fit_STO( 6, 8,11)= 0.295133584690121E+00
c_fit_STO( 7, 8,11)= 0.926348461524726E-01
g_fit_STO( 7, 8,11)= 0.133868235170154E+00
c_fit_STO( 8, 8,11)= 0.141271355023591E-01
g_fit_STO( 8, 8,11)= 0.631318688720054E-01
! chi2=  0.542072906073498E-08

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1, 8,12)= 0.692484290073593E-01
g_fit_STO( 1, 8,12)= 0.846959796395287E+02
c_fit_STO( 2, 8,12)= 0.953110158588053E-01
g_fit_STO( 2, 8,12)= 0.153072293704629E+02
c_fit_STO( 3, 8,12)= 0.131899995559503E+00
g_fit_STO( 3, 8,12)= 0.466153746760789E+01
c_fit_STO( 4, 8,12)= 0.177215060957150E+00
g_fit_STO( 4, 8,12)= 0.170585877564340E+01
c_fit_STO( 5, 8,12)= 0.208526619428896E+00
g_fit_STO( 5, 8,12)= 0.684396226644297E+00
c_fit_STO( 6, 8,12)= 0.184201728516965E+00
g_fit_STO( 6, 8,12)= 0.293448875013278E+00
c_fit_STO( 7, 8,12)= 0.918181980468918E-01
g_fit_STO( 7, 8,12)= 0.132835465641637E+00
c_fit_STO( 8, 8,12)= 0.137724172356456E-01
g_fit_STO( 8, 8,12)= 0.623681920633965E-01
! chi2=  0.497956629267541E-08

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1, 8,13)= 0.388422815200078E-01
g_fit_STO( 1, 8,13)= 0.327222525067302E+03
c_fit_STO( 2, 8,13)= 0.693814922164343E-01
g_fit_STO( 2, 8,13)= 0.444511094506433E+02
c_fit_STO( 3, 8,13)= 0.114167256635150E+00
g_fit_STO( 3, 8,13)= 0.101520918639548E+02
c_fit_STO( 4, 8,13)= 0.170883628688348E+00
g_fit_STO( 4, 8,13)= 0.297291978672760E+01
c_fit_STO( 5, 8,13)= 0.222170740848726E+00
g_fit_STO( 5, 8,13)= 0.101555376605993E+01
c_fit_STO( 6, 8,13)= 0.221330514225058E+00
g_fit_STO( 6, 8,13)= 0.386369915675865E+00
c_fit_STO( 7, 8,13)= 0.127476912866474E+00
g_fit_STO( 7, 8,13)= 0.159244243281345E+00
c_fit_STO( 8, 8,13)= 0.222671833130189E-01
g_fit_STO( 8, 8,13)= 0.691236865355681E-01
! chi2=  0.821951612685088E-09

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1, 8,14)= 0.693163361223335E-01
g_fit_STO( 1, 8,14)= 0.846931276620891E+02
c_fit_STO( 2, 8,14)= 0.954122911090630E-01
g_fit_STO( 2, 8,14)= 0.153106785279944E+02
c_fit_STO( 3, 8,14)= 0.132330460972913E+00
g_fit_STO( 3, 8,14)= 0.465715338297418E+01
c_fit_STO( 4, 8,14)= 0.178038618223853E+00
g_fit_STO( 4, 8,14)= 0.169977931858686E+01
c_fit_STO( 5, 8,14)= 0.209183855895498E+00
g_fit_STO( 5, 8,14)= 0.680137094675435E+00
c_fit_STO( 6, 8,14)= 0.183815282225993E+00
g_fit_STO( 6, 8,14)= 0.290930743741581E+00
c_fit_STO( 7, 8,14)= 0.905880343133265E-01
g_fit_STO( 7, 8,14)= 0.131311728989676E+00
c_fit_STO( 8, 8,14)= 0.132757158228419E-01
g_fit_STO( 8, 8,14)= 0.612549207294788E-01
! chi2=  0.428339165196821E-08

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1, 8,15)= 0.693451245573205E-01
g_fit_STO( 1, 8,15)= 0.846933490270294E+02
c_fit_STO( 2, 8,15)= 0.954589927033277E-01
g_fit_STO( 2, 8,15)= 0.153110667541496E+02
c_fit_STO( 3, 8,15)= 0.132485543087877E+00
g_fit_STO( 3, 8,15)= 0.465530927201208E+01
c_fit_STO( 4, 8,15)= 0.178317010644024E+00
g_fit_STO( 4, 8,15)= 0.169765466121135E+01
c_fit_STO( 5, 8,15)= 0.209404239072837E+00
g_fit_STO( 5, 8,15)= 0.678679603004709E+00
c_fit_STO( 6, 8,15)= 0.183675846048966E+00
g_fit_STO( 6, 8,15)= 0.290054931458167E+00
c_fit_STO( 7, 8,15)= 0.901560120842467E-01
g_fit_STO( 7, 8,15)= 0.130766769470436E+00
c_fit_STO( 8, 8,15)= 0.131056543419370E-01
g_fit_STO( 8, 8,15)= 0.608471454993058E-01
! chi2=  0.400373366470796E-08

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1, 8,16)= 0.693697963816363E-01
g_fit_STO( 1, 8,16)= 0.846928649807612E+02
c_fit_STO( 2, 8,16)= 0.954958432839831E-01
g_fit_STO( 2, 8,16)= 0.153118711889669E+02
c_fit_STO( 3, 8,16)= 0.132622658950664E+00
g_fit_STO( 3, 8,16)= 0.465385385509413E+01
c_fit_STO( 4, 8,16)= 0.178568838557214E+00
g_fit_STO( 4, 8,16)= 0.169578383421488E+01
c_fit_STO( 5, 8,16)= 0.209598391940180E+00
g_fit_STO( 5, 8,16)= 0.677384068556170E+00
c_fit_STO( 6, 8,16)= 0.183544019206098E+00
g_fit_STO( 6, 8,16)= 0.289285456613909E+00
c_fit_STO( 7, 8,16)= 0.897761368161980E-01
g_fit_STO( 7, 8,16)= 0.130295527193849E+00
c_fit_STO( 8, 8,16)= 0.129619936313130E-01
g_fit_STO( 8, 8,16)= 0.604997822802923E-01
! chi2=  0.375864815953168E-08

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1, 9, 2)= 0.479266132085897E-01
g_fit_STO( 1, 9, 2)= 0.153748560961697E+03
c_fit_STO( 2, 9, 2)= 0.700029259529862E-01
g_fit_STO( 2, 9, 2)= 0.281710111581896E+02
c_fit_STO( 3, 9, 2)= 0.106912696587068E+00
g_fit_STO( 3, 9, 2)= 0.787115852295898E+01
c_fit_STO( 4, 9, 2)= 0.138251185498573E+00
g_fit_STO( 4, 9, 2)= 0.279002655482305E+01
c_fit_STO( 5, 9, 2)= 0.172399205823785E+00
g_fit_STO( 5, 9, 2)= 0.117533323724022E+01
c_fit_STO( 6, 9, 2)= 0.229685762484612E+00
g_fit_STO( 6, 9, 2)= 0.498505210449993E+00
c_fit_STO( 7, 9, 2)= 0.214892902180086E+00
g_fit_STO( 7, 9, 2)= 0.184579942620037E+00
c_fit_STO( 8, 9, 2)=-0.125077609474629E-03
g_fit_STO( 8, 9, 2)=-0.156264732131633E+00
c_fit_STO( 9, 9, 2)=-0.453964612708905E-19
g_fit_STO( 9, 9, 2)=-0.176077854695123E+01
! chi2=  0.845071518102453E-08

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1, 9, 3)= 0.491531271435804E-01
g_fit_STO( 1, 9, 3)= 0.153730740132590E+03
c_fit_STO( 2, 9, 3)= 0.709680841833113E-01
g_fit_STO( 2, 9, 3)= 0.281926515585881E+02
c_fit_STO( 3, 9, 3)= 0.109873349297352E+00
g_fit_STO( 3, 9, 3)= 0.791168178959975E+01
c_fit_STO( 4, 9, 3)= 0.151249185202202E+00
g_fit_STO( 4, 9, 3)= 0.270541237582688E+01
c_fit_STO( 5, 9, 3)= 0.186829979051303E+00
g_fit_STO( 5, 9, 3)= 0.107444674520083E+01
c_fit_STO( 6, 9, 3)= 0.201284077783336E+00
g_fit_STO( 6, 9, 3)= 0.466958350014581E+00
c_fit_STO( 7, 9, 3)= 0.159546747139609E+00
g_fit_STO( 7, 9, 3)= 0.211901928034770E+00
c_fit_STO( 8, 9, 3)= 0.543769506552291E-01
g_fit_STO( 8, 9, 3)= 0.986807697317238E-01
c_fit_STO( 9, 9, 3)=-0.374127964895911E-02
g_fit_STO( 9, 9, 3)= 0.424542811239070E-01
! chi2=  0.583339759102229E-08

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1, 9, 4)= 0.497194776775800E-01
g_fit_STO( 1, 9, 4)= 0.153735761471488E+03
c_fit_STO( 2, 9, 4)= 0.714587291989285E-01
g_fit_STO( 2, 9, 4)= 0.281811597815622E+02
c_fit_STO( 3, 9, 4)= 0.110667459968052E+00
g_fit_STO( 3, 9, 4)= 0.793806573810223E+01
c_fit_STO( 4, 9, 4)= 0.154133886603917E+00
g_fit_STO( 4, 9, 4)= 0.269828093259399E+01
c_fit_STO( 5, 9, 4)= 0.190210952049794E+00
g_fit_STO( 5, 9, 4)= 0.105868545564671E+01
c_fit_STO( 6, 9, 4)= 0.198254544309990E+00
g_fit_STO( 6, 9, 4)= 0.458530712490962E+00
c_fit_STO( 7, 9, 4)= 0.150088327915089E+00
g_fit_STO( 7, 9, 4)= 0.211485379764254E+00
c_fit_STO( 8, 9, 4)= 0.548324667762802E-01
g_fit_STO( 8, 9, 4)= 0.101196577478555E+00
c_fit_STO( 9, 9, 4)=-0.854583422183632E-05
g_fit_STO( 9, 9, 4)=-0.600827632059177E-02
! chi2=  0.445914104823020E-08

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1, 9, 5)= 0.499681962928772E-01
g_fit_STO( 1, 9, 5)= 0.153733168469803E+03
c_fit_STO( 2, 9, 5)= 0.718954630323164E-01
g_fit_STO( 2, 9, 5)= 0.281865780168667E+02
c_fit_STO( 3, 9, 5)= 0.111366131278796E+00
g_fit_STO( 3, 9, 5)= 0.792781918205125E+01
c_fit_STO( 4, 9, 5)= 0.154349956039635E+00
g_fit_STO( 4, 9, 5)= 0.269591267147475E+01
c_fit_STO( 5, 9, 5)= 0.185848226262788E+00
g_fit_STO( 5, 9, 5)= 0.106696548621903E+01
c_fit_STO( 6, 9, 5)= 0.178006509708716E+00
g_fit_STO( 6, 9, 5)= 0.482382705696336E+00
c_fit_STO( 7, 9, 5)= 0.121669293465379E+00
g_fit_STO( 7, 9, 5)= 0.254913283752491E+00
c_fit_STO( 8, 9, 5)= 0.793648835593372E-01
g_fit_STO( 8, 9, 5)= 0.148787036960465E+00
c_fit_STO( 9, 9, 5)= 0.267678818334545E-01
g_fit_STO( 9, 9, 5)= 0.822267635453250E-01
! chi2=  0.361008492732274E-08

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1, 9, 6)= 0.501517837457960E-01
g_fit_STO( 1, 9, 6)= 0.153732987529138E+03
c_fit_STO( 2, 9, 6)= 0.721233289898593E-01
g_fit_STO( 2, 9, 6)= 0.281873744930455E+02
c_fit_STO( 3, 9, 6)= 0.111813431047285E+00
g_fit_STO( 3, 9, 6)= 0.792483825410767E+01
c_fit_STO( 4, 9, 6)= 0.154104227620077E+00
g_fit_STO( 4, 9, 6)= 0.269699638525754E+01
c_fit_STO( 5, 9, 6)= 0.184360890466434E+00
g_fit_STO( 5, 9, 6)= 0.107366564356730E+01
c_fit_STO( 6, 9, 6)= 0.184323198338901E+00
g_fit_STO( 6, 9, 6)= 0.482415986995843E+00
c_fit_STO( 7, 9, 6)= 0.140354210487227E+00
g_fit_STO( 7, 9, 6)= 0.237423431946787E+00
c_fit_STO( 8, 9, 6)= 0.679192715784772E-01
g_fit_STO( 8, 9, 6)= 0.125164751363754E+00
c_fit_STO( 9, 9, 6)= 0.140164201787461E-01
g_fit_STO( 9, 9, 6)= 0.694875911487394E-01
! chi2=  0.302663264513475E-08

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1, 9, 7)= 0.502027010103165E-01
g_fit_STO( 1, 9, 7)= 0.153731149887804E+03
c_fit_STO( 2, 9, 7)= 0.724602284217484E-01
g_fit_STO( 2, 9, 7)= 0.281915221608013E+02
c_fit_STO( 3, 9, 7)= 0.111867835294167E+00
g_fit_STO( 3, 9, 7)= 0.791506475318410E+01
c_fit_STO( 4, 9, 7)= 0.154008627310450E+00
g_fit_STO( 4, 9, 7)= 0.270293380019526E+01
c_fit_STO( 5, 9, 7)= 0.186350020136102E+00
g_fit_STO( 5, 9, 7)= 0.107167345106893E+01
c_fit_STO( 6, 9, 7)= 0.186300064092603E+00
g_fit_STO( 6, 9, 7)= 0.476783377870202E+00
c_fit_STO( 7, 9, 7)= 0.140504683199418E+00
g_fit_STO( 7, 9, 7)= 0.232272316184132E+00
c_fit_STO( 8, 9, 7)= 0.659820838314859E-01
g_fit_STO( 8, 9, 7)= 0.120098565179756E+00
c_fit_STO( 9, 9, 7)= 0.114205225503801E-01
g_fit_STO( 9, 9, 7)= 0.645829575429163E-01
! chi2=  0.260987365215933E-08

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1, 9, 8)= 0.503533113211890E-01
g_fit_STO( 1, 9, 8)= 0.153732611747969E+03
c_fit_STO( 2, 9, 8)= 0.724416051469933E-01
g_fit_STO( 2, 9, 8)= 0.281883466706947E+02
c_fit_STO( 3, 9, 8)= 0.112212561665229E+00
g_fit_STO( 3, 9, 8)= 0.792197259923803E+01
c_fit_STO( 4, 9, 8)= 0.154111993969903E+00
g_fit_STO( 4, 9, 8)= 0.269973948537717E+01
c_fit_STO( 5, 9, 8)= 0.185021896677659E+00
g_fit_STO( 5, 9, 8)= 0.107520352085459E+01
c_fit_STO( 6, 9, 8)= 0.188573024289314E+00
g_fit_STO( 6, 9, 8)= 0.477836011306014E+00
c_fit_STO( 7, 9, 8)= 0.144133702788656E+00
g_fit_STO( 7, 9, 8)= 0.228040941396182E+00
c_fit_STO( 8, 9, 8)= 0.631141456322689E-01
g_fit_STO( 8, 9, 8)= 0.114656237665334E+00
c_fit_STO( 9, 9, 8)= 0.911518326693830E-02
g_fit_STO( 9, 9, 8)= 0.604086628011836E-01
! chi2=  0.229140674416782E-08

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1, 9, 9)= 0.503736148836615E-01
g_fit_STO( 1, 9, 9)= 0.153730383216613E+03
c_fit_STO( 2, 9, 9)= 0.726538114399546E-01
g_fit_STO( 2, 9, 9)= 0.281933537410777E+02
c_fit_STO( 3, 9, 9)= 0.112332968911097E+00
g_fit_STO( 3, 9, 9)= 0.791064115904286E+01
c_fit_STO( 4, 9, 9)= 0.153494206672485E+00
g_fit_STO( 4, 9, 9)= 0.270438218489999E+01
c_fit_STO( 5, 9, 9)= 0.184672589772030E+00
g_fit_STO( 5, 9, 9)= 0.107992126858169E+01
c_fit_STO( 6, 9, 9)= 0.189161941082680E+00
g_fit_STO( 6, 9, 9)= 0.479198708210657E+00
c_fit_STO( 7, 9, 9)= 0.144824510062515E+00
g_fit_STO( 7, 9, 9)= 0.227560818304844E+00
c_fit_STO( 8, 9, 9)= 0.628793545248709E-01
g_fit_STO( 8, 9, 9)= 0.113416526845720E+00
c_fit_STO( 9, 9, 9)= 0.864323400269650E-02
g_fit_STO( 9, 9, 9)= 0.587437624968238E-01
! chi2=  0.204480128583373E-08

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1, 9,10)= 0.504476261066923E-01
g_fit_STO( 1, 9,10)= 0.153730986264021E+03
c_fit_STO( 2, 9,10)= 0.726806458533379E-01
g_fit_STO( 2, 9,10)= 0.281918750374211E+02
c_fit_STO( 3, 9,10)= 0.112466284845188E+00
g_fit_STO( 3, 9,10)= 0.791447227457011E+01
c_fit_STO( 4, 9,10)= 0.153992264994601E+00
g_fit_STO( 4, 9,10)= 0.270145908899941E+01
c_fit_STO( 5, 9,10)= 0.185427126316753E+00
g_fit_STO( 5, 9,10)= 0.107625762363072E+01
c_fit_STO( 6, 9,10)= 0.190065042493034E+00
g_fit_STO( 6, 9,10)= 0.475881418458741E+00
c_fit_STO( 7, 9,10)= 0.144613182414163E+00
g_fit_STO( 7, 9,10)= 0.224767069937803E+00
c_fit_STO( 8, 9,10)= 0.613467284187992E-01
g_fit_STO( 8, 9,10)= 0.111302991748792E+00
c_fit_STO( 9, 9,10)= 0.798088386635102E-02
g_fit_STO( 9, 9,10)= 0.570988067086512E-01
! chi2=  0.184508144057501E-08

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1, 9,11)= 0.504799341269026E-01
g_fit_STO( 1, 9,11)= 0.153731864178266E+03
c_fit_STO( 2, 9,11)= 0.727715461470105E-01
g_fit_STO( 2, 9,11)= 0.281901661791791E+02
c_fit_STO( 3, 9,11)= 0.112391288572251E+00
g_fit_STO( 3, 9,11)= 0.791703520846957E+01
c_fit_STO( 4, 9,11)= 0.154088036269884E+00
g_fit_STO( 4, 9,11)= 0.270475812484003E+01
c_fit_STO( 5, 9,11)= 0.185749485118617E+00
g_fit_STO( 5, 9,11)= 0.107578471806097E+01
c_fit_STO( 6, 9,11)= 0.188812468181801E+00
g_fit_STO( 6, 9,11)= 0.476341842160585E+00
c_fit_STO( 7, 9,11)= 0.143911544755442E+00
g_fit_STO( 7, 9,11)= 0.226236668351735E+00
c_fit_STO( 8, 9,11)= 0.624879617226244E-01
g_fit_STO( 8, 9,11)= 0.112167290053600E+00
c_fit_STO( 9, 9,11)= 0.830643824502859E-02
g_fit_STO( 9, 9,11)= 0.569624569000671E-01
! chi2=  0.168097359096930E-08

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1, 9,12)= 0.505312717016304E-01
g_fit_STO( 1, 9,12)= 0.153730840284960E+03
c_fit_STO( 2, 9,12)= 0.727877331557168E-01
g_fit_STO( 2, 9,12)= 0.281922472064823E+02
c_fit_STO( 3, 9,12)= 0.112643389822217E+00
g_fit_STO( 3, 9,12)= 0.791342379945245E+01
c_fit_STO( 4, 9,12)= 0.153971654902271E+00
g_fit_STO( 4, 9,12)= 0.270226542968245E+01
c_fit_STO( 5, 9,12)= 0.186214084900455E+00
g_fit_STO( 5, 9,12)= 0.107598019747681E+01
c_fit_STO( 6, 9,12)= 0.192735351238375E+00
g_fit_STO( 6, 9,12)= 0.471893319031172E+00
c_fit_STO( 7, 9,12)= 0.144903758866476E+00
g_fit_STO( 7, 9,12)= 0.219793702730759E+00
c_fit_STO( 8, 9,12)= 0.583572458699362E-01
g_fit_STO( 8, 9,12)= 0.107443318513022E+00
c_fit_STO( 9, 9,12)= 0.684416859300382E-02
g_fit_STO( 9, 9,12)= 0.544062404981367E-01
! chi2=  0.154498003224961E-08

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1, 9,13)= 0.505560869640932E-01
g_fit_STO( 1, 9,13)= 0.153730647700436E+03
c_fit_STO( 2, 9,13)= 0.728461049722356E-01
g_fit_STO( 2, 9,13)= 0.281927282918219E+02
c_fit_STO( 3, 9,13)= 0.112687895901111E+00
g_fit_STO( 3, 9,13)= 0.791212087004144E+01
c_fit_STO( 4, 9,13)= 0.153905662589566E+00
g_fit_STO( 4, 9,13)= 0.270337773129102E+01
c_fit_STO( 5, 9,13)= 0.186480073853889E+00
g_fit_STO( 5, 9,13)= 0.107622750777727E+01
c_fit_STO( 6, 9,13)= 0.193403183381911E+00
g_fit_STO( 6, 9,13)= 0.470866225809244E+00
c_fit_STO( 7, 9,13)= 0.144877674018295E+00
g_fit_STO( 7, 9,13)= 0.218545112770722E+00
c_fit_STO( 8, 9,13)= 0.576311123923783E-01
g_fit_STO( 8, 9,13)= 0.106453389786085E+00
c_fit_STO( 9, 9,13)= 0.658724746573950E-02
g_fit_STO( 9, 9,13)= 0.536359277896771E-01
! chi2=  0.142906257779060E-08

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1, 9,14)= 0.511706110962377E-01
g_fit_STO( 1, 9,14)= 0.153803431049567E+03
c_fit_STO( 2, 9,14)= 0.715765080269644E-01
g_fit_STO( 2, 9,14)= 0.280443490923773E+02
c_fit_STO( 3, 9,14)= 0.108397744006361E+00
g_fit_STO( 3, 9,14)= 0.818090099629477E+01
c_fit_STO( 4, 9,14)= 0.155590088793801E+00
g_fit_STO( 4, 9,14)= 0.279389531848017E+01
c_fit_STO( 5, 9,14)= 0.193733401368956E+00
g_fit_STO( 5, 9,14)= 0.107640617752023E+01
c_fit_STO( 6, 9,14)= 0.196317213897122E+00
g_fit_STO( 6, 9,14)= 0.459662589156946E+00
c_fit_STO( 7, 9,14)= 0.141317177187076E+00
g_fit_STO( 7, 9,14)= 0.212677006765643E+00
c_fit_STO( 8, 9,14)= 0.548032929281923E-01
g_fit_STO( 8, 9,14)= 0.104322812317003E+00
c_fit_STO( 9, 9,14)= 0.621536081929700E-02
g_fit_STO( 9, 9,14)= 0.528212829444768E-01
! chi2=  0.131390643308769E-08

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1, 9,15)= 0.505960938670062E-01
g_fit_STO( 1, 9,15)= 0.153730368783552E+03
c_fit_STO( 2, 9,15)= 0.729383704446499E-01
g_fit_STO( 2, 9,15)= 0.281934223853524E+02
c_fit_STO( 3, 9,15)= 0.112756291921593E+00
g_fit_STO( 3, 9,15)= 0.791025439917132E+01
c_fit_STO( 4, 9,15)= 0.153826590652945E+00
g_fit_STO( 4, 9,15)= 0.270495798526229E+01
c_fit_STO( 5, 9,15)= 0.186704677229962E+00
g_fit_STO( 5, 9,15)= 0.107679886554676E+01
c_fit_STO( 6, 9,15)= 0.193919735949324E+00
g_fit_STO( 6, 9,15)= 0.470092726749033E+00
c_fit_STO( 7, 9,15)= 0.144815716325103E+00
g_fit_STO( 7, 9,15)= 0.217460066952631E+00
c_fit_STO( 8, 9,15)= 0.570295576017833E-01
g_fit_STO( 8, 9,15)= 0.105472688847113E+00
c_fit_STO( 9, 9,15)= 0.636692308989764E-02
g_fit_STO( 9, 9,15)= 0.527243658577602E-01
! chi2=  0.124288806385338E-08

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1, 9,16)= 0.506103548517612E-01
g_fit_STO( 1, 9,16)= 0.153730279197867E+03
c_fit_STO( 2, 9,16)= 0.729812165893908E-01
g_fit_STO( 2, 9,16)= 0.281936614444480E+02
c_fit_STO( 3, 9,16)= 0.112765196569048E+00
g_fit_STO( 3, 9,16)= 0.790951278154660E+01
c_fit_STO( 4, 9,16)= 0.153834020103205E+00
g_fit_STO( 4, 9,16)= 0.270590647341818E+01
c_fit_STO( 5, 9,16)= 0.187401956440052E+00
g_fit_STO( 5, 9,16)= 0.107569008782665E+01
c_fit_STO( 6, 9,16)= 0.194946389156425E+00
g_fit_STO( 6, 9,16)= 0.467691933052150E+00
c_fit_STO( 7, 9,16)= 0.144502585077661E+00
g_fit_STO( 7, 9,16)= 0.215430255650477E+00
c_fit_STO( 8, 9,16)= 0.558601412605065E-01
g_fit_STO( 8, 9,16)= 0.104195424978035E+00
c_fit_STO( 9, 9,16)= 0.604320307484852E-02
g_fit_STO( 9, 9,16)= 0.519775864800108E-01
! chi2=  0.116731665365522E-08

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1,10, 2)= 0.364333417933986E-01
g_fit_STO( 1,10, 2)= 0.271820015282626E+03
c_fit_STO( 2,10, 2)= 0.517238316866210E-01
g_fit_STO( 2,10, 2)= 0.498372275754285E+02
c_fit_STO( 3,10, 2)= 0.833167620324990E-01
g_fit_STO( 3,10, 2)= 0.140545294480727E+02
c_fit_STO( 4,10, 2)= 0.126307138339173E+00
g_fit_STO( 4,10, 2)= 0.456974409667037E+01
c_fit_STO( 5,10, 2)= 0.173265184652561E+00
g_fit_STO( 5,10, 2)= 0.168274703439585E+01
c_fit_STO( 6,10, 2)= 0.213864010334963E+00
g_fit_STO( 6,10, 2)= 0.673632648078436E+00
c_fit_STO( 7,10, 2)= 0.162464586504181E+00
g_fit_STO( 7,10, 2)= 0.302653170838521E+00
c_fit_STO( 8,10, 2)= 0.138012269335197E+00
g_fit_STO( 8,10, 2)= 0.159862554206575E+00
c_fit_STO( 9,10, 2)=-0.366481684667828E-03
g_fit_STO( 9,10, 2)=-0.112400575246116E+00
c_fit_STO(10,10, 2)=-0.127478179482958E-17
g_fit_STO(10,10, 2)=-0.160846667483305E+01
! chi2=  0.278893222688635E-08

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1,10, 3)= 0.369906858409908E-01
g_fit_STO( 1,10, 3)= 0.271817333371722E+03
c_fit_STO( 2,10, 3)= 0.535359960250411E-01
g_fit_STO( 2,10, 3)= 0.498505861275032E+02
c_fit_STO( 3,10, 3)= 0.836309803599155E-01
g_fit_STO( 3,10, 3)= 0.139898968623201E+02
c_fit_STO( 4,10, 3)= 0.117667501623925E+00
g_fit_STO( 4,10, 3)= 0.478290173530877E+01
c_fit_STO( 5,10, 3)= 0.151761165052138E+00
g_fit_STO( 5,10, 3)= 0.190663549834360E+01
c_fit_STO( 6,10, 3)= 0.180560503085044E+00
g_fit_STO( 6,10, 3)= 0.840188807161769E+00
c_fit_STO( 7,10, 3)= 0.186634758064703E+00
g_fit_STO( 7,10, 3)= 0.390803301656474E+00
c_fit_STO( 8,10, 3)= 0.137865945764589E+00
g_fit_STO( 8,10, 3)= 0.185840920913083E+00
c_fit_STO( 9,10, 3)= 0.398040125431272E-01
g_fit_STO( 9,10, 3)= 0.901295141212720E-01
c_fit_STO(10,10, 3)=-0.383713803879794E-02
g_fit_STO(10,10, 3)= 0.403631933957666E-01
! chi2=  0.186524381792316E-08

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1,10, 4)= 0.373558691274132E-01
g_fit_STO( 1,10, 4)= 0.271817003181971E+03
c_fit_STO( 2,10, 4)= 0.540685954013229E-01
g_fit_STO( 2,10, 4)= 0.498513422588144E+02
c_fit_STO( 3,10, 4)= 0.844056069597326E-01
g_fit_STO( 3,10, 4)= 0.139880537094673E+02
c_fit_STO( 4,10, 4)= 0.118560187864216E+00
g_fit_STO( 4,10, 4)= 0.478438898634978E+01
c_fit_STO( 5,10, 4)= 0.153223037595101E+00
g_fit_STO( 5,10, 4)= 0.190497732173327E+01
c_fit_STO( 6,10, 4)= 0.182204737870604E+00
g_fit_STO( 6,10, 4)= 0.834822416733234E+00
c_fit_STO( 7,10, 4)= 0.183276587364476E+00
g_fit_STO( 7,10, 4)= 0.387352224232799E+00
c_fit_STO( 8,10, 4)= 0.129628890020114E+00
g_fit_STO( 8,10, 4)= 0.187125848236241E+00
c_fit_STO( 9,10, 4)= 0.422055462166468E-01
g_fit_STO( 9,10, 4)= 0.935085528165876E-01
c_fit_STO(10,10, 4)=-0.468717574018206E-03
g_fit_STO(10,10, 4)= 0.398384125923054E-01
! chi2=  0.142879866384062E-08

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1,10, 5)= 0.374449054684824E-01
g_fit_STO( 1,10, 5)= 0.271813865014492E+03
c_fit_STO( 2,10, 5)= 0.546457251668724E-01
g_fit_STO( 2,10, 5)= 0.498593593471393E+02
c_fit_STO( 3,10, 5)= 0.843450248054067E-01
g_fit_STO( 3,10, 5)= 0.139650436012313E+02
c_fit_STO( 4,10, 5)= 0.118958219339375E+00
g_fit_STO( 4,10, 5)= 0.481111697061379E+01
c_fit_STO( 5,10, 5)= 0.158611626653545E+00
g_fit_STO( 5,10, 5)= 0.188722421284838E+01
c_fit_STO( 6,10, 5)= 0.189056287302789E+00
g_fit_STO( 6,10, 5)= 0.806026686518541E+00
c_fit_STO( 7,10, 5)= 0.185076528868293E+00
g_fit_STO( 7,10, 5)= 0.366200365695802E+00
c_fit_STO( 8,10, 5)= 0.121588977808985E+00
g_fit_STO( 8,10, 5)= 0.174391531063755E+00
c_fit_STO( 9,10, 5)= 0.343541461444802E-01
g_fit_STO( 9,10, 5)= 0.869141163849665E-01
c_fit_STO(10,10, 5)= 0.260500630848705E-03
g_fit_STO(10,10, 5)= 0.505960917097343E-01
! chi2=  0.116634136662028E-08

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1,10, 6)= 0.375739290587969E-01
g_fit_STO( 1,10, 6)= 0.271813855835026E+03
c_fit_STO( 2,10, 6)= 0.548366034476913E-01
g_fit_STO( 2,10, 6)= 0.498593823135524E+02
c_fit_STO( 3,10, 6)= 0.846178079447931E-01
g_fit_STO( 3,10, 6)= 0.139649821625838E+02
c_fit_STO( 4,10, 6)= 0.119283457803994E+00
g_fit_STO( 4,10, 6)= 0.481114974445365E+01
c_fit_STO( 5,10, 6)= 0.158844652486337E+00
g_fit_STO( 5,10, 6)= 0.188741987513768E+01
c_fit_STO( 6,10, 6)= 0.188825272721429E+00
g_fit_STO( 6,10, 6)= 0.806247645716323E+00
c_fit_STO( 7,10, 6)= 0.183373811704753E+00
g_fit_STO( 7,10, 6)= 0.366878574043902E+00
c_fit_STO( 8,10, 6)= 0.119287155854262E+00
g_fit_STO( 8,10, 6)= 0.176026558878117E+00
c_fit_STO( 9,10, 6)= 0.353483170770445E-01
g_fit_STO( 9,10, 6)= 0.899156535092783E-01
c_fit_STO(10,10, 6)= 0.229526350047941E-02
g_fit_STO(10,10, 6)= 0.541267450050075E-01
! chi2=  0.977775940904463E-09

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1,10, 7)= 0.377905254182946E-01
g_fit_STO( 1,10, 7)= 0.271819176820010E+03
c_fit_STO( 2,10, 7)= 0.546635621072010E-01
g_fit_STO( 2,10, 7)= 0.498460630190205E+02
c_fit_STO( 3,10, 7)= 0.852145614087121E-01
g_fit_STO( 3,10, 7)= 0.140013305154493E+02
c_fit_STO( 4,10, 7)= 0.120400862750583E+00
g_fit_STO( 4,10, 7)= 0.477810331493838E+01
c_fit_STO( 5,10, 7)= 0.154375489194665E+00
g_fit_STO( 5,10, 7)= 0.189313838158101E+01
c_fit_STO( 6,10, 7)= 0.174050287289156E+00
g_fit_STO( 6,10, 7)= 0.842132918631122E+00
c_fit_STO( 7,10, 7)= 0.165685253965603E+00
g_fit_STO( 7,10, 7)= 0.411126663976126E+00
c_fit_STO( 8,10, 7)= 0.124233103129476E+00
g_fit_STO( 8,10, 7)= 0.213393976136630E+00
c_fit_STO( 9,10, 7)= 0.579718406970210E-01
g_fit_STO( 9,10, 7)= 0.114351171175098E+00
c_fit_STO(10,10, 7)= 0.989818354992049E-02
g_fit_STO(10,10, 7)= 0.630791544954813E-01
! chi2=  0.835370216920314E-09

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1,10, 8)= 0.378321691282427E-01
g_fit_STO( 1,10, 8)= 0.271817496555153E+03
c_fit_STO( 2,10, 8)= 0.548237405533255E-01
g_fit_STO( 2,10, 8)= 0.498501754876706E+02
c_fit_STO( 3,10, 8)= 0.852977885579465E-01
g_fit_STO( 3,10, 8)= 0.139906371964734E+02
c_fit_STO( 4,10, 8)= 0.120223878899364E+00
g_fit_STO( 4,10, 8)= 0.478566204347090E+01
c_fit_STO( 5,10, 8)= 0.154835707222001E+00
g_fit_STO( 5,10, 8)= 0.189612841236587E+01
c_fit_STO( 6,10, 8)= 0.175728903577281E+00
g_fit_STO( 6,10, 8)= 0.839474039143671E+00
c_fit_STO( 7,10, 8)= 0.168026006977224E+00
g_fit_STO( 7,10, 8)= 0.405944985875478E+00
c_fit_STO( 8,10, 8)= 0.124308301907002E+00
g_fit_STO( 8,10, 8)= 0.207866134245094E+00
c_fit_STO( 9,10, 8)= 0.549819520497401E-01
g_fit_STO( 9,10, 8)= 0.109723471829634E+00
c_fit_STO(10,10, 8)= 0.819076854932120E-02
g_fit_STO(10,10, 8)= 0.594401240745767E-01
! chi2=  0.734519163620350E-09

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1,10, 9)= 0.379046428012710E-01
g_fit_STO( 1,10, 9)= 0.271818058578408E+03
c_fit_STO( 2,10, 9)= 0.548476315246198E-01
g_fit_STO( 2,10, 9)= 0.498487310879255E+02
c_fit_STO( 3,10, 9)= 0.854880469543060E-01
g_fit_STO( 3,10, 9)= 0.139948050975355E+02
c_fit_STO( 4,10, 9)= 0.120399894683470E+00
g_fit_STO( 4,10, 9)= 0.478074769451390E+01
c_fit_STO( 5,10, 9)= 0.154148634070917E+00
g_fit_STO( 5,10, 9)= 0.189867962429124E+01
c_fit_STO( 6,10, 9)= 0.175720382060755E+00
g_fit_STO( 6,10, 9)= 0.842638373140349E+00
c_fit_STO( 7,10, 9)= 0.170387258430982E+00
g_fit_STO( 7,10, 9)= 0.405025605575434E+00
c_fit_STO( 8,10, 9)= 0.125239984827278E+00
g_fit_STO( 8,10, 9)= 0.204589916822322E+00
c_fit_STO( 9,10, 9)= 0.529418132918647E-01
g_fit_STO( 9,10, 9)= 0.106561927786185E+00
c_fit_STO(10,10, 9)= 0.715595638893158E-02
g_fit_STO(10,10, 9)= 0.569346483637941E-01
! chi2=  0.654941093445663E-09

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1,10,10)= 0.379547534642940E-01
g_fit_STO( 1,10,10)= 0.271818121723946E+03
c_fit_STO( 2,10,10)= 0.548849779985982E-01
g_fit_STO( 2,10,10)= 0.498485424953649E+02
c_fit_STO( 3,10,10)= 0.856159626165205E-01
g_fit_STO( 3,10,10)= 0.139955018566255E+02
c_fit_STO( 4,10,10)= 0.120456510536697E+00
g_fit_STO( 4,10,10)= 0.477920674882537E+01
c_fit_STO( 5,10,10)= 0.153889146980985E+00
g_fit_STO( 5,10,10)= 0.190066861711027E+01
c_fit_STO( 6,10,10)= 0.176366602182960E+00
g_fit_STO( 6,10,10)= 0.843203115272399E+00
c_fit_STO( 7,10,10)= 0.172127434909406E+00
g_fit_STO( 7,10,10)= 0.402815787752311E+00
c_fit_STO( 8,10,10)= 0.125234027693663E+00
g_fit_STO( 8,10,10)= 0.201633961254656E+00
c_fit_STO( 9,10,10)= 0.512091570945083E-01
g_fit_STO( 9,10,10)= 0.104189366598175E+00
c_fit_STO(10,10,10)= 0.648201575496378E-02
g_fit_STO(10,10,10)= 0.551408238892172E-01
! chi2=  0.591175360261563E-09

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1,10,11)= 0.379935184034003E-01
g_fit_STO( 1,10,11)= 0.271818143517875E+03
c_fit_STO( 2,10,11)= 0.549202193952751E-01
g_fit_STO( 2,10,11)= 0.498484659378551E+02
c_fit_STO( 3,10,11)= 0.857095763202123E-01
g_fit_STO( 3,10,11)= 0.139958374213086E+02
c_fit_STO( 4,10,11)= 0.120511580147529E+00
g_fit_STO( 4,10,11)= 0.477827271528501E+01
c_fit_STO( 5,10,11)= 0.153698581323740E+00
g_fit_STO( 5,10,11)= 0.190204287802948E+01
c_fit_STO( 6,10,11)= 0.176706222394851E+00
g_fit_STO( 6,10,11)= 0.843848169907045E+00
c_fit_STO( 7,10,11)= 0.173408143714906E+00
g_fit_STO( 7,10,11)= 0.401446505609719E+00
c_fit_STO( 8,10,11)= 0.125405802126682E+00
g_fit_STO( 8,10,11)= 0.199464545904921E+00
c_fit_STO( 9,10,11)= 0.498891277013428E-01
g_fit_STO( 9,10,11)= 0.102309732499202E+00
c_fit_STO(10,10,11)= 0.596639768848767E-02
g_fit_STO(10,10,11)= 0.537084219799885E-01
! chi2=  0.538880243380453E-09

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1,10,12)= 0.379186547249306E-01
g_fit_STO( 1,10,12)= 0.271813993708729E+03
c_fit_STO( 2,10,12)= 0.552126516977755E-01
g_fit_STO( 2,10,12)= 0.498590076354136E+02
c_fit_STO( 3,10,12)= 0.853944286901143E-01
g_fit_STO( 3,10,12)= 0.139661999315936E+02
c_fit_STO( 4,10,12)= 0.119831422981300E+00
g_fit_STO( 4,10,12)= 0.480900500229373E+01
c_fit_STO( 5,10,12)= 0.159093819797872E+00
g_fit_STO( 5,10,12)= 0.189109047938155E+01
c_fit_STO( 6,10,12)= 0.187571803311329E+00
g_fit_STO( 6,10,12)= 0.808617297044669E+00
c_fit_STO( 7,10,12)= 0.177498690160842E+00
g_fit_STO( 7,10,12)= 0.371628030973359E+00
c_fit_STO( 8,10,12)= 0.116613323306405E+00
g_fit_STO( 8,10,12)= 0.182199065166767E+00
c_fit_STO( 9,10,12)= 0.407920144845762E-01
g_fit_STO( 9,10,12)= 0.941250001335617E-01
c_fit_STO(10,10,12)= 0.424317857864377E-02
g_fit_STO(10,10,12)= 0.502496188310433E-01
! chi2=  0.497144423821025E-09

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1,10,13)= 0.379425565395577E-01
g_fit_STO( 1,10,13)= 0.271813991086492E+03
c_fit_STO( 2,10,13)= 0.552454301778932E-01
g_fit_STO( 2,10,13)= 0.498590149998920E+02
c_fit_STO( 3,10,13)= 0.854402303386470E-01
g_fit_STO( 3,10,13)= 0.139661735340010E+02
c_fit_STO( 4,10,13)= 0.119889046977148E+00
g_fit_STO( 4,10,13)= 0.480906066837524E+01
c_fit_STO( 5,10,13)= 0.159150385052749E+00
g_fit_STO( 5,10,13)= 0.189098942534569E+01
c_fit_STO( 6,10,13)= 0.187540280800066E+00
g_fit_STO( 6,10,13)= 0.808549057070620E+00
c_fit_STO( 7,10,13)= 0.177399644548739E+00
g_fit_STO( 7,10,13)= 0.371573172520444E+00
c_fit_STO( 8,10,13)= 0.116597529787627E+00
g_fit_STO( 8,10,13)= 0.182036964652188E+00
c_fit_STO( 9,10,13)= 0.407613166231214E-01
g_fit_STO( 9,10,13)= 0.938014052945784E-01
c_fit_STO(10,10,13)= 0.419502238343623E-02
g_fit_STO(10,10,13)= 0.497473951491483E-01
! chi2=  0.459928343044629E-09

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1,10,14)= 0.379634407185959E-01
g_fit_STO( 1,10,14)= 0.271813989578009E+03
c_fit_STO( 2,10,14)= 0.552726939926273E-01
g_fit_STO( 2,10,14)= 0.498590196153873E+02
c_fit_STO( 3,10,14)= 0.854804490706811E-01
g_fit_STO( 3,10,14)= 0.139661541890766E+02
c_fit_STO( 4,10,14)= 0.119935633394027E+00
g_fit_STO( 4,10,14)= 0.480911085673328E+01
c_fit_STO( 5,10,14)= 0.159213928818947E+00
g_fit_STO( 5,10,14)= 0.189087539373958E+01
c_fit_STO( 6,10,14)= 0.187566404306565E+00
g_fit_STO( 6,10,14)= 0.808354480470036E+00
c_fit_STO( 7,10,14)= 0.177356159786474E+00
g_fit_STO( 7,10,14)= 0.371383348829241E+00
c_fit_STO( 8,10,14)= 0.116555046951909E+00
g_fit_STO( 8,10,14)= 0.181791078149001E+00
c_fit_STO( 9,10,14)= 0.406695057713893E-01
g_fit_STO( 9,10,14)= 0.934626432218299E-01
c_fit_STO(10,10,14)= 0.414116928439515E-02
g_fit_STO(10,10,14)= 0.493055864983896E-01
! chi2=  0.427994390195591E-09

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1,10,15)= 0.379833003198874E-01
g_fit_STO( 1,10,15)= 0.271813928588678E+03
c_fit_STO( 2,10,15)= 0.552975423750307E-01
g_fit_STO( 2,10,15)= 0.498591814758633E+02
c_fit_STO( 3,10,15)= 0.855083851658912E-01
g_fit_STO( 3,10,15)= 0.139656530752096E+02
c_fit_STO( 4,10,15)= 0.119970466280773E+00
g_fit_STO( 4,10,15)= 0.480986271762338E+01
c_fit_STO( 5,10,15)= 0.159452897943361E+00
g_fit_STO( 5,10,15)= 0.189000426227313E+01
c_fit_STO( 6,10,15)= 0.187561468010651E+00
g_fit_STO( 6,10,15)= 0.807511917099468E+00
c_fit_STO( 7,10,15)= 0.177479331775249E+00
g_fit_STO( 7,10,15)= 0.370975488972119E+00
c_fit_STO( 8,10,15)= 0.117087437538150E+00
g_fit_STO( 8,10,15)= 0.180890936708305E+00
c_fit_STO( 9,10,15)= 0.400358938339027E-01
g_fit_STO( 9,10,15)= 0.921674064638925E-01
c_fit_STO(10,10,15)= 0.377301467812092E-02
g_fit_STO(10,10,15)= 0.481559116353041E-01
! chi2=  0.400414895047623E-09

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1,10,16)= 0.379984275609729E-01
g_fit_STO( 1,10,16)= 0.271813988201305E+03
c_fit_STO( 2,10,16)= 0.553171639866071E-01
g_fit_STO( 2,10,16)= 0.498590226137749E+02
c_fit_STO( 3,10,16)= 0.855433483973411E-01
g_fit_STO( 3,10,16)= 0.139661477425344E+02
c_fit_STO( 4,10,16)= 0.120017505822961E+00
g_fit_STO( 4,10,16)= 0.480910925744784E+01
c_fit_STO( 5,10,16)= 0.159227498083879E+00
g_fit_STO( 5,10,16)= 0.189092913807730E+01
c_fit_STO( 6,10,16)= 0.187335852482526E+00
g_fit_STO( 6,10,16)= 0.808774043635427E+00
c_fit_STO( 7,10,16)= 0.177072149484402E+00
g_fit_STO( 7,10,16)= 0.371828154465344E+00
c_fit_STO( 8,10,16)= 0.116634539683650E+00
g_fit_STO( 8,10,16)= 0.181961815351716E+00
c_fit_STO( 9,10,16)= 0.408517228548043E-01
g_fit_STO( 9,10,16)= 0.932876435697782E-01
c_fit_STO(10,10,16)= 0.414604714841767E-02
g_fit_STO(10,10,16)= 0.488129955893861E-01
! chi2=  0.375994658938365E-09

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1,11, 2)= 0.273607727269749E-01
g_fit_STO( 1,11, 2)= 0.469213141270166E+03
c_fit_STO( 2,11, 2)= 0.403753800584972E-01
g_fit_STO( 2,11, 2)= 0.860696598462916E+02
c_fit_STO( 3,11, 2)= 0.617780932354256E-01
g_fit_STO( 3,11, 2)= 0.241038670659445E+02
c_fit_STO( 4,11, 2)= 0.904132960985353E-01
g_fit_STO( 4,11, 2)= 0.833432196292942E+01
c_fit_STO( 5,11, 2)= 0.131715294511785E+00
g_fit_STO( 5,11, 2)= 0.316262734230733E+01
c_fit_STO( 6,11, 2)= 0.178633619601403E+00
g_fit_STO( 6,11, 2)= 0.127896627886985E+01
c_fit_STO( 7,11, 2)= 0.236449083511917E+00
g_fit_STO( 7,11, 2)= 0.519472579105876E+00
c_fit_STO( 8,11, 2)= 0.221817610335168E+00
g_fit_STO( 8,11, 2)= 0.187566818229053E+00
c_fit_STO( 9,11, 2)=-0.500654533064752E-04
g_fit_STO( 9,11, 2)=-0.208916771478137E+00
c_fit_STO(10,11, 2)= 0.277528282178491E-07
g_fit_STO(10,11, 2)=-0.506434033011361E+00
c_fit_STO(11,11, 2)=-0.191335769675261E-15
g_fit_STO(11,11, 2)=-0.138617567668132E+01
! chi2=  0.994704745107023E-09

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1,11, 3)= 0.280879512716063E-01
g_fit_STO( 1,11, 3)= 0.469212671472708E+03
c_fit_STO( 2,11, 3)= 0.410053986362845E-01
g_fit_STO( 2,11, 3)= 0.860701841595836E+02
c_fit_STO( 3,11, 3)= 0.637547900997423E-01
g_fit_STO( 3,11, 3)= 0.241080715825976E+02
c_fit_STO( 4,11, 3)= 0.909289306615779E-01
g_fit_STO( 4,11, 3)= 0.830581722009775E+01
c_fit_STO( 5,11, 3)= 0.125645975412673E+00
g_fit_STO( 5,11, 3)= 0.326248321890073E+01
c_fit_STO( 6,11, 3)= 0.160192782949515E+00
g_fit_STO( 6,11, 3)= 0.139463999613079E+01
c_fit_STO( 7,11, 3)= 0.178375561001870E+00
g_fit_STO( 7,11, 3)= 0.644487623673799E+00
c_fit_STO( 8,11, 3)= 0.167050352095851E+00
g_fit_STO( 8,11, 3)= 0.316663012419813E+00
c_fit_STO( 9,11, 3)= 0.110327455054417E+00
g_fit_STO( 9,11, 3)= 0.159414992347027E+00
c_fit_STO(10,11, 3)= 0.255146498526374E-01
g_fit_STO(10,11, 3)= 0.839428909521841E-01
c_fit_STO(11,11, 3)=-0.261501842749168E-02
g_fit_STO(11,11, 3)= 0.338209423087832E-01
! chi2=  0.629119278924405E-09

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1,11, 4)= 0.283705866214372E-01
g_fit_STO( 1,11, 4)= 0.469212732852540E+03
c_fit_STO( 2,11, 4)= 0.414045154931424E-01
g_fit_STO( 2,11, 4)= 0.860700190156572E+02
c_fit_STO( 3,11, 4)= 0.643782105683443E-01
g_fit_STO( 3,11, 4)= 0.241085948292641E+02
c_fit_STO( 4,11, 4)= 0.917282166242007E-01
g_fit_STO( 4,11, 4)= 0.830497270348309E+01
c_fit_STO( 5,11, 4)= 0.126467057481320E+00
g_fit_STO( 5,11, 4)= 0.326339757183598E+01
c_fit_STO( 6,11, 4)= 0.160812991367464E+00
g_fit_STO( 6,11, 4)= 0.139535227559115E+01
c_fit_STO( 7,11, 4)= 0.178284585447593E+00
g_fit_STO( 7,11, 4)= 0.643875802726905E+00
c_fit_STO( 8,11, 4)= 0.162234734206706E+00
g_fit_STO( 8,11, 4)= 0.316925245951045E+00
c_fit_STO( 9,11, 4)= 0.104821077937543E+00
g_fit_STO( 9,11, 4)= 0.162849601629746E+00
c_fit_STO(10,11, 4)= 0.311296482505416E-01
g_fit_STO(10,11, 4)= 0.850262075608474E-01
c_fit_STO(11,11, 4)=-0.147882530098234E-02
g_fit_STO(11,11, 4)= 0.482612546779200E-01
! chi2=  0.481488091336943E-09

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1,11, 5)= 0.285128716754147E-01
g_fit_STO( 1,11, 5)= 0.469212635482527E+03
c_fit_STO( 2,11, 5)= 0.416649443654336E-01
g_fit_STO( 2,11, 5)= 0.860702880119326E+02
c_fit_STO( 3,11, 5)= 0.646483270827210E-01
g_fit_STO( 3,11, 5)= 0.241076865291614E+02
c_fit_STO( 4,11, 5)= 0.923138580927055E-01
g_fit_STO( 4,11, 5)= 0.830679426991733E+01
c_fit_STO( 5,11, 5)= 0.126900713394516E+00
g_fit_STO( 5,11, 5)= 0.326006841212306E+01
c_fit_STO( 6,11, 5)= 0.160348502205334E+00
g_fit_STO( 6,11, 5)= 0.139801226641748E+01
c_fit_STO( 7,11, 5)= 0.179729405347489E+00
g_fit_STO( 7,11, 5)= 0.644286056588878E+00
c_fit_STO( 8,11, 5)= 0.165287232502601E+00
g_fit_STO( 8,11, 5)= 0.312258751717699E+00
c_fit_STO( 9,11, 5)= 0.101878835081592E+00
g_fit_STO( 9,11, 5)= 0.156857664172260E+00
c_fit_STO(10,11, 5)= 0.268880090062041E-01
g_fit_STO(10,11, 5)= 0.814989787986982E-01
c_fit_STO(11,11, 5)=-0.892199135536774E-04
g_fit_STO(11,11, 5)= 0.429820918043915E-01
! chi2=  0.390624031826428E-09

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1,11, 6)= 0.286155888126911E-01
g_fit_STO( 1,11, 6)= 0.469212634115866E+03
c_fit_STO( 2,11, 6)= 0.418052529120604E-01
g_fit_STO( 2,11, 6)= 0.860702917838703E+02
c_fit_STO( 3,11, 6)= 0.648735771455980E-01
g_fit_STO( 3,11, 6)= 0.241076733832704E+02
c_fit_STO( 4,11, 6)= 0.925838877810058E-01
g_fit_STO( 4,11, 6)= 0.830682036434900E+01
c_fit_STO( 5,11, 6)= 0.127240416582325E+00
g_fit_STO( 5,11, 6)= 0.326002616951491E+01
c_fit_STO( 6,11, 6)= 0.160395575141533E+00
g_fit_STO( 6,11, 6)= 0.139807706360855E+01
c_fit_STO( 7,11, 6)= 0.179374593274262E+00
g_fit_STO( 7,11, 6)= 0.644752165311025E+00
c_fit_STO( 8,11, 6)= 0.164875450691838E+00
g_fit_STO( 8,11, 6)= 0.312101647308304E+00
c_fit_STO( 9,11, 6)= 0.100813035751459E+00
g_fit_STO( 9,11, 6)= 0.156332857490707E+00
c_fit_STO(10,11, 6)= 0.267450397577599E-01
g_fit_STO(10,11, 6)= 0.813504545818877E-01
c_fit_STO(11,11, 6)= 0.720390629356336E-03
g_fit_STO(11,11, 6)= 0.473452408353686E-01
! chi2=  0.328705127012219E-09

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1,11, 7)= 0.286842259047866E-01
g_fit_STO( 1,11, 7)= 0.469212631856625E+03
c_fit_STO( 2,11, 7)= 0.419095233994608E-01
g_fit_STO( 2,11, 7)= 0.860702982877909E+02
c_fit_STO( 3,11, 7)= 0.650152578619135E-01
g_fit_STO( 3,11, 7)= 0.241076490701440E+02
c_fit_STO( 4,11, 7)= 0.927994188367530E-01
g_fit_STO( 4,11, 7)= 0.830688224195308E+01
c_fit_STO( 5,11, 7)= 0.127416818931856E+00
g_fit_STO( 5,11, 7)= 0.325985115156131E+01
c_fit_STO( 6,11, 7)= 0.160524783053880E+00
g_fit_STO( 6,11, 7)= 0.139834579232848E+01
c_fit_STO( 7,11, 7)= 0.179529591451974E+00
g_fit_STO( 7,11, 7)= 0.644260440971463E+00
c_fit_STO( 8,11, 7)= 0.163923946166526E+00
g_fit_STO( 8,11, 7)= 0.311646077402172E+00
c_fit_STO( 9,11, 7)= 0.993926688618606E-01
g_fit_STO( 9,11, 7)= 0.156657840504707E+00
c_fit_STO(10,11, 7)= 0.273243654531582E-01
g_fit_STO(10,11, 7)= 0.824449128423240E-01
c_fit_STO(11,11, 7)= 0.149353704234736E-02
g_fit_STO(11,11, 7)= 0.481670133313603E-01
! chi2=  0.282947854139673E-09

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1,11, 8)= 0.287333952584483E-01
g_fit_STO( 1,11, 8)= 0.469212629926848E+03
c_fit_STO( 2,11, 8)= 0.419890209716301E-01
g_fit_STO( 2,11, 8)= 0.860703038496230E+02
c_fit_STO( 3,11, 8)= 0.651136577514927E-01
g_fit_STO( 3,11, 8)= 0.241076282445095E+02
c_fit_STO( 4,11, 8)= 0.929678906509514E-01
g_fit_STO( 4,11, 8)= 0.830693535449999E+01
c_fit_STO( 5,11, 8)= 0.127524354141081E+00
g_fit_STO( 5,11, 8)= 0.325969313095813E+01
c_fit_STO( 6,11, 8)= 0.160674317777131E+00
g_fit_STO( 6,11, 8)= 0.139855478260949E+01
c_fit_STO( 7,11, 8)= 0.179788231790906E+00
g_fit_STO( 7,11, 8)= 0.643547026948668E+00
c_fit_STO( 8,11, 8)= 0.163159861130436E+00
g_fit_STO( 8,11, 8)= 0.311070479725899E+00
c_fit_STO( 9,11, 8)= 0.983640429434359E-01
g_fit_STO( 9,11, 8)= 0.156785011982850E+00
c_fit_STO(10,11, 8)= 0.277926718704556E-01
g_fit_STO(10,11, 8)= 0.829461662625075E-01
c_fit_STO(11,11, 8)= 0.188538270286303E-02
g_fit_STO(11,11, 8)= 0.476864318698059E-01
! chi2=  0.248496583563066E-09

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1,11, 9)= 0.287741954074776E-01
g_fit_STO( 1,11, 9)= 0.469212770020217E+03
c_fit_STO( 2,11, 9)= 0.420417362962923E-01
g_fit_STO( 2,11, 9)= 0.860699270987487E+02
c_fit_STO( 3,11, 9)= 0.652074723780379E-01
g_fit_STO( 3,11, 9)= 0.241088227851921E+02
c_fit_STO( 4,11, 9)= 0.930815496468810E-01
g_fit_STO( 4,11, 9)= 0.830498872537525E+01
c_fit_STO( 5,11, 9)= 0.127408237729259E+00
g_fit_STO( 5,11, 9)= 0.326201102090489E+01
c_fit_STO( 6,11, 9)= 0.161201856741106E+00
g_fit_STO( 6,11, 9)= 0.139917472149798E+01
c_fit_STO( 7,11, 9)= 0.180038798071935E+00
g_fit_STO( 7,11, 9)= 0.641902832031016E+00
c_fit_STO( 8,11, 9)= 0.160086519211063E+00
g_fit_STO( 8,11, 9)= 0.311721533907121E+00
c_fit_STO( 9,11, 9)= 0.963093119942330E-01
g_fit_STO( 9,11, 9)= 0.160294088542455E+00
c_fit_STO(10,11, 9)= 0.306086378517084E-01
g_fit_STO(10,11, 9)= 0.873997558486979E-01
c_fit_STO(11,11, 9)= 0.321910478289301E-02
g_fit_STO(11,11, 9)= 0.504469659446350E-01
! chi2=  0.221389483755105E-09

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1,11,10)= 0.288082725573505E-01
g_fit_STO( 1,11,10)= 0.469212780124908E+03
c_fit_STO( 2,11,10)= 0.420825710035381E-01
g_fit_STO( 2,11,10)= 0.860698994512498E+02
c_fit_STO( 3,11,10)= 0.652805677412240E-01
g_fit_STO( 3,11,10)= 0.241089128461099E+02
c_fit_STO( 4,11,10)= 0.931576096717940E-01
g_fit_STO( 4,11,10)= 0.830482785816867E+01
c_fit_STO( 5,11,10)= 0.127500221318391E+00
g_fit_STO( 5,11,10)= 0.326225250538408E+01
c_fit_STO( 6,11,10)= 0.161207653481489E+00
g_fit_STO( 6,11,10)= 0.139922716748805E+01
c_fit_STO( 7,11,10)= 0.179588739139111E+00
g_fit_STO( 7,11,10)= 0.642399025558532E+00
c_fit_STO( 8,11,10)= 0.159459246866045E+00
g_fit_STO( 8,11,10)= 0.312517867357767E+00
c_fit_STO( 9,11,10)= 0.965004196839202E-01
g_fit_STO( 9,11,10)= 0.160884200603241E+00
c_fit_STO(10,11,10)= 0.311603032196817E-01
g_fit_STO(10,11,10)= 0.873547290573213E-01
c_fit_STO(11,11,10)= 0.322066018870887E-02
g_fit_STO(11,11,10)= 0.494860305660361E-01
! chi2=  0.199899858057072E-09

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1,11,11)= 0.288380877770838E-01
g_fit_STO( 1,11,11)= 0.469212807892918E+03
c_fit_STO( 2,11,11)= 0.421110689019492E-01
g_fit_STO( 2,11,11)= 0.860698235376546E+02
c_fit_STO( 3,11,11)= 0.653494534959637E-01
g_fit_STO( 3,11,11)= 0.241091613899883E+02
c_fit_STO( 4,11,11)= 0.932017409622554E-01
g_fit_STO( 4,11,11)= 0.830437102560962E+01
c_fit_STO( 5,11,11)= 0.127573069486442E+00
g_fit_STO( 5,11,11)= 0.326296649465507E+01
c_fit_STO( 6,11,11)= 0.161371223818404E+00
g_fit_STO( 6,11,11)= 0.139885314495075E+01
c_fit_STO( 7,11,11)= 0.179290652264565E+00
g_fit_STO( 7,11,11)= 0.642254723667091E+00
c_fit_STO( 8,11,11)= 0.158775663221875E+00
g_fit_STO( 8,11,11)= 0.312978250273473E+00
c_fit_STO( 9,11,11)= 0.965534494138392E-01
g_fit_STO( 9,11,11)= 0.161397315907165E+00
c_fit_STO(10,11,11)= 0.316397614629122E-01
g_fit_STO(10,11,11)= 0.874094954131764E-01
c_fit_STO(11,11,11)= 0.325382149614680E-02
g_fit_STO(11,11,11)= 0.488605403113902E-01
! chi2=  0.182271985544239E-09

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1,11,12)= 0.288649227829620E-01
g_fit_STO( 1,11,12)= 0.469212842613603E+03
c_fit_STO( 2,11,12)= 0.421306427913971E-01
g_fit_STO( 2,11,12)= 0.860697285742866E+02
c_fit_STO( 3,11,12)= 0.654142959141990E-01
g_fit_STO( 3,11,12)= 0.241094724631612E+02
c_fit_STO( 4,11,12)= 0.932237569901393E-01
g_fit_STO( 4,11,12)= 0.830379757307686E+01
c_fit_STO( 5,11,12)= 0.127639962562731E+00
g_fit_STO( 5,11,12)= 0.326387151244094E+01
c_fit_STO( 6,11,12)= 0.161558088978627E+00
g_fit_STO( 6,11,12)= 0.139833747433329E+01
c_fit_STO( 7,11,12)= 0.178942054639217E+00
g_fit_STO( 7,11,12)= 0.642099171361900E+00
c_fit_STO( 8,11,12)= 0.158119600281670E+00
g_fit_STO( 8,11,12)= 0.313515542275180E+00
c_fit_STO( 9,11,12)= 0.967115385496743E-01
g_fit_STO( 9,11,12)= 0.161922961165292E+00
c_fit_STO(10,11,12)= 0.320694465508542E-01
g_fit_STO(10,11,12)= 0.874614812620735E-01
c_fit_STO(11,11,12)= 0.327768855311400E-02
g_fit_STO(11,11,12)= 0.483727419982476E-01
! chi2=  0.167555890070159E-09

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1,11,13)= 0.288860534995978E-01
g_fit_STO( 1,11,13)= 0.469212851584830E+03
c_fit_STO( 2,11,13)= 0.421514989357394E-01
g_fit_STO( 2,11,13)= 0.860697038212038E+02
c_fit_STO( 3,11,13)= 0.654591033550049E-01
g_fit_STO( 3,11,13)= 0.241095537403898E+02
c_fit_STO( 4,11,13)= 0.932595373661434E-01
g_fit_STO( 4,11,13)= 0.830364539963299E+01
c_fit_STO( 5,11,13)= 0.127698489133852E+00
g_fit_STO( 5,11,13)= 0.326411655663961E+01
c_fit_STO( 6,11,13)= 0.161593859545222E+00
g_fit_STO( 6,11,13)= 0.139821732645425E+01
c_fit_STO( 7,11,13)= 0.178706892217587E+00
g_fit_STO( 7,11,13)= 0.642258204207741E+00
c_fit_STO( 8,11,13)= 0.157894436335080E+00
g_fit_STO( 8,11,13)= 0.313818738718312E+00
c_fit_STO( 9,11,13)= 0.968796200030772E-01
g_fit_STO( 9,11,13)= 0.162005906404781E+00
c_fit_STO(10,11,13)= 0.321758895448452E-01
g_fit_STO(10,11,13)= 0.872169338623426E-01
c_fit_STO(11,11,13)= 0.324144198511155E-02
g_fit_STO(11,11,13)= 0.478565513769346E-01
! chi2=  0.155107122785967E-09

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1,11,14)= 0.289064080099501E-01
g_fit_STO( 1,11,14)= 0.469212871307202E+03
c_fit_STO( 2,11,14)= 0.421647928192580E-01
g_fit_STO( 2,11,14)= 0.860696495401973E+02
c_fit_STO( 3,11,14)= 0.655063304324955E-01
g_fit_STO( 3,11,14)= 0.241097327280268E+02
c_fit_STO( 4,11,14)= 0.932731337683209E-01
g_fit_STO( 4,11,14)= 0.830330610095536E+01
c_fit_STO( 5,11,14)= 0.127758599123442E+00
g_fit_STO( 5,11,14)= 0.326468004972223E+01
c_fit_STO( 6,11,14)= 0.161696932549341E+00
g_fit_STO( 6,11,14)= 0.139782548249965E+01
c_fit_STO( 7,11,14)= 0.178409652800608E+00
g_fit_STO( 7,11,14)= 0.642285795497427E+00
c_fit_STO( 8,11,14)= 0.157604176534268E+00
g_fit_STO( 8,11,14)= 0.314200609956872E+00
c_fit_STO( 9,11,14)= 0.971158217620689E-01
g_fit_STO( 9,11,14)= 0.162165540193958E+00
c_fit_STO(10,11,14)= 0.323007558195921E-01
g_fit_STO(10,11,14)= 0.870171935191064E-01
c_fit_STO(11,11,14)= 0.320666969213488E-02
g_fit_STO(11,11,14)= 0.474169832119274E-01
! chi2=  0.144425651515301E-09

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1,11,15)= 0.289246712888070E-01
g_fit_STO( 1,11,15)= 0.469212884522861E+03
c_fit_STO( 2,11,15)= 0.421756147501539E-01
g_fit_STO( 2,11,15)= 0.860696128248100E+02
c_fit_STO( 3,11,15)= 0.655480256389656E-01
g_fit_STO( 3,11,15)= 0.241098553787603E+02
c_fit_STO( 4,11,15)= 0.932823332208954E-01
g_fit_STO( 4,11,15)= 0.830306216123636E+01
c_fit_STO( 5,11,15)= 0.127820942234945E+00
g_fit_STO( 5,11,15)= 0.326511943041275E+01
c_fit_STO( 6,11,15)= 0.161771281263573E+00
g_fit_STO( 6,11,15)= 0.139743796759228E+01
c_fit_STO( 7,11,15)= 0.178114968844795E+00
g_fit_STO( 7,11,15)= 0.642389696947151E+00
c_fit_STO( 8,11,15)= 0.157460727145069E+00
g_fit_STO( 8,11,15)= 0.314535373085051E+00
c_fit_STO( 9,11,15)= 0.973937545626746E-01
g_fit_STO( 9,11,15)= 0.162176650438842E+00
c_fit_STO(10,11,15)= 0.323035280368103E-01
g_fit_STO(10,11,15)= 0.867034339103074E-01
c_fit_STO(11,11,15)= 0.314478922902643E-02
g_fit_STO(11,11,15)= 0.469659909551585E-01
! chi2=  0.135170844958300E-09

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1,11,16)= 0.289392260186229E-01
g_fit_STO( 1,11,16)= 0.469212885969830E+03
c_fit_STO( 2,11,16)= 0.421889747746546E-01
g_fit_STO( 2,11,16)= 0.860696084847758E+02
c_fit_STO( 3,11,16)= 0.655759553441942E-01
g_fit_STO( 3,11,16)= 0.241098701329636E+02
c_fit_STO( 4,11,16)= 0.933064228504746E-01
g_fit_STO( 4,11,16)= 0.830302842645715E+01
c_fit_STO( 5,11,16)= 0.127863002511130E+00
g_fit_STO( 5,11,16)= 0.326518806531644E+01
c_fit_STO( 6,11,16)= 0.161780445117580E+00
g_fit_STO( 6,11,16)= 0.139736137538922E+01
c_fit_STO( 7,11,16)= 0.177994237405277E+00
g_fit_STO( 7,11,16)= 0.642495224352818E+00
c_fit_STO( 8,11,16)= 0.157448706459992E+00
g_fit_STO( 8,11,16)= 0.314612315174513E+00
c_fit_STO( 9,11,16)= 0.974986086836420E-01
g_fit_STO( 9,11,16)= 0.162044681669190E+00
c_fit_STO(10,11,16)= 0.322424650215152E-01
g_fit_STO(10,11,16)= 0.864167000192990E-01
c_fit_STO(11,11,16)= 0.310012587279693E-02
g_fit_STO(11,11,16)= 0.466099290874181E-01
! chi2=  0.127080192724916E-09

alpha( 2)= 0.944716504255102E+00
 beta( 2)= 0.100000000000000E+00
c_fit_STO( 1,12, 2)= 0.212562761057255E-01
g_fit_STO( 1,12, 2)= 0.793039630749150E+03
c_fit_STO( 2,12, 2)= 0.306198252788432E-01
g_fit_STO( 2,12, 2)= 0.145472390073800E+03
c_fit_STO( 3,12, 2)= 0.487255553005210E-01
g_fit_STO( 3,12, 2)= 0.407503193332239E+02
c_fit_STO( 4,12, 2)= 0.679545590531481E-01
g_fit_STO( 4,12, 2)= 0.140239311102356E+02
c_fit_STO( 5,12, 2)= 0.994632950225887E-01
g_fit_STO( 5,12, 2)= 0.555642893339299E+01
c_fit_STO( 6,12, 2)= 0.145693086419212E+00
g_fit_STO( 6,12, 2)= 0.224111558110797E+01
c_fit_STO( 7,12, 2)= 0.181423116190035E+00
g_fit_STO( 7,12, 2)= 0.950097470572368E+00
c_fit_STO( 8,12, 2)= 0.205136956321541E+00
g_fit_STO( 8,12, 2)= 0.428328951962141E+00
c_fit_STO( 9,12, 2)= 0.191224736478467E+00
g_fit_STO( 9,12, 2)= 0.175649892324326E+00
c_fit_STO(10,12, 2)=-0.288656249079920E-03
g_fit_STO(10,12, 2)=-0.119727210693168E+00
c_fit_STO(11,12, 2)=-0.279151309763382E-10
g_fit_STO(11,12, 2)=-0.676000816237273E+00
c_fit_STO(12,12, 2)=-0.149440709492466E-17
g_fit_STO(12,12, 2)=-0.160030206461295E+01
! chi2=  0.413348438060518E-09

alpha( 3)= 0.965872774985130E+00
 beta( 3)= 0.100000000000000E+00
c_fit_STO( 1,12, 3)= 0.215832247416248E-01
g_fit_STO( 1,12, 3)= 0.793039459299843E+03
c_fit_STO( 2,12, 3)= 0.316295789238165E-01
g_fit_STO( 2,12, 3)= 0.145473078152385E+03
c_fit_STO( 3,12, 3)= 0.490998626938137E-01
g_fit_STO( 3,12, 3)= 0.407467993373098E+02
c_fit_STO( 4,12, 3)= 0.709200046505098E-01
g_fit_STO( 4,12, 3)= 0.140412153009565E+02
c_fit_STO( 5,12, 3)= 0.988495996902927E-01
g_fit_STO( 5,12, 3)= 0.550918881747263E+01
c_fit_STO( 6,12, 3)= 0.130666267653250E+00
g_fit_STO( 6,12, 3)= 0.236888991291584E+01
c_fit_STO( 7,12, 3)= 0.161189629908686E+00
g_fit_STO( 7,12, 3)= 0.108610702725989E+01
c_fit_STO( 8,12, 3)= 0.170952647879742E+00
g_fit_STO( 8,12, 3)= 0.528053060263405E+00
c_fit_STO( 9,12, 3)= 0.148276940553419E+00
g_fit_STO( 9,12, 3)= 0.272670992654177E+00
c_fit_STO(10,12, 3)= 0.931016146404802E-01
g_fit_STO(10,12, 3)= 0.143449247053137E+00
c_fit_STO(11,12, 3)= 0.173240773638413E-01
g_fit_STO(11,12, 3)= 0.771299437078243E-01
c_fit_STO(12,12, 3)=-0.262445316198037E-02
g_fit_STO(12,12, 3)= 0.332413157853353E-01
! chi2=  0.221185165596165E-09

alpha( 4)= 0.975511624873963E+00
 beta( 4)= 0.100000000000000E+00
c_fit_STO( 1,12, 4)= 0.218130515975567E-01
g_fit_STO( 1,12, 4)= 0.793039439147822E+03
c_fit_STO( 2,12, 4)= 0.319132426597382E-01
g_fit_STO( 2,12, 4)= 0.145473130624258E+03
c_fit_STO( 3,12, 4)= 0.496370035877206E-01
g_fit_STO( 3,12, 4)= 0.407466465691790E+02
c_fit_STO( 4,12, 4)= 0.714652218213560E-01
g_fit_STO( 4,12, 4)= 0.140413562179936E+02
c_fit_STO( 5,12, 4)= 0.999327654401826E-01
g_fit_STO( 5,12, 4)= 0.550954297201870E+01
c_fit_STO( 6,12, 4)= 0.131121132593925E+00
g_fit_STO( 6,12, 4)= 0.236644591571846E+01
c_fit_STO( 7,12, 4)= 0.160770377521425E+00
g_fit_STO( 7,12, 4)= 0.108990374024062E+01
c_fit_STO( 8,12, 4)= 0.174167269810577E+00
g_fit_STO( 8,12, 4)= 0.526265232888283E+00
c_fit_STO( 9,12, 4)= 0.145656027922954E+00
g_fit_STO( 9,12, 4)= 0.267055944516923E+00
c_fit_STO(10,12, 4)= 0.819324504292700E-01
g_fit_STO(10,12, 4)= 0.144479774379223E+00
c_fit_STO(11,12, 4)= 0.237910768930559E-01
g_fit_STO(11,12, 4)= 0.809508330692417E-01
c_fit_STO(12,12, 4)=-0.131611272123263E-02
g_fit_STO(12,12, 4)= 0.456906173062365E-01
! chi2=  0.169737225217881E-09

alpha( 5)= 0.980910434402748E+00
 beta( 5)= 0.100000000000000E+00
c_fit_STO( 1,12, 5)= 0.219263578400078E-01
g_fit_STO( 1,12, 5)= 0.793039437921204E+03
c_fit_STO( 2,12, 5)= 0.321060187625110E-01
g_fit_STO( 2,12, 5)= 0.145473134601179E+03
c_fit_STO( 3,12, 5)= 0.498710236635923E-01
g_fit_STO( 3,12, 5)= 0.407466283168858E+02
c_fit_STO( 4,12, 5)= 0.718987265142471E-01
g_fit_STO( 4,12, 5)= 0.140414251788553E+02
c_fit_STO( 5,12, 5)= 0.100257429760946E+00
g_fit_STO( 5,12, 5)= 0.550927777024094E+01
c_fit_STO( 6,12, 5)= 0.131761148458970E+00
g_fit_STO( 6,12, 5)= 0.236726229680544E+01
c_fit_STO( 7,12, 5)= 0.160789380035816E+00
g_fit_STO( 7,12, 5)= 0.108871307325427E+01
c_fit_STO( 8,12, 5)= 0.171770915897281E+00
g_fit_STO( 8,12, 5)= 0.528108182325493E+00
c_fit_STO( 9,12, 5)= 0.146364682037657E+00
g_fit_STO( 9,12, 5)= 0.268322209794396E+00
c_fit_STO(10,12, 5)= 0.833524871848799E-01
g_fit_STO(10,12, 5)= 0.141669937732336E+00
c_fit_STO(11,12, 5)= 0.209074687727317E-01
g_fit_STO(11,12, 5)= 0.772406908485261E-01
c_fit_STO(12,12, 5)=-0.174243488173077E-03
g_fit_STO(12,12, 5)= 0.401037146400652E-01
! chi2=  0.136715643454771E-09

alpha( 6)= 0.984356832666919E+00
 beta( 6)= 0.100000000000000E+00
c_fit_STO( 1,12, 6)= 0.220012765201499E-01
g_fit_STO( 1,12, 6)= 0.793039436298034E+03
c_fit_STO( 2,12, 6)= 0.322245582104746E-01
g_fit_STO( 2,12, 6)= 0.145473139221844E+03
c_fit_STO( 3,12, 6)= 0.500285228843366E-01
g_fit_STO( 3,12, 6)= 0.407466114724614E+02
c_fit_STO( 4,12, 6)= 0.721583629118042E-01
g_fit_STO( 4,12, 6)= 0.140414666670761E+02
c_fit_STO( 5,12, 6)= 0.100495841783049E+00
g_fit_STO( 5,12, 6)= 0.550916451963638E+01
c_fit_STO( 6,12, 6)= 0.132106586034743E+00
g_fit_STO( 6,12, 6)= 0.236752294231814E+01
c_fit_STO( 7,12, 6)= 0.160768159530767E+00
g_fit_STO( 7,12, 6)= 0.108836635005171E+01
c_fit_STO( 8,12, 6)= 0.170937682576067E+00
g_fit_STO( 8,12, 6)= 0.528710192691113E+00
c_fit_STO( 9,12, 6)= 0.146616702831316E+00
g_fit_STO( 9,12, 6)= 0.268209549546974E+00
c_fit_STO(10,12, 6)= 0.832655400585604E-01
g_fit_STO(10,12, 6)= 0.140082346399405E+00
c_fit_STO(11,12, 6)= 0.200107178777949E-01
g_fit_STO(11,12, 6)= 0.751484084782849E-01
c_fit_STO(12,12, 6)= 0.185198804367570E-03
g_fit_STO(12,12, 6)= 0.418362598923941E-01
! chi2=  0.115254602117819E-09

alpha( 7)= 0.986747894314613E+00
 beta( 7)= 0.100000000000000E+00
c_fit_STO( 1,12, 7)= 0.220556643749937E-01
g_fit_STO( 1,12, 7)= 0.793039435738425E+03
c_fit_STO( 2,12, 7)= 0.323027577499491E-01
g_fit_STO( 2,12, 7)= 0.145473140741207E+03
c_fit_STO( 3,12, 7)= 0.501451219003610E-01
g_fit_STO( 3,12, 7)= 0.407466064108322E+02
c_fit_STO( 4,12, 7)= 0.723238587609512E-01
g_fit_STO( 4,12, 7)= 0.140414755075732E+02
c_fit_STO( 5,12, 7)= 0.100688088562579E+00
g_fit_STO( 5,12, 7)= 0.550915123492193E+01
c_fit_STO( 6,12, 7)= 0.132294594429882E+00
g_fit_STO( 6,12, 7)= 0.236752438752204E+01
c_fit_STO( 7,12, 7)= 0.160743968356892E+00
g_fit_STO( 7,12, 7)= 0.108844399646708E+01
c_fit_STO( 8,12, 7)= 0.170661181704736E+00
g_fit_STO( 8,12, 7)= 0.528820809339505E+00
c_fit_STO( 9,12, 7)= 0.146090406756523E+00
g_fit_STO( 9,12, 7)= 0.268073557359615E+00
c_fit_STO(10,12, 7)= 0.825742789165297E-01
g_fit_STO(10,12, 7)= 0.139939063495536E+00
c_fit_STO(11,12, 7)= 0.202790455415387E-01
g_fit_STO(11,12, 7)= 0.752850101711256E-01
c_fit_STO(12,12, 7)= 0.618802957307556E-03
g_fit_STO(12,12, 7)= 0.434345787712516E-01
! chi2=  0.994013960061617E-10

alpha( 8)= 0.988504266350910E+00
 beta( 8)= 0.100000000000000E+00
c_fit_STO( 1,12, 8)= 0.220967586281247E-01
g_fit_STO( 1,12, 8)= 0.793039435140682E+03
c_fit_STO( 2,12, 8)= 0.323585589028403E-01
g_fit_STO( 2,12, 8)= 0.145473142333472E+03
c_fit_STO( 3,12, 8)= 0.502336795989770E-01
g_fit_STO( 3,12, 8)= 0.407466013944419E+02
c_fit_STO( 4,12, 8)= 0.724396869177374E-01
g_fit_STO( 4,12, 8)= 0.140414824712680E+02
c_fit_STO( 5,12, 8)= 0.100839393000067E+00
g_fit_STO( 5,12, 8)= 0.550915043810252E+01
c_fit_STO( 6,12, 8)= 0.132409432860553E+00
g_fit_STO( 6,12, 8)= 0.236748283164166E+01
c_fit_STO( 7,12, 8)= 0.160716273427507E+00
g_fit_STO( 7,12, 8)= 0.108863475793550E+01
c_fit_STO( 8,12, 8)= 0.170487847758629E+00
g_fit_STO( 8,12, 8)= 0.528879584652504E+00
c_fit_STO( 9,12, 8)= 0.145496091101452E+00
g_fit_STO( 9,12, 8)= 0.268062543205523E+00
c_fit_STO(10,12, 8)= 0.820302407607571E-01
g_fit_STO(10,12, 8)= 0.140119191583470E+00
c_fit_STO(11,12, 8)= 0.207066027777484E-01
g_fit_STO(11,12, 8)= 0.756954341835966E-01
c_fit_STO(12,12, 8)= 0.948065902014289E-03
g_fit_STO(12,12, 8)= 0.436868364899163E-01
! chi2=  0.873609763286178E-10

alpha( 9)= 0.989849200987187E+00
 beta( 9)= 0.100000000000000E+00
c_fit_STO( 1,12, 9)= 0.221290208193080E-01
g_fit_STO( 1,12, 9)= 0.793039434648384E+03
c_fit_STO( 2,12, 9)= 0.324006240320127E-01
g_fit_STO( 2,12, 9)= 0.145473143642487E+03
c_fit_STO( 3,12, 9)= 0.503023840980571E-01
g_fit_STO( 3,12, 9)= 0.407465972876059E+02
c_fit_STO( 4,12, 9)= 0.725267048371542E-01
g_fit_STO( 4,12, 9)= 0.140414880778860E+02
c_fit_STO( 5,12, 9)= 0.100957577205478E+00
g_fit_STO( 5,12, 9)= 0.550915082390418E+01
c_fit_STO( 6,12, 9)= 0.132490178857474E+00
g_fit_STO( 6,12, 9)= 0.236744671927050E+01
c_fit_STO( 7,12, 9)= 0.160690710628948E+00
g_fit_STO( 7,12, 9)= 0.108881439832404E+01
c_fit_STO( 8,12, 9)= 0.170334676382620E+00
g_fit_STO( 8,12, 9)= 0.528944209114403E+00
c_fit_STO( 9,12, 9)= 0.144998517674701E+00
g_fit_STO( 9,12, 9)= 0.268115722018789E+00
c_fit_STO(10,12, 9)= 0.816950268008812E-01
g_fit_STO(10,12, 9)= 0.140305867207138E+00
c_fit_STO(11,12, 9)= 0.210769663888571E-01
g_fit_STO(11,12, 9)= 0.759347781555683E-01
c_fit_STO(12,12, 9)= 0.114914887783532E-02
g_fit_STO(12,12, 9)= 0.434503953098304E-01
! chi2=  0.779699785838617E-10

alpha(10)= 0.990912192506555E+00
 beta(10)= 0.100000000000000E+00
c_fit_STO( 1,12,10)= 0.221642764771145E-01
g_fit_STO( 1,12,10)= 0.793039452298382E+03
c_fit_STO( 2,12,10)= 0.324138741408554E-01
g_fit_STO( 2,12,10)= 0.145473094238917E+03
c_fit_STO( 3,12,10)= 0.503953848844883E-01
g_fit_STO( 3,12,10)= 0.407467681410499E+02
c_fit_STO( 4,12,10)= 0.725215508863754E-01
g_fit_STO( 4,12,10)= 0.140411131588551E+02
c_fit_STO( 5,12,10)= 0.101180274485143E+00
g_fit_STO( 5,12,10)= 0.550996901614608E+01
c_fit_STO( 6,12,10)= 0.132427229312865E+00
g_fit_STO( 6,12,10)= 0.236613149733886E+01
c_fit_STO( 7,12,10)= 0.160592435703491E+00
g_fit_STO( 7,12,10)= 0.108996004071773E+01
c_fit_STO( 8,12,10)= 0.171532545172168E+00
g_fit_STO( 8,12,10)= 0.527889750125052E+00
c_fit_STO( 9,12,10)= 0.144281230386591E+00
g_fit_STO( 9,12,10)= 0.266591443056990E+00
c_fit_STO(10,12,10)= 0.799274660742122E-01
g_fit_STO(10,12,10)= 0.140424130514161E+00
c_fit_STO(11,12,10)= 0.216707908367687E-01
g_fit_STO(11,12,10)= 0.773720995004080E-01
c_fit_STO(12,12,10)= 0.163880489799507E-02
g_fit_STO(12,12,10)= 0.449677838298326E-01
! chi2=  0.703699903727969E-10

alpha(11)= 0.991773530886516E+00
 beta(11)= 0.100000000000000E+00
c_fit_STO( 1,12,11)= 0.221772515114463E-01
g_fit_STO( 1,12,11)= 0.793039434077405E+03
c_fit_STO( 2,12,11)= 0.324596455853854E-01
g_fit_STO( 2,12,11)= 0.145473145152825E+03
c_fit_STO( 3,12,11)= 0.504017463679498E-01
g_fit_STO( 3,12,11)= 0.407465925527315E+02
c_fit_STO( 4,12,11)= 0.726494964977211E-01
g_fit_STO( 4,12,11)= 0.140414942357698E+02
c_fit_STO( 5,12,11)= 0.101128203158766E+00
g_fit_STO( 5,12,11)= 0.550915198719066E+01
c_fit_STO( 6,12,11)= 0.132601047893624E+00
g_fit_STO( 6,12,11)= 0.236739413764722E+01
c_fit_STO( 7,12,11)= 0.160693148454110E+00
g_fit_STO( 7,12,11)= 0.108900401528626E+01
c_fit_STO( 8,12,11)= 0.170185694433107E+00
g_fit_STO( 8,12,11)= 0.528855277961283E+00
c_fit_STO( 9,12,11)= 0.144284432460086E+00
g_fit_STO( 9,12,11)= 0.268047620907645E+00
c_fit_STO(10,12,11)= 0.812329485696124E-01
g_fit_STO(10,12,11)= 0.140460182068636E+00
c_fit_STO(11,12,11)= 0.215577869661399E-01
g_fit_STO(11,12,11)= 0.760809034481258E-01
c_fit_STO(12,12,11)= 0.136568598290250E-02
g_fit_STO(12,12,11)= 0.428592412257919E-01
! chi2=  0.642868417749994E-10

alpha(12)= 0.992485652009950E+00
 beta(12)= 0.100000000000000E+00
c_fit_STO( 1,12,12)= 0.221960096335998E-01
g_fit_STO( 1,12,12)= 0.793039433919574E+03
c_fit_STO( 2,12,12)= 0.324811214802636E-01
g_fit_STO( 2,12,12)= 0.145473145563962E+03
c_fit_STO( 3,12,12)= 0.504387071502016E-01
g_fit_STO( 3,12,12)= 0.407465912825249E+02
c_fit_STO( 4,12,12)= 0.726947241628507E-01
g_fit_STO( 4,12,12)= 0.140414956295602E+02
c_fit_STO( 5,12,12)= 0.101191181958481E+00
g_fit_STO( 5,12,12)= 0.550915368029333E+01
c_fit_STO( 6,12,12)= 0.132641572270050E+00
g_fit_STO( 6,12,12)= 0.236737483702410E+01
c_fit_STO( 7,12,12)= 0.160693600210710E+00
g_fit_STO( 7,12,12)= 0.108907103117982E+01
c_fit_STO( 8,12,12)= 0.170125952088823E+00
g_fit_STO( 8,12,12)= 0.528823320388228E+00
c_fit_STO( 9,12,12)= 0.144029858909936E+00
g_fit_STO( 9,12,12)= 0.268023030453615E+00
c_fit_STO(10,12,12)= 0.810898056681066E-01
g_fit_STO(10,12,12)= 0.140498588349246E+00
c_fit_STO(11,12,12)= 0.217216079911739E-01
g_fit_STO(11,12,12)= 0.760847107750947E-01
c_fit_STO(12,12,12)= 0.142827097899753E-02
g_fit_STO(12,12,12)= 0.426041329612960E-01
! chi2=  0.591541122589250E-10

alpha(13)= 0.993084255002000E+00
 beta(13)= 0.100000000000000E+00
c_fit_STO( 1,12,13)= 0.222205854204712E-01
g_fit_STO( 1,12,13)= 0.793039452147022E+03
c_fit_STO( 2,12,13)= 0.324812147863049E-01
g_fit_STO( 2,12,13)= 0.145473094609751E+03
c_fit_STO( 3,12,13)= 0.505046013587541E-01
g_fit_STO( 3,12,13)= 0.407467664554679E+02
c_fit_STO( 4,12,13)= 0.726661183688669E-01
g_fit_STO( 4,12,13)= 0.140411175025014E+02
c_fit_STO( 5,12,13)= 0.101359414457877E+00
g_fit_STO( 5,12,13)= 0.550996347608445E+01
c_fit_STO( 6,12,13)= 0.132563920248916E+00
g_fit_STO( 6,12,13)= 0.236616018195584E+01
c_fit_STO( 7,12,13)= 0.160433423683075E+00
g_fit_STO( 7,12,13)= 0.109036801856059E+01
c_fit_STO( 8,12,13)= 0.170824802398249E+00
g_fit_STO( 8,12,13)= 0.528689646347139E+00
c_fit_STO( 9,12,13)= 0.143485256607445E+00
g_fit_STO( 9,12,13)= 0.267458979135253E+00
c_fit_STO(10,12,13)= 0.801080405245341E-01
g_fit_STO(10,12,13)= 0.140975855335849E+00
c_fit_STO(11,12,13)= 0.223450520916733E-01
g_fit_STO(11,12,13)= 0.772425734567492E-01
c_fit_STO(12,12,13)= 0.173887473035674E-02
g_fit_STO(12,12,13)= 0.435728166629071E-01
! chi2=  0.547731459035204E-10

alpha(14)= 0.993594489146943E+00
 beta(14)= 0.100000000000000E+00
c_fit_STO( 1,12,14)= 0.222268292591559E-01
g_fit_STO( 1,12,14)= 0.793039433763281E+03
c_fit_STO( 2,12,14)= 0.325141147983730E-01
g_fit_STO( 2,12,14)= 0.145473145953667E+03
c_fit_STO( 3,12,14)= 0.504962023653220E-01
g_fit_STO( 3,12,14)= 0.407465901413797E+02
c_fit_STO( 4,12,14)= 0.727655528101140E-01
g_fit_STO( 4,12,14)= 0.140414961677820E+02
c_fit_STO( 5,12,14)= 0.101287757325907E+00
g_fit_STO( 5,12,14)= 0.550915979685908E+01
c_fit_STO( 6,12,14)= 0.132705861344547E+00
g_fit_STO( 6,12,14)= 0.236734716552515E+01
c_fit_STO( 7,12,14)= 0.160667162578161E+00
g_fit_STO( 7,12,14)= 0.108921018479873E+01
c_fit_STO( 8,12,14)= 0.169966613239694E+00
g_fit_STO( 8,12,14)= 0.528900791482846E+00
c_fit_STO( 9,12,14)= 0.143660028670243E+00
g_fit_STO( 9,12,14)= 0.268094034448170E+00
c_fit_STO(10,12,14)= 0.809560613958561E-01
g_fit_STO(10,12,14)= 0.140580760188760E+00
c_fit_STO(11,12,14)= 0.219718840403128E-01
g_fit_STO(11,12,14)= 0.760451610757649E-01
c_fit_STO(12,12,14)= 0.150831372987378E-02
g_fit_STO(12,12,14)= 0.421770514281929E-01
! chi2=  0.510956743898852E-10

alpha(15)= 0.994034583513567E+00
 beta(15)= 0.100000000000000E+00
c_fit_STO( 1,12,15)= 0.222477455560410E-01
g_fit_STO( 1,12,15)= 0.793039451454492E+03
c_fit_STO( 2,12,15)= 0.325101733121531E-01
g_fit_STO( 2,12,15)= 0.145473096517224E+03
c_fit_STO( 3,12,15)= 0.505517707256979E-01
g_fit_STO( 3,12,15)= 0.407467594991375E+02
c_fit_STO( 4,12,15)= 0.727312868538875E-01
g_fit_STO( 4,12,15)= 0.140411337342015E+02
c_fit_STO( 5,12,15)= 0.101433213333613E+00
g_fit_STO( 5,12,15)= 0.550992406161165E+01
c_fit_STO( 6,12,15)= 0.132629148054552E+00
g_fit_STO( 6,12,15)= 0.236622856033653E+01
c_fit_STO( 7,12,15)= 0.160444475686864E+00
g_fit_STO( 7,12,15)= 0.109035055686535E+01
c_fit_STO( 8,12,15)= 0.170620788118564E+00
g_fit_STO( 8,12,15)= 0.528730620599357E+00
c_fit_STO( 9,12,15)= 0.143176086082243E+00
g_fit_STO( 9,12,15)= 0.267568934904281E+00
c_fit_STO(10,12,15)= 0.801098461666549E-01
g_fit_STO(10,12,15)= 0.141003089880523E+00
c_fit_STO(11,12,15)= 0.225132907163620E-01
g_fit_STO(11,12,15)= 0.770369439411710E-01
c_fit_STO(12,12,15)= 0.175910560789321E-02
g_fit_STO(12,12,15)= 0.429876752077514E-01
! chi2=  0.478480280610739E-10

alpha(16)= 0.994418074640156E+00
 beta(16)= 0.100000000000000E+00
c_fit_STO( 1,12,16)= 0.222532886188265E-01
g_fit_STO( 1,12,16)= 0.793039433500529E+03
c_fit_STO( 2,12,16)= 0.325357379196042E-01
g_fit_STO( 2,12,16)= 0.145473146660011E+03
c_fit_STO( 3,12,16)= 0.505427258339209E-01
g_fit_STO( 3,12,16)= 0.407465877570169E+02
c_fit_STO( 4,12,16)= 0.728120081416362E-01
g_fit_STO( 4,12,16)= 0.140415001044837E+02
c_fit_STO( 5,12,16)= 0.101370230004173E+00
g_fit_STO( 5,12,16)= 0.550915004509140E+01
c_fit_STO( 6,12,16)= 0.132739841566244E+00
g_fit_STO( 6,12,16)= 0.236730710663027E+01
c_fit_STO( 7,12,16)= 0.160865162183687E+00
g_fit_STO( 7,12,16)= 0.108899146512316E+01
c_fit_STO( 8,12,16)= 0.170350959947553E+00
g_fit_STO( 8,12,16)= 0.527924354715341E+00
c_fit_STO( 9,12,16)= 0.143274342659239E+00
g_fit_STO( 9,12,16)= 0.267241817912923E+00
c_fit_STO(10,12,16)= 0.803484853039541E-01
g_fit_STO(10,12,16)= 0.140278133078676E+00
c_fit_STO(11,12,16)= 0.220458377762626E-01
g_fit_STO(11,12,16)= 0.760099194215364E-01
c_fit_STO(12,12,16)= 0.158550173232041E-02
g_fit_STO(12,12,16)= 0.420004905436282E-01
! chi2=  0.450597474738305E-10

end

subroutine build_mapping_ijkl
include 'j.inc'
logical logic
np=(nbasis*(nbasis+1))/2
nint_theo=(np*(np+1))/2
if(nint_theo.gt.nint_max)then
 print*,'expected number of two-electron integrals is',nint_theo
 print*,'declared size of array for two-electron integrals is',nint_max
 print*,' I continue but you shall need to increase nint_max'
endif

kcp=0
do l=1,nbasis
 do k=1,nbasis
  do j=l,nbasis
   do i=k,nbasis
   logic=.true.
   if(i.eq.j.and.l.gt.k)logic=.false.
   if(i.ge.j.and.logic)then
    kcp=kcp+1
    is(kcp)=i
    js(kcp)=j
    ks(kcp)=k
    ls(kcp)=l
    kcp_ijkl(i,j,k,l)=kcp
    endif
   enddo
  enddo
 enddo
enddo
nint=kcp
if(nint.ne.nint_theo)stop 'pb in build_mapping_ijkl'
end

subroutine count_multi_center_integrals
include 'j.inc'

iac_1=0 
iac_2=0 
iac_3=0 
iac_4=0 
do kcp=1,nint
 i=is(kcp)
 k=ks(kcp)
 j=js(kcp)
 l=ls(kcp)
 call compare_nuclei(nucleus_number(i),nucleus_number(j),nucleus_number(k),nucleus_number(l) &
 ,ndiff)
 if(ndiff.eq.1)then
  iac_1=iac_1+1
 endif
 if(ndiff.eq.2)iac_2=iac_2+1
 if(ndiff.eq.3)iac_3=iac_3+1
 if(ndiff.eq.4)iac_4=iac_4+1
enddo
write(*,*)'Number of one-center two-electron integrals =',iac_1
write(*,*)'Number of two-center two-electron integrals =',iac_2
write(*,*)'Number of three-center two-electron integrals =',iac_3
write(*,*)'Number of four-center two-electron integrals =',iac_4
write(*,*)'Total number of (ik|1/r12|jl) =',iac_1+iac_2+iac_3+iac_4
end

subroutine check_basis_type(filename_basis)
include 'j.inc'
character*(128) :: filename_basis

if(filename_basis.eq.'VB0'.or.  &
   filename_basis.eq.'VB1'.or.  &
   filename_basis.eq.'VB2'.or.  &
   filename_basis.eq.'VB3'.or.  &
   filename_basis.eq.'CVB1'.or. &
   filename_basis.eq.'CVB2'.or. &
   filename_basis.eq.'CVB3'    )then
 if(basis_type.ne.'STO')stop 'pb with type of orbitals!'
 return
endif

if(filename_basis.eq.'sto-6g'.or. &
   filename_basis.eq.'cc-pvdz'.or. &
   filename_basis.eq.'cc-pvtz'.or. &
   filename_basis.eq.'cc-pvqz'.or. &
   filename_basis.eq.'cc-pv5z'.or. &
   filename_basis.eq.'cc-pv6z'.or. &
   filename_basis.eq.'cc-cpvdz'.or. &
   filename_basis.eq.'cc-cpvtz'.or. &
   filename_basis.eq.'cc-cpvqz'.or. &
   filename_basis.eq.'cc-cpv5z'.or. &
   filename_basis.eq.'cc-cpv6z'.or. &
   filename_basis.eq.'aug-cc-pvdz'.or. &
   filename_basis.eq.'aug-cc-pvtz'.or. &
   filename_basis.eq.'aug-cc-pvqz'.or. &
   filename_basis.eq.'aug-cc-pv5z'.or. &
   filename_basis.eq.'aug-cc-pv6z'.or. &
   filename_basis.eq.'aug-cc-cpvdz'.or. &
   filename_basis.eq.'aug-cc-cpvtz'.or. &
   filename_basis.eq.'aug-cc-cpvqz'.or. &
   filename_basis.eq.'aug-cc-cpv5z'.or. &
   filename_basis.eq.'aug-cc-cpv6z'    )then
   if(basis_type.ne.'GTO')stop 'pb with type of orbitals!'
 return
endif
write(*,*)'basis set name not recognized'
write(*,*)'Pick among VB0 VB1 VB2 VB3 CVB1 CVB2 CVB3 cc-pvdz cc-pvtz cc-pvqz cc-pv5z cc-pv6z'
write(*,*)'aug-cc-pvdz aug-cc-pvtz aug-cc-pvqz aug-cc-pv5z aug-cc-pv6z'
write(*,*)'aug-cc-pcvdz aug-cc-pcvtz aug-cc-pcvqz aug-cc-pcv5z'
stop
end

!!  I = \int_{-inf}^{inf} dx x^{2n} exp(-g x^2) = sqrt(pi/g)  (2n-1)!!/(2g)**n
double precision function gauss_int(n,g)
implicit double precision (a-h,o-z)
pi=dacos(-1.d0)
gauss_int= dsqrt(pi/g)*dblefact(2*n-1)/(2.d0*g)**n

end

double precision function rnorm_prim(nx,ny,nz,g)
implicit double precision (a-h,o-z)
rnorm_prim=dsqrt( dabs( gauss_int(nx,2.d0*g)*gauss_int(ny,2.d0*g)*gauss_int(nz,2.d0*g)))
end

!! Evaluate rc for a general STO (type='STO')    phi(r)= c(1)* x**nx y **ny z**nz exp(-g(1)r) 
!! or  a general contracted  GTO (type='GTO')    phi(r)=       x**nx y **ny z**nz [ sum_i=1,nc  c(i) exp(-g(i)*r**2)
!!
!! Criterium for choosing rc: rc is such that    
!!  STO:  largest rc such that: 
!!        c(1)* rc**(nx+ny+nz) exp(-g(1)rc) = eps *  rmax**(nx+ny+nz) exp(-g(1)rmax)
!!        where rmax is the value such that r**(nx+ny+nz) exp(-g(1)r) is maximum
!!
!! GTO:  largest rc such that:
!!       rc**(nx+ny+nz)[sum_i=1,nc c(i) exp(-g(i)*rc**2)]= eps*rmax**(nx+ny+nz)[sum_i=1,nc c(i) exp(-g(i)*rmax**2)]

!! out:   r_c(i_orb) and r_infty(i_orb)

subroutine compute_rc_of_orbital(i_orb)
include 'j.inc'

ntot=npower(1,i_orb)+npower(2,i_orb)+npower(3,i_orb)

if(basis_type.eq.'STO')then
 call find_rmax_rc_STO(ntot,n_sto(i_orb),g_slater(i_orb),n_eps,rmax,r_c(i_orb),r_infty(i_orb))
endif  ! STO

if(basis_type.eq.'GTO')then

! estimate of rmax of the  r**ntot u(r) 
 eps=10.d0**(-n_eps)
 r0=dsqrt(dabs(-dlog(eps/100.d0)/g_min(i_orb)))
 nsteps=10000
 dr=r0/nsteps
 phi0_max=0.d0
 do i=1,nsteps
  r0=r0-dr
  finite_range_save=finite_range
  finite_range=.false.
  bid=r0**ntot*dabs(u_orb(i_orb,r0))
finite_range=finite_range_save
!  phi0=dabs(r0**ntot*u_orb(.false.,i_orb,r0))
  phi0=bid
  if(phi0.gt.phi0_max)then
   phi0_max=phi0
   rmax=dabs(r0) 
  endif
 enddo

!! compute r_c of orbital depending on n_eps
 eps=10.d0**(-n_eps)

 phi_target=eps*dabs(phi0_max)
 r0=dsqrt(dabs(-dlog(phi_target/1000.d0)/g_min(i_orb)))
 finite_range_save=finite_range
 finite_range=.false.
 bid=r0**ntot*dabs(u_orb(i_orb,r0))
finite_range=finite_range_save
! diff=phi_target-r0**ntot*dabs(u_orb(.false.,i_orb,r0))
 diff=phi_target-bid
 if(diff.lt.0.d0)stop 'r0 too small'

 dr=r0/nsteps
 do while (diff.gt.0.d0)
  r0=r0-dr
 finite_range_save=finite_range
 finite_range=.false.
 bid=r0**ntot*dabs(u_orb(i_orb,r0))
finite_range=finite_range_save
!  diff=phi_target-r0**ntot*dabs(u_orb(.false.,i_orb,r0))
  diff=phi_target-bid
 enddo
 range_orbital=r0
 finite_range_save=finite_range
 finite_range=.false.
 bid=r0**ntot*dabs(u_orb(i_orb,r0))
finite_range=finite_range_save
! phix=r0**ntot*dabs(u_orb(.false.,i_orb,r0))
 phix=bid
 if((phix-phi_target).gt.0.01d0)then
 finite_range_save=finite_range
 finite_range=.false.
 bid=r0**ntot*dabs(u_orb(i_orb,r0))
finite_range=finite_range_save
!  write(*,*)'phi en r0=',r0**ntot*dabs(u_orb(.false.,i_orb,r0))
  write(*,*)'phi en r0=',bid
  write(*,*)'phi_target =',phi_target
  stop 'pb 1 with GTO in finite range'
 endif
 r_c(i_orb)=r0

!! computation of r_infty where phi is zero
 eps=10.d0**(-16)
 phi_target=eps*dabs(phi0_max)
 r0=dsqrt(dabs(-dlog(phi_target/1000.d0)/g_min(i_orb)))

 dr=r0/nsteps
 finite_range_save=finite_range
 finite_range=.false.
 bid=r0**ntot*dabs(u_orb(i_orb,r0))
finite_range=finite_range_save

! diff=phi_target-r0**ntot*dabs(u_orb(.false.,i_orb,r0))
 diff=phi_target-bid
 if(diff.lt.0.d0)stop 'r0 too small'

 do while (diff.gt.0.d0)
  r0=r0-dr
 finite_range_save=finite_range
 finite_range=.false.
 bid=r0**ntot*dabs(u_orb(i_orb,r0))
finite_range=finite_range_save
!  diff=phi_target-r0**ntot*dabs(u_orb(.false.,i_orb,r0))
   diff=phi_target-bid
 enddo
 range_orbital=r0
 finite_range_save=finite_range
 finite_range=.false.
 bid=r0**ntot*dabs(u_orb(i_orb,r0))
finite_range=finite_range_save

! if(dabs(r0**ntot*u_orb(.false.,i_orb,r0)-phi_target).gt.0.01d0)stop 'pb 2 with GTO in finite range'
 if((bid-phi_target).gt.0.01d0)stop 'pb 2 with GTO in finite range'
 r_infty(i_orb)=r0

endif !GTO
end

subroutine compute_n_star_of_primitive(i_orb)
include 'j.inc'
do mm=1,n_contract(i_orb)
 ntemp=g_contract(mm,i_orb)*r_c(i_orb)**2/dlog(10.d0)
 n_star(mm,i_orb)=max(1,ntemp)
 if(n_star(mm,i_orb).gt.16)n_star(mm,i_orb)=16
 r_c_prim(mm,i_orb)=dsqrt(n_star(mm,i_orb)*dlog(10.d0)/g_contract(mm,i_orb))
 ng(mm,i_orb)=ng_star(n_star(mm,i_orb))
enddo !mm
end

integer function nprox(n,x)
implicit double precision (a-h,o-z)
if(x.lt.n)stop 'pb in nprox'
diff1=x-n
diff2=n+1-x
if(diff1.lt.diff2)then
 nprox=n
else
 nprox=n+1
endif
end

subroutine print_orb
include 'j.inc'
character*(128) :: filename(1000)
character(len=8) fmt,x1
integer unit(1000)
fmt='(I5.5)'
do k=1,nbasis
 write(x1,fmt)k
 unit(k)=200+k
 filename(k)='print_orb/'//trim(x1)
 open(unit=unit(k),file=filename(k))
enddo

npts=1000
do k=1,nbasis
 ntot=npower(1,k)+npower(2,k)+npower(3,k)
 dr=r_infty(k)/npts
 do i=1,npts
  r=(i-1)*dr
  write(unit(k),*)r,r**ntot*u_orb(k,r),r**ntot*u_gauss(k,r)
 enddo
enddo ! nbasis

end

!!!!  u_orb is the radial part of the atomic orbitals for which the 1- and 2-electron integrals are computed
double precision function u_orb(i_orb,r)
include 'j.inc'
u_orb=0.d0
if(finite_range)then
 do mm=1,n_contract(i_orb)
  contrib=c_contract(mm,i_orb)*primitive(r,i_orb,mm,g_contract(mm,i_orb))
  u_orb=u_orb+contrib
 enddo  !mm
else
 do mm=1,n_contract(i_orb)
  if(basis_type.eq.'GTO') &
  contrib=c_contract(mm,i_orb)*dexp(-g_contract(mm,i_orb)*r**2)
  if(basis_type.eq.'STO') &
  contrib=c_contract(mm,i_orb)*r**n_sto(i_orb)*dexp(-g_contract(mm,i_orb)*r   )
  u_orb=u_orb+contrib
 enddo  !mm
endif
return
end
!! derivative of u_orb
double precision function up_orb(i_orb,r)
include 'j.inc'
up_orb=0.d0
if(finite_range)then
 do mm=1,n_contract(i_orb)
  contrib=c_contract(mm,i_orb)*primitive_p(r,i_orb,mm,g_contract(mm,i_orb))
  up_orb=up_orb+contrib
 enddo  !mm
else
 do mm=1,n_contract(i_orb)
  if(basis_type.eq.'GTO') &
  contrib=-2.d0*g_contract(mm,i_orb)*r*c_contract(mm,i_orb)*dexp(-g_contract(mm,i_orb)*r**2)
  if(basis_type.eq.'STO') &
  contrib=(n_sto(i_orb)/r-g_contract(mm,i_orb))*c_contract(mm,i_orb)*r**n_sto(i_orb)*dexp(-g_contract(mm,i_orb)*r   )
  up_orb=up_orb+contrib
 enddo  !mm
endif
return
end

!! u_gauss is the radial part of the end-of-the-day gaussian representation 
double precision function u_gauss(i,r)
include 'j.inc'
u_gauss=0.d0
do mm=1,n_gauss(i)
 u_gauss=u_gauss+(c_gauss(1,mm,i)+c_gauss(2,mm,i)*r**2)*dexp(-g_gauss(mm,i)*r**2)
enddo
end
!! derivative of u_gauss wrt r
double precision function up_gauss(i,r)
include 'j.inc'
up_gauss=0.d0
do mm=1,n_gauss(i)
 gi=g_gauss(mm,i)
 c1=c_gauss(1,mm,i)
 c2=c_gauss(2,mm,i)
 up_gauss=up_gauss+2.d0*r*(c2-gi*c1-gi*c2*r**2)*dexp(-gi*r**2)
enddo
end

double precision function chi2(i_orb)
include 'j.inc'
ntot=npower(1,i_orb)+npower(2,i_orb)+npower(3,i_orb)
npts=1000
dr=r_infty(i_orb)/npts
chi2=0.d0
if(basis_type.eq.'GTO')then
 do k=1,npts
  r=(k-1)*dr
  u1=u_orb(i_orb,r)
  finite_range_save=finite_range
  finite_range=.false.
  u2=u_orb(i_orb,r)
  finite_range=finite_range_save
  chi2=chi2+r*r**ntot*(u1-u2)**2
 enddo !r
endif

if(basis_type.eq.'STO')then
 do k=1,npts
  r=(k-1)*dr
  u1=u_orb  (i_orb,r)
  u2=u_gauss(i_orb,r)
  chi2=chi2+r*r**ntot*(u1-u2)**2
 enddo !r
endif
chi2=chi2/npts

end

subroutine write_STO_in_file_info_basis(i_orb)
include 'j.inc'

if(finite_range)then
 write(21,'(a,i4,3a)')'Orb #',i_orb,'  ',trim(orb_name_full(i_orb)),'  ',trim(ATOM(nucleus_number(i_orb)))
 write(21,'(a,f7.3,a,f7.3)')'STO is about 0 at ',r_infty(i_orb),'  We impose to vanish at rc= ',r_c(i_orb)
else
 write(21,'(a,i4,5a,i4,a)')'Orb #',i_orb,'  ',trim(orb_name_full(i_orb)),'  ',trim(ATOM(nucleus_number(i_orb))), &
 ' STO fitted with',ng0,' gaussians'
endif
write(21,'(a,i2)')'Number of contractions= ',n_contract(i_orb)
rcmax=-1.d0
do mm=1,n_contract(i_orb)
 write(21,'(a5,f10.3,a6,f10.3)')'  c_i= ',c_contract(mm,i_orb),' g_i= ',g_contract(mm,i_orb)
 if(finite_range)then
  write(21,'(a,i3,a,f10.3)')'  n_star of KT primitive=',n_star(mm,i_orb),' KT primitive vanishes at ',r_c_prim(mm,i_orb)
  write(21,'(a,i3,a)')'  Fitted with ',ng(mm,i_orb),' gaussians'
  write(21,*)
  if(r_c_prim(mm,i_orb).gt.rcmax)rcmax=r_c_prim(mm,i_orb)
 endif
enddo
if(finite_range)then
 write(21,'(a,f10.3)')'New r_c for finite-range KT orbital',rcmax
 write(21,'(a,e22.15)')'chi2 (STO orbital and its KT representation)= ',chi2(nbasis)
else
 write(21,'(a,e22.15)')'chi2 (STO orbital and its gaussian representation)= ',chi2(nbasis)
endif
write(21,*)
write(22,'(i4,i3)')i_orb,n_gauss(i_orb)
write(22,'(4e22.15)')(c_gauss(1,mm,i_orb),mm=1,n_gauss(i_orb))
write(22,'(4e22.15)')(c_gauss(2,mm,i_orb),mm=1,n_gauss(i_orb))
write(22,'(4e22.15)')(g_gauss(mm,i_orb)  ,mm=1,n_gauss(i_orb))
write(22,*)
end

subroutine write_GTO_in_file_info_basis(i_orb)
include 'j.inc'

if(finite_range)then
 write(21,'(a,i4,3a)') &
 'Orb #',i_orb,'  ',trim(orb_name_full(i_orb)),'  ',trim(ATOM(nucleus_number(i_orb)))
 write(21,'(a,f7.3,a,f7.3)')'GTO is about 0 at ',r_infty(i_orb),'  We impose to vanish at rc= ',r_c(i_orb)
else
 write(21,'(a,i4,5a)') &
 'Orb #',i_orb,'  ',trim(orb_name_full(i_orb)),'  ',trim(ATOM(nucleus_number(i_orb))),' GTO untouched'
endif
write(21,'(a,i2)')'Number of contractions= ',n_contract(i_orb)
rcmax=-1.d0
do mm=1,n_contract(i_orb)
 write(21,'(a5,f10.3,a6,f10.3)')'  c_i= ',c_contract(mm,i_orb),' g_i= ',g_contract(mm,i_orb)
 if(finite_range)then
  write(21,'(a,i3,a,f10.3)')'  n_star of KT primitive=',n_star(mm,i_orb),' KT primitive vanishes at ',r_c_prim(mm,i_orb)
  write(21,'(a,i3,a)')'  Fitted with ',ng(mm,i_orb),' gaussians'
  write(21,*)
  if(r_c_prim(mm,i_orb).gt.rcmax)rcmax=r_c_prim(mm,i_orb)
 endif
enddo
if(finite_range)then
 write(21,'(a,f10.3)')'New r_c for finite-range KT orbital',rcmax
 write(21,'(a,e22.15)')'chi2 (GTO orbital and its KT representation)= ',chi2(nbasis)
endif
write(21,*)

write(22,'(i4,i3)')i_orb,n_gauss(i_orb)
write(22,'(4e22.15)')(c_gauss(1,mm,i_orb),mm=1,n_gauss(i_orb))
write(22,'(4e22.15)')(c_gauss(2,mm,i_orb),mm=1,n_gauss(i_orb))
write(22,'(4e22.15)')(g_gauss(mm,i_orb)  ,mm=1,n_gauss(i_orb))
write(22,*)

end

!! GTO: fit_kt= exp(-alpha gam r**2 /(1-sqrt(gam)*r/rc)**beta)
!! STO: fit_kt= exp(-alpha gam r/(1-gam*r/rc)**beta)
double precision function fit_kt(n,alp,bet,gam,r,rc)
include 'j.inc'

if(basis_type.eq.'GTO')then
 u=dsqrt(gam)*r
 rc=dsqrt(dabs(-dlog(10.d0**(-n))))
 fit_kt=0.d0
 if(u.lt.rc)then
  arg=alp*u**2/(1.d0-(u/rc)**2)**bet
  if(arg.lt.-dlog(1.d-16))fit_kt=dexp(-arg)
 endif
endif

if(basis_type.eq.'STO')then
 u=gam*r
 rc=dabs(-dlog(10.d0**(-n)))
 fit_kt=0.d0
 if(u.lt.rc)then
  arg=alp*u/(1.d0-u/rc)**bet
  if(arg.lt.-dlog(1.d-16))fit_kt=dexp(-arg)
 endif
endif
end
!! log derivative of fit_kt  (fit_kt)'/fit_kt
double precision function fit_kt_r(n,alp,bet,gam,r,rc)
include 'j.inc'

if(basis_type.eq.'GTO')then
 u=dsqrt(gam)*r
 rc=dsqrt(dabs(-dlog(10.d0**(-n))))
 fit_kt_r=0.d0
 if(u.lt.rc)then
  fit_kt_r=-2.d0*alp*gam*r/(1.d0-(u/rc)**2)**(bet+1.d0)*(1.d0+(bet-1.d0)*(u/rc)**2) &
  *fit_kt(n,alp,bet,gam,r,rc)
 endif
endif

if(basis_type.eq.'STO')then
 u=gam*r
 rc=dabs(-dlog(10.d0**(-n)))
 fit_kt_r=0.d0
 if(u.lt.rc)then
  fit_kt_r=-alp*gam/(1.d0-u/rc)**(bet+1.d0)*(1.d0+(bet-1.d0)*u/rc) &
           *fit_kt(n,alp,bet,gam,r,rc)
 endif
endif

end


subroutine build_c_g_gauss_STO(i_orb)
include 'j.inc'
kcp=0
do mm=1,n_contract(i_orb)
 if(finite_range)then
  do ii=1,ng(mm,i_orb)
   kcp=kcp+1
   c_gauss(1,kcp,i_orb)=c_contract(mm,i_orb)*c_fit_STO(ii,ng(mm,i_orb),n_star(mm,i_orb))
   c_gauss(2,kcp,i_orb)=0.d0
   g_gauss(kcp,i_orb)=g_contract(mm,i_orb)**2*g_fit_STO(ii,ng(mm,i_orb),n_star(mm,i_orb))
  enddo ! ii
 endif
 if(.not.finite_range)then
  if(i_type(i_orb).eq.1)then
   do ii=1,ng(mm,i_orb)
    kcp=kcp+1
    c_gauss(1,kcp,i_orb)=c_contract(mm,i_orb)*c_fit_exp(ii,ng(mm,i_orb))
    c_gauss(2,kcp,i_orb)=0.d0
    g_gauss(kcp,i_orb)=g_contract(mm,i_orb)**2*g_fit_exp(ii,ng(mm,i_orb))
   enddo ! ii
  endif
  if(i_type(i_orb).eq.2)then
   do ii=1,ng(mm,i_orb)
    kcp=kcp+1
    c_gauss(1,kcp,i_orb)=0.d0
    c_gauss(2,kcp,i_orb)=c_contract(mm,i_orb) &
                  *2.d0*c_fit_exp(ii,ng(mm,i_orb))*g_fit_exp(ii,ng(mm,i_orb))*g_contract(mm,i_orb)
    g_gauss(kcp,i_orb)=g_contract(mm,i_orb)**2*g_fit_exp(ii,ng(mm,i_orb))
   enddo ! ii
  endif
  if(i_type(i_orb).eq.3)stop 'i_type=3 not yet coded!'
 endif
enddo !mm
n_gauss(i_orb)=kcp

end


subroutine build_c_g_gauss_GTO(i_orb)
include 'j.inc'
kcp=0
gmin=1.d+06

do mm=1,n_contract(i_orb)
  
 if(finite_range)then
  do ii=1,ng(mm,i_orb)
   kcp=kcp+1
   c_gauss(1,kcp,i_orb)=c_contract(mm,i_orb)*c_fit_GTO(ii,ng(mm,i_orb),n_star(mm,i_orb))
   c_gauss(2,kcp,i_orb)=0.d0
   g_gauss(kcp,i_orb)=g_contract(mm,i_orb)*g_fit_GTO(ii,ng(mm,i_orb),n_star(mm,i_orb))
   if(g_gauss(kcp,i_orb).lt.gmin)gmin=g_gauss(kcp,i_orb)
  enddo ! ii
 endif
 if(.not.finite_range)then
  kcp=kcp+1
  c_gauss(1,kcp,i_orb)=c_contract(mm,i_orb)
  c_gauss(2,kcp,i_orb)=0.d0
  g_gauss(kcp,i_orb)=g_contract(mm,i_orb)
 endif
enddo !mm
n_gauss(i_orb)=kcp
gauss_min(i_orb)=gmin
end

double precision function primitive(r,i_orb,mm,gam)
include 'j.inc'
nkt=n_star(mm,i_orb)
if(basis_type.eq.'GTO')primitive=fit_kt(nkt,alpha(nkt),beta(nkt),gam,r,rc)
if(basis_type.eq.'STO'.and.(n_sto(i_orb).eq.0.or.n_sto(i_orb).eq.2)) &
                     primitive=fit_kt(nkt,alpha(nkt),beta(nkt),gam,r,rc)
if(basis_type.eq.'STO'.and.n_sto(i_orb).eq.1) &
                     primitive=fit_kt_r(nkt,alpha(nkt),beta(nkt),gam,r,rc)
end
!! derivative of primtitive
double precision function primitive_p(r,i_orb,mm,gam)
include 'j.inc'
nkt=n_star(mm,i_orb)
if(basis_type.eq.'GTO')then
primitive_p=fit_kt_r(nkt,alpha(nkt),beta(nkt),gam,r,rc)
endif
if(basis_type.eq.'STO'.and.(n_sto(i_orb).eq.0.or.n_sto(i_orb).eq.2)) &
                     primitive_p=fit_kt_r(nkt,alpha(nkt),beta(nkt),gam,r,rc)
if(basis_type.eq.'STO'.and.n_sto(i_orb).eq.1) &
                     primitive_p=fit_kt_r(nkt,alpha(nkt),beta(nkt),gam,r,rc)
end


subroutine find_rmax_rc_STO(ntot,n_sto,gam,n,rmax,rc,rinf)
implicit double precision (a-h,o-z)
k=ntot+n_sto
rk=dfloat(k)

epsi=epsil(k,n)
x0=(rk+dsqrt(rk**2-4.d0*dlog(0.9d0*epsi)))**2/4.d0
dx=0.001d0
diff=1.d0
do while (diff.ge.0.d0)
 diff=epsi-fx(k,x0)
 x0=x0-dx
enddo
rmax=rk/gam
rc=x0/gam

epsi=epsil(k,16)
x0=(rk+dsqrt(rk**2-4.d0*dlog(0.9d0*epsi)))**2/4.d0
dx=0.001d0
diff=1.d0
do while (diff.ge.0.d0)
 diff=epsi-fx(k,x0)
 x0=x0-dx
enddo
rinf=x0/gam

end

double precision function fx(k,x)
implicit double precision (a-h,o-z)
fx=x**k*dexp(-x)
end

double precision function epsil(k,n)
implicit double precision (a-h,o-z)
rk=dfloat(k)
epsil=rk**rk*dexp(-rk)*10.d0**(-n)
end

