program integrals
include 'j.inc'
!!!!!!!!!!COMMENT IF NOT MPI
  include 'mpif.h'
!!!!!!!!!!COMMENT IF NOT MPI

character*80,charabid
character*80,MOLECULE

logical normalization_OA

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
double precision moy(nint_max), moy2(nint_max), moy2t(nint_max), Kij
!! END BIG ARRAYS nint_max

!! MONTE CARLO PART
integer i_tab_mc(nbasis_max,nbasis_max)
double precision :: d_x(4)

double precision, allocatable :: r12_inv(:),r1(:,:),r2(:,:)
double precision, allocatable :: pi_0(:),pot(:),rkin(:),rkin_G(:)
double precision, allocatable :: rt1(:,:),rt2(:,:)
double precision, allocatable :: ut1(:,:,:),ut2(:,:,:)
double precision, allocatable :: rho(:,:,:),poly(:,:,:)
double precision, allocatable :: weight(:),weight_kin(:)
double precision, allocatable :: rho_G(:,:,:)
double precision, allocatable :: weight_G(:),weight_kin_G(:)
double precision, allocatable :: rjacob(:), rjacob1(:,:,:), rjacob2(:,:,:)
double precision, allocatable :: pi_1(:)

character*(128) :: filename_in
character*(128) :: filename_out_gaus_ijkl
character*(128) :: filename_out_ijkl
character*(128) :: filename_basis

integer mpi_rank
logical MPI
integer ierr, mpi_size, irequest
integer, dimension(MPI_STATUS_SIZE) :: mpi_status
logical mpi_flag
double precision mpi_size_inv

integer*4                      :: seed(33)
integer*4                      :: put(33)

! call getarg(1,filename_in)

allocate(r12_inv(nw),r1(nw,3),r2(nw,3))
allocate(pi_0(nw),pot(nw),rkin(nw),rkin_G(nw))
allocate(rt1(nw,3),rt2(nw,3))
allocate(ut1(3,nw,nbasis_max*nbasis_max),ut2(3,nw,nbasis_max*nbasis_max))
allocate(rho(nw,nbasis_max*nbasis_max,2),poly(nw,nbasis_max*nbasis_max,2))
allocate(weight(nw),weight_kin(nw) )
allocate(rho_G(nw,nbasis_max*nbasis_max,2))
allocate(weight_G(nw),weight_kin_G(nw))

MPI=.false.
mpi_rank=0
ijkl = 0.d0
ijkl2 = 0.d0

!!!!!!!!!!!COMMENT IF NOT MPI
 call mpi_init(ierr)
 call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_rank, ierr)
 call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_size, ierr)
 call sleep(mpi_rank/20)
 MPI=mpi_size > 1
! write(*,*)'mpi_rank=',mpi_rank
if(mpi_rank.eq.0)write(*,*)'mpi_size=',mpi_size
!!!!!!!!!!!COMMENT IF NOT MPI

write(filename_in,'(A4)') 'j_in'

if(mpi_rank.eq.0)write(*,*)'INPUT FILE USED=',filename_in
open(unit=5,file=filename_in)

if(mpi_rank.eq.0)print*,'Simulation number?'
read(5,'(a80)')charabid
read(5,*)num_simulation
if(mpi_rank.eq.0)write(*,*)'num_simulation=',num_simulation

if(mpi_rank.eq.0)then
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
read(5,*)
read(5,'(a80)')charabid
if(mpi_rank.eq.0)write(6,'(a80)')charabid
read(5,*)n_eps
if(mpi_rank.eq.0)write(6,*)n_eps
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
read(5,*)
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
!if(mpi_rank.eq.0)print*,'mono_center_Slater?'
read(5,'(a80)')charabid
read(5,*)
!if(mpi_rank.eq.0)print*,'compute momo-center Slater=',mono_center_Slater
read(5,'(a80)')charabid
read(5,*)normalization_OA
if(mpi_rank.eq.0)print*,'orbital are normalized ?',normalization_OA

!*****************
! END READ INPUT
!*****************

call read_basis(filename_basis)
if(mpi_rank.eq.0)print*,'READING geometry in angstrom'
call read_geometry(MOLECULE)
if(mpi_rank.eq.0)print*,'ENUCL=',enucl

call build_gaussian_expansion_of_orbitals(mpi_rank)
if(i_print_orb.eq.1)call print_orb

if(mpi_rank.eq.0)write(*,*)'**********************************'
if(mpi_rank.eq.0)write(*,*)'Number of basis functions ',nbasis

allocate(rjacob(nw),rjacob1(nw,nbasis,nbasis),rjacob2(nw,nbasis,nbasis))

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
  call one_elect(i,i,Sij,Vij,Kij)
  sqrtSii_gaus(i)=dsqrt(dabs(Sij))
enddo

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

kkk=0
ijkl_gaus = 0.d0
do kcp=mpi_rank+1,nint, mpi_size

 i=is(kcp)
 k=ks(kcp)
 j=js(kcp)
 l=ls(kcp)

 if(dsqrt(precond(i,k)*precond(j,l)).lt.seuil_cauchy)then
  n_zero_cauchy=n_zero_cauchy+1
  ijkl_gaus(kcp)=0.d0
 else
  ijkl_gaus(kcp)=gauss_ijkl(i,k,j,l)
 endif

 k_sort=k_sort+1
 if( nint.ge.10.and.mod(k_sort,nint/10).eq.0)then
  kkk=kkk+1
  call cpu_time(t1)
  write(*,'(a,i3,a,f22.15)')' CPU TIME block',kkk,' of GAUSSIAN TWO-electron',t1-t0
  t0=t1
  k_sort=0
 endif

  rnorm=sqrtSii_gaus(i)*sqrtSii_gaus(k)*sqrtSii_gaus(j)*sqrtSii_gaus(l)
  if(.not.normalization_OA)rnorm=1.d0
  ijkl_gaus(kcp)=ijkl_gaus(kcp)/rnorm

enddo !kcp
!!!!!!!!!!!COMMENT IF NOT MPI
call mpi_allreduce(MPI_IN_PLACE,ijkl_gaus, nint, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!!!!!!!!!!COMMENT IF NOT MPI


if (mpi_rank.eq.0) then
  !! Write gaussian <ij|kl>
  open(unit=14,file='bielec_ao_gaus')
  rewind 14
  do kcp=1,nint
   i=is(kcp)
   k=ks(kcp)
   j=js(kcp)
   l=ls(kcp)
    write(14,'(4(I5,X),D22.15)') i,j,k,l, ijkl_gaus(kcp)
  enddo
  close(14)
endif

if(mpi_rank.eq.0)print*,'*********************************************************************'
if(mpi_rank.eq.0)print*,'IJKL REMOVED DUE TO CAUCHY-SCHWARZ=',n_zero_cauchy,' OVER',nint,'IJKL'
if(mpi_rank.eq.0)print*,'TOTAL NUMBER OF IJKL REMOVED    =',nint_zero+n_zero_cauchy,' LEFT',nint-(nint_zero+n_zero_cauchy)
if(mpi_rank.eq.0)print*,'*******************************************************************************'


!!******************
!! BEGIN MONTE CARLO
!!******************
!!
call cpu_time(t0)
call build_mapping_ijkl

do i=1,nbasis
 do k=1,nbasis
  a_ZV(i,k)=aread
 enddo
enddo

nbl=10

!! Determination of one-center bielectronic
do kcp=1,nint
 i=is(kcp)
 k=ks(kcp)
 j=js(kcp)
 l=ls(kcp)
 call compare_nuclei(nucleus_number(i),nucleus_number(j),nucleus_number(k),nucleus_number(l),ndiff)
 mono_center(kcp)=0
 if(ndiff.eq.1)mono_center(kcp)=1
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

allocate(pi_1(nw))
do kkk=1,nbl

 call cpu_time(t0)

 do kcp=1,nint
  e_S_ijkl(kcp)=0.d0
  e_G_ijkl(kcp)=0.d0
  e_test_ijkl(kcp)=0.d0
 enddo

 k_sort=0
 k_sort2=0

  call cpu_time(t1)
 do kk=1,npts_two_elec/nw

  call draw_configuration(r1,r2)
  do k=1,nbasis
   do i=k,nbasis
    if(i_tab_mc(i,k).eq.1)then
     call compute_r_tilde(2,i,k,r1,r2,rt1,rt2)
     call compute_u_tilde(2,i,k,rt1,rt2,ut1,ut2)
     call compute_densities(2,i,k,ut1,ut2,rho,rho_G,poly)
    endif
   enddo !k
  enddo !i
  call compute_jacobian_and_pi0(r1,r2,rjacob1,rjacob2,pi_0,nbasis)


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

     rjacob(1:nw) = rjacob1(1:nw,i,k)*rjacob2(1:nw,j,l)
     do kw=1,nw
      d_x(1) = ut1(1,kw,ik)-ut2(1,kw,jl)
      d_x(2) = ut1(2,kw,ik)-ut2(2,kw,jl)
      d_x(3) = ut1(3,kw,ik)-ut2(3,kw,jl)
      r12_2 = d_x(1)*d_x(1) + d_x(2)*d_x(2) + d_x(3)*d_x(3)
      factor=rjacob(kw)/pi_0(kw)
      weight  (kw)=factor* rho  (kw,ik,1)*rho  (kw,jl,2)
      weight_G(kw)=factor* rho_G(kw,ik,1)*rho_G(kw,jl,2)
      r12_inv(kw) = real( 1./sqrt(real(r12_2,4)), 8) !simple precision
!      r12_inv(kw) = 1.d0/dsqrt(r12_2)                !double precision
     enddo
     e_S_ijkl(kcp)=e_S_ijkl(kcp) + sum(weight  (1:nw)*r12_inv(1:nw))
     e_G_ijkl(kcp)=e_G_ijkl(kcp) + sum(weight_G(1:nw)*r12_inv(1:nw))
     e_test_ijkl(kcp)=e_test_ijkl(kcp) + sum(r12_inv(1:nw))

   endif
  enddo !kcp


  k_sort=k_sort+1
  k_sort2=k_sort2+nw
  if( npts_two_elec.ge.10.and.mod(k_sort,npts_two_elec/10).eq.0)then
   write(*,*)'mpi_rank= ',mpi_rank,' nsteps=',k_sort2
   k_sort=0
  endif

 enddo !npts_two_elec

 factor = 1.d0 / npts_two_elec
 do kcp=1,nint
  if(mono_center(kcp).eq.0)then
   e_S_ijkl(kcp)=e_S_ijkl(kcp) * factor
   e_G_ijkl(kcp)=e_G_ijkl(kcp) * factor
   e_test_ijkl(kcp)=e_test_ijkl(kcp) * factor
  endif
 enddo

 do kcp=1,nint
  i_value=1
  i=is(kcp)
  k=ks(kcp)
  j=js(kcp)
  l=ls(kcp)
  if(i_value.ne.0)then
   if(mono_center(kcp).eq.0)then
    e_tot_ijkl=e_S_ijkl(kcp)-e_G_ijkl(kcp)+ijkl_gaus(kcp)
    ijkl(kcp)=ijkl(kcp)+e_tot_ijkl
    ijkl2(kcp)=ijkl2(kcp)+e_tot_ijkl*e_tot_ijkl
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
  error_ijkl(kcp)=dsqrt( dabs(ijkl2(kcp)-ijkl(kcp)*ijkl(kcp)))/dsqrt(dfloat(nbl))
 else
  error_ijkl(kcp)=0.d0
 endif
enddo

!!!!!!!!!!COMMENT IF NOT MPI
 mpi_size_inv = 1.d0/dble(mpi_size)
 call MPI_AllReduce(ijkl, moy, nint, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
 moy(:)=moy(:) * mpi_size_inv

 moy2t(:) = (ijkl(:) -  moy(:))**2
 call MPI_reduce(moy2t, moy2, nint, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
 moy2(:)=dsqrt(moy2(:) * mpi_size_inv/ (dble(mpi_size-1)) )
!!!!!!!!!!COMMENT IF NOT MPI

 if(.not.MPI)then
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

 endif

 print*,'ZVMC DONE for mpi_rank=',mpi_rank

!!!!!!!!!!COMMENT IF NOT MPI
 call MPI_BARRIER (MPI_COMM_WORLD, ierr)
 call mpi_finalize(ierr)
!!!!!!!!!!COMMENT IF NOT MPI

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
!!!
r=dsqrt(dx**2+dy**2+dz**2)
phi_orb=dx**npower(1,i_orb)*dy**npower(2,i_orb)*dz**npower(3,i_orb)* u_orb(i_orb,r)
!!!
!r=dsqrt(dx*dx+dy*dy+dz*dz)
!phi_orb=u_orb(i_orb,r)
!do i=1,npower(1,i_orb)
!  phi_orb = phi_orb*dx
!end do
!do i=1,npower(2,i_orb)
!  phi_orb = phi_orb*dy
!end do
!do i=1,npower(3,i_orb)
!  phi_orb = phi_orb*dz
!end do
end

double precision function phi_gaus(i_orb,x)
include 'j.inc'
dimension x(3)
dx=x(1)-center(1,i_orb)
dy=x(2)-center(2,i_orb)
dz=x(3)-center(3,i_orb)
r=dsqrt(dx*dx+dy*dy+dz*dz)
phi_gaus=u_gauss(i_orb,r)
do i=1,npower(1,i_orb)
  phi_gaus = phi_gaus*dx
end do
do i=1,npower(2,i_orb)
  phi_gaus = phi_gaus*dy
end do
do i=1,npower(3,i_orb)
  phi_gaus = phi_gaus*dz
end do
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

subroutine draw_configuration(r1,r2)
include 'j.inc'
dimension r1(nw,3),r2(nw,3)
double precision, parameter :: f = dsqrt(2.d0)
do kw=1,nw
 do ll=1,3
  call random_number(r1(kw,ll))
  r1(kw,ll)=f*dierfc(2.d0-2.d0*r1(kw,ll))
  call random_number(r2(kw,ll))
  r2(kw,ll)=f*dierfc(2.d0-2.d0*r2(kw,ll))
 enddo
enddo
end

subroutine compute_pi0(r1,r2,pi_0)
include 'j.inc'
dimension r1(nw,3),r2(nw,3),pi_0(nw)
double precision, parameter :: pi=dacos(-1.d0)
double precision, parameter :: factor = 1.d0 / (2.d0*pi)**3
do kw=1,nw
 r1_mod_2=r1(kw,1)*r1(kw,1)+r1(kw,2)*r1(kw,2)+r1(kw,3)*r1(kw,3)
 r2_mod_2=r2(kw,1)*r2(kw,1)+r2(kw,2)*r2(kw,2)+r2(kw,3)*r2(kw,3)
 pi_0(kw)=dexp(-0.5d0*(r1_mod_2 + r2_mod_2)) * factor
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

 poly_i=1.d0
 do kk=1,npower(1,i)
  poly_i = poly_i * dx
 enddo
 do kk=1,npower(2,i)
  poly_i = poly_i * dy
 enddo
 do kk=1,npower(3,i)
  poly_i = poly_i * dz
 enddo
 r_i=dsqrt(dx*dx+dy*dy+dz*dz)

 dx=ut1(1,kw,ik)-center(1,k)
 dy=ut1(2,kw,ik)-center(2,k)
 dz=ut1(3,kw,ik)-center(3,k)

 poly_k=1.d0
 do kk=1,npower(1,k)
  poly_k = poly_k * dx
 enddo
 do kk=1,npower(2,k)
  poly_k = poly_k * dy
 enddo
 do kk=1,npower(3,k)
  poly_k = poly_k * dz
 enddo
 r_k=dsqrt(dx*dx+dy*dy+dz*dz)

 poly_ik = poly_i*poly_k
 rho  (kw,ik,1)=poly_ik*u_orb(i,r_i)*u_orb(k,r_k)
 poly (kw,ik,1)=poly_ik
 rho_G(kw,ik,1)=poly_ik*u_gauss(i,r_i)*u_gauss(k,r_k)

enddo
if(i_el.eq.1)return
do kw=1,nw
 rho  (kw,ik,2)=phi_orb (i,ut2(1,kw,ik))*phi_orb (k,ut2(1,kw,ik))
 rho_G(kw,ik,2)=phi_gaus(i,ut2(1,kw,ik))*phi_gaus(k,ut2(1,kw,ik))
enddo
end

subroutine compute_r_tilde(i_el,i,k,r1,r2,rt1,rt2)
include 'j.inc'
integer, intent(in) :: i_el, i, k
double precision, intent(in) :: r1(nw,3),r2(nw,3)
double precision, intent(out) :: rt1(nw,3),rt2(nw,3)
do kw=1,nw
 r1_mod=dsqrt(r1(kw,1)*r1(kw,1)+r1(kw,2)*r1(kw,2)+r1(kw,3)*r1(kw,3))
 call compute_f_fp_jacob(r1_mod,i,k,0,f1,f1p)
 do ll=1,3
  rt1(kw,ll)=f1*r1(kw,ll)
 enddo !!ll
enddo !!kw
if(i_el.eq.1)return
do kw=1,nw
 r2_mod=dsqrt(r2(kw,1)*r2(kw,1)+r2(kw,2)*r2(kw,2)+r2(kw,3)*r2(kw,3))
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
double precision factor
ik = (k-1)*nbasis_max+i
factor = 1.d0/dsqrt(g_min(i)+g_min(k))
do kw=1,nw
 do ll=1,3
  ut1(ll,kw,ik)= factor*rt1(kw,ll) + G_center(ll,i,k)
 enddo !!ll
enddo !!kw
if(i_el.eq.1)return
do kw=1,nw
 do ll=1,3
  ut2(ll,kw,ik)= factor*rt2(kw,ll) + G_center(ll,i,k)
 enddo !!ll
enddo !!kw
end

!TODO
subroutine compute_jacobian(r1,rjacob,n)
include 'j.inc'
dimension r1(nw,3)
dimension rjacob(nw,n,n)
do k=1,n
  do i=k,n
    factor_ik = a_ZV(i,k)/(g_min(i)+g_min(k))
    do kw=1,nw
      f1=factor_ik*dsqrt(r1(kw,1)*r1(kw,1)+r1(kw,2)*r1(kw,2)+r1(kw,3)*r1(kw,3))
      rjacob(kw,i,k)=0.25d0*dabs(f1*f1*f1)
      rjacob(kw,k,i)= rjacob(kw,i,k)
    enddo
  enddo
enddo
end

subroutine compute_jacobian_and_pi0(r1,r2,rjacob1,rjacob2,pi_0,n)
include 'j.inc'
double precision :: r1(nw,3), r2(nw,3),pi_0(nw), rjacob1(nw,n,n), rjacob2(nw,n,n)
double precision :: d1(nw), d2(nw)
double precision, parameter :: factor = 1.d0 / (2.d0*dacos(-1.d0))**3
do kw=1,nw
 r1_mod_2=r1(kw,1)*r1(kw,1)+r1(kw,2)*r1(kw,2)+r1(kw,3)*r1(kw,3)
 r2_mod_2=r2(kw,1)*r2(kw,1)+r2(kw,2)*r2(kw,2)+r2(kw,3)*r2(kw,3)
 d1(kw) = dsqrt(r1_mod_2)
 d2(kw) = dsqrt(r2_mod_2)
 pi_0(kw)= dexp(-0.5d0*(r1_mod_2 + r2_mod_2)) * factor
enddo
do k=1,n
  do i=k,n
    factor_ik = a_ZV(i,k)/(g_min(i)+g_min(k))
    do kw=1,nw
      f1=factor_ik*d1(kw)
      rjacob1(kw,i,k)=0.25d0*dabs(f1*f1*f1)
      rjacob1(kw,k,i)= rjacob1(kw,i,k)

      f2=factor_ik*d2(kw)
      rjacob2(kw,i,k)=0.25d0*dabs(f2*f2*f2)
      rjacob2(kw,k,i)= rjacob2(kw,i,k)
    enddo
  enddo
enddo
end

subroutine compute_f_fp_jacob(r,i,k,ic,f,fp)
include 'j.inc'

 factor = 0.5d0/dsqrt(g_min(i)+g_min(k))
 f=a_ZV(i,k)*r * factor
 if(ic.eq.0)return
 fp=a_ZV(i,k) * factor

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

 call read_fit_SMILES

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

     call build_c_g_gauss_STO(nbasis)
     if(mpi_rank.eq.0)call write_STO_in_file_info_basis(nbasis)

    enddo !m

   endif  !5G

  enddo ! k=1,n_b(i_atom)
 enddo ! i=1,number_nuclei

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
 do mm=1,n_contract(i_orb)
  contrib=c_contract(mm,i_orb)*r**n_sto(i_orb)*dexp(-g_contract(mm,i_orb)*r   )
  u_orb=u_orb+contrib
 enddo  !mm
end
!! derivative of u_orb
double precision function up_orb(i_orb,r)
include 'j.inc'
up_orb=0.d0
 do mm=1,n_contract(i_orb)
  contrib=(n_sto(i_orb)/r-g_contract(mm,i_orb))*c_contract(mm,i_orb)*r**n_sto(i_orb)*dexp(-g_contract(mm,i_orb)*r   )
  up_orb=up_orb+contrib
 enddo  !mm
end

