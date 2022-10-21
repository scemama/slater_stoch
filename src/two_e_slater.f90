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
  double precision               :: d_x(4)

  double precision, allocatable  :: r12_inv(:),r1(:,:),r2(:,:)
  double precision, allocatable  :: pi_0(:),pot(:),rkin(:),rkin_G(:)
  double precision, allocatable  :: rt1(:,:),rt2(:,:)
  double precision, allocatable  :: ut1(:,:,:),ut2(:,:,:)
  double precision, allocatable  :: rho(:,:,:),poly(:,:,:)
  double precision, allocatable  :: weight(:),weight_kin(:)
  double precision, allocatable  :: rho_G(:,:,:)
  double precision, allocatable  :: weight_G(:),weight_kin_G(:)
  double precision, allocatable  :: rjacob(:), rjacob1(:,:,:), rjacob2(:,:,:)
  double precision, allocatable  :: pi_1(:)

  character*(128)                :: filename_in
  character*(128)                :: filename_out_gaus_ijkl
  character*(128)                :: filename_out_ijkl
  character*(128)                :: filename_basis

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
  read(5,*)
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
            call compute_r_tilde(i,k,r1,r2,rt1,rt2)
            call compute_u_tilde(i,k,rt1,rt2,ut1,ut2)
            call compute_densities(i,k,ut1,ut2,rho,rho_G,poly)
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
