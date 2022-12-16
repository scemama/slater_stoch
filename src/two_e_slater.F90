program integrals
  use common_data

#ifdef HAVE_TREXIO
  use trexio
#endif

  implicit none
#ifdef HAVE_MPI
  include 'mpif.h'
#endif

  double precision, parameter    :: seuil_cauchy = 1.d-15

  character(80)                  :: charabid, MOLECULE

  double precision, allocatable  :: precond(:,:)

  double precision, allocatable  :: ijkl_gaus(:)
  double precision, allocatable  :: ijkl(:),ijkl2(:)

  double precision, allocatable  :: e_S_ijkl(:)
  double precision, allocatable  :: e_G_ijkl(:)
  double precision, allocatable  :: error_ijkl(:)
  integer, allocatable           :: mono_center(:)
  double precision, allocatable  :: moy(:), moy2(:), moy2t(:)
  double precision, allocatable  :: G_center(:,:,:)

  double precision               :: Sij, Vij, Kij

  !! MONTE CARLO PART
  integer, allocatable           :: i_tab_mc(:,:)
  double precision               :: d_x(4)

  double precision               :: r12_inv(nw), r1(nw,3), r2(nw,3)
  double precision               :: pi_0(nw),pot(nw),rkin(nw),rkin_G(nw)
  double precision               :: weight(nw),weight_kin(nw)
  double precision               :: weight_G(nw),weight_kin_G(nw)

  double precision, allocatable  :: rt1(:,:),rt2(:,:)
  double precision, allocatable  :: ut1(:,:,:),ut2(:,:,:)
  double precision, allocatable  :: rho(:,:,:),poly(:,:,:)
  double precision, allocatable  :: rho_G(:,:,:)
  double precision, allocatable  :: rjacob(:), rjacob1(:,:,:), rjacob2(:,:,:)
  double precision, allocatable  :: pi_1(:)

  character(128)                 :: filename_in
  character(128)                 :: filename_out_gaus_ijkl
  character(128)                 :: filename_out_ijkl
  character(128)                 :: filename_basis

  integer                        :: mpi_rank, mpi_size

#ifdef HAVE_MPI
  integer                        :: ierr, irequest
  integer                        :: mpi_status(MPI_STATUS_SIZE)
  logical                        :: mpi_flag
  double precision               :: mpi_size_inv
#endif

  integer*4                      :: seed(33)
  integer*4                      :: put(33)

  integer                        :: i_rand, num_simulation, npts_two_elec
  integer                        :: i, j, k, l, k_sort2, kkk, ll, q, rank
  integer*8                      :: kcp
  integer                        :: nint_zero, i_value, n_zero_cauchy, num
  integer                        :: nbl, ndiff, kk, ik, kw, jl, n_ijkl
  double precision               :: t0, t1, factor, r2_mod, r1_mod, r12_2
  double precision               :: e_tot_ijkl, error_max, errmoy_ijkl, dmoy_ijkl
  double precision               :: a_ijkln, diff_ijkl, diff_ijklmax, f

  double precision, external     :: gauss_ijkl
  integer :: xi(8), xj(8), xk(8), xl(8), ii, jj


#ifdef HAVE_TREXIO
  character*(128)   :: trexio_filename
  integer           :: rc
  integer(trexio_t) :: trexio_file
  character(128)    :: err_message
  integer, parameter:: BUFSIZE = 32768
  integer           :: indx(4,BUFSIZE)
  double precision  :: vals(BUFSIZE)
  integer*8         :: icount, offset
#endif


#ifdef HAVE_MPI
  call mpi_init(ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_rank, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_size, ierr)
  if(mpi_rank.eq.0)write(*,*)'mpi_size=',mpi_size
#else
  mpi_rank = 0
  mpi_size = 1
#endif

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

  seed(1) = 1+(num_simulation-1)*10000
#ifdef HAVE_MPI
  seed(1) = mpi_rank+seed(1)
#endif
  do i_rand=2,size(seed)
    seed(i_rand) = i_rand*seed(1)
  enddo
  call random_seed(put=seed)
  call random_number(r1)


  !*****************
  ! BEGIN READ INPUT
  !*****************
  if(mpi_rank.eq.0)print*,'MOLECULE?'
  read(5,'(a80)')charabid
  read(5,*)MOLECULE
  if(mpi_rank.eq.0)write(*,*)MOLECULE
  read(5,'(a80)')charabid
  if(mpi_rank.eq.0)write(6,'(a80)')charabid
  if(mpi_rank.eq.0)print*,'file name for basis set?'
  read(5,*)filename_basis
  if(mpi_rank.eq.0)print*,'*******'
  if(mpi_rank.eq.0)write(*,*)trim(filename_basis)
  if(mpi_rank.eq.0)print*,'*******'
  read(5,'(a80)')charabid
  if(mpi_rank.eq.0)write(6,'(a80)')charabid
  read(5,*)ng0,i_add_2p,i_add_3d,g_thr,i_add_large_g
  if(mpi_rank.eq.0)write(6,*)ng0,i_add_2p,i_add_3d,g_thr,i_add_large_g
  read(5,'(a80)')charabid
  if(mpi_rank.eq.0)write(6,'(a80)')charabid
  read(5,*)npts_two_elec
  if(mpi_rank.eq.0)write(*,*)'Number of MC steps for two-electron integrals=',npts_two_elec, ' times number of blocks'
  if(mod(npts_two_elec,nw).ne.0)stop 'npts_two_elec must be a multiple of nw'
  do i=1,2
    read(5,'(a80)')charabid
    if(mpi_rank.eq.0)write(6,'(a80)')charabid
  enddo
  read(5,*)a_ZV
  if(mpi_rank.eq.0)print*,'a_ZV=',a_ZV

  !*****************
  ! END READ INPUT
  !*****************

  if(mpi_rank.eq.0)print*,'READING geometry in angstrom'
  call read_geometry(MOLECULE)
  if(mpi_rank.eq.0)print*,'ENUCL=',enucl

  call read_basis(filename_basis)
  call build_gaussian_expansion_of_orbitals()

  if(mpi_rank.eq.0)write(*,*)'**********************************'
  if(mpi_rank.eq.0)write(*,*)'Number of basis functions ',nbasis

  allocate(rjacob(nw),rjacob1(nw,nbasis,nbasis),rjacob2(nw,nbasis,nbasis))
  allocate(ut1(nw,3,nbasis*nbasis),ut2(nw,3,nbasis*nbasis))
  allocate(rho(nw,nbasis*nbasis,2),poly(nw,nbasis*nbasis,2))
  allocate(rho_G(nw,nbasis*nbasis,2))


  !! build arrays is(kcp),js(kcp),ks(kcp),ls(kcp) ;  kcp=1,nint and  i=1,nbasis (idem for j,k,l)
  !! nint= total number of two-electron integrals
  call build_mapping_ijkl(nint)

  allocate(ijkl_gaus(nint))
  allocate(ijkl(nint),ijkl2(nint))
  allocate(e_S_ijkl(nint))
  allocate(e_G_ijkl(nint))
  allocate(error_ijkl(nint))
  allocate(mono_center(nint))
  allocate(moy(nint), moy2(nint), moy2t(nint))
  ijkl = 0.d0
  ijkl2 = 0.d0


  if(mpi_rank.eq.0)print*,'**************************************************************************'
  if(mpi_rank.eq.0)print*,'Analytic calculation of APPROXIMATE GAUSSIAN one-electron integrals for ZV'
  if(mpi_rank.eq.0)print*,'**************************************************************************'

  nint_zero=0
  i_value=1
  n_zero_cauchy=0
  kkk=0

  num=nbasis**2
  do i=1,nbasis
    call one_elect(i,i,Sij,Vij,Kij)
  enddo

  if(mpi_rank.eq.0)print*,'*******************************************************************'
  if(mpi_rank.eq.0)print*,'ANALYTIC CALCULATION OF APPROXIMATE GAUSSIAN TWO-ELECTRON INTEGRALS'
  if(mpi_rank.eq.0)print*,'*******************************************************************'

  if (mpi_rank == 0) then
    call count_multi_center_integrals
  end if

  allocate(precond(nbasis,nbasis))
  do k=1,nbasis
    do i=1,nbasis
      precond(i,k)=gauss_ijkl(i,k,i,k)
    enddo
  enddo

  k_sort2=0
  nint_zero=0
  i_value=1

  n_zero_cauchy=0
  nint_zero=0

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
    else
      ijkl_gaus(kcp)=gauss_ijkl(i,k,j,l)
!      ijkl_gaus(kcp)=0.d0
    endif

  enddo !kcp

#ifdef HAVE_MPI
  call mpi_allreduce(MPI_IN_PLACE,n_zero_cauchy, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call mpi_allreduce(MPI_IN_PLACE,ijkl_gaus, nint, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

  call cpu_time(t1)
  if (mpi_rank==0) print *, 'Time for Gaussian integrals : ', t1-t0, ' seconds'

!  if (mpi_rank.eq.0) then
!    !! Write gaussian <ij|kl>
!    open(unit=14,file='bielec_ao_gaus')
!    rewind 14
!    do kcp=1,nint
!      i=is(kcp)
!      k=ks(kcp)
!      j=js(kcp)
!      l=ls(kcp)
!      write(14,'(4(I5,X),D22.15)') i,j,k,l, ijkl_gaus(kcp)
!    enddo
!    close(14)
!  endif

  if(mpi_rank.eq.0)print*,'*********************************************************************'
  if(mpi_rank.eq.0)print*,'IJKL REMOVED DUE TO CAUCHY-SCHWARZ=',n_zero_cauchy,' OVER',nint,'IJKL'
  if(mpi_rank.eq.0)print*,'TOTAL NUMBER OF IJKL REMOVED    =',nint_zero+n_zero_cauchy,' LEFT',nint-(nint_zero+n_zero_cauchy)
  if(mpi_rank.eq.0)print*,'*******************************************************************************'


  !!******************
  !! BEGIN MONTE CARLO
  !!******************
  !!
  call cpu_time(t0)

  nbl=10

  mono_center(:) = 0
  !! Determination of one-center bielectronic
!  do kcp=1,nint
!    i=is(kcp)
!    k=ks(kcp)
!    j=js(kcp)
!    l=ls(kcp)
!    call compare_nuclei(nucleus_number(i),nucleus_number(j),nucleus_number(k),nucleus_number(l),ndiff)
!    if(ndiff.eq.1)mono_center(kcp)=1
!  enddo

  if(mpi_rank.eq.0)then
    print*,'*********************************'
    print*,'ZVMC EXACT TWO-ELECTRON INTEGRALS'
    print*,'*********************************'
    print*,' Number of Monte Carlo steps= ',nbl,' X ',npts_two_elec
    print*,'*******************************************************'
    if(mod(npts_two_elec,nw).ne.0)stop 'npts_two_elec must be a multiple of nw'
  endif

  allocate(G_center(3, nbasis, nbasis))
  do j=1,nbasis
    do i=1,nbasis
      do ll=1,3
        G_center(ll,i,j)= (g_slater(i)*center(ll,i)+g_slater(j)*center(ll,j))/(g_slater(i)+g_slater(j))
      enddo
    enddo
  enddo

  !**************************************************
!  !! Exact calculation of monocenter Slater integrals
!  do kcp=1,nint
!    i=is(kcp)
!    k=ks(kcp)
!    j=js(kcp)
!    l=ls(kcp)
!    if(mono_center(kcp).eq.1)then
!      call compute_int_slater(i,k,j,l,ijkl(kcp))
!    endif
!  enddo
!  !! end calculation ***********************************

!  call cpu_time(t1)
!  if (mpi_rank == 0) print *, 'Time for one-center integrals: ', t1-t0, ' seconds'

!  do kcp=1,nint
!    if(mono_center(kcp).eq.0)then
!      ijkl(kcp)=0.d0
!      ijkl2(kcp)=0.d0
!    endif
!  enddo

  !! Determine which i j are used in computation of W_S and W_G
  !! i_tab_mc(i,k)=1  W_S and W_G can be computed
  allocate(i_tab_mc(nbasis,nbasis))
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
    enddo

    k_sort2=0

    call cpu_time(t1)
    do kk=1,npts_two_elec/nw

      call draw_configuration(r1,r2,nw)
      do k=1,nbasis
        do i=k,nbasis
          if(i_tab_mc(i,k).eq.1)then
            ik = (k-1)*nbasis+i
            factor = 0.5d0 * a_ZV /(g_slater(i)+g_slater(k))
            do kw=1,nw
              r1_mod=dsqrt(r1(kw,1)*r1(kw,1)+r1(kw,2)*r1(kw,2)+r1(kw,3)*r1(kw,3))
              r2_mod=dsqrt(r2(kw,1)*r2(kw,1)+r2(kw,2)*r2(kw,2)+r2(kw,3)*r2(kw,3))
              ut1(kw,1:3,ik)= factor*r1(kw,1:3) * r1_mod + G_center(1:3,i,k)
              ut2(kw,1:3,ik)= factor*r2(kw,1:3) * r2_mod + G_center(1:3,i,k)
            enddo !!kw

            call compute_densities(i,k,ut1,ut2,rho,rho_G,poly)
          endif
        enddo !k
      enddo !i
      call compute_jacobian_and_pi0(r1,r2,rjacob1,rjacob2,pi_0,nbasis)


      do kcp=1,nint

        i=is(kcp)
        k=ks(kcp)
        ik = (k-1)*nbasis+i
        j=js(kcp)
        l=ls(kcp)
        jl = (l-1)*nbasis+j

        if(dsqrt(precond(i,k)*precond(j,l)).lt.seuil_cauchy)then

          n_zero_cauchy=n_zero_cauchy+1

        else

          rjacob(1:nw) = rjacob1(1:nw,i,k)*rjacob2(1:nw,j,l)
          do kw=1,nw
            d_x(1) = ut1(kw,1,ik)-ut2(kw,1,jl)
            d_x(2) = ut1(kw,2,ik)-ut2(kw,2,jl)
            d_x(3) = ut1(kw,3,ik)-ut2(kw,3,jl)
            r12_2 = d_x(1)*d_x(1) + d_x(2)*d_x(2) + d_x(3)*d_x(3)

            factor=rjacob(kw)/pi_0(kw)
            weight  (kw)=factor* rho  (kw,ik,1)*rho  (kw,jl,2)
            weight_G(kw)=factor* rho_G(kw,ik,1)*rho_G(kw,jl,2)
            r12_inv(kw) = real( 1./sqrt(real(r12_2,4)), 8) !simple precision
            !      r12_inv(kw) = 1.d0/dsqrt(r12_2)                !double precision
          enddo
          e_S_ijkl(kcp)=e_S_ijkl(kcp) + sum(weight  (1:nw)*r12_inv(1:nw))
          e_G_ijkl(kcp)=e_G_ijkl(kcp) + sum(weight_G(1:nw)*r12_inv(1:nw))

        endif
      enddo !kcp

      k_sort2=k_sort2+nw
      if( npts_two_elec.ge.10.and.mod(kcp,npts_two_elec/10).eq.0)then
        write(*,*)'mpi_rank= ',mpi_rank,' nsteps=',k_sort2
      endif

    enddo !npts_two_elec

    factor = 1.d0 / dble(npts_two_elec)
    do kcp=1,nint
      if(mono_center(kcp).eq.0)then
        e_S_ijkl(kcp)=e_S_ijkl(kcp) * factor
        e_G_ijkl(kcp)=e_G_ijkl(kcp) * factor
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
          e_tot_ijkl=e_S_ijkl(kcp)-e_G_ijkl(kcp) !+ijkl_gaus(kcp)
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


  if (mpi_size > 1) then

#ifdef HAVE_MPI
    mpi_size_inv = 1.d0/dble(mpi_size)
    call MPI_AllReduce(ijkl, moy, nint, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
    moy(:)=moy(:) * mpi_size_inv

    moy2t(:) = (ijkl(:) -  moy(:))**2
    call MPI_reduce(moy2t, moy2, nint, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
    moy2(:)=dsqrt(moy2(:) * mpi_size_inv/ (dble(mpi_size-1)) )
#else
    stop 'mpi_size should not be >1'
#endif

  else

    moy(:)=ijkl(:)
    moy2(:)=error_ijkl(:)

  endif

  if (mpi_rank == 0) then

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
      if(moy2(kcp).gt.error_max)error_max=moy2(kcp)

      if( dabs(moy(kcp)).gt.1.d-6.and.(moy2(kcp).ne.0.d0))then
        errmoy_ijkl=errmoy_ijkl+moy2(kcp)
        dmoy_ijkl=dmoy_ijkl+1.d0
        diff_ijkl=dabs(moy(kcp))/moy2(kcp)
        if(diff_ijkl.gt.diff_ijklmax)diff_ijklmax=diff_ijkl
        if(diff_ijkl.gt.3.d0)then
          a_ijkln=a_ijkln+diff_ijkl
          n_ijkl=n_ijkl+1
        endif
      endif

    enddo

    write(*,*)'Error mean ijkl=',errmoy_ijkl/dmoy_ijkl,' error_max ',error_max
  end if

  if (mpi_rank == 0) print *, 'Cleaning ERI matrix'

  if (mpi_rank == 0) then
    q=10
    rank=min(20*nbasis, nbasis*nbasis)
    call svd_clean(moy, ijkl_gaus, nint, is, js, ks, ls, nbasis, rank, q, mpi_rank, mpi_size)
    moy(:) = moy(:)+ijkl_gaus(:)
  endif

  if (mpi_rank == 0) then


#ifdef HAVE_TREXIO

    trexio_filename = trim(MOLECULE)//'.h5'
    trexio_file = trexio_open(trexio_filename, 'w', TREXIO_HDF5, rc)
    call trexio_assert(rc, TREXIO_SUCCESS)

    icount = 0_8
    offset = 0_8
    do kcp=1,nint
      if (dabs(moy(kcp)) < 1.d-15) cycle
      icount = icount + 1_8
      vals(icount) = moy(kcp)
      indx(1,icount) = is(kcp)
      indx(2,icount) = js(kcp)
      indx(3,icount) = ks(kcp)
      indx(4,icount) = ls(kcp)
      if (icount == BUFSIZE) then
        rc = trexio_write_ao_2e_int_eri(trexio_file, offset, icount, indx, vals)
        call trexio_assert(rc, TREXIO_SUCCESS)
        offset = offset + icount
        icount = 0_8
      endif
    enddo
    if (icount > 0_8) then
      rc = trexio_write_ao_2e_int_eri(trexio_file, offset, icount, indx, vals)
      call trexio_assert(rc, TREXIO_SUCCESS)
    endif

    rc = trexio_close(trexio_file)
    call trexio_assert(rc, TREXIO_SUCCESS)

#else

    open(unit=104,file=filename_out_ijkl)
    do kcp=1,nint
      i=is(kcp)
      k=ks(kcp)
      j=js(kcp)
      l=ls(kcp)
      write(104,'(4(I5,X),2D22.15)') i,j,k,l, moy(kcp), moy2(kcp)
    enddo
    close(104)

#endif

  endif

  print*,'ZVMC DONE for mpi_rank=',mpi_rank

#ifdef HAVE_MPI
  call MPI_BARRIER (MPI_COMM_WORLD, ierr)
  call mpi_finalize(ierr)
#endif


end

