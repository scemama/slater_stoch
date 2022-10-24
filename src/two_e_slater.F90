program integrals
  include 'j.inc'
#ifdef HAVE_MPI
  include 'mpif.h'
#endif

  double precision, parameter    :: seuil_cauchy = 1.d-15

  character(80)                  :: charabid, MOLECULE

  double precision               :: precond(nbasis_max,nbasis_max)

  !! BEGIN BIG ARRAYS nint_max
  double precision               :: ijkl_gaus(nint_max)
  double precision               :: ijkl(nint_max),ijkl2(nint_max)

  double precision               :: e_S_ijkl(nint_max)
  double precision               :: e_G_ijkl(nint_max)
  double precision               :: e_test_ijkl(nint_max)
  double precision               :: error_ijkl(nint_max)
  double precision               :: mono_center(nint_max)
  double precision               :: moy(nint_max), moy2(nint_max), moy2t(nint_max)
  double precision               :: Sij, Vij, Kij
  !! END BIG ARRAYS nint_max

  !! MONTE CARLO PART
  integer                        :: i_tab_mc(nbasis_max,nbasis_max)
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
  integer                        :: i, j, k, l, k_sort, k_sort2, kkk, kcp, ll
  integer                        :: nint_zero, i_value, n_zero_cauchy, num
  integer                        :: nbl, ndiff, kk, ik, kw, jl, n_ijkl
  double precision               :: aread, t0, t1, factor, r2_mod, r1_mod, r12_2
  double precision               :: e_tot_ijkl, error_max, errmoy_ijkl, dmoy_ijkl
  double precision               :: a_ijkln, diff_ijkl, diff_ijklmax

  double precision, external     :: gauss_ijkl


  allocate(r12_inv(nw),r1(nw,3),r2(nw,3))
  allocate(pi_0(nw),pot(nw),rkin(nw),rkin_G(nw))
  allocate(ut1(3,nw,nbasis_max*nbasis_max),ut2(3,nw,nbasis_max*nbasis_max))
  allocate(rho(nw,nbasis_max*nbasis_max,2),poly(nw,nbasis_max*nbasis_max,2))
  allocate(weight(nw),weight_kin(nw) )
  allocate(rho_G(nw,nbasis_max*nbasis_max,2))
  allocate(weight_G(nw),weight_kin_G(nw))

  mpi_rank = 0
  mpi_size = 1
  ijkl = 0.d0
  ijkl2 = 0.d0

#ifdef HAVE_MPI
  call mpi_init(ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_rank, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_size, ierr)
  call sleep(mpi_rank/20.)
  if(mpi_rank.eq.0)write(*,*)'mpi_size=',mpi_size
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

#ifdef HAVE_MPI
  seed(1) = mpi_rank+1+(num_simulation-1)*10000
#endif
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
  read(5,*)aread
  if(mpi_rank.eq.0)print*,'aread=',aread

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
  nint_zero=0
  i_value=1
  n_zero_cauchy=0
  k_sort=0
  call cpu_time(t0)
  kkk=0

  num=nbasis**2
  do i=1,nbasis
    call one_elect(i,i,Sij,Vij,Kij)
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

  enddo !kcp

#ifdef HAVE_MPI
  call mpi_allreduce(MPI_IN_PLACE,ijkl_gaus, nint, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif


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
    endif
  enddo
  !! end calculation ***********************************


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
            ik = (k-1)*nbasis_max+i
            factor = 0.5d0 * a_ZV(i,k) /(g_min(i)+g_min(k))
            do kw=1,nw
              r1_mod=dsqrt(r1(kw,1)*r1(kw,1)+r1(kw,2)*r1(kw,2)+r1(kw,3)*r1(kw,3))
              r2_mod=dsqrt(r2(kw,1)*r2(kw,1)+r2(kw,2)*r2(kw,2)+r2(kw,3)*r2(kw,3))
              ut1(1:3,kw,ik)= factor*r1(kw,1:3) * r1_mod + G_center(1:3,i,k)
              ut2(1:3,kw,ik)= factor*r2(kw,1:3) * r2_mod + G_center(1:3,i,k)
            enddo !!kw

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

#ifdef HAVE_MPI
  call MPI_BARRIER (MPI_COMM_WORLD, ierr)
  call mpi_finalize(ierr)
#endif

end
