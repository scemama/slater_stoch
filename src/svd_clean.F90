subroutine svd_clean(moy, nint, is, js, ks, ls, nbasis, rank, q)
#ifdef HAVE_MMAP
  use mmap_module
  use iso_c_binding
#endif
  implicit none
#ifdef HAVE_MPI
#include <mpif.h>
#endif

  ! Removes negative eigenvalues of the ERI matrix

  integer*8, intent(in)          :: nint
  integer, intent(in)            :: nbasis, q
  integer, intent(inout)         :: rank
  double precision, intent(inout) :: moy(nint)
  integer, dimension(nint), intent(in) :: is, js, ks, ls

  double precision, allocatable  :: H_diag(:)
  double precision, allocatable  :: V(:,:), E(:)

  integer                        :: NMAX
  integer                        :: i, j, k, l, kk, iter, shift, n, mmax
  integer*8                      :: ii, jj, kcp
  double precision               :: r1, r2
  double precision, allocatable  :: U(:,:), D(:), Vt(:,:)
  double precision, pointer      :: W(:,:), W_work(:,:)
  double precision, allocatable  :: UY(:,:)
  double precision, allocatable  :: P(:,:), Y(:,:), Z(:,:), rd(:,:)
  double precision, parameter    :: dtwo_pi = 2.d0*dacos(-1.d0)


  double precision, external     :: dnrm2, ddot
  logical                        :: converged

  integer                        :: lwork, info, nremoved, npass
  double precision, allocatable  :: moy_work(:)
  integer                        :: w0, w1, wsize

#ifdef HAVE_MPI
  integer                        :: ierr
  integer                        :: mpi_status(MPI_STATUS_SIZE)
  integer                        :: mpi_rank, mpi_size
#endif

#ifdef HAVE_MMAP
  type(c_ptr) :: ptr_W, ptr_W_work
  integer     :: fd_W, fd_W_work
  logical, parameter     :: do_mmap = .True.
#else
  logical, parameter     :: do_mmap = .False.
#endif

  n = nbasis*nbasis
  nremoved = 0

#ifdef HAVE_MPI
    call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_rank, ierr)
    if (ierr /= MPI_SUCCESS) then
       print *, 'MPI error:', __FILE__, ':', __LINE__
       stop -1
    endif
    call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_size, ierr)
    if (ierr /= MPI_SUCCESS) then
       print *, 'MPI error:', __FILE__, ':', __LINE__
       stop -1
    endif
#else
    mpi_rank = 0
    mpi_size = 1
#endif

    w0 = n/mpi_size+1
    w1 = min(n,(mpi_rank+1)*w0)
    w0 = (mpi_rank)*w0+1
    wsize = w1-w0+1

    if (mpi_size == 1) then

#ifdef HAVE_MMAP
      call mmap('W.bin', (/ int(n,8), int(n,8) /), 8, fd_W, .False., ptr_w)
      call c_f_pointer(ptr_W, W, (/ n, n /))

      call mmap('W_work.bin', (/ int(n,8), int(n,8) /), 8, fd_W_work, .False., ptr_w_work)
      call c_f_pointer(ptr_W_work, W_work, (/ n, n /))
#else
      allocate(W(n,n), W_work(n,n))
#endif

    else

      allocate(W(n,w0:w1), W_work(n,w0:w1))

    endif

    W_work(:,:) = 0.d0
    do kcp=1,nint
      i = is(kcp) ; j = js(kcp) ; k = ks(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) W_work(ii,jj) = moy(kcp)

      i = ks(kcp) ; j = js(kcp) ; k = is(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) W_work(ii,jj) = moy(kcp)

      i = ks(kcp) ; j = ls(kcp) ; k = is(kcp) ; l = js(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) W_work(ii,jj) = moy(kcp)

      i = is(kcp) ; j = ls(kcp) ; k = ks(kcp) ; l = js(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) W_work(ii,jj) = moy(kcp)

      i = js(kcp) ; j = is(kcp) ; k = ls(kcp) ; l = ks(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) W_work(ii,jj) = moy(kcp)

      i = js(kcp) ; j = ks(kcp) ; k = ls(kcp) ; l = is(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) W_work(ii,jj) = moy(kcp)

      i = ls(kcp) ; j = ks(kcp) ; k = js(kcp) ; l = is(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) W_work(ii,jj) = moy(kcp)

      i = ls(kcp) ; j = is(kcp) ; k = js(kcp) ; l = ks(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) W_work(ii,jj) = moy(kcp)
    enddo

    W(:,:) = W_work(:,:)

   allocate(U(n,rank), P(n,rank), Z(n,rank))
   allocate(D(rank))
   allocate(Y(rank,n), Vt(rank,n))
   allocate(UY(rank,rank))

    do npass=1,n,rank/2
      if (mpi_rank == 0) then
        print *, ''
        print *, 'Vectors ', npass, '--', min(n,npass+rank/2-1), ' / ', n
      endif
#ifdef HAVE_MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
#endif

      ! ---
      ! P is a normal random matrix (n,rank)
      if (mpi_rank == 0) then
        allocate(rd(n,2))
        do k=1,rank
          call random_number(rd)
          do i=1,n
            r1 = dsqrt(-2.d0*dlog(rd(i,1)))
            r2 = dtwo_pi*rd(i,2)
            P(i,k) = r1*dcos(r2)
          enddo
          r1 = dnrm2(n,P(1:n,k),1)
          call dscal(n,1.d0/r1,P(1:n,k),1)
        enddo
        deallocate(rd)
      endif
#ifdef HAVE_MPI
      call MPI_BCAST(P, n*rank, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
#endif


      ! Z(n,rank) = W_work(n,n).P(n,rank)
      if (wsize > 0) then
        call dgemm('N', 'N', n, rank, wsize, 1.d0, &
          W_work(1,w0), size(W_work,1), &
          P(w0,1), size(P,1), 0.d0, Z, size(Z,1))
      else
        Z = 0.d0
      endif

      ! QR factorization of Z
#ifdef HAVE_MPI
      call MPI_Reduce(Z, U, n*rank, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
      if (mpi_rank == 0) then
        Z(:,:) = U(:,:)
#endif
        call ortho_qr(Z,size(Z,1),n,rank)
#ifdef HAVE_MPI
      endif
      call MPI_BCAST(Z, n*rank, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
#endif

      ! Power iterations
      do i=1,q
        ! P(n,rank) = W_work^t(n,n).Z(n,rank)
        P=0.d0
        if (wsize > 0) then
          call dgemm('T', 'N', wsize, rank, n, 1.d0, &
              W_work(1,w0), size(W_work,1), Z, size(Z,1), 0.d0, &
              P(w0,1), size(P,1))
        endif
#ifdef HAVE_MPI
        call MPI_AllReduce(MPI_IN_PLACE, P, n*rank, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
        if (ierr /= MPI_SUCCESS) then
           print *, 'MPI error:', __FILE__, ':', __LINE__
           stop -1
        endif
#endif

        ! Z(n,rank) = W_work(n,n).P(n,rank)
        if (wsize > 0) then
          call dgemm('N', 'N', n, rank, wsize, 1.d0, &
              W_work(1,w0), size(W_work,1), P(w0,1), size(P,1), 0.d0, &
              Z, size(Z,1))
        else
          Z = 0.d0
        endif

#ifdef HAVE_MPI
        call MPI_Reduce(Z, U, n*rank, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
        if (ierr /= MPI_SUCCESS) then
           print *, 'MPI error:', __FILE__, ':', __LINE__
           stop -1
        endif
        if (mpi_rank == 0) then
          Z(:,:) = U(:,:)
#endif
          call ortho_qr(Z,size(Z,1),n,rank)
#ifdef HAVE_MPI
        endif
        call MPI_BCAST(Z, n*rank, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
           print *, 'MPI error:', __FILE__, ':', __LINE__
           stop -1
        endif
      enddo
#endif


      ! Y(rank,n) = Zt(rank,n).W_work(n,n)
      Y = 0.d0
      if (wsize > 0) then
        call dgemm('T', 'N', rank, wsize, n, 1.d0, &
          Z, size(Z,1), W_work(1,w0), size(W_work,1), 0.d0, &
          Y(1,w0), size(Y,1))
      endif

#ifdef HAVE_MPI
      call MPI_Reduce(Y, Vt, n*rank, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif

      if (mpi_rank == 0) then
        Y(:,:) = Vt(:,:)
#endif

        ! SVD of Y
        call svd(Y,size(Y,1),UY,size(UY,1),D,Vt,size(Vt,1),rank,n)

        ! U(n,rank) = Z(n,rank).UY(rank,rank)
        call dgemm('N','N',n,rank,rank,1.d0,Z,size(Z,1),UY,size(UY,1),0.d0,U,size(U,1))

        ! ---

        do kk=1,rank/2
          r1 = ddot(n, U(1,kk), 1, Vt(kk,1), size(Vt,1))
          if (r1 < 0.d0) then
            D(kk) = -D(kk)
          endif
          if (dabs(r1) < 0.9d0) then
            D(kk) = 0.d0
          endif
        end do

        print *, 'Smallest found singular value    : ', D(rank)

#ifdef HAVE_MPI
      endif

      call MPI_BCAST(D, rank, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
      call MPI_BCAST(U, n*rank, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
      call MPI_BCAST(Vt, n*rank, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         print *, 'MPI error:', __FILE__, ':', __LINE__
         stop -1
      endif
#endif

      ! Remove positive eigenvalues from W_work
      do kk=1,rank/2
        do jj=w0,w1
          do ii=1,n
            W_work(ii,jj) = W_work(ii,jj) - U(ii,kk) * dabs(D(kk)) * Vt(kk,jj)
          enddo
        enddo
      enddo

      ! Remove negative eigenvalues from W
      do kk=1,rank/2
        if (D(kk) >= 0.d0) cycle
        do jj=w0,w1
          do ii=1,n
            W(ii,jj) = W(ii,jj) + U(ii,kk) * D(kk) * Vt(kk,jj)
          enddo
        enddo
        nremoved = nremoved+1
      enddo

      if (mpi_rank == 0) then
        print *, 'Removed ', nremoved, ' components'
        print *, D(1:3)
      endif

    end do

    deallocate(D, Vt, U, P, Z, UY, Y)


    moy(:) = 0.d0
    do kcp=1,nint
      i = is(kcp) ; j = js(kcp) ; k = ks(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) moy(kcp) = moy(kcp) + W(ii,jj)*0.125d0

      i = ks(kcp) ; j = js(kcp) ; k = is(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) moy(kcp) = moy(kcp) + W(ii,jj)*0.125d0

      i = ks(kcp) ; j = ls(kcp) ; k = is(kcp) ; l = js(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) moy(kcp) = moy(kcp) + W(ii,jj)*0.125d0

      i = is(kcp) ; j = ls(kcp) ; k = ks(kcp) ; l = js(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) moy(kcp) = moy(kcp) + W(ii,jj)*0.125d0

      i = js(kcp) ; j = is(kcp) ; k = ls(kcp) ; l = ks(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) moy(kcp) = moy(kcp) + W(ii,jj)*0.125d0

      i = js(kcp) ; j = ks(kcp) ; k = ls(kcp) ; l = is(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) moy(kcp) = moy(kcp) + W(ii,jj)*0.125d0

      i = ls(kcp) ; j = ks(kcp) ; k = js(kcp) ; l = is(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) moy(kcp) = moy(kcp) + W(ii,jj)*0.125d0

      i = ls(kcp) ; j = is(kcp) ; k = js(kcp) ; l = ks(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if ((jj >= w0).and.(jj <= w1)) moy(kcp) = moy(kcp) + W(ii,jj)*0.125d0
    enddo

#ifdef HAVE_MPI
    call MPI_AllReduce(MPI_IN_PLACE, moy, int(nint,4), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
    if (ierr /= MPI_SUCCESS) then
       print *, 'MPI error:', __FILE__, ':', __LINE__
       stop -1
    endif
#endif


    if (mpi_size == 1) then
#ifdef HAVE_MMAP
      call munmap( (/ int(n,8), int(n,8) /), 8, fd_W, ptr_W )
      call munmap( (/ int(n,8), int(n,8) /), 8, fd_W_work, ptr_W_work )
      call system('rm -f -- W.bin W_work.bin')
#else
      deallocate(W, W_work)
#endif
    endif

end subroutine svd_clean
