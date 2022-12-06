subroutine davidson_clean(moy, nint, is, js, ks, ls, nbasis)
  implicit none
  ! Removes negative eigenvalues of the ERI matrix

  integer, intent(in) :: nint, nbasis
  double precision, intent(inout) :: moy(nint)
  integer, dimension(nint), intent(in) :: is, js, ks, ls

  double precision, allocatable :: H_diag(:)
  double precision, allocatable :: V(:,:), E(:)

  integer :: NMAX
  integer :: i, j, k, l, kcp, ii, jj, iter, shift, n, m
  double precision :: r1, r2
  double precision, allocatable :: W(:,:), U(:,:)
  double precision, allocatable :: small_h(:,:), small_s(:,:), small_v(:,:)

  double precision, external :: dnrm2
  logical :: converged

  integer :: lwork, info
  double precision, allocatable :: work(:), residual_norm(:), E_tmp(:)

  n = nbasis*nbasis

  if (n < 10) then

    ! Lapack
    ! ======

    allocate (W(n,n),U(n,n), E_tmp(n))
    W = 0.d0

    do kcp=1,nint
      i = is(kcp) ; j = js(kcp) ; k = ks(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      W(ii,jj) = moy(kcp)

      i = ks(kcp) ; j = js(kcp) ; k = is(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      W(ii,jj) = moy(kcp)

      i = ks(kcp) ; j = ls(kcp) ; k = is(kcp) ; l = js(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      W(ii,jj) = moy(kcp)

      i = is(kcp) ; j = ls(kcp) ; k = ks(kcp) ; l = js(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      W(ii,jj) = moy(kcp)

      i = js(kcp) ; j = is(kcp) ; k = ls(kcp) ; l = ks(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      W(ii,jj) = moy(kcp)

      i = js(kcp) ; j = ks(kcp) ; k = ls(kcp) ; l = is(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      W(ii,jj) = moy(kcp)

      i = ls(kcp) ; j = ks(kcp) ; k = js(kcp) ; l = is(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      W(ii,jj) = moy(kcp)

      i = ls(kcp) ; j = is(kcp) ; k = js(kcp) ; l = ks(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      W(ii,jj) = moy(kcp)
    enddo

    U = W

    lwork = -1
    allocate(work(1))
    call dsyev('V', 'U', n, &
         U, size(U,1), &
         E_tmp, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)

    allocate(work(lwork))
    call dsyev('V', 'U', n, &
         U, size(U,1), &
         E_tmp, work, lwork, info)
    deallocate(work)

    if (info /= 0) then
       print *, 'info=', info
       stop 'DSYGV Diagonalization failed'
    endif

    do l=1,n
      if (E_tmp(l) >= 0.d0) exit
    end do
    l=l-1

    do k=1,l
     do j=1,n
       do i=1,n
          W(i,j) = W(i,j) - E_tmp(k) * U(i,k) * U(j,k)
        end do
      end do
    end do


    do kcp=1,nint
      i = is(kcp) ; j = js(kcp) ; k = ks(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      moy(kcp) = W(ii,jj)
    enddo

    deallocate(W,U)



  else



    ! Davidson
    ! ========

    m=4
    allocate (H_diag(n), E(m), V(n, m))

    ! Diagonal elements
    ! -----------------

    H_diag = 0.d0
    do kcp=1,nint
      i = is(kcp) ; j = js(kcp) ; k = ks(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if (ii == jj) H_diag(ii) = moy(kcp)

      i = ks(kcp) ; j = js(kcp) ; k = is(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if (ii == jj) H_diag(ii) = moy(kcp)

      i = ks(kcp) ; j = ls(kcp) ; k = is(kcp) ; l = js(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if (ii == jj) H_diag(ii) = moy(kcp)

      i = is(kcp) ; j = ls(kcp) ; k = ks(kcp) ; l = js(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if (ii == jj) H_diag(ii) = moy(kcp)

      i = js(kcp) ; j = is(kcp) ; k = ls(kcp) ; l = ks(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if (ii == jj) H_diag(ii) = moy(kcp)

      i = js(kcp) ; j = ks(kcp) ; k = ls(kcp) ; l = is(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if (ii == jj) H_diag(ii) = moy(kcp)

      i = ls(kcp) ; j = ks(kcp) ; k = js(kcp) ; l = is(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if (ii == jj) H_diag(ii) = moy(kcp)

      i = ls(kcp) ; j = is(kcp) ; k = js(kcp) ; l = ks(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      if (ii == jj) H_diag(ii) = moy(kcp)
    enddo

    NMAX=min(10, (n/m))

    allocate(U(n,NMAX*m), W(n,NMAX*m), residual_norm(m) )
    allocate(small_h(NMAX*m, NMAX*m), small_s(NMAX*m, NMAX*m), small_v(NMAX*m,NMAX*m), E_tmp(NMAX*m) )

    ! Initialize with normalized gaussian random vectors
    do j=1,m
       do i=1,n
          call random_number(r1)
          call random_number(r2)
          r1 = dsqrt(-2.d0*dlog(r1))
          r2 = 2.d0*dacos(-1.d0) * r2
          U(i,j) = r1*dcos(r2)
       enddo
       r1 = dnrm2(n,U(1:n,j),1)
       call dscal(n,1.d0/r1,U(1:n,j),1)
    enddo

    converged = .False.
    do while (.not. converged)
      do iter=1,NMAX-1
         shift = m*(iter-1)

         ! Compute |W_k> = \sum_i |i><i|H|u_k>
         ! ------------------------------------

!         call dgemm('N', 'N', n, m, n, 1.d0, &
!              H, size(H,1), &
!              U(1,shift+1), size(U,1), 0.d0, &
!              W(1,shift+1), size(W,1) )


         ! Compute h_kl = <u_k | u_l>
         ! --------------------------

         call dgemm('T', 'N', shift+m, shift+m, n, 1.d0, &
              U(1,1), size(U,1), &
              U(1,1), size(U,1), 0.d0, &
              small_s, size(small_s,1) )

         ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
         ! ------------------------------------------

         call dgemm('T', 'N', shift+m, shift+m, n, 1.d0, &
              U(1,1), size(U,1), &
              W(1,1), size(W,1), 0.d0, &
              small_h, size(small_h,1) )


         ! Diagonalize small H matrix
         ! --------------------------

         small_v = small_h

         lwork = -1
         allocate(work(1))
         call dsygv(1, 'V', 'U', shift+m, &
              small_v, size(small_v,1), &
              small_s, size(small_s,1), &
              E_tmp, work, lwork, info)
         lwork = int(work(1))
         deallocate(work)

         allocate(work(lwork))
         call dsygv(1, 'V', 'U', shift+m, &
              small_v, size(small_v,1), &
              small_s, size(small_s,1), &
              E_tmp, work, lwork, info)
         deallocate(work)

         if (info /= 0) then
            print *, 'info=', info
            stop 'DSYGV Diagonalization failed'
         endif


         ! Compute eigenvalue for each eigenvector
         ! ---------------------------------------

         call dgemm('N', 'N', shift+m, shift+m, shift+m, 1.d0, &
              small_h, size(small_h,1), &
              small_v, size(small_v,1), 0.d0, &
              small_s, size(small_s,1) )

         call dgemm('T', 'N', shift+m, shift+m, shift+m, 1.d0, &
              small_v, size(small_v,1), &
              small_s, size(small_s,1), 0.d0, &
              small_h, size(small_h,1) )

         do i=1,m
            E(i) = small_h(i,i)
         enddo


         ! Express eigenvectors in original basis
         ! --------------------------------------

         call dgemm('N', 'N', n, m, shift+m, 1.d0, &
              U, size(U,1), &
              small_v, size(small_v,1), 0.d0, &
              U(1,shift+m+1), size(U,1) )

         call dgemm('N', 'N', n, m, shift+m, 1.d0, &
              W, size(W,1), &
              small_v, size(small_v,1), 0.d0, &
              W(1,shift+m+1), size(W,1) )


         ! Compute residual vector and davidson step
         ! -----------------------------------------

         do j=1,m
            do i=1,n
               U(i,shift+m+j) = (E(j) * U(i,shift+m+j) - W(i,shift+m+j)) &
                    /max(H_diag(i) - E(j), 1.d-2)
            enddo
            residual_norm(j) = dnrm2(n, U(1:n,shift+m+j), 1)
            call dscal(n,1.d0/residual_norm(j),U(1:n,shift+m+j),1)
            print *, iter, j, E(j), residual_norm(j)
         enddo
         print *, '---'

         converged = (maxval(residual_norm(1:m)) < 1.d-08)

         if (converged) exit
      end do

      ! Re-contract U and update S and W
      ! --------------------------------

      call dgemm('N', 'N', n, m, shift+m, 1.d0,  &
          U, size(U,1), small_v, size(small_v,1), 0.d0, &
          V, size(V,1))

      U(1:n,1:m) = V(1:n,1:m)

    end do


  endif

end subroutine davidson_clean
