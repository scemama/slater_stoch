subroutine svd_clean(moy, nint, is, js, ks, ls, nbasis)
  implicit none
  ! Removes negative eigenvalues of the ERI matrix

  integer*8, intent(in) :: nint
  integer, intent(in) :: nbasis
  double precision, intent(inout) :: moy(nint)
  integer, dimension(nint), intent(in) :: is, js, ks, ls

  double precision, allocatable :: H_diag(:)
  double precision, allocatable :: V(:,:), E(:)

  integer :: NMAX
  integer :: i, j, k, l, kk, iter, shift, n, m, mmax
  integer*8 :: ii, jj, kcp
  double precision :: r1, r2, hij
  double precision, allocatable :: W(:,:), U(:,:), D(:), Vt(:,:)
  double precision, allocatable :: Wt(:,:), Ut(:,:)
  double precision, allocatable :: small_h(:,:), small_s(:,:), small_v(:,:)

  double precision, external :: dnrm2, ddot
  logical :: converged

  integer :: lwork, info, rank
  double precision, allocatable :: work(:), residual_norm(:), E_tmp(:)
  integer :: xi(8), xj(8), xk(8), xl(8)

  n = nbasis*nbasis

  if (nbasis < 400) then

    ! Lapack
    ! ======
    rank = min(n, 1000)

    allocate (W(n,n))
    W(:,:) = 0.d0

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

    rank = n/2
    allocate(U(n,rank))
    allocate(D(n))
    allocate(Vt(n,n))
    do while (rank<=n)
      D=1.d0
      deallocate(U)
      allocate(U(n,rank))
      call randomized_svd(W, size(W,1), U, size(U,1), D, Vt, size(Vt,1), &
        n, n, 10, rank)
      if (D(rank) > 1.d-6) then
        print *, 'rank=', rank, D(rank)
        rank = min(n,rank+rank/2)
      else
        exit
      endif
    end do
    print *, 'rank=', rank
!    call svd(W,size(W,1),U,size(U,1),D,Vt,size(Vt,1),n,n)

    W(:,:) = 0.d0
    do k=1,rank
     if (ddot(n, U(1,k), 1, Vt(k,1), size(Vt,1)) < 0.d0) cycle
     do j=1,n
       do i=1,n
          W(i,j) = W(i,j) + D(k) * U(i,k) * Vt(k,j)
        end do
      end do
    end do

    deallocate(D, Vt, U)

    do kcp=1,nint
      i = is(kcp) ; j = js(kcp) ; k = ks(kcp) ; l = ls(kcp)
      ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
      moy(kcp) = W(ii,jj)
    enddo

    deallocate(W)
    return


  else

    return


    ! Davidson
    ! ========


    mmax=nbasis
    m=mmax
    NMAX=nbasis
    allocate (E(m), V(n, m))
    allocate(U(n,NMAX*m), W(n,NMAX*m), residual_norm(m) )
    allocate(Ut(m,n), Wt(m,n))
    allocate(small_h(NMAX*m, NMAX*m), small_s(NMAX*m, NMAX*m), small_v(NMAX*m,NMAX*m), E_tmp(NMAX*m) )

    do while (m>0)
      m=mmax

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

           do jj=1,n
             do ii=1,m
               Ut(ii,jj) = U(jj,shift+ii)
             enddo
           enddo


           Wt = 0.d0
           do kcp=1,nint

             i = is(kcp) ; j = js(kcp) ; k = ks(kcp) ; l = ls(kcp)
             xi(1) = i ; xj(1) = j ; xk(1) = k ; xl(1) = l
             xi(2) = i ; xj(2) = l ; xk(2) = k ; xl(2) = j
             xi(3) = k ; xj(3) = j ; xk(3) = i ; xl(3) = l
             xi(4) = k ; xj(4) = l ; xk(4) = i ; xl(4) = j
             xi(5) = j ; xj(5) = i ; xk(5) = l ; xl(5) = k
             xi(6) = j ; xj(6) = k ; xk(6) = l ; xl(6) = i
             xi(7) = l ; xj(7) = i ; xk(7) = j ; xl(7) = k
             xi(8) = l ; xj(8) = k ; xk(8) = j ; xl(8) = i

             do i=2,8
               do j=1,i-1
                 if ( (xi(i) == xi(j)).and. &
                      (xj(i) == xj(j)).and. &
                      (xk(i) == xk(j)).and. &
                      (xl(i) == xl(j)) ) then
                     xi(i) = 0
                     exit
                 endif
               enddo
             enddo

             do i=1,8
               if (xi(i) == 0) cycle
               ii = xi(i) + (xj(i)-1)*nbasis ; jj = xk(i) + (xl(i)-1)*nbasis
               Wt(1:m,jj) = Wt(1:m,jj) +  Ut(1:m,ii) * moy(kcp)
             enddo

           enddo

           do ii=1,m
             do jj=1,n
               W(jj,shift+ii) = Wt(ii,jj)
             enddo
           enddo


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
                 U(i,shift+m+j) = E(j) * U(i,shift+m+j) - W(i,shift+m+j)
              enddo
              residual_norm(j) = dnrm2(n, U(1:n,shift+m+j), 1)
              call dscal(n,1.d0/residual_norm(j),U(1:n,shift+m+j),1)
           enddo
           print *, iter, real(E(1)), real(residual_norm(1)), real(E(m)), real(residual_norm(m))

           converged = (residual_norm(m/2+1) < 1.d-08).or.(E(m)<0.d0)

           if (converged) exit
        end do

        ! Re-contract U and update S and W
        ! --------------------------------

        call dgemm('N', 'N', n, m, shift+m, 1.d0,  &
            U, size(U,1), small_v, size(small_v,1), 0.d0, &
            V, size(V,1))

        U(1:n,1:m) = V(1:n,1:m)

      end do

      do kk=1,m
        if (E(kk) >= 0.d0) exit
      enddo
      m = kk-1
      do kcp=1,nint
        i = is(kcp) ; j = js(kcp) ; k = ks(kcp) ; l = ls(kcp)
        ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
        do kk=1,m
          moy(kcp) = moy(kcp) - E(kk) * U(ii,kk) * U(jj,kk)
        enddo
      enddo

    enddo

  endif

end subroutine svd_clean
