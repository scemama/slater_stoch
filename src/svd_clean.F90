subroutine svd_clean(moy, nint, is, js, ks, ls, nbasis, rank, q)
  implicit none
  ! Removes negative eigenvalues of the ERI matrix

  integer*8, intent(in)          :: nint
  integer, intent(in)            :: nbasis, rank, q
  double precision, intent(inout) :: moy(nint)
  integer, dimension(nint), intent(in) :: is, js, ks, ls

  double precision, allocatable  :: H_diag(:)
  double precision, allocatable  :: V(:,:), E(:)

  integer                        :: NMAX
  integer                        :: i, j, k, l, kk, iter, shift, n, m, mmax
  integer*8                      :: ii, jj, kcp
  double precision               :: r1, r2, hij
  double precision, allocatable  :: W(:,:), U(:,:), D(:), Vt(:,:)
  double precision, allocatable  :: Wt(:,:), Ut(:,:)
  double precision, allocatable  :: Pt(:,:), Zt(:,:), UY(:,:)
  double precision, allocatable  :: P(:,:), Z(:,:)
  double precision, parameter    :: dtwo_pi = 2.d0*dacos(-1.d0)


  double precision, external     :: dnrm2, ddot
  logical                        :: converged

  integer                        :: lwork, info
  double precision, allocatable  :: work(:), residual_norm(:), E_tmp(:)
  integer                        :: xi(8), xj(8), xk(8), xl(8)

  n = nbasis*nbasis

  if (nbasis < 4) then

    ! Lapack
    ! ======

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

    allocate(U(n,rank))
    allocate(D(rank))
    allocate(Vt(rank,n))
    D=1.d0
    call randomized_svd(W, size(W,1), U, size(U,1), D, Vt, size(Vt,1),&
        n, n, q, rank)
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

    ! P is a normal random matrix (n,r)
    allocate(Pt(rank,n), Zt(rank,n))

    call random_number(Pt)
    call random_number(Zt)
    do i=1,n
      do k=1,rank
        r1 = Pt(k,i)
        r2 = Zt(k,i)
        r1 = dsqrt(-2.d0*dlog(r1))
        r2 = dtwo_pi*r2
        Pt(k,i) = r1*dcos(r2)
      enddo
    enddo

    ! Z(m,r) = A(m,n).P(n,r)
    Zt = 0.d0
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
          if ( (xi(i) == xi(j)).and.                                 &
                (xj(i) == xj(j)).and.                                &
                (xk(i) == xk(j)).and.                                &
                (xl(i) == xl(j)) ) then
            xi(i) = 0
            exit
          endif
        enddo
      enddo

      do i=1,8
        if (xi(i) == 0) cycle
        ii = xi(i) + (xj(i)-1)*nbasis ; jj = xk(i) + (xl(i)-1)*nbasis
        Zt(1:rank,jj) = Zt(1:rank,jj) +  Pt(1:rank,ii) * moy(kcp)
      enddo

    enddo

    allocate(Z(n,rank))
    do k=1,rank
      do i=1,n
        Z(i,k) = Zt(k,i)
      enddo
    enddo
    call ortho_qr(Z,size(Z,1),n,rank)
    do i=1,n
      do k=1,rank
        Zt(k,i) = Z(i,k)
      enddo
    enddo
    deallocate(Z)

    ! Power iterations
    do iter=1,q
      ! P(n,r) = At(n,m).Z(m,r)
      Pt = 0.d0
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
            if ( (xi(i) == xi(j)).and.                               &
                  (xj(i) == xj(j)).and.                              &
                  (xk(i) == xk(j)).and.                              &
                  (xl(i) == xl(j)) ) then
              xi(i) = 0
              exit
            endif
          enddo
        enddo

        do i=1,8
          if (xi(i) == 0) cycle
          ii = xi(i) + (xj(i)-1)*nbasis ; jj = xk(i) + (xl(i)-1)*nbasis
          Pt(1:rank,jj) = Pt(1:rank,jj) +  Zt(1:rank,ii) * moy(kcp)
        enddo

      enddo

      Zt = 0.d0
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
            if ( (xi(i) == xi(j)).and.                               &
                  (xj(i) == xj(j)).and.                              &
                  (xk(i) == xk(j)).and.                              &
                  (xl(i) == xl(j)) ) then
              xi(i) = 0
              exit
            endif
          enddo
        enddo

        do i=1,8
          if (xi(i) == 0) cycle
          ii = xi(i) + (xj(i)-1)*nbasis ; jj = xk(i) + (xl(i)-1)*nbasis
          Zt(1:rank,jj) = Zt(1:rank,jj) +  Pt(1:rank,ii) * moy(kcp)
        enddo

      enddo
      allocate(Z(n,rank))
      do k=1,rank
        do i=1,n
          Z(i,k) = Zt(k,i)
        enddo
      enddo
      call ortho_qr(Z,size(Z,1),n,rank)
      do i=1,n
        do k=1,rank
          Zt(k,i) = Z(i,k)
        enddo
      enddo
      deallocate(Z)

    enddo

    ! Pt(r,n) = Zt(r,m).A(m,n)
    Pt = 0.d0
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
          if ( (xi(i) == xi(j)).and.                                 &
                (xj(i) == xj(j)).and.                                &
                (xk(i) == xk(j)).and.                                &
                (xl(i) == xl(j)) ) then
            xi(i) = 0
            exit
          endif
        enddo
      enddo

      do i=1,8
        if (xi(i) == 0) cycle
        ii = xi(i) + (xj(i)-1)*nbasis ; jj = xk(i) + (xl(i)-1)*nbasis
        Pt(1:rank,jj) = Pt(1:rank,jj) +  Zt(1:rank,ii) * moy(kcp)
      enddo

    enddo

    ! SVD of Pt
    allocate(UY(rank,rank), D(rank), Vt(rank,n))
    call svd(Pt,size(Pt,1),UY,size(UY,1),D,Vt,size(Vt,1),rank,n)
    deallocate(Pt)

    ! U(m,r) = Z(m,r).UY(r,r)
    allocate(U(n,rank))
    call dgemm('T','N',n,rank,rank,1.d0,Zt,size(Zt,1),UY,size(UY,1),0.d0,U,size(U,1))
    deallocate(UY,Zt)

    moy(:) = 0.d0
    do kk=1,rank
      if (ddot(n, U(1,kk), 1, Vt(kk,1), size(Vt,1)) < 0.d0) cycle
      do kcp=1, nint
        i = is(kcp) ; j = js(kcp) ; k = ks(kcp) ; l = ls(kcp)
        ii = i + (j-1)*nbasis ; jj = k + (l-1)*nbasis
        moy(kcp) = moy(kcp) + D(kk) * U(ii,kk) * Vt(kk,jj)
      enddo
    enddo

    deallocate(D, Vt, U)

  endif

end subroutine svd_clean
