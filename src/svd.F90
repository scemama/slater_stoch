subroutine svd(A, LDA, U, LDU, D, Vt, LDVt, m, n)
  implicit none
  ! Compute A = U.D.Vt
  ! LDx : leftmost dimension of x
  ! Dimension of A is m x n

  integer, intent(in)             :: LDA, LDU, LDVt, m, n
  double precision, intent(in)    :: A(LDA,n)
  double precision, intent(out)   :: U(LDU,min(m,n))
  double precision,intent(out)    :: Vt(LDVt,n)
  double precision,intent(out)    :: D(min(m,n))
  double precision,allocatable    :: work(:)
  integer                         :: info, lwork, i, j, k

  double precision,allocatable    :: A_tmp(:,:)
  allocate (A_tmp(LDA,n))
  do k=1,n
    do i=1,m
      A_tmp(i,k) = A(i,k)
    enddo
  enddo

  ! Find optimal size for temp arrays
  allocate(work(1))
  lwork = -1
  call dgesvd('S','S', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  lwork = max(int(work(1)), 10*MIN(M,N))
  deallocate(work)

  allocate(work(lwork))
  call dgesvd('S','S', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  deallocate(A_tmp,work)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

  do j=1,min(m,n)
    do i=1,m
      if (dabs(U(i,j)) < 1.d-14)  U(i,j) = 0.d0
    enddo
  enddo

  do j=1,n
    do i=1,n
      if (dabs(Vt(i,j)) < 1.d-14) Vt(i,j) = 0.d0
    enddo
  enddo

end

