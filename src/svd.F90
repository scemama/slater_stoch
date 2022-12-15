subroutine svd(A, LDA, U, LDU, D, Vt, LDVt, m, n)
  implicit none
  ! Compute A = U.D.Vt
  ! LDx : leftmost dimension of x
  ! Dimension of A is m x n

  integer, intent(in)             :: LDA, LDU, LDVt, m, n
  double precision, intent(in)    :: A(LDA,n)
  double precision, intent(out)   :: U(LDU,n)
  double precision, intent(out)   :: Vt(LDVt,n)
  double precision, intent(out)   :: D(n)
  double precision, allocatable   :: work(:)
  integer                         :: info, lwork, i, j, k

  double precision,allocatable    :: A_tmp(:,:)
  allocate (A_tmp(m,n))
  do k=1,n
    do i=1,m
      A_tmp(i,k) = A(i,k)
    enddo
  enddo

  ! Find optimal size for temp arrays
  allocate(work(1))
  lwork = -1
  call dgesvd('S','S', m, n, A_tmp, size(A_tmp,1),                             &
      D, U, size(U,1), Vt, size(Vt,1), work, lwork, info)
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  lwork = max(int(work(1)), 10*MIN(M,N))
  deallocate(work)

  allocate(work(lwork))
  call dgesvd('S','S', m, n, A_tmp, size(A_tmp,1),                             &
      D, U, size(U,1), Vt, size(Vt,1), work, lwork, info)
  deallocate(A_tmp,work)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

end

