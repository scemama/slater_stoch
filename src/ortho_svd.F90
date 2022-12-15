subroutine ortho_svd(A,LDA,m,n)
  implicit none
  ! Orthogonalization via SVD
  ! A : matrix to orthogonalize
  ! LDA : leftmost dimension of A
  ! m : Number of rows of A
  ! n : Number of columns of A
  integer, intent(in)            :: m,n, LDA
  double precision, intent(inout) :: A(LDA,n)
  double precision, allocatable :: U(:,:), D(:), Vt(:,:)
  integer :: i,j

  if (m < n) then
    call ortho_qr(A,LDA,m,n)
  endif
  allocate(U(m,m), D(max(m,n)), Vt(n,n))
  call SVD(A,LDA,U,size(U,1),D,Vt,size(Vt,1),m,n)
  do j=1,n
    do i=1,m
      A(i,j) = U(i,j)
    enddo
  enddo
  deallocate(U,D, Vt)

end

