subroutine ortho_qr(A,LDA,m,n)
  implicit none
  ! Orthogonalization using Q.R factorization
  !
  ! A : matrix to orthogonalize
  !
  ! LDA : leftmost dimension of A
  !
  ! m : Number of rows of A
  !
  ! n : Number of columns of A
  !
  integer, intent(in)            :: m,n, LDA
  double precision, intent(inout) :: A(LDA,n)

  integer                        :: LWORK, INFO
  double precision, allocatable  :: TAU(:), WORK(:)

  allocate (TAU(min(m,n)), WORK(1))

  LWORK=-1
  call dgeqrf( m, n, A, LDA, TAU, WORK, LWORK, INFO )
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  LWORK=max(n,int(WORK(1)))

  deallocate(WORK)
  allocate(WORK(LWORK))
  call dgeqrf(m, n, A, LDA, TAU, WORK, LWORK, INFO )

  LWORK=-1
  call dorgqr(m, n, n, A, LDA, TAU, WORK, LWORK, INFO)
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  LWORK=max(n,int(WORK(1)))
  deallocate(WORK)
  allocate(WORK(LWORK))
  call dorgqr(m, n, n, A, LDA, TAU, WORK, LWORK, INFO)

  deallocate(WORK,TAU)
end

