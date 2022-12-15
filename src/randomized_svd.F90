subroutine randomized_svd(A,LDA,U,LDU,D,Vt,LDVt,m,n,q,r)
  implicit none
! Randomized SVD: rank r, q power iterations
! 1. Sample column space of A with P: Z = A.P where P is random with r+p columns.
! 2. Power iterations : Z <- X . (Xt.Z)
! 3. Z = Q.R
! 4. Compute SVD on projected Qt.X = U' . S. Vt
! 5. U = Q U'

  integer, intent(in)             :: LDA, LDU, LDVt, m, n, q, r
  double precision, intent(in)    :: A(LDA,n)
  double precision, intent(out)   :: U(LDU,r)
  double precision,intent(out)    :: Vt(LDVt,r)
  double precision,intent(out)    :: D(r)
  integer                         :: i, j, k

  double precision, parameter     :: dtwo_pi = 2.d0*dacos(-1.d0)
  double precision,allocatable    :: Z(:,:), P(:,:), Y(:,:), UY(:,:)
  double precision :: r1,r2
  allocate(P(n,r), Z(m,r))

  ! P is a normal random matrix (n,r)
  do k=1,r
    do i=1,n
      call random_number(r1)
      call random_number(r2)
      r1 = dsqrt(-2.d0*dlog(r1))
      r2 = dtwo_pi*r2
      P(i,k) = r1*dcos(r2)
    enddo
  enddo

  ! Z(m,r) = A(m,n).P(n,r)
  call dgemm('N','N',m,r,n,1.d0,A,size(A,1),P,size(P,1),0.d0,Z,size(Z,1))

  ! Power iterations
  do k=1,q
    ! P(n,r) = At(n,m).Z(m,r)
    call dgemm('T','N',n,r,m,1.d0,A,size(A,1),Z,size(Z,1),0.d0,P,size(P,1))
    ! Z(m,r) = A(m,n).P(n,r)
    call dgemm('N','N',m,r,n,1.d0,A,size(A,1),P,size(P,1),0.d0,Z,size(Z,1))
  enddo

  deallocate(P)

  ! QR factorization of Z
  call ortho_svd(Z,size(Z,1),m,r)

  allocate(Y(r,n), UY(r,r))
  ! Y(r,n) = Zt(r,m).A(m,n)
  call dgemm('T','N',r,n,m,1.d0,Z,size(Z,1),A,size(A,1),0.d0,Y,size(Y,1))

  ! SVD of Y
  call svd(Y,size(Y,1),UY,size(UY,1),D,Vt,size(Vt,1),r,n)
  deallocate(Y)

  ! U(m,r) = Z(m,r).UY(r,r)
  call dgemm('N','N',m,r,r,1.d0,Z,size(Z,1),UY,size(UY,1),0.d0,U,size(U,1))
  deallocate(UY,Z)
end


