double precision function norm_sto(gam,n)
  implicit none
  double precision :: gam,g
  double precision, external :: fact
  integer :: n, p
  p=n+n
  g=2.d0*gam
  norm_sto=dsqrt( g**(2*n+1)/fact(p))
end

