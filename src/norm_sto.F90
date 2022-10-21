double precision function norm_sto(gam,n)
  implicit none
  double precision gam,g,fact
  integer n,p
  p=2*n
  g=2.d0*gam
  norm_sto=dsqrt( g**(p+1)/fact(p))
end

