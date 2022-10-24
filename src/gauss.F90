double precision function gauss(sigma)
  implicit none
  double precision, intent(in) :: sigma
  double precision :: dpi, r1, r2
  dpi=2.d0*dacos(-1.d0)
  call random_number(r1)
  call random_number(r2)
  gauss=sigma*dsqrt(-2.d0*dlog(r1)*dcos(dpi*r2))
end
