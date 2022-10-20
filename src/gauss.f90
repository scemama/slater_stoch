double precision function gauss(sigma)
  implicit double precision(a-h,o-z)
  dpi=2.d0*dacos(-1.d0)
  call random_number(r1)
  call random_number(r2)
  gauss=sigma*dsqrt(-2.d0*dlog(r1)*dcos(dpi*r2))
end
