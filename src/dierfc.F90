! inverse of error function in double precision
!
double precision function dierfc(y)
  implicit none
  double precision :: qa, qb, qc, qd
  double precision :: q0, q1, q2, q3, q4
  double precision :: pa, pb
  double precision :: p0, p1, p2, p3, p4, p5, p6, p7, p8, p9
  double precision :: p10, p11, p12, p13, p14, p15, p16, p17, p18, p19
  double precision :: p20, p21, p22
  double precision :: x, y, z, s, t, u, w
  parameter (                                                        &
      qa = 9.16461398268964d-01,                                     &
      qb = 2.31729200323405d-01,                                     &
      qc = 4.88826640273108d-01,                                     &
      qd = 1.24610454613712d-01,                                     &
      q0 = 4.99999303439796d-01,                                     &
      q1 = 1.16065025341614d-01,                                     &
      q2 = 1.50689047360223d-01,                                     &
      q3 = 2.69999308670029d-01,                                     &
      q4 = -7.28846765585675d-02)
  parameter (                                                        &
      pa = 3.97886080735226000d+00,                                  &
      pb = 1.20782237635245222d-01,                                  &
      p0 = 2.44044510593190935d-01,                                  &
      p1 = 4.34397492331430115d-01,                                  &
      p2 = 6.86265948274097816d-01,                                  &
      p3 = 9.56464974744799006d-01,                                  &
      p4 = 1.16374581931560831d+00,                                  &
      p5 = 1.21448730779995237d+00,                                  &
      p6 = 1.05375024970847138d+00,                                  &
      p7 = 7.13657635868730364d-01,                                  &
      p8 = 3.16847638520135944d-01,                                  &
      p9 = 1.47297938331485121d-02,                                  &
      p10 = -1.05872177941595488d-01,                                &
      p11 = -7.43424357241784861d-02)
  parameter (                                                        &
      p12 = 2.20995927012179067d-03,                                 &
      p13 = 3.46494207789099922d-02,                                 &
      p14 = 1.42961988697898018d-02,                                 &
      p15 = -1.18598117047771104d-02,                                &
      p16 = -1.12749169332504870d-02,                                &
      p17 = 3.39721910367775861d-03,                                 &
      p18 = 6.85649426074558612d-03,                                 &
      p19 = -7.71708358954120939d-04,                                &
      p20 = -3.51287146129100025d-03,                                &
      p21 = 1.05739299623423047d-04,                                 &
      p22 = 1.12648096188977922d-03)
  z = y
  if (y .gt. 1.d0) z = 2.d0 - y
  w = qa - dlog(z)
  u = dsqrt(w)
  s = (qc + dlog(u)) / w
  t = 1.d0 / (u + qb)
  x = u * (1 - s * (0.5d0 + s * qd)) -                               &
      ((((q4 * t + q3) * t + q2) * t + q1) * t + q0) * t
  t = pa / (pa + x)
  u = t - 0.5d0
  s = (((((((((p22 * u + p21) * u + p20) * u +                       &
      p19) * u + p18) * u + p17) * u + p16) * u +                    &
      p15) * u + p14) * u + p13) * u + p12
  s = ((((((((((((s * u + p11) * u + p10) * u +                      &
      p9) * u + p8) * u + p7) * u + p6) * u + p5) * u +              &
      p4) * u + p3) * u + p2) * u + p1) * u + p0) * t -              &
      z * dexp(x * x - pb)
  x = x + s * (1.d0 + x * s)
  if (y .gt. 1.d0) x = -x
  dierfc = x
end
