subroutine draw_configuration(r1,r2,nw)
  double precision, intent(out) :: r1(nw,3),r2(nw,3)
  double precision, parameter   :: f = dsqrt(2.d0)
  double precision, external    :: dierfc
  integer :: kw, ll
  double precision :: x(nw,3,2)

  call random_number(x)
  x = 2.d0 - 2.d0*x

#ifdef HAVE_VDERFINVC

  call vderfcinv(3*nw, x(1,1,1), r1)
  call vderfcinv(3*nw, x(1,1,2), r2)
  r1 = r1*f
  r2 = r2*f

#else

  do ll=1,3
    do kw=1,nw
      r1(kw,ll)=f*dierfc(x(kw,ll,1))
      r2(kw,ll)=f*dierfc(x(kw,ll,2))
    enddo
  enddo

#endif

end

