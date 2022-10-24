subroutine draw_configuration(r1,r2,nw)
  double precision, intent(out) :: r1(nw,3),r2(nw,3)
  double precision, parameter   :: f = dsqrt(2.d0)
  double precision, external    :: dierfc
  integer :: kw, ll
  call random_number(r1)
  call random_number(r2)
  do ll=1,3
    do kw=1,nw
      r1(kw,ll)=f*dierfc(2.d0-2.d0*r1(kw,ll))
      r2(kw,ll)=f*dierfc(2.d0-2.d0*r2(kw,ll))
    enddo
  enddo
end

