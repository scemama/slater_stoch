subroutine draw_configuration(r1,r2)
  include 'j.inc'
  dimension r1(nw,3),r2(nw,3)
  double precision, parameter    :: f = dsqrt(2.d0)
  do kw=1,nw
    do ll=1,3
      call random_number(r1(kw,ll))
      r1(kw,ll)=f*dierfc(2.d0-2.d0*r1(kw,ll))
      call random_number(r2(kw,ll))
      r2(kw,ll)=f*dierfc(2.d0-2.d0*r2(kw,ll))
    enddo
  enddo
end

