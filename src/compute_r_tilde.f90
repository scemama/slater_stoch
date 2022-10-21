subroutine compute_r_tilde(i,k,r1,r2,rt1,rt2)
  include 'j.inc'
  integer, intent(in)            :: i, k
  double precision, intent(in)   :: r1(nw,3),r2(nw,3)
  double precision, intent(out)  :: rt1(nw,3),rt2(nw,3)
  factor = 0.5d0 * a_ZV(i,k)/dsqrt(g_min(i)+g_min(k))
  do kw=1,nw
    r1_mod=dsqrt(r1(kw,1)*r1(kw,1)+r1(kw,2)*r1(kw,2)+r1(kw,3)*r1(kw,3))
    r2_mod=dsqrt(r2(kw,1)*r2(kw,1)+r2(kw,2)*r2(kw,2)+r2(kw,3)*r2(kw,3))
    do ll=1,3
      rt1(kw,ll)=r1(kw,ll) * r1_mod * factor
      rt2(kw,ll)=r2(kw,ll) * r2_mod * factor
    enddo !!ll
  enddo !!kw
end

