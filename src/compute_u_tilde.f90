subroutine compute_u_tilde(i,k,rt1,rt2,ut1,ut2)
  include 'j.inc'
  dimension rt1(nw,3),rt2(nw,3)
  dimension ut1(3,nw,nbasis_max*nbasis_max),ut2(3,nw,nbasis_max*nbasis_max)
  double precision factor
  ik = (k-1)*nbasis_max+i
  factor = 1.d0/dsqrt(g_min(i)+g_min(k))
  do kw=1,nw
    do ll=1,3
      ut1(ll,kw,ik)= factor*rt1(kw,ll) + G_center(ll,i,k)
      ut2(ll,kw,ik)= factor*rt2(kw,ll) + G_center(ll,i,k)
    enddo !!ll
  enddo !!kw
end
