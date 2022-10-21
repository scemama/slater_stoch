subroutine compute_jacobian_and_pi0(r1,r2,rjacob1,rjacob2,pi_0,n)
  include 'j.inc'
  double precision               :: r1(nw,3), r2(nw,3),pi_0(nw), rjacob1(nw,n,n), rjacob2(nw,n,n)
  double precision               :: d1(nw), d2(nw)
  double precision, parameter    :: factor = 1.d0 / (2.d0*dacos(-1.d0))**3
  do kw=1,nw
    r1_mod_2=r1(kw,1)*r1(kw,1)+r1(kw,2)*r1(kw,2)+r1(kw,3)*r1(kw,3)
    r2_mod_2=r2(kw,1)*r2(kw,1)+r2(kw,2)*r2(kw,2)+r2(kw,3)*r2(kw,3)
    d1(kw) = dsqrt(r1_mod_2)
    d2(kw) = dsqrt(r2_mod_2)
    pi_0(kw)= dexp(-0.5d0*(r1_mod_2 + r2_mod_2)) * factor
  enddo
  do k=1,n
    do i=k,n
      factor_ik = a_ZV(i,k)/(g_min(i)+g_min(k))
      do kw=1,nw
        f1=factor_ik*d1(kw)
        rjacob1(kw,i,k)=0.25d0*dabs(f1*f1*f1)
        rjacob1(kw,k,i)= rjacob1(kw,i,k)

        f2=factor_ik*d2(kw)
        rjacob2(kw,i,k)=0.25d0*dabs(f2*f2*f2)
        rjacob2(kw,k,i)= rjacob2(kw,i,k)
      enddo
    enddo
  enddo
end
