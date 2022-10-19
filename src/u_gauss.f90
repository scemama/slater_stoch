!! u_gauss is the radial part of the end-of-the-day gaussian representation
double precision function u_gauss(i,r)
include 'j.inc'
u_gauss=0.d0
do mm=1,n_gauss(i)
 u_gauss=u_gauss+(c_gauss(1,mm,i)+c_gauss(2,mm,i)*r**2)*dexp(-g_gauss(mm,i)*r**2)
enddo
end
