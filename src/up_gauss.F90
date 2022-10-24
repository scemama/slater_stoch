!! derivative of u_gauss wrt r
double precision function up_gauss(i,r)
  use common_data
  up_gauss=0.d0
  do mm=1,n_gauss(i)
    gi=g_gauss(mm,i)
    c1=c_gauss(1,mm,i)
    c2=c_gauss(2,mm,i)
    up_gauss=up_gauss+2.d0*r*(c2-gi*c1-gi*c2*r**2)*dexp(-gi*r*r)
  enddo
end


