double precision function chi2(i_orb)
  include 'j.inc'
  ntot=npower(1,i_orb)+npower(2,i_orb)+npower(3,i_orb)
  npts=1000
  dr=r_infty(i_orb)/npts
  chi2=0.d0
  
  do k=1,npts
    r=(k-1)*dr
    u1=u_orb  (i_orb,r)
    u2=u_gauss(i_orb,r)
    chi2=chi2+r*r**ntot*(u1-u2)**2
  enddo !r
  chi2=chi2/npts
  
end
