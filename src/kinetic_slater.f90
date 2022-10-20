double precision function kinetic_slater(nslat,lslat,mslat,gamma)
  implicit none
  integer nslat(4),lslat(4),mslat(4),i,nslatemp(4)
  integer l1,l2,m1,m2,n2,n1,n12
  double precision gamma(4),int,norm,norm_sto,g2,overlap_slater,g1,g12,fact
  double precision t1,t2,t3,t4
  
  kinetic_slater=0.d0
  
  l1=lslat(1)
  l2=lslat(2)
  if(l1.ne.l2)return
  m1=mslat(1)
  m2=mslat(2)
  if(m1.ne.m2)return
  g1=gamma(1)
  g2=gamma(2)
  g12=g1+g2
  n1=nslat(1)
  n2=nslat(2)
  n12=n1+n2
  
  norm=1.d0
  do i=1,2
    norm=norm*norm_sto(gamma(i),nslat(i))
  enddo
  
  t1=n2*(n2-1)*fact(n12-2)/g12**(n12-1)
  t2=2.d0*n2*g2*fact(n12-1)/g12**n12
  t3=g2**2*fact(n12)/g12**(n12+1)
  t4=l2*(l2+1)*fact(n12-2)/g12**(n12-1)
  
  kinetic_slater=-0.5d0*norm*(t1-t2+t3-t4)
  
end

