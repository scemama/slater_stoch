!! Computation of overlap (n1l1m1 | n2l2m2)=  N_n1*N_n2 !(n1+n2)!(g1+g2)**(n1+n2+1)   l1=l1 m1=m2
!!
double precision function overlap_slater(nslat,lslat,mslat,gamma)
  implicit none
  integer nslat(4),lslat(4),mslat(4),i
  integer l1,l2,m1,m2,n12
  double precision gamma(4),int,norm,norm_sto,g12,fact
  
  
  overlap_slater=0.d0
  
  l1=lslat(1)
  l2=lslat(2)
  if(l1.ne.l2)return
  m1=mslat(1)
  m2=mslat(2)
  if(m1.ne.m2)return
  
  norm=1.d0
  do i=1,2
    norm=norm*norm_sto(gamma(i),nslat(i))
  enddo
  
  n12=nslat(1)+nslat(2)
  g12=gamma(1)+gamma(2)
  
  overlap_slater=norm*fact(n12)/g12**(n12+1)
  
end

