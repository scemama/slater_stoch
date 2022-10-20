subroutine compute_lmin_lmax(l1,l2,l3,l4,m1,m2,m3,m4,lmin,lmax)
  implicit none
  integer l1,l2,l3,l4,m1,m2,m3,m4,lmin,lmax
  integer mx12,l12min,mx34,l34min
  
  lmax=min(l1+l2,l3+l4)
  
  mx12=max(iabs(l1-l2),iabs(m1-m2))
  if(mod(lmax+mx12,2).eq.0)then
    l12min=mx12
  else
    l12min=mx12+1
  endif
  
  mx34=max(iabs(l3-l4),iabs(m3-m4))
  if(mod(lmax+mx34,2).eq.0)then
    l34min=mx34
  else
    l34min=mx34+1
  endif
  
  lmin=max(l12min,l34min)
  
end

