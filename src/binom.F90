double precision function binom(i,j)
 implicit none
 integer :: i,j,k
 double precision :: fact,prod
 if(i.ge.0)then
  binom=fact(i)/(fact(j)*fact(i-j))
  return
 else
  if(j.eq.0)then
   binom=1.d0
   return
  endif
  prod=1.d0
  do k=0,j-1
   prod=prod*(i-k)
  enddo
  binom=prod/fact(j)
 endif
end

