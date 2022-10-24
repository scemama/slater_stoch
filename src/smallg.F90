double precision function smallg(n,m,g,gp,l)
  implicit none
  integer n,m,l,i
  double precision g,gp,fact
  smallg=0.d0
  do i=0,n+l
    smallg=smallg+fact(n+l)*fact(n+m-1-i)/fact(n+l-i)/g**(i+1)/gp**(n+m-i)
  enddo
end



