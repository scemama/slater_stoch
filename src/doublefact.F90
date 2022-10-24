double precision function doublefact(n)
  implicit none
  integer, intent(in) :: n
  double precision :: d
  integer :: i
  doublefact=1.d0
  if(n.le.2)return
  d=1.d0
  do i=n,1,-2
    d=d*dble(i)
  enddo
  doublefact=d
end

