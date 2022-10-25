!! u_gauss is the radial part of the end-of-the-day gaussian representation
double precision function u_gauss(i,r)
  use common_data
  integer, intent(in) :: i
  double precision, intent(in) :: r
  integer :: mm
  double precision :: r2
  r2 = r*r
  u_gauss=0.d0
  do mm=1,n_gauss(i)
    u_gauss=u_gauss+(c_gauss(mm,i,1)+c_gauss(mm,i,2)*r2)*dexp(-g_gauss(mm,i)*r2)
  enddo
end
