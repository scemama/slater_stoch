!!!!  u_orb is the radial part of the atomic orbitals for which the 1- and 2-electron integrals are computed
double precision function u_orb(i_orb,r)
  use common_data
  integer, intent(in) :: i_orb
  double precision, intent(in) :: r
  integer :: mm
  double precision :: r2
  u_orb=0.d0
  do mm=1,n_contract(i_orb)
    u_orb=u_orb+c_contract(mm,i_orb)*r**n_sto(i_orb)*dexp(-g_contract(mm,i_orb)*r)
  enddo
end
