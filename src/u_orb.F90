!!!!  u_orb is the radial part of the atomic orbitals for which the 1- and 2-electron integrals are computed
double precision function u_orb(i_orb,r)
  use common_data
  integer, intent(in) :: i_orb
  double precision, intent(in) :: r
  integer :: mm
  double precision :: r2
  u_orb=0.d0
  select case (n_sto(i_orb))
    case (0)
      do mm=1,n_contract(i_orb)
        u_orb=u_orb+c_contract(mm,i_orb)*dexp(-g_contract(mm,i_orb)*r)
      enddo
    case (1)
      do mm=1,n_contract(i_orb)
        u_orb=u_orb+c_contract(mm,i_orb)*r*dexp(-g_contract(mm,i_orb)*r)
      enddo
    case (2)
      r2 = r*r
      do mm=1,n_contract(i_orb)
        u_orb=u_orb+c_contract(mm,i_orb)*r2*dexp(-g_contract(mm,i_orb)*r)
      enddo
    case (3)
      r2 = r*r
      do mm=1,n_contract(i_orb)
        u_orb=u_orb+c_contract(mm,i_orb)*r2*r*dexp(-g_contract(mm,i_orb)*r)
      enddo
    case (4)
      do mm=1,n_contract(i_orb)
        u_orb=u_orb+c_contract(mm,i_orb)*r2*r2*dexp(-g_contract(mm,i_orb)*r)
      enddo
    case default
      do mm=1,n_contract(i_orb)
        u_orb=u_orb+c_contract(mm,i_orb)*r**n_sto(i_orb)*dexp(-g_contract(mm,i_orb)*r)
      enddo
  end select
end
