!!!!  u_orb is the radial part of the atomic orbitals for which the 1- and 2-electron integrals are computed
double precision function u_orb(i_orb,r)
  include 'j.inc'
  u_orb=0.d0
  do mm=1,n_contract(i_orb)
    if(basis_type.eq.'GTO') then
      contrib=c_contract(mm,i_orb)*dexp(-g_contract(mm,i_orb)*r**2)
    endif
    if(basis_type.eq.'STO') then
      contrib=c_contract(mm,i_orb)*r**n_sto(i_orb)*dexp(-g_contract(mm,i_orb)*r   )
      u_orb=u_orb+contrib
    endif
  enddo  !mm
end
