subroutine allocate_basis()
  use common_data
  implicit none

  allocate(nucleus_number(nbasis))
  allocate(n_STO(nbasis))
  allocate(npower(3,nbasis), i_type(nbasis), n_contract(nbasis))
  allocate(n_gauss(nbasis))
  allocate(ng(n_contract_max, nbasis))
  allocate(center(3,nbasis))
  allocate(orb_name(nbasis))
  allocate(orb_name_full(nbasis))

  allocate(a_ZV(nbasis,nbasis))
  allocate(c_gauss(n_gauss_max*ng_max,nbasis,2))
  allocate(g_gauss(n_gauss_max*ng_max,nbasis))
  allocate(g_slater(nbasis))
  allocate(g_contract(n_contract_max,nbasis))
  allocate(c_contract(n_contract_max,nbasis))

end
