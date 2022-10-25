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

end
