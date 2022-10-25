subroutine allocate_basis()
  use common_data
  implicit none

  allocate(nucleus_number(nbasis_max))
  allocate(n_STO(nbasis_max))
  allocate(npower(3,nbasis_max), i_type(nbasis_max), n_contract(nbasis_max))
  allocate(n_gauss(nbasis_max))
  allocate(ng(n_contract_max, nbasis_max))
  allocate(center(3,nbasis_max))
  allocate(orb_b(nbasis_max, n_atoms_max))
  allocate(orb_name(nbasis_max))
  allocate(orb_name_full(nbasis_max))


end
