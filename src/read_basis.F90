subroutine read_basis(filename_basis)
  use common_data
  character(128), intent(in)    :: filename_basis
  character(80) ATOM_READ
  character(80) orb

  integer :: i, i_atom
  double precision :: gam

  integer, external :: number_atom
  !! read basis set for the 36 first atoms
  !! n_b(i) number of basis (S,P,D,F, or G) centered on nucleus i=1,36
  !! orb_b(k,i) name of the basis function, that is S,P,D,F or G for k=1,n_b(i)
  !! n_cont_b(k,i) number of contracted primitives for the radial part
  !! for the moment STO have only one primtive
  !! u_STO(r)=r**n_sto* exp(-gi r)
  !! u_GTO(r) = sum_i c_i exp(-gi r**2)
  !! g_i and c_i  gamma_b(k,m,i),coef_b(k,m,i) m=1,n_cont_b(k,i)

  open(1,file=filename_basis)
  do i=1,n_atoms_max
    n_b(i)=0
  enddo
  do i=1,n_atoms_max
    read(1,*,end=1000)ATOM_READ,orb,gam
    i_atom=number_atom(ATOM_READ)
    n_b(i_atom)=n_b(i_atom)+1
    if(n_b(i_atom).gt.nbasis_max)stop 'in read_basis: increase nbasis_max'
    gamma_b(n_b(i_atom),1,i_atom)=gam
    orb_b(n_b(i_atom),i_atom)=orb
    n_cont_b(n_b(i_atom),i_atom)=1
    coef_b(n_b(i_atom),1,i_atom)=1.d0
  enddo
1000  continue

end
