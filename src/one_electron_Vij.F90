double precision function one_electron_Vij(na,ra,gamA,nb,rb,gamB)
  include 'j.inc'

  integer, intent(in) :: na(3),nb(3)
  double precision, intent(in) :: ra(3),rb(3), gamA, gamB

  integer :: dim_int, kn
  double precision :: ac

  double precision, external :: NAI_pol_mult

  dim_int = (max(maxval(na), maxval(nb)) *24 + 4)*8

  ac=0.d0
  do kn=1,number_nuclei
    ac=ac-charge(kn)*NAI_pol_mult(ra,rb,na,nb,gamA,gamB,centers_nuclei(1,kn),dim_int)
  enddo

  one_electron_Vij=ac

end
