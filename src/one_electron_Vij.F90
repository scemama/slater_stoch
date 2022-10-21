double precision function one_electron_Vij(na,ra,gamA,nb,rb,gamB)
  include 'j.inc'
  dimension na(3),ra(3),nb(3),rb(3)
  
  double precision NAI_pol_mult
  
  integer dim_int_tab(4)
  integer dim_int
  
  dim_int_tab(1) = maxval(na)
  dim_int_tab(2) = maxval(nb)
  dim_int_tab(3) = 0
  dim_int_tab(4) = 0
  dim_int = (maxval(dim_int_tab) *24 + 4)*8
  
  ac=0.d0
  do kn=1,number_nuclei
    ac=ac-charge(kn)*NAI_pol_mult(ra,rb,na,nb,gamA,gamB,centers_nuclei(1,kn),dim_int)
  enddo
  
  one_electron_Vij=ac
  
end
