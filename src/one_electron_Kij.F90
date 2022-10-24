double precision function one_electron_Kij(na,ra,gamA,nb,rb,gamB)
  implicit none
  integer, intent(in)            :: na(3), nb(3)
  double precision, intent(in)   :: ra(3), rb(3), gamA, gamB
  double precision               :: rp(3)
  integer                        :: na2(3)
  double precision :: sab, term1, sab_2, term2, term3, sab_12
  integer :: l, li

  double precision, external :: one_electron_Sij


  sab=one_electron_Sij(na,ra,gamA,nb,rb,gamB)
  term1=gamA*(3.d0+2.d0*(na(1)+na(2)+na(3)))*sab

  sab_2=0.d0
  do l=1,3
    do li=1,3
      if(li.eq.l)then
        na2(li)=na(li)+2
      else
        na2(li)=na(li)
      endif
    enddo
    sab_2=sab_2+one_electron_Sij(na2,ra,gamA,nb,rb,gamB)
  enddo

  term2=-2.d0*gamA**2*sab_2

  sab_12=0.d0
  do l=1,3
    do li=1,3
      if(li.eq.l.and.na(li).ge.2)then
        na2(li)=na(li)-2
      else
        na2(li)=na(li)
      endif
    enddo
    sab_12=sab_12+na(l)*(na(l)-1)*one_electron_Sij(na2,ra,gamA,nb,rb,gamB)
  enddo

  term3=-0.5d0*sab_12

  one_electron_Kij=term1+term2+term3

end
