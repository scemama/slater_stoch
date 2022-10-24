! F is defined as:

! F(k1,g1,k2,g2,l)= int_0^inf r1^k1 exp(-g1 r1)
!                                             [ 1/r1**(l+1) int_0^r1 r2**(k2+l) exp(-g2 r2)
!                                             +   r1**l     int_r1^inf  r2**(k2-l-1) exp(-g2 r2) ]
!
! k1 integer ge 0
! k2 integer ge 0
! k1 and k2 ge (l+1)
!
! We can show:
!
!  F= (k1-l-1)! (k2+l)!/g1**(k1-l)/g2**(k2+l+1)  + \sum_{i=0}^{k2-l-1} (k2-l-1)!(k1+k2-1-i)!/(k2-l-1-i)!/g2**(i+1)/(g1+g2)**(k1+k2-i)
!                                                 - sum_{i=0}^{k2+l  } (k2+l  )!(k1+k2-1-i)!/(k2+l-i  )!/g2**(i+1)/(g1+g2)**(k1+k2-i)
!
! This formula can been checked with bigF_num
!
! Here F must be computed with  k1=n12  k2=n34 g1 -->g12 et g2 --->g13
!
!
double precision function bigF(n12,g12,n34,g34,l)
  implicit none
  integer :: n12,n34,l,l12,l34
  double precision :: g12,g34,g,first_term,smallg,fact

  l12=-l
  l34=l+1
  g=g12+g34

  first_term= fact(n34+l34-1)/g34**(n34+l34)*fact(n12+l12-1)/g12**(n12+l12)

  bigF=first_term+ smallg(n34,n12,g34,g,l12-1)-smallg(n34,n12,g34,g,l34-1)
end

