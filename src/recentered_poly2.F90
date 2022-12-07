subroutine recentered_poly2(P_new,x_A,x_P,a,P_new2,x_B,x_Q,b)
  implicit none
!  BEGIN_DOC
  ! Recenter two polynomials
!  END_DOC
  integer, intent(in)            :: a,b
  double precision, intent(in)   :: x_A,x_P,x_B,x_Q
  double precision, intent(out)  :: P_new(0:a),P_new2(0:b)
  double precision               :: pows_a(-2:a+b+4), pows_b(-2:a+b+4)
  double precision               :: binom_func
  double precision               :: binom_transp
  integer                        :: i,j,k,l, minab, maxab

!#!/usr/bin/env python
!from math import exp, log
!
!def logfact(k):
!  logfact = 0.
!  for i in range(1,k+1):
!     logfact += log(float(i))
!  return logfact
!
!def binom_transp(j,i):
!   return exp( logfact(i)-logfact(j)-logfact(i-j) )
!
!def main():
!    print("  select case(a)")
!    for a in range(0,5):
!        print("    case(%d)"%(a))
!        for i in range(2,a+1):
!            print(f"      pows_a({i}) = pows_a({i-1})*pows_a(1)")
!        print(f"      P_new (0) =  pows_a({a})")
!        for i in range(1,a+1):
!            binom = "%fd0"%(binom_transp(a-i,a))
!            print(f"      P_new ({i}) =  {binom} * pows_a({a-i})")
!    print("  end select")
!    print("  select case(b)")
!    for b in range(0,5):
!        print("    case(%d)"%(b))
!        for i in range(2,b+1):
!            print(f"      pows_b({i}) = pows_b({i-1})*pows_b(1)")
!        print(f"      P_new2(0) =  pows_b({b})")
!        for i in range(1,b+1):
!            binom = "%fd0"%(binom_transp(b-i,b))
!            print(f"      P_new2({i}) =  {binom} * pows_b({b-i})")
!    print("  end select")




  if ((a<0).or.(b<0) ) return
  pows_a(0) = 1.d0
  pows_a(1) = (x_P - x_A)
  pows_b(0) = 1.d0
  pows_b(1) = (x_Q - x_B)

  select case(a)
    case(0)
      P_new (0) =  pows_a(0)
    case(1)
      P_new (0) =  pows_a(1)
      P_new (1) =  pows_a(0)
    case(2)
      pows_a(2) = pows_a(1)*pows_a(1)
      P_new (0) =  pows_a(2)
      P_new (1) =  2.000000d0 * pows_a(1)
      P_new (2) =  pows_a(0)
    case(3)
      pows_a(2) = pows_a(1)*pows_a(1)
      pows_a(3) = pows_a(2)*pows_a(1)
      P_new (0) =  pows_a(3)
      P_new (1) =  3.000000d0 * pows_a(2)
      P_new (2) =  3.000000d0 * pows_a(1)
      P_new (3) =  pows_a(0)
    case(4)
      pows_a(2) = pows_a(1)*pows_a(1)
      pows_a(3) = pows_a(2)*pows_a(1)
      pows_a(4) = pows_a(3)*pows_a(1)
      P_new (0) =  pows_a(4)
      P_new (1) =  4.000000d0 * pows_a(3)
      P_new (2) =  6.000000d0 * pows_a(2)
      P_new (3) =  4.000000d0 * pows_a(1)
      P_new (4) =  pows_a(0)
    case default
      do i =  2,a
        pows_a(i) = pows_a(i-1)*pows_a(1)
      enddo
      P_new (0) =  pows_a(a)
      do i =  1,a
        P_new(i) =  binom_func(a,a-i) * pows_a(a-i)
      enddo
  end select
  select case(b)
    case(0)
      P_new2(0) =  pows_b(0)
    case(1)
      P_new2(0) =  pows_b(1)
      P_new2(1) =  pows_b(0)
    case(2)
      pows_b(2) = pows_b(1)*pows_b(1)
      P_new2(0) =  pows_b(2)
      P_new2(1) =  2.000000d0 * pows_b(1)
      P_new2(2) =  pows_b(0)
    case(3)
      pows_b(2) = pows_b(1)*pows_b(1)
      pows_b(3) = pows_b(2)*pows_b(1)
      P_new2(0) =  pows_b(3)
      P_new2(1) =  3.000000d0 * pows_b(2)
      P_new2(2) =  3.000000d0 * pows_b(1)
      P_new2(3) =  pows_b(0)
    case(4)
      pows_b(2) = pows_b(1)*pows_b(1)
      pows_b(3) = pows_b(2)*pows_b(1)
      pows_b(4) = pows_b(3)*pows_b(1)
      P_new2(0) =  pows_b(4)
      P_new2(1) =  4.000000d0 * pows_b(3)
      P_new2(2) =  6.000000d0 * pows_b(2)
      P_new2(3) =  4.000000d0 * pows_b(1)
      P_new2(4) =  pows_b(0)
    case default
      do i =  2,b
        pows_b(i) = pows_b(i-1)*pows_b(1)
      enddo
      P_new2(0) =  pows_b(b)
      do i =  1,b
        P_new2(i) =  binom_func(b,b-i) * pows_b(b-i)
      enddo
  end select

end

