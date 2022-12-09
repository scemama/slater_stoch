subroutine give_explicit_poly_and_gaussian(P_new,P_center,p,fact_k,iorder,alpha,beta,a,b,A_center,B_center,dim)
!  BEGIN_DOC
! Transforms the product of
! (x-x_A)^a(1) (x-x_B)^b(1) (x-x_A)^a(2) (y-y_B)^b(2) (z-z_A)^a(3) (z-z_B)^b(3)
! exp(-(r-A)^2 alpha) exp(-(r-B)^2 beta)
! into
! fact_k *[sum (l_x = 0,i_order(1)) P_new(l_x,1)*(x-P_center(1))^l_x ] exp (- p (x-P_center(1))^2 )
!        *[sum (l_y = 0,i_order(2)) P_new(l_y,2)*(y-P_center(2))^l_y ] exp (- p (y-P_center(2))^2 )
!        *[sum (l_z = 0,i_order(3)) P_new(l_z,3)*(z-P_center(3))^l_z ] exp (- p (z-P_center(3))^2 )
!  END_DOC
  implicit none

integer, parameter :: max_dim = 511

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15

  integer, intent(in)            :: dim
  integer, intent(in)            :: a(3),b(3)         ! powers : (x-xa)**a_x = (x-A(1))**a(1)
  double precision, intent(in)   :: alpha, beta       ! exponents
  double precision, intent(in)   :: A_center(3)       ! A center
  double precision, intent(in)   :: B_center (3)      ! B center
  double precision, intent(out)  :: P_center(3)       ! new center
  double precision, intent(out)  :: p                 ! new exponent
  double precision, intent(out)  :: fact_k            ! constant factor
  double precision, intent(out)  :: P_new(0:max_dim,3)! polynomial
  integer, intent(out)           :: iorder(3)         ! i_order(i) = order of the polynomials

  double precision               :: P_a(0:max_dim,3), P_b(0:max_dim,3)
  integer                        :: n_new,i,j

  iorder(1) = 0
  iorder(2) = 0
  iorder(3) = 0
  P_new(0,1) = 0.d0
  P_new(0,2) = 0.d0
  P_new(0,3) = 0.d0

  call gaussian_product(alpha,A_center,beta,B_center,fact_k,p,P_center)
  if (fact_k < thresh) then
    fact_k = 0.d0
    return
  endif

  !DIR$ FORCEINLINE
  call recentered_poly(P_a(0,1),A_center(1),P_center(1),a(1))
  !DIR$ FORCEINLINE
  call recentered_poly(P_b(0,1),B_center(1),P_center(1),b(1))
  iorder(1) = a(1) + b(1)
  do i=0,iorder(1)
    P_new(i,1) = 0.d0
  enddo
  n_new=0
  call multiply_poly(P_a(0,1),a(1),P_b(0,1),b(1),P_new(0,1),n_new)

  !DIR$ FORCEINLINE
  call recentered_poly(P_a(0,2),A_center(2),P_center(2),a(2))
  !DIR$ FORCEINLINE
  call recentered_poly(P_b(0,2),B_center(2),P_center(2),b(2))
  iorder(2) = a(2) + b(2)
  do i=0,iorder(2)
    P_new(i,2) = 0.d0
  enddo
  n_new=0
  call multiply_poly(P_a(0,2),a(2),P_b(0,2),b(2),P_new(0,2),n_new)

  !DIR$ FORCEINLINE
  call recentered_poly(P_a(0,3),A_center(3),P_center(3),a(3))
  !DIR$ FORCEINLINE
  call recentered_poly(P_b(0,3),B_center(3),P_center(3),b(3))
  iorder(3) = a(3) + b(3)
  do i=0,iorder(3)
    P_new(i,3) = 0.d0
  enddo
  n_new=0
  call multiply_poly(P_a(0,3),a(3),P_b(0,3),b(3),P_new(0,3),n_new)

end

subroutine recentered_poly(P_new,x_A,x_P,a)
  implicit none
!  BEGIN_DOC
  ! Recenter two polynomials
!  END_DOC
  integer, intent(in)            :: a
  double precision, intent(in)   :: x_A,x_P
  double precision, intent(out)  :: P_new(0:a)
  double precision               :: pows_a(-2:100)
  double precision               :: binom_func
  integer                        :: i,j,k,l

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
!    for a in range(0,9):
!        print("    case(%d)"%(a))
!        for i in range(2,a+1):
!            print(f"      pows_a({i}) = pows_a({i-1})*pows_a(1)")
!        print(f"      P_new (0) =  pows_a({a})")
!        for i in range(1,a+1):
!            binom = "%fd0"%(binom_transp(a-i,a))
!            print(f"      P_new ({i}) =  {binom} * pows_a({a-i})")
!    print("  end select")


  if ((a<0)) return
  pows_a(0) = 1.d0
  pows_a(1) = (x_P - x_A)

  select case(a)
    case(0)
      P_new (0) =  pows_a(0)
    case(1)
      P_new (0) =  pows_a(1)
      P_new (1) =  pows_a(0)
    case(2)
      pows_a(2) = pows_a(1)*pows_a(1)
      P_new (0) =  pows_a(2)
      P_new (1) =  2.d0 * pows_a(1)
      P_new (2) =  pows_a(0)
    case(3)
      pows_a(2) = pows_a(1)*pows_a(1)
      pows_a(3) = pows_a(2)*pows_a(1)
      P_new (0) =  pows_a(3)
      P_new (1) =  3.d0 * pows_a(2)
      P_new (2) =  3.d0 * pows_a(1)
      P_new (3) =  pows_a(0)
    case(4)
      pows_a(2) = pows_a(1)*pows_a(1)
      pows_a(3) = pows_a(2)*pows_a(1)
      pows_a(4) = pows_a(2)*pows_a(2)
      P_new (0) =  pows_a(4)
      P_new (1) =  4.d0 * pows_a(3)
      P_new (2) =  6.d0 * pows_a(2)
      P_new (3) =  4.d0 * pows_a(1)
      P_new (4) =  pows_a(0)
    case(5)
      pows_a(2) = pows_a(1)*pows_a(1)
      pows_a(3) = pows_a(2)*pows_a(1)
      pows_a(4) = pows_a(3)*pows_a(1)
      pows_a(5) = pows_a(4)*pows_a(1)
      P_new (0) =  pows_a(5)
      P_new (1) =  5.000000d0 * pows_a(4)
      P_new (2) =  10.000000d0 * pows_a(3)
      P_new (3) =  10.000000d0 * pows_a(2)
      P_new (4) =  5.000000d0 * pows_a(1)
      P_new (5) =  1.000000d0 * pows_a(0)
    case(6)
      pows_a(2) = pows_a(1)*pows_a(1)
      pows_a(3) = pows_a(2)*pows_a(1)
      pows_a(4) = pows_a(3)*pows_a(1)
      pows_a(5) = pows_a(4)*pows_a(1)
      pows_a(6) = pows_a(5)*pows_a(1)
      P_new (0) =  pows_a(6)
      P_new (1) =  6.000000d0 * pows_a(5)
      P_new (2) =  15.000000d0 * pows_a(4)
      P_new (3) =  20.000000d0 * pows_a(3)
      P_new (4) =  15.000000d0 * pows_a(2)
      P_new (5) =  6.000000d0 * pows_a(1)
      P_new (6) =  1.000000d0 * pows_a(0)
    case(7)
      pows_a(2) = pows_a(1)*pows_a(1)
      pows_a(3) = pows_a(2)*pows_a(1)
      pows_a(4) = pows_a(3)*pows_a(1)
      pows_a(5) = pows_a(4)*pows_a(1)
      pows_a(6) = pows_a(5)*pows_a(1)
      pows_a(7) = pows_a(6)*pows_a(1)
      P_new (0) =  pows_a(7)
      P_new (1) =  7.000000d0 * pows_a(6)
      P_new (2) =  21.000000d0 * pows_a(5)
      P_new (3) =  35.000000d0 * pows_a(4)
      P_new (4) =  35.000000d0 * pows_a(3)
      P_new (5) =  21.000000d0 * pows_a(2)
      P_new (6) =  7.000000d0 * pows_a(1)
      P_new (7) =  1.000000d0 * pows_a(0)
    case(8)
      pows_a(2) = pows_a(1)*pows_a(1)
      pows_a(3) = pows_a(2)*pows_a(1)
      pows_a(4) = pows_a(3)*pows_a(1)
      pows_a(5) = pows_a(4)*pows_a(1)
      pows_a(6) = pows_a(5)*pows_a(1)
      pows_a(7) = pows_a(6)*pows_a(1)
      pows_a(8) = pows_a(7)*pows_a(1)
      P_new (0) =  pows_a(8)
      P_new (1) =  8.000000d0 * pows_a(7)
      P_new (2) =  28.000000d0 * pows_a(6)
      P_new (3) =  56.000000d0 * pows_a(5)
      P_new (4) =  70.000000d0 * pows_a(4)
      P_new (5) =  56.000000d0 * pows_a(3)
      P_new (6) =  28.000000d0 * pows_a(2)
      P_new (7) =  8.000000d0 * pows_a(1)
      P_new (8) =  1.000000d0 * pows_a(0)
    case default
      do i =  2,a
        pows_a(i) = pows_a(i-1)*pows_a(1)
      enddo
      P_new (0) =  pows_a(a)
      do i =  1,a
        P_new(i) =  binom_func(a,a-i) * pows_a(a-i)
      enddo
  end select

end

