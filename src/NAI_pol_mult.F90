double precision function NAI_pol_mult(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in)
  ! function that calculate the folowing integral :
  !       int{dr} of (x-A_x)^ax (x-B_X)^bx exp(-alpha (x-A_x)^2 - beta (x-B_x)^2 ) 1/(r-R_c)

  implicit none
  double precision,intent(in)    :: C_center(3),A_center(3),B_center(3),alpha,beta
  integer                        :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,n_pt
  double precision               :: P_center(3)
  double precision               :: d(0:n_pt_in),pouet,coeff,rho,dist,const,pouet_2,p,p_inv,factor
  double precision               :: I_n_special_exact,integrate_bourrin
  double precision               :: V_e_n,const_factor,dist_integral,two_pi
  double precision               :: accu,epsilo,rint
  integer                        :: n_pt_in,n_pt_out,lmax

  two_pi = 2.d0 * dacos(-1.d0)
  p = alpha + beta
  p_inv = 1.d0/p
  rho = alpha * beta * p_inv
  dist = 0.d0
  dist_integral = 0.d0
  do i = 1, 3
    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
    dist = dist+ (A_center(i) - B_center(i))*(A_center(i) - B_center(i))
    dist_integral =  dist_integral +(P_center(i) - C_center(i))*(P_center(i) - C_center(i))
  enddo
  const_factor = dist*rho
  const = p * dist_integral
  factor = dexp(-const_factor)
  coeff = two_pi * factor * p_inv
  lmax = 20

  do i = 0, n_pt_in
    d(i) = 0.d0
  enddo
  n_pt =  2 * ( (power_A(1) + power_B(1)) +(power_A(2) + power_B(2)) +(power_A(3) + power_B(3)) )
  if (n_pt == 0) then
    epsilo = 1.d0
    pouet = rint(0,const)
    NAI_pol_mult = coeff * pouet
    return
  endif

  call give_polynom_mult_center_mono_elec(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out)
  if(n_pt_out<0)then
    NAI_pol_mult = 0.d0
    return
  endif
  accu = 0.d0

  ! 1/r1 standard attraction integral
  epsilo = 1.d0
  ! sum of integrals of type : int {t,[0,1]}  exp-(rho.(P-Q)^2 * t^2) * t^i
  do i =0 ,n_pt_out,2
    if(abs(d(i)).lt.1.d-9)cycle
    accu = accu + d(i) * rint(i/2,const)
  enddo
  NAI_pol_mult = accu * coeff

end
