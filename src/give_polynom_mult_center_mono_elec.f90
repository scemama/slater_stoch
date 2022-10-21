subroutine give_polynom_mult_center_mono_elec(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out)
  !!!! subroutine that returns the explicit polynom in term of the "t" variable of the following polynomw :: 
  !!!!         I_x1(a_x, d_x,p,q) * I_x1(a_y, d_y,p,q) * I_x1(a_z, d_z,p,q)
  !!!! it is for the nuclear electron atraction
  implicit none
  integer, intent(in)            :: n_pt_in
  integer,intent(out)            :: n_pt_out
  double precision, intent(in)   :: A_center(3), B_center(3),C_center(3)
  double precision, intent(in)   :: alpha,beta
  integer, intent(in)            :: power_A(3), power_B(3)
  integer                        :: a_x,b_x,a_y,b_y,a_z,b_z
  double precision               :: d(0:n_pt_in)
  double precision               :: d1(0:n_pt_in)
  double precision               :: d2(0:n_pt_in)
  double precision               :: d3(0:n_pt_in)
  double precision               :: accu,  pq_inv, p10_1, p10_2, p01_1, p01_2
  double precision               :: p,P_center(3),rho,p_inv,p_inv_2
  double precision               :: R1x(0:2), B01(0:2), R1xp(0:2),R2x(0:2)
  integer                        :: n_pt1,n_pt2,n_pt3,dim,i
  integer                        :: n_pt_tmp
  !print*,'n_pt_in = ',n_pt_in
  accu = 0.d0
  !COMPTEUR irp_rdtsc1 = irp_rdtsc()
  p = alpha+beta
  p_inv = 1.d0/p
  p_inv_2 = 0.5d0/p
  do i =1, 3
    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
  enddo
  ! print*,'passed the P_center'
  
  R1x(0)  = (P_center(1) - A_center(1))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(1) - C_center(1))
  ! R1x = (P_x - A_x) - (P_x - C_x) t^2
  R1xp(0)  = (P_center(1) - B_center(1))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(1) - C_center(1))
  !R1xp = (P_x - B_x) - (P_x - C_x) t^2
  R2x(0)  =  p_inv_2
  R2x(1)  = 0.d0
  R2x(2)  = -p_inv_2
  !R2x  = 0.5 / p - 0.5/p t^2
  do i = 0,n_pt_in
    d(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d1(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d2(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d3(i) = 0.d0
  enddo
  n_pt1 = n_pt_in
  n_pt2 = n_pt_in
  n_pt3 = n_pt_in
  a_x = power_A(1)
  b_x = power_B(1)
  call I_x1_pol_mult_mono_elec(a_x,b_x,R1x,R1xp,R2x,d1,n_pt1,n_pt_in)
  ! print*,'passed the first I_x1'
  if(n_pt1<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif
  
  
  R1x(0)  = (P_center(2) - A_center(2))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(2) - C_center(2))
  ! R1x = (P_x - A_x) - (P_x - C_x) t^2
  R1xp(0)  = (P_center(2) - B_center(2))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(2) - C_center(2))
  !R1xp = (P_x - B_x) - (P_x - C_x) t^2
  a_y = power_A(2)
  b_y = power_B(2)
  call I_x1_pol_mult_mono_elec(a_y,b_y,R1x,R1xp,R2x,d2,n_pt2,n_pt_in)
  ! print*,'passed the second I_x1'
  if(n_pt2<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif
  
  
  R1x(0)  = (P_center(3) - A_center(3))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(3) - C_center(3))
  ! R1x = (P_x - A_x) - (P_x - C_x) t^2
  R1xp(0)  = (P_center(3) - B_center(3))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(3) - C_center(3))
  !R2x  = 0.5 / p - 0.5/p t^2
  a_z = power_A(3)
  b_z = power_B(3)
  
  ! print*,'a_z = ',a_z
  ! print*,'b_z = ',b_z
  call I_x1_pol_mult_mono_elec(a_z,b_z,R1x,R1xp,R2x,d3,n_pt3,n_pt_in)
  ! print*,'passed the third I_x1'
  if(n_pt3<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif
  n_pt_tmp = 0
  call multiply_poly_2(d1,n_pt1,d2,n_pt2,d,n_pt_tmp,n_pt_in)
  do i = 0,n_pt_tmp
    d1(i) = 0.d0
  enddo
  n_pt_out = 0
  call multiply_poly_2(d ,n_pt_tmp ,d3,n_pt3,d1,n_pt_out,n_pt_in)
  do i = 0, n_pt_out
    d(i) = d1(i)
  enddo
  
end

