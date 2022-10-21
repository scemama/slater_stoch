double precision function bielec_integral(alpha,beta,delta,gama,A_center,B_center,&
      C_center,D_center,power_A,power_B,power_C,power_D)

  implicit none
  !DOC
  !  integral of the AO basis <ik|jl> =  (ij|kl)
  !                           <AC|BD> =  (AB|CD)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  !DOC

  integer, parameter             :: max_dim = 511
  integer,intent(in)             :: power_A(3),power_B(3),power_C(3),power_D(3)
  double precision, intent(in)   :: alpha,beta,delta,gama
  double precision, intent(in)   :: A_center(3), B_center(3), C_center(3), D_center(3)
  integer                        :: i,j,k,l
  integer                        :: p,q,r,s
  integer                        :: num_i,num_j,num_k,num_l,dim1
  double precision               :: integral
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: ao_bielec_integral_schwartz_accel
  double precision               :: p_inv,q_inv
  double precision               :: general_primitive_integral

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision, parameter    :: pi =  dacos(-1.d0)
  double precision, parameter    :: sqpi =  dsqrt(dacos(-1.d0))
  double precision, parameter    :: pi_5_2 =  34.9868366552d0
  double precision, parameter    :: dfour_pi =  4.d0*dacos(-1.d0)
  double precision, parameter    :: dtwo_pi =  2.d0*dacos(-1.d0)
  double precision, parameter    :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
  double precision, parameter    :: thresh = 1.d-15
  double precision, parameter    :: cx_lda = -0.73855876638202234d0
  double precision, parameter    :: c_2_4_3 = 2.5198420997897464d0
  double precision, parameter    :: cst_lda = -0.93052573634909996d0
  double precision, parameter    :: c_4_3 = 1.3333333333333333d0
  double precision, parameter    :: c_1_3 = 0.3333333333333333d0

  dim1 = 511

  call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
      alpha,beta,power_A,power_B,A_center,B_center,dim1)
  p_inv = 1.d0/pp

  call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
      delta,gama,power_C,power_D,C_center,D_center,dim1)
  q_inv = 1.d0/qq

  bielec_integral = general_primitive_integral(dim1,                 &
      P_new,P_center,fact_p,pp,p_inv,iorder_p,                       &
      Q_new,Q_center,fact_q,qq,q_inv,iorder_q)

end
