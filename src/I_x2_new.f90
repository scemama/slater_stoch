recursive subroutine I_x2_new(c,B_10,B_01,B_00,res,n_pt)
  implicit none
  !EGIN_DOC
  !  recursive function involved in the bielectronic integral
  !ND_DOC
  integer, intent(in)            :: c, n_pt
  double precision, intent(in)   :: B_10(1000),B_01(1000),B_00(1000)
  double precision, intent(out)  :: res(1000)
  integer                        :: i
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter             :: max_dim = 511
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

  if(c==1)then
    do i=1,n_pt
      res(i) = 0.d0
    enddo
  elseif(c==0) then
    do i=1,n_pt
      res(i) = 1.d0
    enddo
  else
    call I_x1_new(0,c-2,B_10,B_01,B_00,res,n_pt)
    do i=1,n_pt
      res(i) =  (c-1) * B_01(i) * res(i)
    enddo
  endif
end


