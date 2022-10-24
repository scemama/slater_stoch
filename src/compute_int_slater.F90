subroutine compute_int_slater(i,k,j,l,ijkl_gaus)
  use common_data
  integer, intent(in) :: i, j, k, l
  double precision, intent(out) :: ijkl_gaus

  double precision :: gam(4)
  integer :: n_i(3),n_k(3),n_j(3),n_l(3),nix(4),niy(4),niz(4),power(4)
  character(80) :: orb

  double precision, external :: ijkl_slater_xnynzn

  gam(1)=g_slater(i)
  gam(2)=g_slater(k)
  gam(3)=g_slater(j)
  gam(4)=g_slater(l)

  nix(1)=npower(1,i)
  nix(2)=npower(1,k)
  nix(3)=npower(1,j)
  nix(4)=npower(1,l)

  niy(1)=npower(2,i)
  niy(2)=npower(2,k)
  niy(3)=npower(2,j)
  niy(4)=npower(2,l)

  niz(1)=npower(3,i)
  niz(2)=npower(3,k)
  niz(3)=npower(3,j)
  niz(4)=npower(3,l)

  orb=orb_name(i)
  if(orb.ne.'1S'.and.orb.ne.'2S'.and.orb.ne.'3S'.and.orb.ne.'2P'.and.orb.ne.'3P'.and.orb.ne.'3D')then
    print*,'orb=',orb
    stop 'pb in compute_in_slater'
  endif

  power(1)=0
  if(orb.eq.'2S')power(1)=1
  if(orb.eq.'3S')power(1)=2
  if(orb.eq.'3P')power(1)=1

  orb=orb_name(k)
  power(2)=0
  if(orb.eq.'2S')power(2)=1
  if(orb.eq.'3S')power(2)=2
  if(orb.eq.'3P')power(2)=1

  orb=orb_name(j)
  power(3)=0
  if(orb.eq.'2S')power(3)=1
  if(orb.eq.'3S')power(3)=2
  if(orb.eq.'3P')power(3)=1

  orb=orb_name(l)
  power(4)=0
  if(orb.eq.'2S')power(4)=1
  if(orb.eq.'3S')power(4)=2
  if(orb.eq.'3P')power(4)=1

  ijkl_gaus=ijkl_slater_xnynzn(nix,niy,niz,power,gam)
end
