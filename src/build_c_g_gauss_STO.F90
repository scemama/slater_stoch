subroutine build_c_g_gauss_STO(i_orb)
  use common_data
  integer, intent(in) :: i_orb
  integer :: kcp, mm, ii
  kcp=0
  do mm=1,n_contract(i_orb)
    if(i_type(i_orb).eq.1)then
      do ii=1,ng(mm,i_orb)
        kcp=kcp+1
        c_gauss(1,kcp,i_orb)=c_contract(mm,i_orb)*c_fit_exp(ii,ng(mm,i_orb))
        c_gauss(2,kcp,i_orb)=0.d0
        g_gauss(kcp,i_orb)=g_contract(mm,i_orb)**2*g_fit_exp(ii,ng(mm,i_orb))
      enddo ! ii
    endif
    if(i_type(i_orb).eq.2)then
      do ii=1,ng(mm,i_orb)
        kcp=kcp+1
        c_gauss(1,kcp,i_orb)=0.d0
        c_gauss(2,kcp,i_orb)=c_contract(mm,i_orb)                    &
            *2.d0*c_fit_exp(ii,ng(mm,i_orb))*g_fit_exp(ii,ng(mm,i_orb))*g_contract(mm,i_orb)
        g_gauss(kcp,i_orb)=g_contract(mm,i_orb)**2*g_fit_exp(ii,ng(mm,i_orb))
      enddo ! ii
    endif
    if(i_type(i_orb).eq.3)stop 'i_type=3 not yet coded!'
  enddo !mm
  n_gauss(i_orb)=kcp

end


