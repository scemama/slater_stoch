!!cc Calculation of  I = <phiA|phiB> = /int dxdydz  phiA(r) phiB(r)
!!cc *************************************************************
!!cc  r=(x,y,z)
!!cc  A=[ra(1),ra(2),ra(3]
!!cc  B=[rb(1),rb(2),rb(3)]
!!cc phiA(r)=[x-ra(1)]**na(1)*[y-ra(2)]**na(2)*[z-ra(3)]**na(3)*exp[-gamA*(x-ra(1))**2-gamA(y-ra(2))**2-gamA(z-ra(3))**2)]
!!cc phiB(r)=[x-rb(1)]**nb(1)*[y-rb(2)]**nb(2)*[z-rb(3)]**nb(3)*exp[-gamB*(x-rb(1))**2-gamB(y-rb(2))**2-gamB(z-rb(3))**2)]

!!cc  mu=gamA*gamB/(gamA+gamB)
!!cc  rp(l)=(gamA*ra(l)+gamB*rb(l))/(gamA+gamB)

!!c  I = exp[-mu*(ra-rb)**2] * Ix*Iy*Iz
!!cc    where:
!!cc  Ix = /int dx [x-ra(1)]**na(1)*[x-rb(1)]**nb(1)*exp[-(gamA+gamB)*(x-rp(1))**2)]
!!cc  Iy = /int dy [y-ra(2)]**na(2)*[y-rb(2)]**nb(2)*exp[-(gamA+gamB)*(y-rp(2))**2)]
!!cc  Iz = /int dz [z-ra(3)]**na(3)*[z-rb(3)]**nb(3)*exp[-(gamA+gamB)*(z-rp(3))**2)]
!!cc
!!cc  one_electron_sij = I

double precision function one_electron_Sij(na,ra,gamA,nb,rb,gamB)
  implicit none
  integer, intent(in)            :: na(3),nb(3)
  double precision, intent(in)   :: ra(3),rb(3), gamA, gamB
  double precision :: rp(3)
  double precision, parameter :: pi=dacos(-1.d0)
  double precision :: gamtot, rmu, AmB, arg
  integer :: l

  double precision, external :: one_electron_I

  gamtot=gamA+gamB
  rmu=gamA*gamB/gamtot

  do l=1,3
    rp(l)=(gamA*ra(l)+gamB*rb(l))/gamtot
  enddo
  AmB=(ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2

  one_electron_sij=0.d0
  arg=rmu*AmB
  if(arg.gt.-dlog(1.d-15))return

  one_electron_sij=1.d0
  do l=1,3
    one_electron_sij=one_electron_sij                                &
        * one_electron_I(na(l),ra(l),nb(l),rb(l),gamtot,rp(l))

  enddo
  one_electron_sij=dexp(-arg)*one_electron_sij
end
