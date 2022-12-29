!! exp(-r**2) is replaced by KT(n,alpha(n),beta(n),r,rc)
!!                          =exp(-alpha(n) r**2 / ( 1-(r/rc(n))**2)**beta(n))
!!
!! where n is such that exp(-rc(n)**2)=10^(-n)
!!
!! KT is expanded as   KT = \sum_{i=1}^ng  c_fit_GTO(i,ng,n)  exp(-g_fit_GTO(i,ng,n) r**2)
!!

subroutine build_mapping_ijkl(nint_out)
  use common_data
  integer*8, intent(out) :: nint_out
  logical :: logic
  integer*8 :: np, nint_theo
  integer*8 :: kcp
  integer :: i, j, k, l
  np=(nbasis*(nbasis+1))/2
  nint_theo=(np*(np+1))/2

  allocate(is(nint_theo))
  allocate(js(nint_theo))
  allocate(ks(nint_theo))
  allocate(ls(nint_theo))

  kcp=0
  do l=1,nbasis
    do k=1,nbasis
      do j=l,nbasis
        do i=k,nbasis
          logic=.true.
          if(i.eq.j.and.l.gt.k)logic=.false.
          if(i.ge.j.and.logic)then
            kcp=kcp+1
            is(kcp)=i
            js(kcp)=j
            ks(kcp)=k
            ls(kcp)=l
          endif
        enddo
      enddo
    enddo
  enddo
  nint_out=kcp
  if(nint_out.ne.nint_theo)stop 'pb in build_mapping_ijkl'
end

