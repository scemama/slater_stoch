!! exp(-r**2) is replaced by KT(n,alpha(n),beta(n),r,rc)
!!                          =exp(-alpha(n) r**2 / ( 1-(r/rc(n))**2)**beta(n))
!!
!! where n is such that exp(-rc(n)**2)=10^(-n)
!!
!! KT is expanded as   KT = \sum_{i=1}^ng  c_fit_GTO(i,ng,n)  exp(-g_fit_GTO(i,ng,n) r**2)
!!

subroutine build_mapping_ijkl
  include 'j.inc'
  logical :: logic
  integer :: np, nint_theo
  integer :: kcp, i, j, k, l
  np=(nbasis*(nbasis+1))/2
  nint_theo=(np*(np+1))/2
  if(nint_theo.gt.nint_max)then
    print*,'expected number of two-electron integrals is',nint_theo
    print*,'declared size of array for two-electron integrals is',nint_max
    print*,' I continue but you shall need to increase nint_max'
  endif
  
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
            kcp_ijkl(i,j,k,l)=kcp
          endif
        enddo
      enddo
    enddo
  enddo
  nint=kcp
  if(nint.ne.nint_theo)stop 'pb in build_mapping_ijkl'
end

