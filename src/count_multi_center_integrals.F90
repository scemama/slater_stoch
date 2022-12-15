subroutine count_multi_center_integrals
  use common_data
  integer :: iac_1, iac_2, iac_3, iac_4
  integer*8 :: kcp
  integer :: i, j, k, l, ndiff

  iac_1=0
  iac_2=0
  iac_3=0
  iac_4=0
  do kcp=1,nint
    i=is(kcp)
    k=ks(kcp)
    j=js(kcp)
    l=ls(kcp)
    call compare_nuclei(nucleus_number(i),nucleus_number(j),nucleus_number(k),nucleus_number(l)&
        ,ndiff)
    if(ndiff.eq.1)then
      iac_1=iac_1+1
    endif
    if(ndiff.eq.2)iac_2=iac_2+1
    if(ndiff.eq.3)iac_3=iac_3+1
    if(ndiff.eq.4)iac_4=iac_4+1
  enddo
  write(*,*)'Number of one-center two-electron integrals =',iac_1
  write(*,*)'Number of two-center two-electron integrals =',iac_2
  write(*,*)'Number of three-center two-electron integrals =',iac_3
  write(*,*)'Number of four-center two-electron integrals =',iac_4
  write(*,*)'Total number of (ik|1/r12|jl) =',iac_1+iac_2+iac_3+iac_4
end

