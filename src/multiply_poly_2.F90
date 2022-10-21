subroutine multiply_poly_2(b,nb,c,nc,d,nd,dim_int)
  implicit none
  ! D(t) += B(t)*C(t)
  integer, intent(in)            :: dim_int
  integer, intent(in)            :: nb, nc
  integer, intent(inout)         :: nd
  double precision, intent(in)   :: b(0:dim_int), c(0:dim_int)
  double precision, intent(inout) :: d(0:dim_int)
  
  integer                        :: ndtmp
  integer                        :: ib, ic, id
  if(nc==-1.or.nb==-1)then
    return
  endif
  ndtmp = nb+nc
  do ib=0,nb
    !print*,ib
    if ( b(ib)==0.d0) then
      cycle
    endif
    do ic = 0,nc
      d(ib+ic) = d(ib+ic) + c(ic) * b(ib)
    enddo
  enddo
  ndtmp = max(ndtmp,nd)
  !do while (d(ndtmp) == 0.d0)
  !   ndtmp -= 1
  !   if(ndtmp.lt.0)then
  !    exit
  !   endif
  !enddo
  do ndtmp = max(ndtmp,nd),0,-1
    if (d(ndtmp) /= 0.d0) then
      exit
    endif
  enddo
  nd = ndtmp
  
end

