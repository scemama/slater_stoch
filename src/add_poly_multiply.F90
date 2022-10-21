subroutine add_poly_multiply(b,nb,cst,d,nd)
  implicit none
!  BEGIN_DOC
  ! Add a polynomial multiplied by a constant
  ! D(t) += cst * B(t)
!  END_DOC

    integer                        :: ib, ic, id

  integer, intent(in)            :: nb
  integer, intent(inout)         :: nd
  double precision, intent(in)   :: b(0:nb),cst
  double precision, intent(inout) :: d(0:max(nb,nd))

  nd = max(nd,nb)
  if (nd /= -1) then
    do ib=0,nb
      d(ib) = d(ib) + cst*b(ib)
    enddo
    do while ( d(nd) == 0.d0 )
      nd = nd - 1
      if (nd < 0) then
        exit
      endif
    enddo
  endif

end

