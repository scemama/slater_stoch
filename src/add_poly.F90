subroutine add_poly(b,nb,c,nc,d,nd)
  implicit none
!  BEGIN_DOC
  ! Add two polynomials
  ! D(t) += B(t)+C(t)
!  END_DOC

  integer                        :: ib, ic, id

  integer, intent(inout)         :: nb, nc
  integer, intent(out)           :: nd
  double precision, intent(in)   :: b(0:nb), c(0:nc)
  double precision, intent(out)  :: d(0:nb+nc)

  nd = nb+nc
  do ib=0,max(nb,nc)
    d(ib) = d(ib) + c(ib) + b(ib)
  enddo
  do while ( (d(nd) == 0.d0).and.(nd>=0) )
    nd = nd -1
    if (nd < 0) then
      exit
    endif
  enddo

end

