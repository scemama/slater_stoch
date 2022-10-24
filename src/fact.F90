double precision function fact(n)
  implicit none
  !EGIN_DOC
  ! n!
  !ND_DOC
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1
  integer                        :: i
  double precision :: logfact

  if (n<=memomax) then
    if (n<2) then
      fact = 1.d0
    else
      fact = memo(n)
    endif
    return
  endif

  memo(1) = 1.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)*dble(i)
  enddo
  memomax = min(n,100)
  fact = dexp(logfact(n))
end function

double precision function fact_inv(i)
 implicit none
 integer :: i
 double precision :: fact
 fact_inv = 1.d0/fact(i)
end

