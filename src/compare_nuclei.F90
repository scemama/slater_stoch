subroutine compare_nuclei(ni,nj,nk,nl,ndiff)
  implicit none
  integer, intent(in)            :: ni, nj, nk, nl
  integer, intent(out)           :: ndiff

  ! 1111
  if ( (ni == nj) &
  .and.(nj == nk) &
  .and.(nk == nl) ) then
    ndiff=1
    return
  endif

  ! 2111
  if ( (ni /= nj) &
  .and.(nj == nk) &
  .and.(nk == nl) ) then
    ndiff=2
    return
  endif

  ! 1211
  if ( (ni /= nj) &
  .and.(ni == nk) &
  .and.(nk == nl) ) then
    ndiff=2
    return
  endif

  ! 1121
  if ( (ni == nj) &
  .and.(nj /= nk) &
  .and.(nl == ni) ) then
    ndiff=2
    return
  endif

  ! 1112
  if ( (ni == nj) &
  .and.(nj == nk) &
  .and.(nk /= nl) ) then
    ndiff=2
    return
  endif

  ! 1212
  if ( (ni /= nj) &
  .and.(nk /= nl) &
  .and.(ni == nk) &
  .and.(nj == nl) ) then
    ndiff=2
    return
  endif

  ! 2211
  if ( (ni == nj) &
  .and.(nk == nl) &
  .and.(ni /= nk) ) then
    ndiff=2
    return
  endif

  ! 1221
  if ( (ni /= nj) &
  .and.(nj == nk) &
  .and.(nk /= nl) &
  .and.(ni == nl) ) then
    ndiff=2
    return
  endif


  !11xy
  if ( (ni == nj) &
  .and.(ni /= nk) &
  .and.(ni /= nl) &
  .and.(nk /= nl) ) then
    ndiff=3
    return
  endif

  !1x1y
  if ( (ni == nk) &
  .and.(ni /= nj) &
  .and.(ni /= nl) &
  .and.(nj /= nl) ) then
    ndiff=3
    return
  endif

  !1xy1
  if ( (ni == nl) &
  .and.(ni /= nj) &
  .and.(ni /= nk) &
  .and.(nj == nk) ) then
    ndiff=3
    return
  endif

  !x1y1
  if ( (nj == nl) &
  .and.(nj /= ni) &
  .and.(nj /= nk) &
  .and.(ni /= nk) ) then
    ndiff=3
    return
  endif

  !x11y
  if ( (nj == nk) &
  .and.(nj /= ni) &
  .and.(nj /= nl) &
  .and.(ni /= nl) ) then
    ndiff=3
    return
  endif

  !xy11
  if ( (nk == nl) &
  .and.(nk /= ni) &
  .and.(nk /= nj) &
  .and.(ni /= nj) ) then
    ndiff=3
    return
  endif

  !1234
  if ( (ni /= nj) &
  .and.(nj /= nk) &
  .and.(nk /= nl) &
  .and.(nl /= ni) &
  .and.(nl /= nj) &
  .and.(ni /= nk) ) then
    ndiff=4
    return
  endif

end

