subroutine compare_nuclei(ni,nj,nk,nl,ndiff)
  implicit none
  integer, intent(in)            :: ni, nj, nk, nl
  integer, intent(out)           :: ndiff

  ! 1111
  if( (ni.eq.nj).and.(nj.eq.nk).and.(nk.eq.nl) ) then
    ndiff=1
    return
  endif

  ! 2111
  if( (ni.ne.nj).and.(nj.eq.nk).and.(nk.eq.nl) ) then
    ndiff=2
    return
  endif

  ! 1211
  if( (ni.ne.nj).and.(nk.eq.ni).and.(nk.eq.nl) ) then
    ndiff=2
    return
  endif

  ! 1121
  if( (ni.eq.nj).and.(nj.ne.nk).and.(nl.eq.ni) ) then
    ndiff=2
    return
  endif

  ! 1112
  if( (ni.eq.nj).and.(nj.eq.nk).and.(nk.ne.nl) ) then
    ndiff=2
    return
  endif

  ! 1212
  if((ni.ne.nj).and.(nk.ne.nl).and.(ni.eq.nk).and.(nj.eq.nl) ) then
    ndiff=2
    return
  endif

  ! 2211
  if((ni.eq.nj).and.(nk.eq.nl).and.(ni.ne.nk) ) then
    ndiff=2
    return
  endif

  ! 1221
  if((ni.ne.nj).and.(nj.eq.nk).and.(nk.ne.nl).and.(ni.eq.nl) ) then
    ndiff=2
    return
  endif


  !11xy
  if( (ni.eq.nj).and.(ni.ne.nk).and.(ni.ne.nl).and.(nk.ne.nl) ) then
    ndiff=3
    return
  endif

  !1x1y
  if( (ni.eq.nk).and.(ni.ne.nj).and.(ni.ne.nl).and.(nj.ne.nl) ) then
    ndiff=3
    return
  endif

  !1xy1
  if( (ni.eq.nl).and.(ni.ne.nj).and.(ni.ne.nk).and.(nj.eq.nk) ) then
    ndiff=3
    return
  endif

  !x1y1
  if( (nj.eq.nl).and.(nj.ne.ni).and.(nj.ne.nk).and.(ni.ne.nk) ) then
    ndiff=3
    return
  endif

  !x11y
  if( (nj.eq.nk).and.(nj.ne.ni).and.(nj.ne.nl).and.(ni.ne.nl) ) then
    ndiff=3
    return
  endif

  !xy11
  if( (nk.eq.nl).and.(nk.ne.ni).and.(nk.ne.nj).and.(ni.ne.nj) ) then
    ndiff=3
    return
  endif

  !1234
  if( (ni.ne.nj).and.(nj.ne.nk).and.(nk.ne.nl).and.(nl.ne.ni).and.(nl.ne.nj).and.(ni.ne.nk) ) then
    ndiff=4
    return
  endif

end

