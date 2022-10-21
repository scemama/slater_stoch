integer function number_atom(ATOM)
  character*80 ATOM
  character*80 HYDROGEN
  character*80 HELIUM
  character*80 LITHIUM
  character*80 BERYLLIUM
  character*80 BORON
  character*80 CARBON
  character*80 NITROGEN
  character*80 OXYGEN
  character*80 FLUORINE
  character*80 NEON
  character*80 SODIUM
  character*80 MAGNESIUM
  character*80 ALUMINUM
  character*80 SILICON
  character*80 PHOSPHORUS
  character*80 SULFUR
  character*80 CHLORINE
  character*80 ARGON

  if(ATOM.eq.'H')then
    number_atom=1
    return
  endif
  if(ATOM.eq.'He')then
    number_atom=2
    return
  endif
  if(ATOM.eq.'Li')then
    number_atom=3
    return
  endif
  if(ATOM.eq.'Be')then
    number_atom=4
    return
  endif
  if(ATOM.eq.'B')then
    number_atom=5
    return
  endif
  if(ATOM.eq.'C')then
    number_atom=6
    return
  endif
  if(ATOM.eq.'N')then
    number_atom=7
    return
  endif
  if(ATOM.eq.'O')then
    number_atom=8
    return
  endif
  if(ATOM.eq.'F')then
    number_atom=9
    return
  endif
  if(ATOM.eq.'Ne')then
    number_atom=10
    return
  endif
  if(ATOM.eq.'Na')then
    number_atom=11
    return
  endif
  if(ATOM.eq.'Mg')then
    number_atom=12
    return
  endif
  if(ATOM.eq.'Al')then
    number_atom=13
    return
  endif
  if(ATOM.eq.'Si')then
    number_atom=14
    return
  endif
  if(ATOM.eq.'P')then
    number_atom=15
    return
  endif
  if(ATOM.eq.'S')then
    number_atom=16
    return
  endif
  if(ATOM.eq.'Cl')then
    number_atom=17
    return
  endif
  if(ATOM.eq.'Ar')then
    i_atom=18
    return
  endif
  stop 'ATOM not defined'
end

