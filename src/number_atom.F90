integer function number_atom(atom)
  implicit none
  character(2), intent(in)       :: atom
  character(2), parameter        :: a(36) = (/                         &
      'H ' ,                                    'He',                  &
      'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',                  &
      'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar',                  &
      'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr' /)
  integer :: i
  do i=1,size(a)
    if (atom == a(i)) then
      number_atom = i
      return
    end if
  end do
  print *, ':'//trim(atom)//':'
  stop 'Atom not defined'
end

