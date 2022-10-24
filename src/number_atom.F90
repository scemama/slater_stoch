integer function number_atom(ATOM)
  implicit none
  character(80), intent(in)      :: ATOM
  character(2), parameter        :: a(18) = (/                       &
      'H' ,                                   'He',                  &
      'Li', 'Be', 'B' , 'C' , 'N', 'O', 'F' , 'Ne',                  &
      'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'  /)
  integer :: i
  do i=1,size(a)
    if (ATOM == a(i)) then
      number_atom = i
      return
    end if
  end do
  stop 'ATOM not defined'
end

