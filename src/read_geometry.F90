subroutine read_geometry(MOLECULE)
  include 'j.inc'
  character(80), intent(in) :: MOLECULE
  character(80) charabid,filename

  double precision, parameter :: conversion_ang_to_ua=1.d0/0.52917721092d0
  integer :: i, j, l
  double precision :: rij

  integer, external :: number_atom

  filename=trim(MOLECULE)//'.xyz'

  open(2,file=filename)
  read(2,*)number_nuclei
  if(number_nuclei.gt.number_nuclei_max)stop 'increase number_nuclei_max'
  read(2,'(a80)')charabid

  do i=1,number_nuclei
    read(2,*)ATOM(i),centers_nuclei(1,i),centers_nuclei(2,i),centers_nuclei(3,i)
    charge(i)=number_atom(ATOM(i))
    do l=1,3
      centers_nuclei(l,i)=conversion_ang_to_ua*centers_nuclei(l,i)
    enddo
  enddo

  !! Calculation of the nuclear energy'
  !*******************************
  enucl=0.d0
  do i=1,number_nuclei
    do j=i+1,number_nuclei

      rij=dsqrt( (centers_nuclei(1,i)-centers_nuclei(1,j))**2        &
          +(centers_nuclei(2,i)-centers_nuclei(2,j))**2              &
          +(centers_nuclei(3,i)-centers_nuclei(3,j))**2 )
      enucl=enucl+charge(i)*charge(j)/rij
    enddo
  enddo
end

