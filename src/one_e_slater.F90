program one_e_integrals
  use common_data

  character(80) :: charabid
  character(80) :: MOLECULE

  double precision :: S_ij_STO_ex(nbasis_max,nbasis_max)
  double precision :: V_ij_STO_ex(nbasis_max,nbasis_max)
  double precision :: K_ij_STO_ex(nbasis_max,nbasis_max)

  double precision :: Sij, Vij, Kij

  character(128) :: filename_in
  character(128) :: filename_out_s_ex
  character(128) :: filename_out_v_ex
  character(128) :: filename_out_k_ex
  character(128) :: filename_basis

  integer :: i, k

  write(filename_in,'(A4)') 'j_in'

  write(*,*)'INPUT FILE USED=',filename_in
  open(unit=5,file=filename_in)

  print*,'Simulation number?'
  read(5,'(a80)')charabid
  read(5,*)

  write(filename_out_s_ex,'(A13)')'overlap_ao_ex'
  write(filename_out_v_ex,'(A13)')'nuclear_ao_ex'
  write(filename_out_k_ex,'(A13)')'kinetic_ao_ex'
  !*****************
  ! BEGIN READ INPUT
  !*****************
  print*,'MOLECULE?'
  read(5,'(a80)')charabid
  read(5,*)MOLECULE
  write(*,*)MOLECULE
  read(5,'(a80)')charabid
  write(6,'(a80)')charabid
  print*,'file name for basis set?'
  read(5,*)filename_basis
  print*,'*******'
  write(*,*)trim(filename_basis)
  print*,'*******'
  read(5,'(a80)')charabid
  write(6,'(a80)')charabid
  read(5,*)ng0,i_add_2p,i_add_3d,g_thr,i_add_large_g
  write(6,*)ng0,i_add_2p,i_add_3d,g_thr,i_add_large_g

  !*****************
  ! END READ INPUT
  !*****************

  call read_basis(filename_basis)

  print*,'READING geometry in angstrom'
  call read_geometry(MOLECULE)
  print*,'ENUCL=',enucl

  ng0=20
  call build_gaussian_expansion_of_orbitals()

  write(*,*)'**********************************'
  write(*,*)'Number of basis functions ',nbasis

  do i=1,nbasis
    do k=1,nbasis
      call one_elect(i,k,Sij,Vij,Kij)
      S_ij_STO_ex(i,k)=Sij
      V_ij_STO_ex(i,k)=Vij
      K_ij_STO_ex(i,k)=Kij
    enddo
  enddo

  open(unit=104,file=filename_out_s_ex)
  open(unit=105,file=filename_out_v_ex)
  open(unit=106,file=filename_out_k_ex)
  do i=1,nbasis
    do k=1,nbasis
      write(104,'(2(I5,X),D22.15)') i,k,S_ij_STO_ex(i,k)
      write(105,'(2(I5,X),D22.15)') i,k,V_ij_STO_ex(i,k)
      write(106,'(2(I5,X),D22.15)') i,k,K_ij_STO_ex(i,k)
    enddo
  enddo
  close(104)
  close(105)
  close(106)

end

