program one_e_integrals
  use common_data

#ifdef HAVE_TREXIO
  use trexio
#endif

  implicit none

  character(80) :: charabid
  character(80) :: MOLECULE

  double precision, allocatable :: S_ij_STO_ex(:,:)
  double precision, allocatable :: V_ij_STO_ex(:,:)
  double precision, allocatable :: K_ij_STO_ex(:,:)

  double precision :: Sij, Vij, Kij

  character(128) :: filename_in
  character(128) :: filename_out_s_ex
  character(128) :: filename_out_v_ex
  character(128) :: filename_out_k_ex
  character(128) :: filename_basis

  integer :: i, k

#ifdef HAVE_TREXIO
  character*(128)   :: trexio_filename
  integer           :: rc
  integer(trexio_t) :: trexio_file
  character(128)    :: err_message
#endif


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

  print*,'READING geometry in angstrom'
  call read_geometry(MOLECULE)
  print*,'ENUCL=',enucl

  call read_basis(filename_basis)

  ng0=20
  call build_gaussian_expansion_of_orbitals()

  write(*,*)'**********************************'
  write(*,*)'Number of basis functions ',nbasis

  allocate(S_ij_STO_ex(nbasis,nbasis))
  allocate(V_ij_STO_ex(nbasis,nbasis))
  allocate(K_ij_STO_ex(nbasis,nbasis))

  do i=1,nbasis
    do k=1,nbasis
      call one_elect(i,k,Sij,Vij,Kij)
      S_ij_STO_ex(i,k)=Sij
      V_ij_STO_ex(i,k)=Vij
      K_ij_STO_ex(i,k)=Kij
    enddo
  enddo

#ifdef HAVE_TREXIO

  trexio_filename = trim(MOLECULE)//'.h5'
  trexio_file = trexio_open(trexio_filename, 'w', TREXIO_HDF5, rc)
  call trexio_assert(rc, TREXIO_SUCCESS)

  call trexio_write_geometry(trexio_file)
  call trexio_write_basis(trexio_file)

  ! Write Integrals
  ! ---------------




  rc = trexio_write_ao_1e_int_kinetic(trexio_file, K_ij_STO_ex)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_ao_1e_int_potential_n_e(trexio_file, V_ij_STO_ex)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_ao_1e_int_overlap(trexio_file, S_ij_STO_ex)
  call trexio_assert(rc, TREXIO_SUCCESS)

#else

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

#endif
end

