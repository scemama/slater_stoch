program one_e_integrals
  include 'j.inc'

  character*80,charabid
  character*80,MOLECULE

  logical normalization_OA

  ! Arrays for one-electron integrals
  double precision S_ij_STO_ex(nbasis_max,nbasis_max),               &
      V_ij_STO_ex(nbasis_max,nbasis_max),                            &
      K_ij_STO_ex(nbasis_max,nbasis_max)
  double precision sqrtSii_STO_ex(nbasis_max), Kij

  !! MONTE CARLO PART
  character*(128)                :: filename_in
  character*(128)                :: filename_out_s_ex
  character*(128)                :: filename_out_v_ex
  character*(128)                :: filename_out_k_ex
  character*(128)                :: filename_basis

  write(filename_in,'(A4)') 'j_in'

  write(*,*)'INPUT FILE USED=',filename_in
  open(unit=5,file=filename_in)

  print*,'Simulation number?'
  read(5,'(a80)')charabid
  read(5,*)
  !write(*,*)'num_simulation=',num_simulation

  write(filename_out_s_ex,'(A13)')'overlap_ao_ex'
  write(filename_out_v_ex,'(A13)')'nuclear_ao_ex'
  write(filename_out_k_ex,'(A13)')'kinetic_ao_ex'
  basis_type='STO'
  !*****************
  ! BEGIN READ INPUT
  !*****************
  print*,'MOLECULE?'
  read(5,'(a80)')charabid
  read(5,*)MOLECULE
  write(*,*)MOLECULE
  do i=1,3
    read(5,'(a80)')charabid
    write(6,'(a80)')charabid
  enddo
  read(5,*)
  write(*,*)'basis_type=',basis_type
  read(5,'(a80)')charabid
  write(6,'(a80)')charabid
  print*,'file name for basis set?'
  read(5,*)filename_basis
  print*,'*******'
  write(*,*)trim(filename_basis)
  print*,'*******'
  do i=1,3
    read(5,'(a80)')charabid
    write(6,'(a80)')charabid
  enddo
  read(5,*)
  read(5,'(a80)')charabid
  write(6,'(a80)')charabid
  read(5,*)
  !if(mpi_rank.eq.0)write(6,*)n_eps
  do i=1,5
    read(5,'(a80)')charabid
    write(6,'(a80)')charabid
  enddo
  read(5,*)ng0,i_add_2p,i_add_3d,g_thr,i_add_large_g
  write(6,*)ng0,i_add_2p,i_add_3d,g_thr,i_add_large_g
  do i=1,3
    read(5,'(a80)')charabid
    write(6,'(a80)')charabid
  enddo
  read(5,*)level
  do i=1,2
    read(5,'(a80)')charabid
    write(6,'(a80)')charabid
  enddo
  read(5,*)
  do i=1,2
    read(5,'(a80)')charabid
    write(6,'(a80)')charabid
  enddo
  read(5,*)
  !if(mpi_rank.eq.0)write(*,*)'Number of MC steps for one-electron integrals=',npts_one_elec, ' times number of blocks'
  read(5,'(a80)')charabid
  write(6,'(a80)')charabid
  read(5,*)
  do i=1,2
    read(5,'(a80)')charabid
    write(6,'(a80)')charabid
  enddo
  read(5,*)
  read(5,'(a80)')charabid
  read(5,*)
  !if(mpi_rank.eq.0)print*,'nocc=',nocc
  !if(mpi_rank.eq.0)print*,'Number iter for SCF calculations?'
  read(5,'(a80)')charabid
  read(5,*)
  !if(mpi_rank.eq.0)print*,'niter=',niter_SCF
  !! Determination of one-center bielectronic
  !if(mpi_rank.eq.0)print*,'mono_center_Slater?'
  read(5,'(a80)')charabid
  read(5,*)
  !if(mpi_rank.eq.0)print*,'compute momo-center Slater=',mono_center_Slater
  read(5,'(a80)')charabid
  read(5,*)normalization_OA
  print*,'orbital are normalized ?',normalization_OA

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
    sqrtSii_STO_ex(i)=dsqrt(dabs(S_ij_STO_ex(i,i)))
  enddo

  do i=1,nbasis
    do k=1,nbasis
      rnorm=sqrtSii_STO_ex(i)*sqrtSii_STO_ex(k)
      if(.not.normalization_OA)rnorm=1.d0
      S_ij_STO_ex(i,k)=S_ij_STO_ex(i,k)/rnorm
      V_ij_STO_ex(i,k)=V_ij_STO_ex(i,k)/rnorm
      K_ij_STO_ex(i,k)=K_ij_STO_ex(i,k)/rnorm
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

