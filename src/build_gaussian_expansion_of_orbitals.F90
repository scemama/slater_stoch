subroutine build_gaussian_expansion_of_orbitals()
  use common_data
  character(80) :: orb
  integer :: i, j, k, l, m, mm
  integer :: i_atom

  integer, external :: number_atom

  open(21,file='info_basis')
  open(22,file='orbital_coefficients_gaussian_expansion')
  open(23,file='info_simplified')


  call read_fit_SMILES

  nbasis=0

  do i=1,number_nuclei

    i_atom=number_atom(ATOM(i))

    do k=1,n_b(i_atom)

      orb=orb_b(k,i_atom)

      if(orb.eq.'1S')then
        nbasis=nbasis+1
        if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
        orb_name(nbasis)='1S'
        orb_name_full(nbasis)='1S'
        n_STO(nbasis)=0
        i_type(nbasis)=1
        do l=1,3
          npower(l,nbasis)=0
          center(l,nbasis)=centers_nuclei(l,i)
        enddo
        nucleus_number(nbasis)=i

        ng(1,nbasis)=ng0
        if(gamma_b(k,1,i_atom).gt.g_thr)then
          ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
        endif
        n_contract(nbasis)=n_cont_b(k,i_atom)
        do mm=1,n_contract(nbasis)
          g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
          c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
        enddo

        if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
        g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
        g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

        call build_c_g_gauss_STO(nbasis)
        call write_STO_in_file_info_basis(nbasis)

      endif  !1S

      if(orb.eq.'2S')then

        nbasis=nbasis+1
        if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
        orb_name(nbasis)='2S'
        orb_name_full(nbasis)='2S'
        n_STO(nbasis)=1
        i_type(nbasis)=2

        do l=1,3
          npower(l,nbasis)=0
          center(l,nbasis)=centers_nuclei(l,i)
        enddo
        nucleus_number(nbasis)=i

        ng(1,nbasis)=ng0
        if(gamma_b(k,1,i_atom).gt.g_thr)then
          ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
        endif
        n_contract(nbasis)=n_cont_b(k,i_atom)
        do mm=1,n_contract(nbasis)
          g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
          c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
        enddo

        if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
        g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
        g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

        call build_c_g_gauss_STO(nbasis)
        call write_STO_in_file_info_basis(nbasis)

      endif !2S

      if(orb.eq.'3S')then

        nbasis=nbasis+1
        if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
        orb_name(nbasis)='3S'
        orb_name_full(nbasis)='3S'
        n_STO(nbasis)=2
        i_type(nbasis)=2

        do l=1,3
          npower(l,nbasis)=0
          center(l,nbasis)=centers_nuclei(l,i)
        enddo
        nucleus_number(nbasis)=i

        ng(1,nbasis)=ng0
        if(gamma_b(k,1,i_atom).gt.g_thr)then
          ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
        endif
        n_contract(nbasis)=n_cont_b(k,i_atom)
        do mm=1,n_contract(nbasis)
          g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
          c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
        enddo

        if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
        g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
        g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

        call build_c_g_gauss_STO(nbasis)
        call write_STO_in_file_info_basis(nbasis)

      endif !3S

      if(orb.eq.'2P')then

        do m=1,3
          nbasis=nbasis+1
          if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
          n_STO(nbasis)=0
          i_type(nbasis)=1
          orb_name(nbasis)='2P'
          do l=1,3
            center(l,nbasis)=centers_nuclei(l,i)
            if(l.eq.m)then
              npower(l,nbasis)=1
            else
              npower(l,nbasis)=0
            endif
            nucleus_number(nbasis)=i
          enddo
          if(m.eq.1)orb_name_full(nbasis)='2P_X'
          if(m.eq.2)orb_name_full(nbasis)='2P_Y'
          if(m.eq.3)orb_name_full(nbasis)='2P_Z'

          ng(1,nbasis)=ng0+i_add_2p
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
          endif
          n_contract(nbasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(nbasis)
            g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
          g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
          g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

          call build_c_g_gauss_STO(nbasis)
          call write_STO_in_file_info_basis(nbasis)

        enddo ! m

      endif !2P

      if(orb.eq.'3P')then

        do m=1,3
          nbasis=nbasis+1
          if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
          if(m.eq.1)orb_name_full(nbasis)='3P_X'
          if(m.eq.2)orb_name_full(nbasis)='3P_Y'
          if(m.eq.3)orb_name_full(nbasis)='3P_Z'
          n_STO(nbasis)=1
          i_type(nbasis)=2
          orb_name(nbasis)='3P'
          do l=1,3
            center(l,nbasis)=centers_nuclei(l,i)
            if(l.eq.m)then
              npower(l,nbasis)=1
            else
              npower(l,nbasis)=0
            endif
          enddo
          nucleus_number(nbasis)=i

          ng(1,nbasis)=ng0+i_add_2p
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
          endif
          n_contract(nbasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(nbasis)
            g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
          g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
          g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

          call build_c_g_gauss_STO(nbasis)
          call write_STO_in_file_info_basis(nbasis)

        enddo ! m

      endif !3P

      if(orb.eq.'3D')then

        do m=1,6
          nbasis=nbasis+1
          if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
          n_STO(nbasis)=0
          i_type(nbasis)=1
          orb_name(nbasis)='3D'
          do l=1,3
            center(l,nbasis)=centers_nuclei(l,i)
          enddo
          if(m.eq.1)then
            npower(1,nbasis)=2
            npower(2,nbasis)=0
            npower(3,nbasis)=0
            orb_name_full(nbasis)='3D_XX'
          endif
          if(m.eq.2)then
            npower(1,nbasis)=1
            npower(2,nbasis)=1
            npower(3,nbasis)=0
            orb_name_full(nbasis)='3D_XY'
          endif
          if(m.eq.3)then
            npower(1,nbasis)=1
            npower(2,nbasis)=0
            npower(3,nbasis)=1
            orb_name_full(nbasis)='3D_XZ'
          endif
          if(m.eq.4)then
            npower(1,nbasis)=0
            npower(2,nbasis)=2
            npower(3,nbasis)=0
            orb_name_full(nbasis)='3D_YY'
          endif
          if(m.eq.5)then
            npower(1,nbasis)=0
            npower(2,nbasis)=1
            npower(3,nbasis)=1
            orb_name_full(nbasis)='3D_YZ'
          endif
          if(m.eq.6)then
            npower(1,nbasis)=0
            npower(2,nbasis)=0
            npower(3,nbasis)=2
            orb_name_full(nbasis)='3D_ZZ'
          endif
          nucleus_number(nbasis)=i

          ng(1,nbasis)=ng0+i_add_3d
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
          endif
          n_contract(nbasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(nbasis)
            g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
          g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
          g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

          call build_c_g_gauss_STO(nbasis)
          call write_STO_in_file_info_basis(nbasis)

        enddo ! m

      endif  !3D

      if(orb.eq.'4F')then

        do m=1,10
          nbasis=nbasis+1
          if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
          n_STO(nbasis)=0
          i_type(nbasis)=1
          orb_name(nbasis)='4F'
          do l=1,3
            center(l,nbasis)=centers_nuclei(l,i)
          enddo
          if(m.eq.1)then
            npower(1,nbasis)=3
            npower(2,nbasis)=0
            npower(3,nbasis)=0
            orb_name_full(nbasis)='4F_XXX'
          endif
          if(m.eq.2)then
            npower(1,nbasis)=2
            npower(2,nbasis)=1
            npower(3,nbasis)=0
            orb_name_full(nbasis)='4F_XXY'
          endif
          if(m.eq.3)then
            npower(1,nbasis)=2
            npower(2,nbasis)=0
            npower(3,nbasis)=1
            orb_name_full(nbasis)='4F_XXZ'
          endif
          if(m.eq.4)then
            npower(1,nbasis)=1
            npower(2,nbasis)=2
            npower(3,nbasis)=0
            orb_name_full(nbasis)='4F_XYY'
          endif
          if(m.eq.5)then
            npower(1,nbasis)=1
            npower(2,nbasis)=1
            npower(3,nbasis)=1
            orb_name_full(nbasis)='4F_XYZ'
          endif
          if(m.eq.6)then
            npower(1,nbasis)=1
            npower(2,nbasis)=0
            npower(3,nbasis)=2
            orb_name_full(nbasis)='4F_XZZ'
          endif
          if(m.eq.7)then
            npower(1,nbasis)=0
            npower(2,nbasis)=3
            npower(3,nbasis)=0
            orb_name_full(nbasis)='4F_YYY'
          endif
          if(m.eq.8)then
            npower(1,nbasis)=0
            npower(2,nbasis)=2
            npower(3,nbasis)=1
            orb_name_full(nbasis)='4F_YYZ'
          endif
          if(m.eq.9)then
            npower(1,nbasis)=0
            npower(2,nbasis)=1
            npower(3,nbasis)=2
            orb_name_full(nbasis)='4F_YZZ'
          endif
          if(m.eq.10)then
            npower(1,nbasis)=0
            npower(2,nbasis)=0
            npower(3,nbasis)=3
            orb_name_full(nbasis)='4F_ZZZ'
          endif
          nucleus_number(nbasis)=i

          ng(1,nbasis)=ng0+i_add_3d
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
          endif
          n_contract(nbasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(nbasis)
            g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
          g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
          g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

          call build_c_g_gauss_STO(nbasis)
          call write_STO_in_file_info_basis(nbasis)

        enddo ! m

      endif  !4F

      if(orb.eq.'5G')then

        do m=1,15
          nbasis=nbasis+1
          if(nbasis.gt.nbasis_max)stop 'increase nbasis_max'
          n_STO(nbasis)=0
          i_type(nbasis)=1
          orb_name(nbasis)='5G'
          do l=1,3
            center(l,nbasis)=centers_nuclei(l,i)
          enddo
          if(m.eq.1)then
            npower(1,nbasis)=4
            npower(2,nbasis)=0
            npower(3,nbasis)=0
            orb_name_full(nbasis)='5G_XXXX'
          endif
          if(m.eq.2)then
            npower(1,nbasis)=3
            npower(2,nbasis)=1
            npower(3,nbasis)=0
            orb_name_full(nbasis)='5G_XXXY'
          endif
          if(m.eq.3)then
            npower(1,nbasis)=2
            npower(2,nbasis)=2
            npower(3,nbasis)=0
            orb_name_full(nbasis)='5G_XXYY'
          endif
          if(m.eq.4)then
            npower(1,nbasis)=1
            npower(2,nbasis)=3
            npower(3,nbasis)=0
            orb_name_full(nbasis)='5G_XYYY'
          endif
          if(m.eq.5)then
            npower(1,nbasis)=0
            npower(2,nbasis)=4
            npower(3,nbasis)=0
            orb_name_full(nbasis)='5G_YYYY'
          endif
          if(m.eq.6)then
            npower(1,nbasis)=3
            npower(2,nbasis)=0
            npower(3,nbasis)=1
            orb_name_full(nbasis)='5G_XXXZ'
          endif
          if(m.eq.7)then
            npower(1,nbasis)=2
            npower(2,nbasis)=1
            npower(3,nbasis)=1
            orb_name_full(nbasis)='5G_XXYZ'
          endif
          if(m.eq.8)then
            npower(1,nbasis)=1
            npower(2,nbasis)=2
            npower(3,nbasis)=1
            orb_name_full(nbasis)='5G_XYYZ'
          endif
          if(m.eq.9)then
            npower(1,nbasis)=0
            npower(2,nbasis)=3
            npower(3,nbasis)=1
            orb_name_full(nbasis)='5G_YYYZ'
          endif
          if(m.eq.10)then
            npower(1,nbasis)=2
            npower(2,nbasis)=0
            npower(3,nbasis)=2
            orb_name_full(nbasis)='5G_XXZZ'
          endif
          if(m.eq.11)then
            npower(1,nbasis)=1
            npower(2,nbasis)=1
            npower(3,nbasis)=2
            orb_name_full(nbasis)='5G_XYZZ'
          endif
          if(m.eq.12)then
            npower(1,nbasis)=0
            npower(2,nbasis)=2
            npower(3,nbasis)=2
            orb_name_full(nbasis)='5G_YYZZ'
          endif
          if(m.eq.13)then
            npower(1,nbasis)=1
            npower(2,nbasis)=0
            npower(3,nbasis)=3
            orb_name_full(nbasis)='5G_XZZZ'
          endif
          if(m.eq.14)then
            npower(1,nbasis)=0
            npower(2,nbasis)=1
            npower(3,nbasis)=3
            orb_name_full(nbasis)='5G_YZZZ'
          endif
          if(m.eq.15)then
            npower(1,nbasis)=0
            npower(2,nbasis)=0
            npower(3,nbasis)=4
            orb_name_full(nbasis)='5G_ZZZZ'
          endif
          nucleus_number(nbasis)=i

          ng(1,nbasis)=ng0+i_add_3d
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,nbasis)=ng(1,nbasis)+i_add_large_g
          endif
          n_contract(nbasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(nbasis)
            g_contract(mm,nbasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,nbasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(nbasis).gt.1)stop 'no contraction for STO'
          g_min(nbasis)=g_contract(n_contract(nbasis),nbasis)
          g_slater(nbasis)=g_contract(n_contract(nbasis),nbasis)

          call build_c_g_gauss_STO(nbasis)
          call write_STO_in_file_info_basis(nbasis)

        enddo !m

      endif  !5G

    enddo ! k=1,n_b(i_atom)
  enddo ! i=1,number_nuclei


  do i=1,nbasis
    do j=1,nbasis
      dist_ij(i,j)=                                                  &
          dsqrt((center(1,i)-center(1,j))**2+(center(2,i)-center(2,j))**2+(center(3,i)-center(3,j))**2)
    enddo
  enddo

end
