subroutine build_gaussian_expansion_of_orbitals()
  use common_data
  character(80) :: orb
  integer :: i, j, k, l, m, mm
  integer :: i_atom, ibasis

  integer, external :: number_atom

  open(21,file='info_basis')
  open(22,file='orbital_coefficients_gaussian_expansion')
  open(23,file='info_simplified')

  call read_fit_SMILES

  ibasis = 0
  do i=1,number_nuclei

    i_atom=number_atom(ATOM(i))

    do k=1,n_b(i_atom)

      orb=orb_b(k,i_atom)

      if(orb.eq.'1S')then
        ibasis=ibasis+1
        orb_name(ibasis)='1S'
        orb_name_full(ibasis)='1S'
        n_STO(ibasis)=0
        i_type(ibasis)=1
        do l=1,3
          npower(l,ibasis)=0
          center(l,ibasis)=centers_nuclei(l,i)
        enddo
        nucleus_number(ibasis)=i

        ng(1,ibasis)=ng0
        if(gamma_b(k,1,i_atom).gt.g_thr)then
          ng(1,ibasis)=ng(1,ibasis)+i_add_large_g
        endif
        n_contract(ibasis)=n_cont_b(k,i_atom)
        do mm=1,n_contract(ibasis)
          g_contract(mm,ibasis)=gamma_b(k,mm,i_atom)
          c_contract(mm,ibasis)=coef_b(k,mm,i_atom)
        enddo

        if(n_contract(ibasis).gt.1)stop 'no contraction for STO'
        g_slater(ibasis)=g_contract(n_contract(ibasis),ibasis)

        call build_c_g_gauss_STO(ibasis)
        call write_STO_in_file_info_basis(ibasis)

      endif  !1S

      if(orb.eq.'2S')then

        ibasis=ibasis+1
        orb_name(ibasis)='2S'
        orb_name_full(ibasis)='2S'
        n_STO(ibasis)=1
        i_type(ibasis)=2

        do l=1,3
          npower(l,ibasis)=0
          center(l,ibasis)=centers_nuclei(l,i)
        enddo
        nucleus_number(ibasis)=i

        ng(1,ibasis)=ng0
        if(gamma_b(k,1,i_atom).gt.g_thr)then
          ng(1,ibasis)=ng(1,ibasis)+i_add_large_g
        endif
        n_contract(ibasis)=n_cont_b(k,i_atom)
        do mm=1,n_contract(ibasis)
          g_contract(mm,ibasis)=gamma_b(k,mm,i_atom)
          c_contract(mm,ibasis)=coef_b(k,mm,i_atom)
        enddo

        if(n_contract(ibasis).gt.1)stop 'no contraction for STO'
        g_slater(ibasis)=g_contract(n_contract(ibasis),ibasis)

        call build_c_g_gauss_STO(ibasis)
        call write_STO_in_file_info_basis(ibasis)

      endif !2S

      if(orb.eq.'3S')then

        ibasis=ibasis+1
        orb_name(ibasis)='3S'
        orb_name_full(ibasis)='3S'
        n_STO(ibasis)=2
        i_type(ibasis)=2

        do l=1,3
          npower(l,ibasis)=0
          center(l,ibasis)=centers_nuclei(l,i)
        enddo
        nucleus_number(ibasis)=i

        ng(1,ibasis)=ng0
        if(gamma_b(k,1,i_atom).gt.g_thr)then
          ng(1,ibasis)=ng(1,ibasis)+i_add_large_g
        endif
        n_contract(ibasis)=n_cont_b(k,i_atom)
        do mm=1,n_contract(ibasis)
          g_contract(mm,ibasis)=gamma_b(k,mm,i_atom)
          c_contract(mm,ibasis)=coef_b(k,mm,i_atom)
        enddo

        if(n_contract(ibasis).gt.1)stop 'no contraction for STO'
        g_slater(ibasis)=g_contract(n_contract(ibasis),ibasis)

        call build_c_g_gauss_STO(ibasis)
        call write_STO_in_file_info_basis(ibasis)

      endif !3S

      if(orb.eq.'2P')then

        do m=1,3
          ibasis=ibasis+1
          n_STO(ibasis)=0
          i_type(ibasis)=1
          orb_name(ibasis)='2P'
          do l=1,3
            center(l,ibasis)=centers_nuclei(l,i)
            if(l.eq.m)then
              npower(l,ibasis)=1
            else
              npower(l,ibasis)=0
            endif
            nucleus_number(ibasis)=i
          enddo
          if(m.eq.1)orb_name_full(ibasis)='2P_X'
          if(m.eq.2)orb_name_full(ibasis)='2P_Y'
          if(m.eq.3)orb_name_full(ibasis)='2P_Z'

          ng(1,ibasis)=ng0+i_add_2p
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,ibasis)=ng(1,ibasis)+i_add_large_g
          endif
          n_contract(ibasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(ibasis)
            g_contract(mm,ibasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,ibasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(ibasis).gt.1)stop 'no contraction for STO'
          g_slater(ibasis)=g_contract(n_contract(ibasis),ibasis)

          call build_c_g_gauss_STO(ibasis)
          call write_STO_in_file_info_basis(ibasis)

        enddo ! m

      endif !2P

      if(orb.eq.'3P')then

        do m=1,3
          ibasis=ibasis+1
          if(m.eq.1)orb_name_full(ibasis)='3P_X'
          if(m.eq.2)orb_name_full(ibasis)='3P_Y'
          if(m.eq.3)orb_name_full(ibasis)='3P_Z'
          n_STO(ibasis)=1
          i_type(ibasis)=2
          orb_name(ibasis)='3P'
          do l=1,3
            center(l,ibasis)=centers_nuclei(l,i)
            if(l.eq.m)then
              npower(l,ibasis)=1
            else
              npower(l,ibasis)=0
            endif
          enddo
          nucleus_number(ibasis)=i

          ng(1,ibasis)=ng0+i_add_2p
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,ibasis)=ng(1,ibasis)+i_add_large_g
          endif
          n_contract(ibasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(ibasis)
            g_contract(mm,ibasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,ibasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(ibasis).gt.1)stop 'no contraction for STO'
          g_slater(ibasis)=g_contract(n_contract(ibasis),ibasis)

          call build_c_g_gauss_STO(ibasis)
          call write_STO_in_file_info_basis(ibasis)

        enddo ! m

      endif !3P

      if(orb.eq.'3D')then

        do m=1,6
          ibasis=ibasis+1
          n_STO(ibasis)=0
          i_type(ibasis)=1
          orb_name(ibasis)='3D'
          do l=1,3
            center(l,ibasis)=centers_nuclei(l,i)
          enddo
          if(m.eq.1)then
            npower(1,ibasis)=2
            npower(2,ibasis)=0
            npower(3,ibasis)=0
            orb_name_full(ibasis)='3D_XX'
          endif
          if(m.eq.2)then
            npower(1,ibasis)=1
            npower(2,ibasis)=1
            npower(3,ibasis)=0
            orb_name_full(ibasis)='3D_XY'
          endif
          if(m.eq.3)then
            npower(1,ibasis)=1
            npower(2,ibasis)=0
            npower(3,ibasis)=1
            orb_name_full(ibasis)='3D_XZ'
          endif
          if(m.eq.4)then
            npower(1,ibasis)=0
            npower(2,ibasis)=2
            npower(3,ibasis)=0
            orb_name_full(ibasis)='3D_YY'
          endif
          if(m.eq.5)then
            npower(1,ibasis)=0
            npower(2,ibasis)=1
            npower(3,ibasis)=1
            orb_name_full(ibasis)='3D_YZ'
          endif
          if(m.eq.6)then
            npower(1,ibasis)=0
            npower(2,ibasis)=0
            npower(3,ibasis)=2
            orb_name_full(ibasis)='3D_ZZ'
          endif
          nucleus_number(ibasis)=i

          ng(1,ibasis)=ng0+i_add_3d
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,ibasis)=ng(1,ibasis)+i_add_large_g
          endif
          n_contract(ibasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(ibasis)
            g_contract(mm,ibasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,ibasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(ibasis).gt.1)stop 'no contraction for STO'
          g_slater(ibasis)=g_contract(n_contract(ibasis),ibasis)

          call build_c_g_gauss_STO(ibasis)
          call write_STO_in_file_info_basis(ibasis)

        enddo ! m

      endif  !3D

      if(orb.eq.'4F')then

        do m=1,10
          ibasis=ibasis+1
          n_STO(ibasis)=0
          i_type(ibasis)=1
          orb_name(ibasis)='4F'
          do l=1,3
            center(l,ibasis)=centers_nuclei(l,i)
          enddo
          if(m.eq.1)then
            npower(1,ibasis)=3
            npower(2,ibasis)=0
            npower(3,ibasis)=0
            orb_name_full(ibasis)='4F_XXX'
          endif
          if(m.eq.2)then
            npower(1,ibasis)=2
            npower(2,ibasis)=1
            npower(3,ibasis)=0
            orb_name_full(ibasis)='4F_XXY'
          endif
          if(m.eq.3)then
            npower(1,ibasis)=2
            npower(2,ibasis)=0
            npower(3,ibasis)=1
            orb_name_full(ibasis)='4F_XXZ'
          endif
          if(m.eq.4)then
            npower(1,ibasis)=1
            npower(2,ibasis)=2
            npower(3,ibasis)=0
            orb_name_full(ibasis)='4F_XYY'
          endif
          if(m.eq.5)then
            npower(1,ibasis)=1
            npower(2,ibasis)=1
            npower(3,ibasis)=1
            orb_name_full(ibasis)='4F_XYZ'
          endif
          if(m.eq.6)then
            npower(1,ibasis)=1
            npower(2,ibasis)=0
            npower(3,ibasis)=2
            orb_name_full(ibasis)='4F_XZZ'
          endif
          if(m.eq.7)then
            npower(1,ibasis)=0
            npower(2,ibasis)=3
            npower(3,ibasis)=0
            orb_name_full(ibasis)='4F_YYY'
          endif
          if(m.eq.8)then
            npower(1,ibasis)=0
            npower(2,ibasis)=2
            npower(3,ibasis)=1
            orb_name_full(ibasis)='4F_YYZ'
          endif
          if(m.eq.9)then
            npower(1,ibasis)=0
            npower(2,ibasis)=1
            npower(3,ibasis)=2
            orb_name_full(ibasis)='4F_YZZ'
          endif
          if(m.eq.10)then
            npower(1,ibasis)=0
            npower(2,ibasis)=0
            npower(3,ibasis)=3
            orb_name_full(ibasis)='4F_ZZZ'
          endif
          nucleus_number(ibasis)=i

          ng(1,ibasis)=ng0+i_add_3d
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,ibasis)=ng(1,ibasis)+i_add_large_g
          endif
          n_contract(ibasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(ibasis)
            g_contract(mm,ibasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,ibasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(ibasis).gt.1)stop 'no contraction for STO'
          g_slater(ibasis)=g_contract(n_contract(ibasis),ibasis)

          call build_c_g_gauss_STO(ibasis)
          call write_STO_in_file_info_basis(ibasis)

        enddo ! m

      endif  !4F

      if(orb.eq.'5G')then

        do m=1,15
          ibasis=ibasis+1
          n_STO(ibasis)=0
          i_type(ibasis)=1
          orb_name(ibasis)='5G'
          do l=1,3
            center(l,ibasis)=centers_nuclei(l,i)
          enddo
          if(m.eq.1)then
            npower(1,ibasis)=4
            npower(2,ibasis)=0
            npower(3,ibasis)=0
            orb_name_full(ibasis)='5G_XXXX'
          endif
          if(m.eq.2)then
            npower(1,ibasis)=3
            npower(2,ibasis)=1
            npower(3,ibasis)=0
            orb_name_full(ibasis)='5G_XXXY'
          endif
          if(m.eq.3)then
            npower(1,ibasis)=2
            npower(2,ibasis)=2
            npower(3,ibasis)=0
            orb_name_full(ibasis)='5G_XXYY'
          endif
          if(m.eq.4)then
            npower(1,ibasis)=1
            npower(2,ibasis)=3
            npower(3,ibasis)=0
            orb_name_full(ibasis)='5G_XYYY'
          endif
          if(m.eq.5)then
            npower(1,ibasis)=0
            npower(2,ibasis)=4
            npower(3,ibasis)=0
            orb_name_full(ibasis)='5G_YYYY'
          endif
          if(m.eq.6)then
            npower(1,ibasis)=3
            npower(2,ibasis)=0
            npower(3,ibasis)=1
            orb_name_full(ibasis)='5G_XXXZ'
          endif
          if(m.eq.7)then
            npower(1,ibasis)=2
            npower(2,ibasis)=1
            npower(3,ibasis)=1
            orb_name_full(ibasis)='5G_XXYZ'
          endif
          if(m.eq.8)then
            npower(1,ibasis)=1
            npower(2,ibasis)=2
            npower(3,ibasis)=1
            orb_name_full(ibasis)='5G_XYYZ'
          endif
          if(m.eq.9)then
            npower(1,ibasis)=0
            npower(2,ibasis)=3
            npower(3,ibasis)=1
            orb_name_full(ibasis)='5G_YYYZ'
          endif
          if(m.eq.10)then
            npower(1,ibasis)=2
            npower(2,ibasis)=0
            npower(3,ibasis)=2
            orb_name_full(ibasis)='5G_XXZZ'
          endif
          if(m.eq.11)then
            npower(1,ibasis)=1
            npower(2,ibasis)=1
            npower(3,ibasis)=2
            orb_name_full(ibasis)='5G_XYZZ'
          endif
          if(m.eq.12)then
            npower(1,ibasis)=0
            npower(2,ibasis)=2
            npower(3,ibasis)=2
            orb_name_full(ibasis)='5G_YYZZ'
          endif
          if(m.eq.13)then
            npower(1,ibasis)=1
            npower(2,ibasis)=0
            npower(3,ibasis)=3
            orb_name_full(ibasis)='5G_XZZZ'
          endif
          if(m.eq.14)then
            npower(1,ibasis)=0
            npower(2,ibasis)=1
            npower(3,ibasis)=3
            orb_name_full(ibasis)='5G_YZZZ'
          endif
          if(m.eq.15)then
            npower(1,ibasis)=0
            npower(2,ibasis)=0
            npower(3,ibasis)=4
            orb_name_full(ibasis)='5G_ZZZZ'
          endif
          nucleus_number(ibasis)=i

          ng(1,ibasis)=ng0+i_add_3d
          if(gamma_b(k,1,i_atom).gt.g_thr)then
            ng(1,ibasis)=ng(1,ibasis)+i_add_large_g
          endif
          n_contract(ibasis)=n_cont_b(k,i_atom)
          do mm=1,n_contract(ibasis)
            g_contract(mm,ibasis)=gamma_b(k,mm,i_atom)
            c_contract(mm,ibasis)=coef_b(k,mm,i_atom)
          enddo

          if(n_contract(ibasis).gt.1)stop 'no contraction for STO'
          g_slater(ibasis)=g_contract(n_contract(ibasis),ibasis)

          call build_c_g_gauss_STO(ibasis)
          call write_STO_in_file_info_basis(ibasis)

        enddo !m

      endif  !5G

    enddo ! k=1,n_b(i_atom)
  enddo ! i=1,number_nuclei

  nbasis = ibasis

end
