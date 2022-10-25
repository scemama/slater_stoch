!! Build the gaussian expansion of each atomic orbitals
!!*****************************************************
!!
!! For i=1,nbasis
!!
!!  phi_i(r) centered at [ center(1,i),center(2,i),center(3,i) ]
!!  phi_i(r) = x^npower(1,i) y^npower(2,i) z^npower(3,i)  u_i(r)
!!
!!  u_i(r)= sum_{m=1,n_gauss(i)} [ c_gauss(m,i,1) + c_gauss(m,i,2) r^2 ] exp[-g_gauss(m,i) r^2]
!!
!!  i_type(i)=1  :     c_gauss(...,2)=0 and only first part is computed
!!  i_type(i)=2  :     c_gauss(...,1)=0 and only second part is computed
!!  i_type(i)=3  :     c_gauss(...,1) and c_gauss(...,2) computed


!! ROUTINE GAUSS_IJKL
!********************

! A_x=center(1,i)
! A_y=center(2,i)
! A_z=center(3,i)

! nA_x=npower(1,i)
! nA_y=npower(2,i)
! nA_z=npower(3,i)

! idem for B,C,D for orbitals i,k,j,l
!
! The following quantity is computed:

! gauss_ijkl_new = int d1 d2 (x1-A_x)**nA_x * (y1-A_y)**nA_y * (z1-A_z)**nA_z *  u_i(r1)
!                            (x1-B_x)**nB_x * (y1-B_y)**nB_y * (z1-B_z)**nB_z *  u_k(r1)
! 1/r_12
!                            (x2-C_x)**nC_x * (y2-C_y)**nC_y * (z2-C_z)**nC_z * u_j(r2)
!                            (x2-D_x)**nD_x * (y2-D_y)**nD_y * (z2-D_z)**nD_z * u_l(r2)
!
! with the radial parts given by
!
!  u_i(r1) = sum_{i1=1,ngauss(1)}  [ c0(i1,1)+ c2(i1,1)*(r1-A)**2] exp[-gg(i1,1)*(r1-A)**2]
!  u_k(r1) = sum_{i2=1,ngauss(2)}  [ c0(i2,1)+ c2(i2,1)*(r1-A)**2] exp[-gg(i2,1)*(r1-B)**2]
!
!  u_j(r2) = sum_{i3=1,ngauss(3)}  [ c0(i3,1)+ c2(i3,1)*(r2-C)**2] exp[-gg(i3,1)*(r2-C)**2]
!  u_l(r2) = sum_{i4=1,ngauss(4)}  [ c0(i4,1)+ c2(i4,1)*(r2-D)**2] exp[-gg(i4,1)*(r2-D)**2]

! To save computational cost:
!
! it=1 the r**2-component is not calculated (c2=0)
! it=2 only the r**2-component is calculated (c1=0)
! it=3 both r**0 and r**2-component are calculated

double precision function gauss_ijkl(i,k,j,l)
  use common_data
  integer, intent(in) :: i, j, k, l

  integer :: n_orb(4),nc(4),n_c(3,4,4)
  integer :: i_o, kk, m, ii, ll, i1, i2, i3, i4, ii1, ii2, ii3, ii4

  double precision :: d(n_gauss_max,4,4)
  double precision :: rint

  double precision, external :: bielec_integral

  n_orb(1)=i
  n_orb(2)=k
  n_orb(3)=j
  n_orb(4)=l

  do kk=1,4

    i_o=n_orb(kk)

    if(i_type(i_o).eq.1)then
      nc(kk)=1
      do m=1,n_gauss(i_o)
        do ii=1,nc(kk)
          d(m,ii,kk)=c_gauss(m,i_o,1)
        enddo
      enddo
      do ll=1,3
        do ii=1,nc(kk)
          n_c(ll,ii,kk)=npower(ll,i_o)
        enddo
      enddo
    endif

    if(i_type(i_o).eq.2)then
      nc(kk)=3
      do m=1,n_gauss(i_o)
        do ii=1,nc(kk)
          d(m,ii,kk)=c_gauss(m,i_o,2)
        enddo
      enddo
      do ll=1,3
        do ii=1,nc(kk)
          if(ll.eq.ii)then
            n_c(ll,ii,kk)=npower(ll,i_o)+2
          else
            n_c(ll,ii,kk)=npower(ll,i_o)
          endif
        enddo
      enddo
    endif

    if(i_type(i_o).eq.3)then
      nc(kk)=4
      do m=1,n_gauss(i_o)
        do ii=1,nc(kk)
          if(ii.eq.1)d(m,ii,kk)=c_gauss(m,i_o,1)
          if(ii.gt.1)d(m,ii,kk)=c_gauss(m,i_o,2)
        enddo
      enddo
      do ll=1,3
        do ii=1,nc(kk)
          if(ll.eq.(ii-1))then
            n_c(ll,ii,kk)=npower(ll,i_o)+2
          else
            n_c(ll,ii,kk)=npower(ll,i_o)
          endif
        enddo
      enddo
    endif

  enddo ! i_o=1,4

  rint=0.d0

  do i1=1,n_gauss(n_orb(1))
    do i2=1,n_gauss(n_orb(2))
      do i3=1,n_gauss(n_orb(3))
        do i4=1,n_gauss(n_orb(4))

          do ii1=1,nc(1)
            do ii2=1,nc(2)
              do ii3=1,nc(3)
                do ii4=1,nc(4)

                  rint=rint+                                         &
                      d(i1,ii1,1)*d(i2,ii2,2)*d(i3,ii3,3)*d(i4,ii4,4)*&
                      bielec_integral(g_gauss(i1,i),g_gauss(i2,k),g_gauss(i3,j),g_gauss(i4,l),&
                      center(1,i),center(1,k),center(1,j),center(1,l),&
                      n_c(1,ii1,1),n_c(1,ii2,2),n_c(1,ii3,3),n_c(1,ii4,4))

                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  gauss_ijkl=rint

end

