subroutine one_elect(i,k,Sij,Vij,Kij)
  use common_data

  integer, intent(in)            :: i, k
  double precision, intent(out)  :: Sij, Vij, Kij

  integer :: n_orb(4),nc(4),n_c(3,4,4)
  double precision :: d(n_gauss_max,4,4)
  integer :: ii, ii1, ii2, kk, ll, i_o, m, i1, i2

  double precision, external :: one_electron_Sij
  double precision, external :: one_electron_Vij
  double precision, external :: one_electron_Kij

  n_orb(1)=i
  n_orb(2)=k

  do kk=1,2

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

  enddo ! i_o=1,2

  Sij=0.d0
  Vij=0.d0
  Kij=0.d0

  do i1=1,n_gauss(n_orb(1))
    do i2=1,n_gauss(n_orb(2))
      do ii1=1,nc(1)
        do ii2=1,nc(2)

          !    write(37,*)'in one-lect=',d(i1,ii1,1)*d(i2,ii2,2)

          Sij=Sij+                                                   &
              d(i1,ii1,1)*d(i2,ii2,2)*                               &
              one_electron_Sij(n_c(1,ii1,1),center(1,i),g_gauss(i1,i),n_c(1,ii2,2),center(1,k),g_gauss(i2,k))

          Vij=Vij+                                                   &
              d(i1,ii1,1)*d(i2,ii2,2)*                               &
              one_electron_Vij(n_c(1,ii1,1),center(1,i),g_gauss(i1,i),n_c(1,ii2,2),center(1,k),g_gauss(i2,k))

          Kij=Kij+                                                   &
              d(i1,ii1,1)*d(i2,ii2,2)*                               &
              one_electron_Kij(n_c(1,ii1,1),center(1,i),g_gauss(i1,i),n_c(1,ii2,2),center(1,k),g_gauss(i2,k))

        enddo
      enddo
    enddo
  enddo

end

