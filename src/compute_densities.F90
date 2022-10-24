subroutine compute_densities(i,k,ut1,ut2,rho,rho_G,poly)
  use common_data
  integer, intent(in) :: i, k
  double precision, intent(in)  :: ut1(3,nw,nbasis_max*nbasis_max)
  double precision, intent(in)  :: ut2(3,nw,nbasis_max*nbasis_max)
  double precision, intent(out) :: rho  (nw,nbasis_max*nbasis_max,2)
  double precision, intent(out) :: poly  (nw,nbasis_max*nbasis_max,2)
  double precision, intent(out) :: rho_G(nw,nbasis_max*nbasis_max,2)

  integer :: ik, kw, kk

  double precision :: dx, dy, dz, poly_i, poly_k, r_i, r_k, poly_ik, r_ik

  double precision, external :: u_orb, u_gauss

  ik = (k-1)*nbasis_max+i
  do kw=1,nw
    dx=ut1(1,kw,ik)-center(1,i)
    dy=ut1(2,kw,ik)-center(2,i)
    dz=ut1(3,kw,ik)-center(3,i)

    poly_i=1.d0
    do kk=1,npower(1,i)
      poly_i = poly_i * dx
    enddo
    do kk=1,npower(2,i)
      poly_i = poly_i * dy
    enddo
    do kk=1,npower(3,i)
      poly_i = poly_i * dz
    enddo
    r_i=dsqrt(dx*dx+dy*dy+dz*dz)

    dx=ut1(1,kw,ik)-center(1,k)
    dy=ut1(2,kw,ik)-center(2,k)
    dz=ut1(3,kw,ik)-center(3,k)

    poly_k=1.d0
    do kk=1,npower(1,k)
      poly_k = poly_k * dx
    enddo
    do kk=1,npower(2,k)
      poly_k = poly_k * dy
    enddo
    do kk=1,npower(3,k)
      poly_k = poly_k * dz
    enddo
    r_k=dsqrt(dx*dx+dy*dy+dz*dz)

    poly_ik = poly_i*poly_k
    rho  (kw,ik,1)=poly_ik*u_orb(i,r_i)*u_orb(k,r_k)
    poly (kw,ik,1)=poly_ik
    rho_G(kw,ik,1)=poly_ik*u_gauss(i,r_i)*u_gauss(k,r_k)


    dx=ut2(1,kw,ik)-center(1,i)
    dy=ut2(2,kw,ik)-center(2,i)
    dz=ut2(3,kw,ik)-center(3,i)

    poly_i=1.d0
    do kk=1,npower(1,i)
      poly_i = poly_i * dx
    enddo
    do kk=1,npower(2,i)
      poly_i = poly_i * dy
    enddo
    do kk=1,npower(3,i)
      poly_i = poly_i * dz
    enddo
    r_i=dsqrt(dx*dx+dy*dy+dz*dz)

    dx=ut2(1,kw,ik)-center(1,k)
    dy=ut2(2,kw,ik)-center(2,k)
    dz=ut2(3,kw,ik)-center(3,k)

    poly_k=1.d0
    do kk=1,npower(1,k)
      poly_k = poly_k * dx
    enddo
    do kk=1,npower(2,k)
      poly_k = poly_k * dy
    enddo
    do kk=1,npower(3,k)
      poly_k = poly_k * dz
    enddo
    r_k=dsqrt(dx*dx+dy*dy+dz*dz)

    poly_ik = poly_i*poly_k
    rho  (kw,ik,2)=poly_ik*u_orb(i,r_i)*u_orb(k,r_k)
    poly (kw,ik,2)=poly_ik
    rho_G(kw,ik,2)=poly_ik*u_gauss(i,r_i)*u_gauss(k,r_k)

  enddo

end
