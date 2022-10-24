double precision function ijkl_slater_xnynzn(nix,niy,niz,power,gamma)
  implicit none
  
  integer nix(4),niy(4),niz(4)
  integer power(4)
  double precision gamma(4)
  
  integer max_ni_max
  parameter(max_ni_max=2)
  
  integer nslat(4),lslat(4),mslat(4)
  integer lslat_max(4)
  double precision sij,tij,vij,vijkl
  double precision overlap_slater,kinetic_slater,nuclear_slater,ijkl_slater
  
  integer nx,ny,nz,l,m,max_ni,k,ntot
  integer m1,m2,m3,m4,l1,l2,l3,l4
  double complex imag
  double complex c(0:max_ni_max,0:max_ni_max,0:max_ni_max,0:3*max_ni_max,-3*max_ni_max:3*max_ni_max)
  double complex rint,coef
  double precision norm,norm_sto
  double precision pi,rac4pi,rac4pi_3,racpi_15,c0
  
  !sij=overlap_slater(nslat,lslat,mslat,gamma)
  !tij=kinetic_slater(nslat,lslat,mslat,gamma)
  !vij=nuclear_slater(nslat,lslat,mslat,gamma)
  !write(*,*)sij,tij,vij
  !vijkl=ijkl_slater(nslat,lslat,mslat,gamma)/sij**2
  !write(*,*)vijkl
  
  max_ni=2
  
  do k=1,4
    ntot=nix(k)+niy(k)+niz(k)
    if(ntot.gt.max_ni)stop ' x^nx x^nx x^nx Slater no yet implemented'
  enddo
  
  pi=dacos(-1.d0)
  rac4pi=dsqrt(4.d0*pi)
  rac4pi_3=dsqrt(4.d0*pi/3.d0)
  racpi_15=dsqrt(pi/15.d0)
  
  imag=(0.d0,1.d0)
  
  do nx=0,max_ni
    do ny=0,max_ni
      do nz=0,max_ni
        do l=0,nx+ny+nz
          do m=-l,l
            c(nx,ny,nz,l,m)=0.d0
          enddo
        enddo
      enddo
    enddo
  enddo
  
  !! Decomposition of x^nx y^ny z^nz over spherical harmonics (complex)
  !!
  !!******************************************************************
  !! x^nx y^ny z^nz = \sum_{l=0,lmax} {m=-l,l}  c(nx,ny,nz,l,m) Y(l,m)
  !! lmax=nx+ny+nz
  !!*****************************************************************
  
  ! 1
  c(0,0,0,0,0) = rac4pi
  
  !x
  c(1,0,0,1,-1)=  rac4pi_3*1.d0/dsqrt(2.d0)
  c(1,0,0,1, 1)= -rac4pi_3*1.d0/dsqrt(2.d0)
  
  !y
  c(0,1,0,1,-1)=  imag*rac4pi_3*1.d0/dsqrt(2.d0)
  c(0,1,0,1, 1)=  imag*rac4pi_3*1.d0/dsqrt(2.d0)
  
  !z
  c(0,0,1,1, 0)=  rac4pi_3
  
  !xy
  c(1,1,0,2,-2)=  imag*racpi_15*dsqrt(2.d0)
  c(1,1,0,2, 2)= -imag*racpi_15*dsqrt(2.d0)
  
  !yz
  c(0,1,1,2,-1)=  imag*racpi_15*dsqrt(2.d0)
  c(0,1,1,2, 1)=  imag*racpi_15*dsqrt(2.d0)
  
  !xz
  c(1,0,1,2,-1)=  racpi_15*dsqrt(2.d0)
  c(1,0,1,2, 1)= -racpi_15*dsqrt(2.d0)
  
  !!x2
  c(2,0,0,0, 0)=  rac4pi/3.d0
  c(2,0,0,2,-2)=  racpi_15*dsqrt(2.d0)
  c(2,0,0,2, 2)=  racpi_15*dsqrt(2.d0)
  c(2,0,0,2, 0)= -racpi_15*2.d0/dsqrt(3.d0)
  
  !!y2
  c(0,2,0,0, 0)=  rac4pi/3.d0
  c(0,2,0,2,-2)= -racpi_15*dsqrt(2.d0)
  c(0,2,0,2, 2)= -racpi_15*dsqrt(2.d0)
  c(0,2,0,2, 0)= -racpi_15*2.d0/dsqrt(3.d0)
  
  !!z2
  c(0,0,2,0, 0)=  rac4pi/3.d0
  c(0,0,2,2, 0)=  racpi_15*4.d0/dsqrt(3.d0)
  
  do k=1,4
    nslat(k)=1+nix(k)+niy(k)+niz(k)+power(k)
    lslat_max(k)=nix(k)+niy(k)+niz(k)
  enddo
  
  rint=(0.d0,0.d0)
  
  do l1=0,lslat_max(1)
    do m1=-l1,l1
      do l2=0,lslat_max(2)
        do m2=-l2,l2
          do l3=0,lslat_max(3)
            do m3=-l3,l3
              do l4=0,lslat_max(4)
                do m4=-l4,l4
                  
                  lslat(1)=l1
                  lslat(2)=l2
                  lslat(3)=l3
                  lslat(4)=l4
                  
                  mslat(1)=-m1
                  mslat(2)=m2
                  mslat(3)=-m3
                  mslat(4)=m4
                  
                  coef=c(nix(1),niy(1),niz(1),l1,m1)*c(nix(2),niy(2),niz(2),l2,m2) *&
                      c(nix(3),niy(3),niz(3),l3,m3)*c(nix(4),niy(4),niz(4),l4,m4)
                  
                  if(abs(coef).gt.1.d-20)then
                    rint=rint + coef*(-1.d0)**(m1+m3)*ijkl_slater(nslat,lslat,mslat,gamma)
                  endif
                  
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  
  norm=1.d0
  do k=1,4
    norm=norm*norm_sto(gamma(k),nslat(k))
  enddo
  
  rint=rint/norm
  
  ijkl_slater_xnynzn=dreal(rint)
  
end

