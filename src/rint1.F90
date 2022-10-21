double precision function rint1(n,rho)
  implicit none
  !  BEGIN_DOC
  ! Standard version of rint
  !  END_DOC
  integer, intent(in)            :: n
  double precision, intent(in)   :: rho
  !!!!!!  double precision, parameter :: eps=1.d-13
  double precision, parameter    :: eps=1.d-15
  double precision               :: rho_tmp, diff
  double precision               :: fact_inv
  integer                        :: k
  double precision, parameter, dimension(100) :: inv_int = (/ 1.d0, 1.d0/2.d0, 1.d0/3.d0,&
      1.d0/4.d0, 1.d0/5.d0, 1.d0/6.d0, 1.d0/7.d0, 1.d0/8.d0, 1.d0/9.d0, 1.d0/10.d0,&
      1.d0/11.d0, 1.d0/12.d0, 1.d0/13.d0, 1.d0/14.d0, 1.d0/15.d0, 1.d0/16.d0, 1.d0/17.d0,&
      1.d0/18.d0, 1.d0/19.d0, 1.d0/20.d0, 1.d0/21.d0, 1.d0/22.d0, 1.d0/23.d0, 1.d0/24.d0,&
      1.d0/25.d0, 1.d0/26.d0, 1.d0/27.d0, 1.d0/28.d0, 1.d0/29.d0, 1.d0/30.d0, 1.d0/31.d0,&
      1.d0/32.d0, 1.d0/33.d0, 1.d0/34.d0, 1.d0/35.d0, 1.d0/36.d0, 1.d0/37.d0, 1.d0/38.d0,&
      1.d0/39.d0, 1.d0/40.d0, 1.d0/41.d0, 1.d0/42.d0, 1.d0/43.d0, 1.d0/44.d0, 1.d0/45.d0,&
      1.d0/46.d0, 1.d0/47.d0, 1.d0/48.d0, 1.d0/49.d0, 1.d0/50.d0, 1.d0/51.d0, 1.d0/52.d0,&
      1.d0/53.d0, 1.d0/54.d0, 1.d0/55.d0, 1.d0/56.d0, 1.d0/57.d0, 1.d0/58.d0, 1.d0/59.d0,&
      1.d0/60.d0, 1.d0/61.d0, 1.d0/62.d0, 1.d0/63.d0, 1.d0/64.d0, 1.d0/65.d0, 1.d0/66.d0,&
      1.d0/67.d0, 1.d0/68.d0, 1.d0/69.d0, 1.d0/70.d0, 1.d0/71.d0, 1.d0/72.d0, 1.d0/73.d0,&
      1.d0/74.d0, 1.d0/75.d0, 1.d0/76.d0, 1.d0/77.d0, 1.d0/78.d0, 1.d0/79.d0, 1.d0/80.d0,&
      1.d0/81.d0, 1.d0/82.d0, 1.d0/83.d0, 1.d0/84.d0, 1.d0/85.d0, 1.d0/86.d0, 1.d0/87.d0,&
      1.d0/88.d0, 1.d0/89.d0, 1.d0/90.d0, 1.d0/91.d0, 1.d0/92.d0, 1.d0/93.d0, 1.d0/94.d0,&
      1.d0/95.d0, 1.d0/96.d0, 1.d0/97.d0, 1.d0/98.d0, 1.d0/99.d0, 1.d0/100.d0 /)
  
  if (n<50) then
    rint1=inv_int(n+n+1)
    rho_tmp = 1.d0
    do k=1,20
      rho_tmp = -rho_tmp*rho
      diff=rho_tmp*fact_inv(k)*inv_int(ishft(k+n,1)+1)
      
      rint1=rint1+diff
      if (dabs(diff) > eps) then
        cycle
      endif
      return
    enddo
  else
    rint1=1.d0/dfloat(n+n+1)
    rho_tmp = 1.d0
    do k=1,20
      rho_tmp = -rho_tmp*rho
      diff=rho_tmp*fact_inv(k)/dfloat(shiftl(k+n,1)+1)
      rint1=rint1+diff
      if (dabs(diff) > eps) then
        cycle
      endif
      return
    enddo
  endif
  
  write(*,*)'n rhod=',n,rho
  write(*,*)'diff=',diff,' pb in rint1 k too large!'
  stop 1
end

