double precision function ijkl_slater_cartesian(nslat,lslat,mslat,gamma)
  implicit none
  integer nslat(4),lslat(4),mslat(4)
  double precision gamma(4),ijkl_slater,rac2
  integer sig(4),epsk1,epsk2,epsk3,epsk4,eps(4),i,mslat_k(4)
  double complex int,coef(-1:1,-1:1,-1:1),imag
  parameter(rac2=1.41421356237310d0)
  
  do i=1,4
    sig(i)=(-1)**mslat(i)
    eps(i)=mslat(i)*sig(i)
  enddo
  
  coef(-1,-1,-1)=(0.d0, 1.d0)
  coef(-1,-1,+1)=(0.d0, 1.d0)
  
  coef(-1, 0,-1)=(0.d0, 0.d0)
  coef(-1, 0,+1)=(0.d0, 0.d0)
  
  coef(-1,+1,-1)=(0.d0, 1.d0)
  coef(-1,+1,+1)=(0.d0,-1.d0)
  
  coef( 0,-1,-1)=(0.d0, 0.d0)
  coef( 0,-1,+1)=(0.d0, 0.d0)
  
  coef( 0, 0,-1)=(rac2, 0.d0)
  coef( 0, 0,+1)=(rac2, 0.d0)
  
  coef( 0,+1,-1)=(0.d0, 0.d0)
  coef( 0,+1,+1)=(0.d0, 0.d0)
  
  coef(+1,-1,-1)=(1.d0, 0.d0)
  coef(+1,-1,+1)=(1.d0, 0.d0)
  
  coef(+1, 0,-1)=(0.d0, 0.d0)
  coef(+1, 0,+1)=(0.d0, 0.d0)
  
  coef(+1,+1,-1)=(0.d0,-1.d0)
  coef(+1,+1,+1)=(0.d0, 1.d0)
  
  int=(0.d0,0.d0)
  
  do epsk1=-1,1
    do epsk2=-1,1
      do epsk3=-1,1
        do epsk4=-1,1
          
          mslat_k(1)=epsk1*iabs(mslat(1))
          mslat_k(2)=epsk2*iabs(mslat(2))
          mslat_k(3)=epsk3*iabs(mslat(3))
          mslat_k(4)=epsk4*iabs(mslat(4))
          
          int=int +conjg(coef(eps(1),sig(1),epsk1))*coef(eps(2),sig(2),epsk2)&
              *conjg(coef(eps(3),sig(3),epsk3))*coef(eps(4),sig(4),epsk4)*&
              ijkl_slater(nslat,lslat,mslat_k,gamma)
        enddo
      enddo
    enddo
  enddo
  
  print*,'int=',int
  
  ijkl_slater_cartesian=int
end


