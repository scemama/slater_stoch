recursive subroutine I_x2_pol_mult_mono_elec(c,R1x,R1xp,R2x,d,nd,dim)
 implicit none
 integer , intent(in) :: dim
 double precision :: d(0:dim)
 integer,intent(inout) :: nd
 integer, intent(in):: c
 double precision, intent(in) :: R1x(0:2),R1xp(0:2),R2x(0:2)
 integer :: i,ix,nx,ny
 double precision :: X(0:dim),Y(0:dim)
!print*,'X2,c = ',c
!print*,'nd_in = ',nd

 if(c==0) then
   nd = 0
   d(0) = 1.d0
!  print*,'nd  IX2 = ',nd
   return
 elseif ((nd<0).or.(c<0))then
   nd = -1
    return
 else
     do ix=0,dim
       X(ix) = 0.d0
       Y(ix) = 0.d0
     enddo
     nx = 0
     call I_x1_pol_mult_mono_elec(0,c-2,R1x,R1xp,R2x,X,nx,dim)
!      print*,'nx 0,c-2 = ',nx
       do ix=0,nx
         X(ix) = X(ix)*(c-1)
       enddo
       call multiply_poly_2(X,nx,R2x,2,d,nd,dim)
!      print*,'nd = ',nd
       ny = 0
       do ix=0,dim
        Y(ix) = 0.d0
       enddo

       call I_x1_pol_mult_mono_elec(0,c-1,R1x,R1xp,R2x,Y,ny,dim)
!      print*,'ny = ',ny
!      do ix=0,ny
!        print*,'Y(ix) = ',Y(ix)
!      enddo
       if(ny.ge.0)then
        call multiply_poly_2(Y,ny,R1xp,2,d,nd,dim)
       endif
 endif
end

