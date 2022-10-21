!! Here gam=gamA+gamB
!! one_electron_I= Ix
      double precision function one_electron_I(na,ra,nb,rb,gam,rp)
      implicit double precision(a-h,o-z)
      pi=dacos(-1.d0)
      one_electron_I=0.d0

      rint=0.d0
      do kx=0,na
       do lx=0,nb
        iskip=0
        kbig=na+nb-kx-lx
        if(mod(kbig,2).eq.1)iskip=1
        if(iskip.eq.0)then
         kbig=kbig/2
         f3=doublefact(2*kbig-1)
         f2=(1.d0/(2.d0*gam))**kbig
         rint=rint+f2*f3*binom(na,kx)*binom(nb,lx)*(rp-ra)**kx*(rp-rb)**lx
        endif
       enddo
      enddo

      one_electron_I=dsqrt(pi/gam)*rint
      end
