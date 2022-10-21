double precision function doublefact(n)
        implicit double precision(a-h,o-z)
        doublefact=1.d0
        if(n.le.2)return
        d=1.d0
        do i=n,1,-2
          d=d*dfloat(i)
        enddo
        doublefact=d
end

