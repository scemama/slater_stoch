! Computation of I= (nlm_1 nlm_2 | nlm_3 nlm_4)= /int d vector_1 d vector_2  (n1 l1 m1)*(r1)(n2 l2 m2)(r1) 1/r12 (n3 l3 m3)*(r2) (n4 l4 m4)(r2)
!!
!! where (nlm)*(vector r) = N_n r^{n-1} exp(-g r) Y_l^m(theta,phi)   N_n=sqrt((2g)**(2n+1)/(2n)!)
!! with Y_l^m(theta,phi) = i^(m+|m|) ([(2l+1)*(l-|m|)!]/[4pi*(l+|m|)!])^1/2 P_l^|m|(cos(theta))  exp(i m phi)
!!
!! l=0,1,2,...,n-1
!! m=-l...,l
!!
!! Here P_l^|m|(cos(theta)) = "Associated Legendre Polynomials wikipedia" computed with routine:  SUBROUTINE LPMN(MM,M,N,X,PM)
!!
!! In input of LPMN:  n=l (n=0,1,...)  m=0,1,...,n x=cos(theta) 0 (< or =)  x (< ! or =) 1
!!
!! the routine computes:   PM(m,n) for n=0,...,N (number N in input) and  m=0,..,n
!!
!! Exemples (see 'Associated Legendre Polynomilas wikipedia')
!!   P_{0}^{0}(x)=1
!!   P_{1}^{-1}(x)=-1/2 P_{1}^{1}(x)
!!   P_{1}^{0}(x)=x
!!   P_{1}^{1}(x)=-(1-x^2)^{1/2}
!!   P_{2}^{-2}(x)=1/24 P_{2}^{2}(x)
!!   P_{2}^{-1}(x)=-1/6 P_{2}^{1}(x)
!!   P_{2}^{0}(x)=1/2 (3x^{2}-1)
!!   P_{2}^{1}(x)=-3x(1-x^2)^{1/2}
!!   P_{2}^{2}(x)=3(1-x^2)
!!
!!
!!*************
!!  Conditions: l1+l2+l3+l4 even and  m2-m1=m3-m4  otherwise I=0
!!*************
!!
!! I= N_n1*N_n2*N_n3*N_n4  \sum_{lmin}^{lmax},2  4*pi/(2l+1)  [G_{l1m1 lm}^{l2m2} G_{l1m1 lm}^{l2m2}] * bigF(n1+n2,g1+g2,n3+n4,g3+g4,l)
!!
!!  Gaunt coefficient:  G_{l1m1 l2m2}^lm =  int_dOmega  Ylm^*  Yl1m1  Yl2m2       m=m1+m2 (0 otherwise)
!!
!!  lmin=max(l12_min, l34_min)
!!  lmax=min(l12_max, l34_max)
!!
!!  lij_max=li+lj
!!  lij_min=  max(|li-lj|,|mi-mj|)     if [lij_max + max(|li-lj|,|mi-mj|)] = even
!!  lij_min=  max(|li-lj|,|mi-mj|)+1   if [lij_max + max(|li-lj|,|mi-mj|)] = odd
!!
!!
!! bifF defined as:
!!
!! bigF(k1,g1,k2,g2,l)= int_0^inf r1^k1 exp(-g1 r1) [1/r1**(l+1) int_0^r1 r2**(k2+l) exp(-g2 r2) + r1**l int_r1^inf r2**(k2-l-1) exp(-g2 r2)]
!!
!! k1 integer larger or equal to 0 , k2 integer larger or equal to 0, k1 and k2 larger and equal to (l+1)
!!
!! We can show that
!!
!! bigF=  [(k1-l-1)!(k2+l)!/g1**(k1-l)/g2**(k2+l+1)]
!!      + \sum_{i=0}^{k2-l-1} [(k2-l-1)!(k1+k2-1-i)!/(k2-l-1-i)!/g2**(i+1)/(g1+g2)**(k1+k2-i)]
!!      - \sum_{i=0}^{k2+l  } [(k2+l)!(k1+k2-1-i)!/(k2+l-i)!/g2**(i+1)/(g1+g2)**(k1+k2-i)]
!!
!! This formula can been checked with bigF_num
!!
!!
double precision function ijkl_slater(nslat,lslat,mslat,gamma)
  implicit none
  integer nslat(4),lslat(4),mslat(4),i,l,lmin,lmax
  integer l1,l2,l3,l4,m1,m2,m3,m4,n12,n34,m
  double precision gamma(4),pi,int,norm,norm_sto,gaunt,gamm12,gamm34,bigF
  
  ijkl_slater=0.d0
  
  l1=lslat(1)
  l2=lslat(2)
  l3=lslat(3)
  l4=lslat(4)
  if(mod(l1+l2+l3+l4,2).eq.1)return
  
  m1=mslat(1)
  m2=mslat(2)
  m3=mslat(3)
  m4=mslat(4)
  if((m2-m1).ne.(m3-m4))then
    return
  else
    m=m2-m1
  endif
  
  pi=dacos(-1.d0)
  
  norm=1.d0
  do i=1,4
    norm=norm*norm_sto(gamma(i),nslat(i))
  enddo
  
  call compute_lmin_lmax(l1,l2,l3,l4,m1,m2,m3,m4,lmin,lmax)
  
  n12=nslat(1)+nslat(2)
  n34=nslat(3)+nslat(4)
  gamm12=gamma(1)+gamma(2)
  gamm34=gamma(3)+gamma(4)
  
  int=0.d0
  do l=lmin,lmax,2
    int=int+gaunt(l2,l,l1,m2,m,m1)*gaunt(l3,l,l4,m3,m,m4)*bigF(n12,gamm12,n34,gamm34,l)/(2.d0*l+1.d0)
  enddo
  ijkl_slater=4.d0*pi*int*norm
  
end

