!    IMPORTANT:    call as follows:  gaunt(l,l1,l2,m,m1,m2) = /int d_Omega Y*_l^m  Y_l1^m1  Y_l2^m2
!
!      Y_l^m = (-1)**[(m+|m|)/2] C_lm  P_l^|m| ( cos[theta] ) exp(i m phi)
!
!    C_lm =sqrt(  (2*l+1)/(4pi) (l-|m|)!/(l+|m|)!  )
!
!     m=m1+m2
!
!     Be careful in the routine  l1 --> m  l2---> l1  l3 --> l2
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: gaunt
! !INTERFACE:
Real (8) Function gaunt (l1, l2, l3, m1, m2, m3)
! !INPUT/OUTPUT PARAMETERS:
! l1, l2, l3 : angular momentum quantum numbers (in,integer)
! m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
! Returns the Gaunt coefficient given by
! $$ \langle Y^{l_1}_{m_1}|Y^{l_2}_{m_2}|Y^{l_3}_{m_3} \rangle
! = (-1)^{m_1}\left[\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi} \right]
! ^{\frac{1}{2}}
! \begin{pmatrix} l_1 & l_2 & l_3 \\ 0 & 0 & 0 \end{pmatrix}
! \begin{pmatrix} l_1 & l_2 & l_3 \\ -m_1 & m_2 & m_3 \end{pmatrix}. $$
! Suitable for $l_i$ less than 50.
!
! !REVISION HISTORY:
! Created November 2002 (JKD)
!EOP
!BOC
      Implicit None
      ! arguments
      Integer                        :: l1
      Integer                        :: l2
      Integer                        :: l3
      Integer                        :: m1
      Integer                        :: m2
      Integer                        :: m3
      ! local variables
      Integer                        :: j, j1, j2, j3, jh
      Real (8)                       :: t1
      ! real constant 1/sqrt(4*pi)
      Real (8), Parameter            :: c1 = 0.28209479177387814347d0
      ! external functions
      Real (8)                       :: wigner3j, factr, factnm
      External wigner3j, factr, factnm
      If ((l1 .Lt. 0) .Or. (l2 .Lt. 0) .Or. (l3 .Lt. 0) .Or. (Abs(m1)&
                                                                     &.Gt. l1) .Or. (Abs(m2) .Gt. l2) .Or. (Abs(m3) .Gt. l3)) Then
        Write (*,*)
        Write (*, '("Error(gaunt): non-physical arguments :")')
        Write (*, '("l1 = ", I8, " l2 = ", I8, " l3 = ", I8)') l1, l2,&
                                                                     &l3
        Write (*, '("m1 = ", I8, " m2 = ", I8, " m3 = ", I8)') m1, m2,&
                                                                     &m3
        Write (*,*)
        Stop
      End If
      If ((l1 .Gt. 50) .Or. (l2 .Gt. 50) .Or. (l3 .Gt. 50)) Then
        Write (*,*)
        Write (*, '("Error(gaunt): angular momenta out of range : ", 3&
                                                                     &I8)') l1, l2, l3
        Write (*,*)
        Stop
      End If
      If (m1-m2-m3 .Ne. 0) Then
        gaunt = 0.d0
        Return
      End If
      j1 = l2 - l1 + l3
      j2 = l1 - l2 + l3
      j3 = l1 + l2 - l3
      If ((j1 .Lt. 0) .Or. (j2 .Lt. 0) .Or. (j3 .Lt. 0)) Then
        gaunt = 0.d0
        Return
      End If
      j = l1 + l2 + l3
      If (Mod(j, 2) .Ne. 0) Then
        gaunt = 0.d0
        Return
      End If
      jh = j / 2
      t1 = Sqrt (dble((2*l1+1)*(2*l2+1)*(2*l3+1))*factr(j1, j+1)*factnm(j2, 1)*factnm(j3, 1))
      t1 = t1 * factr (jh, jh-l1) / (factnm(jh-l2, 1)*factnm(jh-l3, 1))
      gaunt = t1 * c1 * wigner3j (l1, l2, l3,-m1, m2, m3)
      If (Mod(m1+jh, 2) .Ne. 0) gaunt = - gaunt
      Return
End Function
!EOC

!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: factr
! !INTERFACE:
Real (8) Function factr (n, d)
! !INPUT/OUTPUT PARAMETERS:
! n : numerator (in,integer)
! d : denominator (in,integer)
! !DESCRIPTION:
! Returns the ratio $n!/d!$ for $n,d\ge 0$. Performs no under- or overflow
! checking.
!
! !REVISION HISTORY:
! Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
      ! arguments
      Integer, Intent (In)           :: n
      Integer, Intent (In)           :: d
      ! local variables
      Integer                        :: i
      If (n .Lt. 0) Then
        Write (*,*)
        Write (*, '("Error(factr): n < 0 : ", I8)') n
        Write (*,*)
        Stop
      End If
      If (d .Lt. 0) Then
        Write (*,*)
        Write (*, '("Error(factr): d < 0 : ", I8)') d
        Write (*,*)
        Stop
      End If
      If (n .Lt. d) Then
        factr = 1.d0 / dble (n+1)
        Do i = n + 2, d
          factr = factr / dble (i)
        End Do
      Else If (n .Eq. d) Then
        factr = 1.d0
      Else
        factr = dble (d+1)
        Do i = d + 2, n
          factr = factr * dble (i)
        End Do
      End If
      Return
End Function
!EOC

!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: factnm
! !INTERFACE:
Real (8) Function factnm (n, m)
! !INPUT/OUTPUT PARAMETERS:
! n : input (in,integer)
! m : order of multifactorial (in,integer)
! !DESCRIPTION:
! Returns the multifactorial
! $$ n\underbrace{!!\,...\,!}_{m\,{\rm times}}=\prod_{i\ge 0,\,n-im>0}n-im $$
! for $n,\,m \ge 0$. $n$ should be less than 150.
!
! !REVISION HISTORY:
! Created January 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Integer, Intent (In) :: m
! local variables
      Integer :: i, j
      Real (8) :: f1 (24), f2 (38)
      Data f1 / 1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0, 5040.d0, &
     & 40320.d0, 362880.d0, 3628800.d0, 39916800.d0, 479001600.d0, &
     & 6227020800.d0, 87178291200.d0, 1307674368000.d0, &
     & 20922789888000.d0, 355687428096000.d0, 6402373705728000.d0, &
     & 121645100408832000.d0, 2432902008176640000.d0, &
     & 51090942171709440000.d0, 1124000727777607680000.d0, &
     & 25852016738884976640000.d0, 620448401733239439360000.d0 /
      Data f2 / 1.d0, 2.d0, 3.d0, 8.d0, 15.d0, 48.d0, 105.d0, 384.d0, &
     & 945.d0, 3840.d0, 10395.d0, 46080.d0, 135135.d0, 645120.d0, &
     & 2027025.d0, 10321920.d0, 34459425.d0, 185794560.d0, &
     & 654729075.d0, 3715891200.d0, 13749310575.d0, 81749606400.d0, &
     & 316234143225.d0, 1961990553600.d0, 7905853580625.d0, &
     & 51011754393600.d0, 213458046676875.d0, 1428329123020800.d0, &
     & 6190283353629375.d0, 42849873690624000.d0, &
     & 191898783962510625.d0, 1371195958099968000.d0, &
     & 6332659870762850625.d0, 46620662575398912000.d0, &
     & 221643095476699771875.d0, 1678343852714360832000.d0, &
     & 8200794532637891559375.d0, 63777066403145711616000.d0 /
! fast return if possible
If (n .Eq. 0) Then
  factnm = 1.d0
  Return
End If
If (m .Eq. 1) Then
  If ((n .Ge. 1) .And. (n .Le. 24)) Then
    factnm = f1 (n)
    Return
  End If
End If
If (m .Eq. 2) Then
  If ((n .Ge. 1) .And. (n .Le. 38)) Then
    factnm = f2 (n)
    Return
  End If
End If
If (n .Lt. 0) Then
  Write (*,*)
  Write (*, '("Error(factnm): n < 0 : ", I8)') n
  Write (*,*)
  Stop
End If
If (m .Le. 0) Then
  Write (*,*)
  Write (*, '("Error(factnm): m <= 0 : ", I8)') m
  Write (*,*)
  Stop
End If
If (n .Gt. 150) Then
  Write (*,*)
  Write (*, '("Error(factnm): n out of range : ", I8)') n
  Write (*,*)
  Stop
End If
If (m .Eq. 1) Then
  factnm = f1 (24)
  Do i = 25, n
    factnm = factnm * dble (i)
  End Do
Else
  j = n / m
  If (Mod(n, m) .Eq. 0) j = j - 1
  factnm = dble (n)
  Do i = 1, j
    factnm = factnm * dble (n-i*m)
  End Do
End If
Return
End Function
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: wigner3j
! !INTERFACE:
real (8) Function wigner3j (j1, j2, j3, m1, m2, m3)
! !INPUT/OUTPUT PARAMETERS:
! j1, j2, j3 : angular momentum quantum numbers (in,integer)
! m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
! Returns the Wigner $3j$-symbol. There are many equivalent definitions for
! the $3j$-symbols, the following provides high accuracy for $j\le 50$
! \begin{align*}
! &\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix}=(-1)^
! {j1+j2+m3}\\
! &\times\sqrt{\frac{(j_1+m_1)!(j_2+m_2)!(j_3+m_3)!(j_3-m_3)!(j_1-m_1)!
! (j_2-m_2)!}{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!(1+j_1+j_2+j_3)!}}
! \times\sum_{\max(0,j_2-j_3-m_1,j_1-j_3+m_2)}^
! {\min(j_1+j_2-j_3,j_1-m_1,j_2+m_2)}\\
! &(-1)^k\frac{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!}
! {(j_3-j_1-m_2+k)!(j_3-j_2+m_1+k)!(j_1+j_2-j_3-k)!k!(j_1-m_1-k)!
! (j_2+m_2-k)}.
! \end{align*}
!
! !REVISION HISTORY:
! Created November 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: j1
      Integer, Intent (In) :: j2
      Integer, Intent (In) :: j3
      Integer, Intent (In) :: m1
      Integer, Intent (In) :: m2
      Integer, Intent (In) :: m3
! local variables
      Integer :: k, k1, k2, l1, l2, l3, n1, n2
      Real (8) :: sgn, sum, t1
! external functions
      Real (8) :: factnm, factr
      External factnm, factr
! check input variables
      If ((j1 .Lt. 0) .Or. (j2 .Lt. 0) .Or. (j3 .Lt. 0) .Or. (Abs(m1)&
        .Gt. j1) .Or. (Abs(m2) .Gt. j2) .Or. (Abs(m3) .Gt. j3)) Then
        Write (*,*)
        Write (*, '("Error(wigner3j): invalid arguments :")')
        Write (*, '("j1 = ", I8, " j2 = ", I8, " j3 = ", I8)') j1, j2, j3
        Write (*, '("m1 = ", I8, " m2 = ", I8, " m3 = ", I8)') m1, m2, m3
        Write (*,*)
        Stop
      End If
      If ((j1 .Eq. 0) .And. (j2 .Eq. 0) .And. (j3 .Eq. 0)) Then
        wigner3j = 1.d0
        Return
      End If
      If ((j1 .Gt. 50) .Or. (j2 .Gt. 50) .Or. (j3 .Gt. 50)) Then
        Write (*,*)
        Write (*, '("Error(wigner3j): angular momenta out of range : ", 3I8)') j1, j2, j3
        Write (*,*)
        Stop
      End If
      l1 = j2 - j1 + j3
      l2 = j1 - j2 + j3
      l3 = j1 + j2 - j3
      If ((m1+m2+m3 .Ne. 0) .Or. (l1 .Lt. 0) .Or. (l2 .Lt. 0) .Or. (l3&
                                                                     &.Lt. 0)) Then
        wigner3j = 0.d0
        Return
      End If
      n1 = j1 - m1
      n2 = j2 + m2
      k1 = Max (0, j2-j3-m1, j1-j3+m2)
      k2 = Min (l3, n1, n2)
      If (Mod(k1+j1+j2+m3, 2) .Ne. 0) Then
        sgn = - 1.d0
      Else
        sgn = 1.d0
      End If
      sum = 0.d0
      Do k = k1, k2
        t1 = sgn * factr (l1, l1-n2+k) * factr (l2, l2-n1+k) * factr (l3, l3-k)
        sum = sum + t1 / (factnm(k, 1)*factnm(n1-k, 1)*factnm(n2-k, 1))
        sgn = - sgn
      End Do
      t1 = factr (j1+m1, l1) * factr (j2+m2, l2) * factr (j3+m3, l3)
      t1 = t1 * factr (j3-m3, 1+j1+j2+j3) * factnm (j1-m1, 1) * factnm (j2-m2, 1)
      wigner3j = sum * Sqrt (t1)
      Return
End Function
      !EOC

