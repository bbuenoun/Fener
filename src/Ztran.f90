!*==ZTRAN.spg  processed by SPAG 6.72Dc at 05:59 on  5 Nov 2016
      SUBROUTINE ZTRAN(Esp,Con,Den,Ces,Res,Ncap,U,Aoz,Aiz,Afz,Bz,Ncoe,  &
                     & Dt)
!
      IMPLICIT NONE
!*--ZTRAN6
!*** Start of declarations inserted by SPAG
      REAL Dt
      INTEGER*2 i , id , ii , ip , ivar , j , mm3 , mm4 , Ncap , Ncoe , &
              & nnn , nraiz , nreg
!*** End of declarations inserted by SPAG
      DOUBLE PRECISION rr(20) , beta(20) , Esp(20) , Con(20) , Den(20) ,&
                     & Ces(20) , Res(20) , raiz(80) , func(80,4) ,      &
                     & der(80,4) , dn(80,3) , ff(2,2) , Bz(80) , pol(80)&
                     & , Aoz(80) , Aiz(80) , Afz(80) , U , c1x , c1y ,  &
                     & c1z , t , pol6(80) , sumaa , sumab , sumac ,     &
                     & sumad
      DIMENSION ip(5)

!*** ztran(esp,cond,dens,csubp,res,ncap,umuro,a,c,b,d,ns,deltat)
!*** deltat - Paso de tiempo (s)
!*** ncap - Numero de capas


      U = 0.
      DO i = 1 , Ncap
         IF ( Res(i).EQ.0 ) THEN
            rr(i) = Esp(i)/Con(i)
            beta(i) = Esp(i)*DSQRT(Den(i)*Ces(i)/Con(i))
         ELSE
            rr(i) = Res(i)
            beta(i) = 0.
         ENDIF
         U = U + rr(i)
      ENDDO
      U = 1./U
      Ncoe = 0
      CALL POLOS(rr,beta,raiz,der,func,Ncap,Dt,nraiz,ip)
      nreg = 2
      DO ii = 1 , 4
         nreg = nreg + 1
         nreg = nreg + 1
      ENDDO
      DO ii = 1 , 4
         nreg = nreg + 1
         nreg = nreg + 1
      ENDDO
      IF ( nraiz.EQ.0 ) THEN
         Aoz(1) = U
         Afz(1) = U
         Aiz(1) = U
         Bz(1) = 1
         Ncoe = 1
         RETURN
      ENDIF
      DO i = 1 , 80
         Bz(i) = 0.
         pol(i) = 0.
         pol6(i) = 0.
         Aoz(i) = 0.
         Aiz(i) = 0.
         Afz(i) = 0.
      ENDDO
      nnn = 0
      Bz(1) = 1.
      pol(1) = 1.
      DO i = 1 , nraiz
         pol(2) = -DEXP(-raiz(i)*Dt)
         IF ( DABS(pol(2)).LT.1.D-16 ) GOTO 100
         CALL POLIN(Bz,pol,nnn,1)
      ENDDO
 100  id = nnn + 1
      CALL ORIGE(rr,beta,Ncap,ff)
      c1y = -U**2*ff(1,2)
      c1x = c1y + ff(2,2)*U
      c1z = c1y + ff(1,1)*U
      DO i = 1 , nraiz
         dn(i,2) = 1./(raiz(i)**2*der(i,2))
         dn(i,1) = dn(i,2)*func(i,4)
         dn(i,3) = dn(i,2)*func(i,1)
      ENDDO
      DO i = 1 , id
         pol(i) = Bz(i)
      ENDDO
      Aoz(1) = 1./Dt
      Aoz(2) = -2./Dt
      Aoz(3) = 1./Dt
      CALL POLIN(pol,Aoz,nnn,2)
      ivar = nnn + 1
      Aoz(1) = 0.
      Aoz(2) = 0.
      Aoz(3) = 0.
      DO i = 1 , 80
         t = i*Dt
         DO j = 1 , nraiz
            IF ( raiz(j)*t.GE.40. ) GOTO 150
            Aoz(i) = Aoz(i) + DEXP(-raiz(j)*t)*dn(j,1)
            Aiz(i) = Aiz(i) + DEXP(-raiz(j)*t)*dn(j,3)
            Afz(i) = Afz(i) + DEXP(-raiz(j)*t)*dn(j,2)
            IF ( j.GT.10 .AND.                                          &
               & (DABS(DEXP(-raiz(j)*t)*dn(j,2)).LT.1.D-16) ) GOTO 150
         ENDDO
 150     Aoz(i) = Aoz(i) + U*t + c1x
         Aiz(i) = Aiz(i) + U*t + c1z
         Afz(i) = Afz(i) + U*t + c1y
         IF ( i.GT.10 .AND. (DABS(Afz(i)).LT.1.D-16) ) GOTO 200
      ENDDO
      i = 81
 200  mm3 = i - 1
99001 FORMAT (10X,D15.7)
      CALL POLIN(Aoz,pol,mm3,nnn)
      mm3 = i - 1
      CALL POLIN(Aiz,pol,mm3,nnn)
      mm4 = i - 1
      CALL POLIN(Afz,pol,mm4,nnn)
      DO i = 1 , 80
         IF ( DABS(Bz(i)).LT.1.D-6 ) THEN
            Ncoe = i - 1
            GOTO 300
         ENDIF
      ENDDO
 300  sumaa = 0.
      sumab = 0.
      sumac = 0.
      sumad = 0.
      DO i = 1 , Ncoe + 1
         sumaa = sumaa + Aoz(i)
         sumab = sumab + Afz(i)
         sumac = sumac + Aiz(i)
         sumad = sumad + Bz(i)
      ENDDO
      Aiz(Ncoe+1) = Aiz(Ncoe+1) + sumad*U - sumac
      Afz(Ncoe+1) = Afz(Ncoe+1) + sumad*U - sumab
      Aoz(Ncoe+1) = Aoz(Ncoe+1) + sumad*U - sumaa
      END
!*==POLIN.spg  processed by SPAG 6.72Dc at 05:59 on  5 Nov 2016
      SUBROUTINE POLIN(A,B,M,N)
!     SUBRUTINA PARA CALCULO DEL PRODUCTO DE DOS POLINOMIOS
      IMPLICIT NONE
!*--POLIN134
!*** Start of declarations inserted by SPAG
      INTEGER*2 i , inn , j , k , l , M , N , nn
!*** End of declarations inserted by SPAG
      DOUBLE PRECISION A(80) , B(80) , c(80)
      k = M + N + 1
      IF ( k.GT.80 ) k = 80
      DO i = 1 , 80
         c(i) = 0.
      ENDDO
      DO i = 1 , k
         l = M + 1
         IF ( i.LT.l ) l = i
         j = i - N
         IF ( j.LT.1 ) j = 1
         DO nn = j , l
            inn = i - nn + 1
            c(i) = c(i) + A(nn)*B(inn)
         ENDDO
         IF ( i.GE.5 .AND. DABS(c(i)).LT..5D-6 ) GOTO 100
      ENDDO
      i = k + 1
 100  M = i - 1
      DO i = 1 , 80
         A(i) = c(i)
      ENDDO
      END
!*==POLOS.spg  processed by SPAG 6.72Dc at 05:59 on  5 Nov 2016
      SUBROUTINE POLOS(Rr,Beta,Raiz,Der,Func,M,Dt,Nraiz,Ip)
      IMPLICIT NONE
!*--POLOS164
!*** Start of declarations inserted by SPAG
      REAL Dt
      INTEGER*2 i , ia , ic , icon , Ip , j , jj , k , kkll , kx , ky , &
              & l , last , M , nn , Nraiz
!*** End of declarations inserted by SPAG
      DOUBLE PRECISION Rr(20) , Beta(20) , Raiz(80) , Der(80,4) ,       &
                     & Func(80,4) , f(2,2) , ff(2,2) , r1 , r2 , r3 ,   &
                     & f1 , f2 , fp1 , fp2 , rtemp
      DIMENSION icon(12) , Ip(5)
      DATA icon/1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12/
      kkll = 0
      ia = 0
      DO i = 1 , 80
         DO j = 1 , 4
            Func(i,j) = 0.
            Der(i,j) = 0.
         ENDDO
         Raiz(i) = 0.
      ENDDO
      last = 0
      r1 = 30./Dt
      r3 = 100./Dt
      Nraiz = 0
      DO l = 1 , 200
         Nraiz = Nraiz + 1
         IF ( Nraiz.GT.4 ) ia = 1
         IF ( Nraiz.EQ.1 ) GOTO 150
         CALL MATRI(Rr,Beta,r1,M,f,ff,2)
         f1 = f(1,2)
         fp1 = ff(1,2)
         ic = 0
         DO i = 1 , 20
99001       FORMAT (5X,I3)
            r2 = r1 + (r3-r1)/20.*i
            CALL MATRI(Rr,Beta,r2,M,f,ff,2)
            f2 = f(1,2)
            fp2 = ff(1,2)
            IF ( (f1.GT.0. .AND. f2.LE.0.) .OR.                         &
               & (f1.LE.0. .AND. f2.GT.0.) ) GOTO 100
            IF ( .NOT.((fp1.GT.0. .AND. fp2.GT.0.) .OR.                 &
               & (fp1.LE.0. .AND. fp2.LE.0.)) ) THEN
               ic = ic + 1
               IF ( ic.EQ.2 ) GOTO 100
            ENDIF
            f1 = f2
            fp1 = fp2
         ENDDO
         Nraiz = Nraiz - 1
         last = last - 1
 50      IF ( last.EQ.0 ) RETURN
         r1 = .0001/Dt
         IF ( last.NE.1 ) r1 = Raiz(last-1) + 0.00001/Dt
         r3 = Raiz(last) - .00001/Dt
         GOTO 400
 100     r3 = r2
 150     CALL MATRI(Rr,Beta,r1,M,f,ff,1)
         f1 = f(1,2)
         DO i = 1 , 25
            nn = 10*i
            DO j = 1 , nn
               r2 = r1 + j*(r3-r1)/nn
               CALL MATRI(Rr,Beta,r2,M,f,ff,1)
               f2 = f(1,2)
               IF ( (f1.GT.0. .AND. f2.GT.0.) .OR.                      &
                  & (f1.LE.0. .AND. f2.LE.0.) ) GOTO 200
 160           rtemp = (r1+r2)/2.
               CALL MATRI(Rr,Beta,rtemp,M,f,ff,1)
               IF ( f(1,2).LT.0 ) THEN
                  IF ( f1.GT.0. ) GOTO 170
               ELSEIF ( f(1,2).EQ.0 ) THEN
                  GOTO 190
               ELSEIF ( f1.LE.0. ) THEN
                  GOTO 170
               ENDIF
               f1 = f(1,2)
               r1 = rtemp
               GOTO 180
 170           f2 = f(1,2)
               r2 = rtemp
               kkll = kkll + 1
 180           IF ( DABS((r1-r2)/r1).GT.1.D-14 ) GOTO 160
 190           CALL MATRI(Rr,Beta,r2,M,f,ff,2)
               GOTO 250
 200        ENDDO
         ENDDO
         Nraiz = 0
         RETURN
 250     DO i = 1 , Nraiz
            j = i + 1
            IF ( r2.LE.Raiz(i) ) GOTO 300
         ENDDO
 300     last = last + 1
         IF ( Nraiz.NE.1 ) THEN
            jj = Nraiz + 1
            DO i = j , Nraiz
               jj = jj - 1
               Raiz(jj) = Raiz(jj-1)
               DO k = 1 , 4
                  Der(jj,k) = Der(jj-1,k)
                  Func(jj,k) = Func(jj-1,k)
               ENDDO
            ENDDO
         ENDIF
         j = j - 1
         Raiz(j) = r2
         DO k = 1 , 4
            kx = (k+1)/2
            ky = k/kx
            Der(j,k) = ff(kx,ky)
            Func(j,k) = f(kx,ky)
         ENDDO
         GOTO 50
 400  ENDDO
      END
!*==MATRI.spg  processed by SPAG 6.72Dc at 05:59 on  5 Nov 2016
      SUBROUTINE MATRI(Rr,Beta,W,M,F,Ff,Icont)
!     SUBRUTINA PARA CALCULO DE LA MATRIZ DE TRANSMISION
      IMPLICIT NONE
!*--MATRI283
!*** Start of declarations inserted by SPAG
      INTEGER*2 i , Icont , j , k , l , M , n
!*** End of declarations inserted by SPAG
      DOUBLE PRECISION Rr(20) , Beta(20) , f1(20,2,2) , f2(20,2,2) ,    &
                     & f3(20,2,2) , F(2,2) , Ff(2,2) , p , sq , W
      DO i = 1 , M
         sq = DSQRT(W)
         p = sq*Beta(i)
         IF ( p.NE.0. .AND. Icont.NE.3 ) THEN
            f1(i,1,1) = DCOS(p)
            f1(i,1,2) = Rr(i)/p*DSIN(p)
            f1(i,2,1) = -p/Rr(i)*DSIN(p)
            f1(i,2,2) = f1(i,1,1)
            IF ( Icont.NE.1 ) THEN
               f2(i,1,1) = Beta(i)*DSIN(Beta(i)*sq)/(2.*sq)
               f2(i,1,2) = -Rr(i)*DCOS(Beta(i)*sq)/(2.*W) + Rr(i)       &
                         & *DSIN(Beta(i)*sq)/(Beta(i)*2.*sq**3)
               f2(i,2,1) = Beta(i)**2*DCOS(Beta(i)*sq)/(2.*Rr(i))       &
                         & + DSIN(Beta(i)*sq)*Beta(i)/(2.*sq*Rr(i))
               f2(i,2,2) = f2(i,1,1)
            ENDIF
         ELSE
            f1(i,1,1) = 1.
            f1(i,1,2) = Rr(i)
            f1(i,2,1) = 0.
            f1(i,2,2) = 1.
            IF ( Icont.NE.1 ) THEN
               f2(i,1,1) = Beta(i)**2/2.
               f2(i,1,2) = Beta(i)**2*Rr(i)/6.
               f2(i,2,1) = Beta(i)**2/Rr(i)
               f2(i,2,2) = f2(i,1,1)
            ENDIF
         ENDIF
      ENDDO
      IF ( M.GT.1 ) THEN
         DO k = 1 , 2
            DO l = 1 , 2
               Ff(k,l) = 0.
            ENDDO
         ENDDO
         IF ( Icont.NE.1 ) THEN
            DO i = 1 , M
               DO j = 1 , M
                  DO k = 1 , 2
                     DO l = 1 , 2
                        f3(j,k,l) = f1(j,k,l)
                        IF ( i.EQ.j ) f3(j,k,l) = f2(j,k,l)
                     ENDDO
                  ENDDO
                  IF ( j.NE.1 ) THEN
                     DO k = 1 , 2
                        DO l = 1 , 2
                           F(k,l) = 0.
                           DO n = 1 , 2
                              F(k,l) = F(k,l) + f3(j-1,k,n)*f3(j,n,l)
                           ENDDO
                           f3(j,k,l) = F(k,l)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               DO k = 1 , 2
                  DO l = 1 , 2
                     Ff(k,l) = Ff(k,l) + F(k,l)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         DO i = 2 , M
            DO k = 1 , 2
               DO l = 1 , 2
                  F(k,l) = 0.
               ENDDO
            ENDDO
            DO k = 1 , 2
               DO l = 1 , 2
                  DO n = 1 , 2
                     F(k,l) = F(k,l) + f1(i-1,k,n)*f1(i,n,l)
                  ENDDO
               ENDDO
            ENDDO
            DO k = 1 , 2
               DO l = 1 , 2
                  f1(i,k,l) = F(k,l)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO k = 1 , 2
            DO l = 1 , 2
               IF ( Icont.NE.1 ) Ff(k,l) = f2(1,k,l)
               F(k,l) = f1(1,k,l)
            ENDDO
         ENDDO
         RETURN
      ENDIF
      END
!*==ORIGE.spg  processed by SPAG 6.72Dc at 05:59 on  5 Nov 2016
      SUBROUTINE ORIGE(Rr,Beta,M,Mp)
      IMPLICIT NONE
!*--ORIGE384
!*** Start of declarations inserted by SPAG
      INTEGER*2 i , j , k , l , ll , lll , M
!*** End of declarations inserted by SPAG
!     SUBRUTINA PARA CALCULO DEL ORIGEN
      DOUBLE PRECISION Rr(20) , Beta(20) , Mp(2,2) , a(20,2,2) ,        &
                     & b(20,2,2) , d(2,2) , f(2,2) , temp(2,2) , p , r
      DO i = 1 , 2
         DO j = 1 , 2
            Mp(i,j) = 0.
         ENDDO
      ENDDO
      DO i = 1 , M
         p = Beta(i)*Beta(i)
         r = Rr(i)
         a(i,1,1) = 1.
         a(i,1,2) = r
         a(i,2,1) = 0.
         a(i,2,2) = 1.
         b(i,1,1) = p/2.
         b(i,1,2) = r*p/6.
         b(i,2,1) = p/r
         b(i,2,2) = p/2.
      ENDDO
      IF ( M.NE.1 ) THEN
         DO i = 1 , M
            DO k = 1 , M
               IF ( i.EQ.k ) THEN
                  IF ( i.NE.1 ) THEN
                     d(1,1) = a(1,1,1)
                     d(1,2) = a(1,1,2)
                     d(2,1) = a(1,2,1)
                     d(2,2) = a(1,2,2)
                  ELSE
                     d(1,1) = b(1,1,1)
                     d(1,2) = b(1,1,2)
                     d(2,1) = b(1,2,1)
                     d(2,2) = b(1,2,2)
                  ENDIF
               ENDIF
               DO j = 2 , M
                  IF ( i.EQ.k ) THEN
                     IF ( i.EQ.j ) THEN
                        f(1,1) = b(j,1,1)
                        f(1,2) = b(j,1,2)
                        f(2,1) = b(j,2,1)
                        f(2,2) = b(j,2,2)
                     ELSE
                        f(1,1) = a(j,1,1)
                        f(1,2) = a(j,1,2)
                        f(2,1) = a(j,2,1)
                        f(2,2) = a(j,2,2)
                     ENDIF
                     DO l = 1 , 2
                        DO ll = 1 , 2
                           temp(l,ll) = 0.
                           DO lll = 1 , 2
                              temp(l,ll) = temp(l,ll) + d(l,lll)        &
                               & *f(lll,ll)
                           ENDDO
                        ENDDO
                     ENDDO
                     DO l = 1 , 2
                        DO ll = 1 , 2
                           d(l,ll) = temp(l,ll)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            DO l = 1 , 2
               DO ll = 1 , 2
                  Mp(l,ll) = Mp(l,ll) + temp(l,ll)
               ENDDO
            ENDDO
         ENDDO
         GOTO 99999
      ENDIF
      DO i = 1 , 2
         DO j = 1 , 2
            Mp(i,j) = b(1,i,j)
         ENDDO
      ENDDO
      RETURN
99999 END
