C***WARNING***WARNING***WARNING***WARNING***WARNING***WARNING***
C
C This version of the hydrocarbon potential is currently
C (8/3/95) unpublished. There may still be some small changes
C forthcoming, particularly with respect to interstitial
C hydrogen in diamond. For updates, e-mail D. Brenner at
C dwb@ripley.mte.ncsu.edu.
C
C input coordinate file coord.d is written over; remember to
C keep a backup copy!!!!!
C
C***WARNING***WARNING***WARNING***WARNING***WARNING***WARNING***

      SUBROUTINE CAGUTS
C
C CALCULATE TWO-BODY FORCES AND NEIGHBOR LIST for hydrocarbons
C
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'
c
      DIMENSION RPP(nlmax,3),RR(3),RI(3)
C
      do i=1,np
           eatom(i) = 0.0d0
      enddo
C
      IF(LCHK.EQ.1) THEN
C
C Set up neighbor list
C
           K=0
           DO 302 I=1,NP
                NABORS(I)=K+1
                DO 299 L=1,3
                     RI(L)=R0(I,L)
299             CONTINUE
                KI=KTYPE(I)
C
c cuts out all but C,H,Si, and Ge
C
                if(ki.ge.5) go to 302
C
                DO 301 J=1,NP
C
                     IF(I.EQ.J) GO TO 301
C
                     KJ=KTYPE(J)
C
c cuts out all but C,H,Si, and Ge
C
                     if(kj.ge.5) go to 301
                     RLIS=RLIST(KI,KJ)
C
                     RSQ=0.0D0
                     DO 298 L=1,3
                          RR(L)=RI(L)-R0(J,L)
                          RR(L)=RR(L) -
     &                          CUBE(L)*ANINT(RR(L)/CUBE(L))
                          RSQ=RSQ+RR(L)*RR(L)
                          IF(RSQ.GT.RLIS) GO TO 301
298                  CONTINUE
C
405                  CONTINUE
                     K=K+1
                     LIST(K)=J
                     IVCT2B(K)=I
                     JVCT2B(K)=J
C
301             CONTINUE
302        CONTINUE
C
           NABORS(NP+1)=K+1
           KEND=K
           if(kend.gt.nlmax) then
                 write(*,*) 'kend exceeds nlmax'
                 write(*,*) 'kend,nlmax = ',kend,nlmax
                 write(*,*) 'increase nlmax and recompile' 
                 include 'close.inc' 
                 stop 
           endif
c           write(*,*) 'kend= ',kend
      ENDIF
c
      DO 320 K=1,KEND
C
           I=IVCT2B(K)
           J=JVCT2B(K)
           KI=KTYPE(I)
           KJ=KTYPE(J)
C
           LCHECK(K)=0
           RSQ=0.0D0
           DO L=1,3
                RR(L)=R0(I,L)-R0(J,L)
                RR(L)=RR(L) - CUBE(L)*ANINT(RR(L)/CUBE(L))
                RSQ=RSQ+RR(L)*RR(L)
                COR(K,L)=RR(L)
           ENDDO
c
           IF(RSQ.GT.RMAX(KI,KJ)) GOTO 320
           if((kj.le.2).and.(ki.le.2)) LCHECK(K)=1
           if((kj.ge.3).and.(ki.ge.3)) LCHECK(K)=2
           RC=SQRT(RSQ)
           rt = rc/ddtab(ki,kj)

           it=0
           if(int(rt)+1 .gt. ntab - 1) then
              it=ntab-1
           else
              it=int(rt)+1
           end if

           RCOR(K)=RC
C
           WW(K)=TABFC(ki,kj,it)
     &            +(TABFC(ki,kj,it+1)-TABFC(ki,kj,it))*(rt-it+1)
           DWW(K)=TABDFC(ki,kj,it)
     &           +(TABDFC(ki,kj,it+1)-TABDFC(ki,kj,it))*(rt-it+1)
C
          EXX1(K) = atable(ki,kj,it)
     &              +(atable(ki,kj,it+1)-atable(ki,kj,it))*(rt-it+1)
          DEXX1(K) = datable(ki,kj,it) +
     &        (datable(ki,kj,it+1)-datable(ki,kj,it))*(rt - it +1)
C
          IF(I.GE.J) GO TO 320
C
          vv = rtable(ki,kj,it)
     &              +(rtable(ki,kj,it+1)-rtable(ki,kj,it))*(rt-it+1)
          rp = drtable(ki,kj,it)
     &              +(drtable(ki,kj,it+1)-drtable(ki,kj,it))*(rt-it+1)
          tote = tote + vv
 
          eatom(i) = eatom(i) + vv/2.0d0
          eatom(j) = eatom(j) + vv/2.0d0


          DO 318 L=1,3
               RPP(K,L)=RP*RR(L)
318       CONTINUE
320   CONTINUE
C
      DO 321 K=1,KEND
           if(lcheck(k).eq.0) go to 321
           I=IVCT2B(K)
           J=JVCT2B(K)
           IF(I.GE.J) GO TO 321
            DO 322 L=1,3
                RNP(I,L)=RNP(I,L) + RPP(K,L)
                RNP(J,L)=RNP(J,L) - RPP(K,L)
322        CONTINUE
321   CONTINUE
C
C
      if(noa(1)+noa(2).ne.0) call pibond
c      if(noa(3)+noa(4).ne.0) call sili_germ
C
      RETURN
C
      END
c

