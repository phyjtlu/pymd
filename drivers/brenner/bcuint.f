C
      SUBROUTINE BCUINT(KL,KI,XX1,XX2,NH,NC,ANSY,ANSY1,ANSY2)
C
C Bicubic spline
C
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'
c
      ANSY=0.0d0
      ANSY1=0.0d0
      ANSY2=0.0d0
C
      IF((KI.EQ.0).OR.(NH.EQ.0).OR.(NC.EQ.0)) THEN
           WRITE(6,*) KI,KL,NH,NC
           include 'close.inc'
           STOP
      ENDIF
C
      DO 32 J=1,16
           X=CLM(KI,NH,NC,J)*
     &       (XX1**IN2(J,1))*(XX2**IN2(J,2))
           ANSY=ANSY+X
           ANSY1=ANSY1+X*IN2(J,1)/XX1
           ANSY2=ANSY2+X*IN2(J,2)/XX2
32    CONTINUE
C
      RETURN
      END

