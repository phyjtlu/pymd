      SUBROUTINE PARAM
C
C Parameters for hydrocarbons
C
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'
C
      data spgc/
     &  0.2817216000000E+00, 0.1062912000000E+01, 0.2136736000000E+01,
     &  0.2533952000000E+01, 0.1554736000000E+01, 0.3863296000000E+00,
     &  0.2817216000000E+00, 0.1062912000000E+01, 0.2136736000000E+01,
     &  0.2533952000000E+01, 0.1554736000000E+01, 0.3863296000000E+00,
     &  0.6900668660000E+00, 0.5460691360000E+01, 0.2301345680000E+02,
     &  0.5491519344000E+02, 0.6862037040000E+02, 0.3470897779200E+02,
     &  0.3754490870000E+00, 0.1407252749388E+01, 0.2255103926323E+01,
     &  0.2028902219952E+01, 0.1426981217906E+01, 0.5063107994308E+00,
     &  0.2718558000000E+00, 0.4892727456293E+00,-0.4328199017473E+00,
     & -0.5616795197048E+00, 0.1270874966906E+01,-0.3750409108350E-01/
C
      DATA IGC/16*4,2*3,2*2,5*1/
C
      DATA SPGH/270.467795364007301,1549.701314596994564
     &,3781.927258631323866,4582.337619544424228,2721.538161662818368,
     &630.658598136730774,16.956325544514659,-21.059084522755980,
     &-102.394184748124742,-210.527926707779059,-229.759473570467513,
     &-94.968528666251945,19.065031149937783,2.017732531534021,
     &-2.566444502991983,3.291353893907436,-2.653536801884563,
     &0.837650930130006/
C
      DATA IGH/18*3,4*2,3*1/
C
      DATA XH/200*0.0/,XH1/200*0.0/,XH2/200*0.0/
C
      DATA ATT/3.20D0/,XQM/3.70D0/

C
C Zero bicubic spline coefficients
C
      DO 125 I=1,2
           DO 124 L=1,10
                DO 123 M=1,10
                     DO 122 J=1,16
                          CLM(I,L,M,J)=0.0
122                  CONTINUE
123             CONTINUE
124        CONTINUE
125   CONTINUE
C
C bicubic spline
C
      READ(15,*) I2D
C
C read integer values of bicubic spline

C
  119 CONTINUE
      READ(15,*) I,J,K,XHH
      IF(I.LE.0) GO TO 121
      XH(I,J,K)=XHH
      GO TO 119
  121 CONTINUE
C
C
      XH1(2,3,1)=(XH(2,4,1)-XH(2,2,1))/2.0D0
      XH1(2,2,2)=(XH(2,3,2)-XH(2,1,2))/2.0D0
C
      XH2(2,2,2)=(XH(2,2,3)-XH(2,2,1))/2.0D0
      XH2(2,1,3)=(XH(2,1,4)-XH(2,1,2))/2.0D0
C
C Read bicubic spline coefficients
C
  126 CONTINUE
      READ(15,*,END=127) I,L,M
      READ(15,*) (CLM(I,L,M,J),J=1,16)
      GO TO 126
  127 CONTINUE
C
      IC=0
C
      DO 9 I=1,4
           DO 8 J=1,4
                     IC=IC+1
                     IN2(IC,1)=I-1
                     IN2(IC,2)=J-1
8               CONTINUE
9         CONTINUE
C
C Read tricubic spline coefficients
C
      IC=0
C
      DO 7 I=1,4
           DO 6 J=1,4
                DO 5 K=1,4
                     IC=IC+1
                     IN3(IC,1)=I-1
                     IN3(IC,2)=J-1
                     IN3(IC,3)=K-1
5               CONTINUE
6         CONTINUE
7     CONTINUE
C
      READ(14,*) I3D
  129 CONTINUE
      READ(14,*,END=130) L,M,N
      READ(14,*) (CLMN(1,L,M,N,I),I=1,64)
      GO TO 129
  130 CONTINUE
C
      READ(17,*) I3D
  229 CONTINUE


      READ(17,*,END=230) L,M,N
      READ(17,*) (CLMN(3,L,M,N,I),I=1,64)
      GO TO 229
  230 CONTINUE
C
      READ(18,*) I3D
  239 CONTINUE
      READ(18,*,END=240) L,M,N
      READ(18,*) (CLMN(2,L,M,N,I),I=1,64)
      GO TO 239
  240 CONTINUE
C
      DO 134 L=1,10
           DO 133 M=1,10
                DO 132 I=1,64
                     CLMN(1,L,M,10,I)=CLMN(1,L,M,9,I)
                     CLMN(2,L,M,10,I)=CLMN(2,L,M,9,I)
                     DO 131 N=6,10
                          CLMN(3,L,M,N,I)=CLMN(3,L,M,5,I)
131                  CONTINUE
132             CONTINUE
133        CONTINUE
134   CONTINUE
C
C Read tricubic spline coefficients for torsional potential
C
      READ(16,*) ITD
  135 CONTINUE
      READ(16,*,END=136) L,M,N
      READ(16,*) (TLMN(L,M,N,I),I=1,64)
      GO TO 135
  136 CONTINUE
C
      DO 140 L=1,10
           DO 139 M=1,10
                DO 138 N=4,10
                     DO 137 I=1,64
                          TLMN(L,M,N,I)=TLMN(L,M,3,I)
137                  CONTINUE
138             CONTINUE
139        CONTINUE
140   CONTINUE
C
      IF((ITD.NE.I2D).OR.(ITD.NE.I3D)) THEN
            WRITE(6,*) 'INCOMPATABLE POTENTIAL TYPES'
            include 'close.inc'
            STOP
      ENDIF

C
      PQ=PI/(XQM-ATT)

      do i=1,4
          do j=1,4
               AD(i,j)  = 0.0d0
               AXL(i,j) = 0.0d0
               BD(i,j)  = 0.0d0



               BXL(i,j) = 0.0d0
               CD(i,j)  = 0.0d0
               CXL(i,j) = 0.0d0
               DD(i,j)  = 0.0d0
               DXL(i,j) = 0.0d0
               ED(i,j)  = 0.0d0
               RB1(i,j) = 0.0d0
               RB2(i,j) = 0.0d0
               RMAX(i,j) = 0.0d0
               PID(i,j) = 1.0d0
               CHI(i,j) = 1.0d0
               do k = 1,4
                xdb(i,j,k) = 0.0d0
               enddo
          enddo
      enddo

C
C*** IMPORTANT***************
C                           *
C TO INCLUDE DIHEDRAL TERMS *
C SET NDIHED=2, OTHERWISE   *
C SET NDIHED=10             *
C                           *
C****************************
C
      NDIHED=2
C
C CARBON
C
      AD(1,1)=12388.79197798375D0
      AXL(1,1)=4.720452312717397D0
      BD(1,1)=17.56740646508968D0
      BXL(1,1)=1.433213249951261D0
      CD(1,1)=30.71493208065162D0
      CXL(1,1)=1.382691250599169D0
      DD(1,1)=10953.54416216992D0
      DXL(1,1)=4.746539060659529D0
      ED(1,1)=0.3134602960832605d0
      RB1(1,1)=1.7d0
      RB2(1,1)=2.0d0
      RMAX(1,1)=RB2(1,1)
      PID(1,1)=PI/(RB2(1,1)-RB1(1,1))
C
C Hydrogen
C
      AD(2,2)=29.6325931D0
      AXL(2,2)=1.715892169856421D0
      BD(2,2)=0.0D0
      BXL(2,2)=1.0D0
      CD(2,2)=0.0D0
      CXL(2,2)=1.0D0
      DD(2,2)=32.81735574722296D0
      DXL(2,2)=3.536298648376465D0
      ED(2,2)=0.3704714870452888d0
c 
      RB1(2,2)=1.10d0


      RB2(2,2)=1.70d0
      RMAX(2,2)=RB2(2,2)
      PID(2,2)=PI/(RB2(2,2)-RB1(2,2))
C
c CARBON-HYDROGEN
C
       AD(2,1)=32.35518665873256
       AXL(2,1)=1.434458059249837
       DD(2,1)=149.9409872288120
       DXL(2,1)= 4.102549828548784
       ED(2,1)=0.3407757282257080
c
      BD(2,1)=0.0D0
      BXL(2,1)=1.0D0
      CD(2,1)=0.0D0
      CXL(2,1)=1.0D0
C
      AD(1,2)=AD(2,1)
      AXL(1,2)=AXL(2,1)
      BD(1,2)=BD(2,1)
      BXL(1,2)=BXL(2,1)
      CD(1,2)=CD(2,1)
      CXL(1,2)=CXL(2,1)
      DD(1,2)=DD(2,1)
      DXL(1,2)=DXL(2,1)
      ED(1,2)=ED(2,1)
C
      RB1(2,1)=1.3d0
      RB2(2,1)=1.8d0
      RMAX(2,1)=RB2(2,1)
      PID(2,1)=PI/(RB2(2,1)-RB1(2,1))
      PIDT=PI/0.30D0
C
      RB1(1,2)=RB1(2,1)
      RB2(1,2)=RB2(2,1)
      RMAX(1,2)=RB2(1,2)
      PID(1,2)=PI/(RB2(1,2)-RB1(1,2))
C
      DO 12 I=1,2
           DO 11 J=1,2
                DO 10 K=1,2
                     XDB(I,J,K)=0.0d0
                     REG(I,J,K)=1.0d0
10              CONTINUE
11         CONTINUE
12    CONTINUE 
C
      XXDB=4.0D0
      RHH=0.7415886997d0
      RCH=1.09d0
      XDB(2,2,2)=4.0D0
C
      XDB(2,1,2)=4.0D0
      XDB(2,2,1)=4.0D0
C
      XDB(2,1,1)=4.0D0
      XDB(1,2,1)=0.0D0
      XDB(1,2,2)=0.0d0
C
      REG(2,1,2)=EXP(XDB(2,1,2)*(RHH-RCH))
      REG(2,2,1)=EXP(XDB(2,2,1)*(RCH-RHH))
C
C TERSOFF-III SILICON
C
      DXL(3,3)=2.4799d0
      AXL(3,3)=1.7322d0
      DD(3,3)=1830.8d0
      AD(3,3)=471.18d0
      XTN2(3)=0.78734d0
      XTN1(3)=1/(2.0*XTN2(3))
      ADB(3)=1.0999D-6
      CDB(3)=1.0039D+5
      CDB2(3) = CDB(3)*cdb(3)
      DDB(3)=16.218d0
      DDB2(3) = DDB(3)*DDB(3)
      HDB(3)=0.59826d0
      RB1(3,3)=2.7d0
      RB2(3,3)=3.0d0
      RMAX(3,3)=RB2(3,3)
      PID(3,3)=PI/(RB2(3,3)-RB1(3,3))
      CHI(3,3) = 1.0d0
C
C TERSOFF GERMANIUM
C
      DXL(4,4) = 2.4451d0
      AXL(4,4) = 1.7047d0
      DD(4,4) = 1769.0d0
      AD(4,4) = 0.50d0 * 419.23d0
      XTN2(4) = 0.75627d0
      XTN1(4) = 1/(2.0*XTN2(4))
      ADB(4) = 9.0166d-07
      CDB(4) = 1.0643d+5
      CDB2(4) = CDB(4)*CDB(4)
      DDB(4) = 15.652D0
      DDB2(4) = DDB(4)*DDB(4)
      HDB(4) = -0.43884d0
      RB1(4,4)=2.7d0
      RB2(4,4)=3.0d0
      RMAX(4,4)=RB2(4,4)
      PID(4,4)=PI/(RB2(4,4)-RB1(4,4))
      CHI(4,4) = 1.0d0
C
C Mixed SILICON-GERMANIUM
C
      DXL(4,3) = (DXL(4,4)+DXL(3,3))/2.0d0
      DXL(3,4) = DXL(4,3)
      AXL(4,3) = (AXL(4,4)+AXL(3,3))/2.0d0
      AXL(3,4) = AXL(4,3)
      DD(4,3) = sqrt(DD(4,4)*DD(3,3))
      DD(3,4) = DD(4,3)
      AD(4,3) = sqrt(AD(3,3)*AD(4,4))
      AD(3,4) = AD(4,3)
      RB1(4,3) = sqrt(RB1(3,3)*RB1(4,4))
      RB1(3,4) = RB1(4,3)
      RB2(4,3) = sqrt(RB2(3,3)*RB2(4,4))
      RB2(3,4) = RB2(4,3)
      RMAX(4,3) = RB2(4,3)
      RMAX(3,4) = RMAX(4,3)
      PID(4,3) = PI/(RB2(4,3)-RB1(4,3))
      PID(3,4) = PID(4,3)
      CHI(4,3) = 1.00061d0
      CHI(3,4) = CHI(4,3)
c
      SIGMA=1.0d0
      EPSI=11605.0d0
C
      DO 61 I=1,4
           DO 60 J=1,4
                RLIST(I,J)=(RMAX(I,J)+RLL)**2
                RMAX(I,J)=RMAX(I,J)**2
   60      CONTINUE
   61 CONTINUE
      RETURN
110   FORMAT(4I5)
120   FORMAT(4E20.6)
      END

