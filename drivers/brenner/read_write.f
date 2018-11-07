
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_data
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'

      character*24 xyzfilename
      character*2 CATOM
      character*1 svar

C
C READ INPUT DATA
C



      READ(13,*) KVC,MAXKB,KFLAG,nxmol
      READ(13,*) PSEED,RLL,TEM
      READ(13,*) IPOT

c 20 1 5 1  / # steps, # steps between data, thermostat, xmol writes
c 0.6955 0.500 3000.  /Random # seed, neighbor list , temperature(K)
c  1                       / =1 REBO (C,H,Si,Ge), =2 tight-binding for C


      KFLAG=-1
      
c     write out max 50 steps and min. 10 steps
      MAXKB=NINT(KVC/50.)
      if(MAXKB.lt.1) MAXKB=1

      
      TEM=0.

      ktmax = 4
      ilj = 0

1     continue

      READ(13,*,end=2) natom, xma, epst, sigt
c      write(*,*) 'hola',natom, xma, epst, sigt
      if(natom.le.0) go to 1
      ilj = 1
      if(kt(natom).eq.0) then
           ktmax = ktmax + 1
           if(ktmax.gt.ntypes) then
                write(*,*) 'Maximum ntypes of ',ntypes,' exceeded'
                write(*,*) 'Change NTYPES in common_n.inc and recompile'
                include 'close.inc'
                stop
           endif
           kt(natom) = ktmax
           kt2(ktmax) = natom
      endif
      xmass(kt(natom)) = xma
      sig(kt(natom),kt(natom)) = sigt
      eps(kt(natom),kt(natom)) = epst
      go to 1

2     continue


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
 7913 stop 'Error: xyz file not found'
  100 FORMAT(4I6)
  200 FORMAT(4F12.6)
  300 FORMAT(3E20.11)
  350 FORMAT(2I5,3E20.11,I3)
  360 FORMAT(I5,3E20.11)
 1800 FORMAT(20A2)

      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_data1
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
c     WRITE(9,500) HEAD
c     WRITE(9,600) KVC*MAXKB,MAXKB
c     WRITE(9,700) 1,1,1,ECONV
c     WRITE(9,900) TEM,PSEED
c     call flush(9)

      WRITE(6,500) HEAD
      WRITE(6,600) KVC*MAXKB,MAXKB
      WRITE(6,*) ' '
      WRITE(6,*) 
     .     ' Step       time(fs)   Energy/atom(eV)   Energy(eV)'

c      WRITE(6,700) 1,1,1,ECONV
c      WRITE(6,900) TEM,PSEED

      return

  500 FORMAT(/,'* CLASSICAL DYNAMICS SIMULATION OF ',20A2)
  600 FORMAT('TOTAL MD STEPS = ',I6,' DATA WRITTEN EVERY ',I4,
C    &' STEPS WITH ',F7.4,' fs/STEP')
     &' STEPS')
  700 FORMAT(/,'UNITS OF LENGTH, MASS, TIME AND ENERGY:',I2,' A ',I2,
     &' AMU ',I2,' fs ',F8.4,' eV')
  800 FORMAT(/,'NEIGHBOR LIST PARAMETERS: ',F12.3,' A ')
  900 FORMAT(/,'LANGEVIN PARAMETERS: ',F12.3,' k PSEED= ',F9.6,/)
      end 


ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_data2
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'

c
C  CALCULATE KINETIC ENERGY
C
      XX=0.0d0
c      DO J=1,3
c           DO I=1,NP
c                XX=XX+(R1(I,J)**2)*XMASS(KTYPE(I))
c           enddo
c      enddo

c      EETOT=TOTE+XX/(4.0d0*DELTSQ)*ECONV
C
      EETOT=TOTE


c      WRITE(9,1500) LSTEP,EETOT/FLOAT(NP),
c     &   ENPR*XX/6.0d0/DELTSQ*ECONV,DELTA,TIME

      
      WRITE(6,1501) LSTEP,TIME,EETOT/FLOAT(NP),EETOT
      ETOT(LSTEP)=EETOT

c      call flush(9)
c      IF((KFLAG.NE.6).AND.(KFLAG.NE.8)) THEN
c         WRITE(6,*) 'Total Energy (eV): ',TOTE
c      ENDIF

      return


 1500 FORMAT(I6,F14.8,F10.3,2F9.3)
 1501 FORMAT(I6,1X,F14.2,5X,F10.4,5X,F10.4)
      end 


      subroutine write_data3  
      IMPLICIT REAL*8(A-H,O-Z)
      include 'common_files.inc'
c     REWIND 11
      WRITE(93,1800) HEAD
      WRITE(93,100) NP,IDUM,NRA,NLA
      WRITE(93,300) TTIME,DELTA
      WRITE(93,300) (CUBE(N),N=1,3)
      write(93,*) "OUT.xyz"

C
c      DO 11 I=1,NP
c           WRITE(11,350) I,KT2(KTYPE(I)),(R0(I,N),N=1,3),itr(i)
c11    CONTINUE
C
c      DO 12 I=1,NP
c           WRITE(11,360) I,((R1(I,N)/DELTA),N=1,3)
c12    CONTINUE
C
c      DO 13 I=1,NP
c           WRITE(11,360) I,((R2(I,N)),N=1,3)
c13    CONTINUE
C
c      DO 14 I=1,NP
c           WRITE(11,360) I,((R3(I,N)),N=1,3)
c14    CONTINUE
C
c      DO 15 I=1,NP
c           WRITE(11,360) I,((R4(I,N)),N=1,3)
c15    CONTINUE

      return 




  100 FORMAT(4I6)
  200 FORMAT(4F12.6)
  300 FORMAT(3E20.11)
  350 FORMAT(2I5,3E20.11,I3)
  360 FORMAT(I5,3E20.11)
  500 FORMAT('*CLASSICAL DYNAMICS SIMULATION OF ',20A2)
  600 FORMAT(/,'TOTAL STEPS= ',I6,' DATA WRITTEN EVERY ',I4,
C    &' STEPS WITH ',F7.4,' fs/STEP')
     &' STEPS')
  700 FORMAT(/,'UNITS OF LENGTH, MASS, TIME AND ENERGY:',I2,' A ',I2,
     &' AMU ',I2,' fs ',F8.4,' eV')
  800 FORMAT(/,'NEIGHBOR LIST PARAMETERS: ',F12.3,' A ')
  900 FORMAT(/,'LANGEVIN PARAMETERS: ',F12.3,' k PSEED= ',F9.6,/)
 1200 FORMAT('NEIGHBOR LIST UPDATES: ',/)
 1300 FORMAT(10I4,/)
 1400 FORMAT(8X,'ENERGY(eV)',5X,'T',6X,'TSTEP(fs)  TIME(fs)',/)
 1500 FORMAT(I6,F14.8,F10.3,2F9.3)
 1800 FORMAT(20A2)
      end 

