      subroutine setin

      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'

      DATA KT/2,4*0,1,7*0,3,17*0,4,68*0/
      DATA KT2/6,1,14,32,96*0/
c
      do i=1,ntypes
          xmass(i) = 0.0d0
          noa(i) = 0
          do j=1,ntypes
               sig(i,j) = 0.0d0
               eps(i,j) = 0.0d0
          enddo
      enddo
      xmass(1) = 12.0d0
      xmass(2) = 1.0d0
      xmass(3) = 28.0d0
      xmass(4) = 72.0d0
      PI=ACOS(-1.D0)
      BOLZ=1.380662d0
      EPSI=11604.5D0
      AVO=6.02205
      ECONV=(1.0D0/(AVO*BOLZ*EPSI))*1.0D+07
      return 
      end 

      subroutine setpp

      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'

      if(ipot.eq.1) then
           CALL PARAM
           call mtable
      endif
      if(ilj.ne.0) call ljparam
      call ljcset
      return
      end 

      subroutine setmd 
      IMPLICIT REAL*8(A-H,O-Z)
c
      include 'common_files.inc'

            IF(NMA.ne.0) THEN
                  ENPR=EPSI/FLOAT(NMA)
            ELSE
                 ENPR=0.0D0
            ENDIF
C
            xkea = 0.0d0
            time = 0.0d0
            lchk = 1
c            write(9,1400)
            return 
 1400 FORMAT(8X,'ENERGY(eV)',5X,'T',6X,'TSTEP(fs)  TIME(fs)',/)
            end 
