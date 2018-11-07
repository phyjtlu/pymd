C##############################################################C
C                                                              C
C                                                              C
C  THESE PROGRAMS PERFORM MOLECULAR DYNAMICS SIMULATIONS WITH  C
C  BOND-ORDER POTENTIALS FOR HYDROCARBON, SILICON AND          C
C  GERMANIUM; TIGHT BINDING FOR CARBON; AND LENNARD-JONES      C
C  WITH INPUTTED PARAMETERS. THE LATTER POTENTIALS ARE         C 
C  TRUNCATED AT SMALL DISTANCES FOR C-C, H-H, AND C-H PAIRS    C 
C                                                              C
C   Units: mass = AMU's, length = Angstroms, time = fs         C
C          energy = eV's                                       C
C                                                              C
C   DOCUMENTATION CAN BE FOUND IN:                             C
C      /MD/Documentation                                       C 
C                                                              C
C##############################################################C
c


      subroutine brennerf(npfromC,CellfromC,AnfromC,RfromC,
     .                   Ea2C,F2C,Etot2C) 


      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'common_files.inc'

      integer,intent(in) :: npfromC
      integer,intent(in) ::  AnfromC(npfromC)
      real*8,intent(in) ::  CellfromC(3)
      real*8,intent(in) ::  RfromC(3*npfromC)
      real*8,intent(out) ::  Ea2C(npfromC)
      real*8,intent(out) ::  Etot2C
      real*8,intent(out) ::  F2C(3*npfromC)


      real*8 epsibug
      parameter(epsibug=0.000000000001)
      logical firsttime 
      save firsttime
      data firsttime / .true. /

c 
c************************* 
c set up and input data  * 
c*************************
c 

c      firsttime=.true.
c open input/output files 
      if(firsttime) then
c     some LJ stuff:
         open(13,file='./datafiles/input.d',status='old')
c     
c     data files for REBO potential 
C     
         open(14,file='./datafiles/inter3d_iv_new.d',
     &        status='old')
         open(15,file='./datafiles/inter2d_iv.d',
     &status='old')
         open(16,file='./datafiles/inter3dtors.d',
     &        status='old')
         open(17,file='./datafiles/inter3d_h.d',
     &        status='old')
         open(18,file='./datafiles/inter3d_ch.d',
     &        status='old')
      end if


      np=npfromC
      
      if(np.gt.npmax) then
           write(*,*) 'np= ',np,' greater than npmax= ',npmax
           write(*,*) 'increase npmax and recompile'
           include 'close.inc'
           stop
      endif


c debug:
c      write(*,*) firsttime
c      write(*,*) "npfromC: ",npfromC
c      write(*,*) "CellfromC: ",CellfromC
c      write(*,*) "RfromC: "
c      do i=0,npfromC-1
c         write(*,*) i+1,AnfromC(i+1),(RfromC(i*3+j),j=1,3)
c      end do
c      write(*,*) 'hej'




c initialize 
      if(firsttime) call setin 
C read input data 

      cube(1)=CellfromC(1)
      cube(2)=CellfromC(2)
      cube(3)=CellfromC(3)

      nma = 0
      nta = 0
      DO I=1,NP
         itr(i)=1
         NATOM=AnfromC(i)         
         do N=1,3
c there is a bug somewhere which cause caguts to return NANs when
c R0 has exact zero elements ... uh..uhh..
            R0(I,N)=RfromC(3*(I-1)+N)+epsibug

         
            if(R0(I,N).gt.cube(n)/2.0d0) R0(I,N)=R0(I,N)-cube(n)
         enddo
         KTYPE(I)=KT(NATOM)
         if(ktype(i).eq.0) then
            write(*,*) 'unknown atom type for atom ',k
            include 'close.inc'
            stop
         endif
         noa(ktype(i)) = noa(ktype(i)) + 1
         if(itr(i).ne.2) then
            nma = nma + 1
            mlist(nma) = i
            if(itr(i).eq.1) then
               nta = nta + 1
               nlist(nta) = i
            endif
         endif
      ENDDO

C
C ESTABLISH INITIAL POSITIONS FOR NEIGHBOR LIST UPDATE
C
      DO 8 I=1,3
         DO 7 J=1,NP
            R0L(J,I)=R0(J,I)
 7       CONTINUE
 8    CONTINUE
C
      VOL=CUBE(1)*CUBE(2)*CUBE(3)
C
      DO 6 I=1,3
         CUBE2(I)=CUBE(I)/2.0D0
 6    CONTINUE
C     

      TTCONV=2.0d0/3.0d0/FLOAT(NP)

      if(firsttime) call read_data 
C setup potential parameters  

      if(firsttime) call setpp 
c setup predictor-corrector coefficients
c write out data headers 
c      call write_data1 
c     calculate energy and forces 

      LCHK=1                    ! setup neighbor list
      CALL MODEL 

     
c forces
            DO J=1,3
               DO I=1,NP
c                  write(*,*) j,i,rnp(i,j),  F2c((i-1)*3 + j)
                  F2c((i-1)*3 + j)=rnp(i,j)
               ENDDO
            ENDDO


c energy
            esum=0.0
            etot2c=TOTE
            DO I=1,NP
               ea2c(i)=eatom(i)
               esum=esum+eatom(i)
c                   write(*,*) "f:eatom i ",i,eatom(i)
            ENDDO
c            write(*,*) "f:",esum,TOTE
       



      if(firsttime) then
         close(13)
         close(14)
         close(15)
         close(16)
         close(17)
         close(18)
      end if



      firsttime=.false.
      return
      END


C
C add included subroutines
C
      include 'subroutines.inc' 
