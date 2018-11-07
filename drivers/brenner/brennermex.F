#include "fintrf.h"
C
#if 0
C     
C     brennerf.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     brenner.f
C     Inputs: cell = vector with unitcell dimensions, 
C     R = list of x1,y1,z1,x2,y2,z2,..., eain = (initial energy for each 
C     atom) list with NA zeros, Fin = initial forces (vector with 3*NA 0's)
C     Outputs: F (force) list (F1x,F1y,F1z,F2x,F2y,...),ea vector (energy),
C     etot = total energy.
      
C     This is a MEX-file for MATLAB.
C     Tue Gunst, 2010.

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C
      mwpointer plhs(*), prhs(*)
      mwpointer mxCreateDoubleMatrix, mxGetPr
      mwpointer cell_pr, an_pr, R_pr, eain_pr, Fin_pr, etot_pr
      mwpointer F_pr, ea_pr

C-----------------------------------------------------------------------
C

      integer nlhs, nrhs, mxIsNumeric
      mwsize mxGetM, mxGetN
      mwsize m, n, size
      integer NA, np, an(1000)
      real*8  cell(3), R(3000), eain(1000), Fin(3000)
      real*8  anapprox(1000), etot, F(3000), ea(1000)

C     Check for proper number of arguments. 
      if(nrhs .lt. 4) then
         call mexErrMsgTxt('Four inputs required.')
      elseif(nlhs .lt. 3) then
         call mexErrMsgTxt('Three outputs required.')
      endif

C     Get the size of the input array (m rows, n collums).
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      size = m*n
      np = mxGetN(prhs(3))


C     Column * row should be smaller than 3000
      if(size.gt.3000) then
         call mexErrMsgTxt('Row * column must be <= 3000.')
      endif
      
C     Check to insure the array is numeric (not strings).
      if(mxIsNumeric(prhs(4)) .eq. 0) then
         call mexErrMsgTxt('Input must be a numeric array.')
      endif

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(m,n,0)
      plhs(2) = mxCreateDoubleMatrix(1,np,0)
      plhs(3) = mxCreateDoubleMatrix(1,1,0)

C      NA_pr = mxGetPr(prhs(1))
      cell_pr = mxGetPr(prhs(1))
C      an_pr = mxGetPr(prhs(2))
      R_pr = mxGetPr(prhs(2))
      eain_pr = mxGetPr(prhs(3))
      Fin_pr = mxGetPr(prhs(4))

      F_pr = mxGetPr(plhs(1))
      ea_pr = mxGetPr(plhs(2))
      etot_pr = mxGetPr(plhs(3))

C     Load the data into Fortran arrays.
C      call mxCopyPtrToInteger1(NA_pr,NA,1)
      NA=int(np)
      call mxCopyPtrToReal8(cell_pr,cell,3)
C      call mxCopyPtrToInteger1(an_pr,anapprox,np)
      an(:) = 6
C      an = int(anapprox)
      call mxCopyPtrToReal8(R_pr,R,size)
      call mxCopyPtrToReal8(eain_pr,eain,np)
      call mxCopyPtrToReal8(Fin_pr,Fin,size)

C     Call the computational subroutine. 
      call brennerf(NA,cell,an,R, eain,Fin,etotTest)
C      write(56,*) "Your favorite integer is: %i", np(1)

C     Load the data into the outputs to MATLAB
      call mxCopyReal8ToPtr(Fin,F_pr,size)
      call mxCopyReal8ToPtr(eain,ea_pr,np)
      call mxCopyReal8ToPtr(etotTest,etot_pr,1)

      return
      end

      subroutine brennerf(npfromC,CellfromC,AnfromC,RfromC,
     .                   Ea2C,F2C,Etot2C) 


      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'common_files.inc'

      integer AnfromC(npfromC)
      real*8 CellfromC(3)
      real*8 RfromC(3*npfromC)
      real*8 Ea2C(npfromC)
      real*8 Etot2C
      real*8 F2C(3*npfromC)


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
      write(55,*) firsttime
      write(55,*) "npfromC: ",npfromC
      write(55,*) "CellfromC: ",CellfromC
      write(55,*) "RfromC: "
      do i=0,npfromC-1
         write(55,*) i+1,AnfromC(i+1),(RfromC(i*3+j),j=1,3)
      end do
      write(55,*) 'hej'
      close(55)

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
