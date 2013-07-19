!-------------------------------------------------------
!
!     USE THIS FILE TO COMBINE THE SNAP_X.XXX FILES
!     INTO GLOB.XXX FILES.
!
!     SET THE NUMBER OF PROCESSORS
!     USED AND THE SIZE OF THE ARRAYS
!
!     INPUT: nfile = number of files to be read
!            startfile = file to start with
!            fileincrement = the increment between files
!
!-------------------------------------------------------
!
!---- MODULE FOR VARIABLES
      MODULE var_init
      IMPLICIT NONE
      SAVE
!--- inputs
      INTEGER :: global_NX,global_NY,global_NZ
      INTEGER :: locNX,locNY,locNZ,numprocs
      INTEGER :: myid,proc_coords(3),proc_dims(3)
      DOUBLE PRECISION :: ZEIT,DELT
      INTEGER :: num_iter,cartN,NCOSYS
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: UPPERBND
      logical :: restart
      integer :: tmpflag
      DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793D0
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XXA,XXB,DXA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XYA,XYB,DYA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XZA,XZB,DZA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: V,G,H
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: RH,T,ER
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: PG,CS
      DOUBLE PRECISION :: TinitBottom,Tinittop,rhinittop,rhinitbottom
      DOUBLE PRECISION :: hydrotime,Maccreated
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: XKAPR,flimx
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: EdepFor,EdepAdv
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: EdepRad
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: xkapA,xkapP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: STELLARINPUT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: delta_ER
!--- outputs
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: global_XXA,global_XXB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: global_DXA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: global_XYA,global_XYB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: global_DYA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: global_XZA,global_XZB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: global_DZA

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_V,global_G
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_H
      
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_RH
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_T
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_ER
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_PG
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_CS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_xkapR
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_flimx
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_EdepFor
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_EdepAdv
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_EdepRad
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_xkapA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: 
     $     global_stellarinput
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_delta_ER
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_xkapP


      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: var1,var2,var3
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: threeDUPPERBND
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: global_upperbnd

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: cart_T,cart_RH
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: cart_VX,cart_VY
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: cart_VZ
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: cart_EdepFor
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: cart_EdepAdv
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: cart_EdepRad
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: cartvar1,cartvar2,cartvar3
      END MODULE var_init

      MODULE var_alloc
      USE var_init
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatearrays(startfile)
          IMPLICIT NONE
          integer :: startfile,test
          character fnum*3
          character filename*12
          print *,'restart',restart
          if(restart) then
             write(fnum,'(i3.3)') startfile
             filename = 'RSRT_000.'//fnum
c             print *,'filename test: ',filename,fnum
             OPEN (12,FILE=filename,STATUS='OLD',
     $            FORM='UNFORMATTED')
          else
             write(fnum,'(i3.3)') startfile
             filename = 'SNAP_000.'//fnum
c             print *,'filename test: ',filename,' fnum=',fnum
             OPEN (12,FILE=filename,STATUS='OLD',
     $            FORM='UNFORMATTED')
          endif
          rewind(12)
          READ (12) global_NX,global_NY,global_NZ
          READ (12) locNX,locNY,locNZ,numprocs
          CLOSE(12)
          print *, 'Global array sizes',global_NX,global_NY,global_NZ
          print *, 'Local  array sizes',locNX,locNY,locNZ
          print *, 'Number of processors',numprocs

          test  = INT(numprocs**(1.d0/2.d0))
          if(mod(global_NY,test).ne.0.or.
     %         mod(global_NZ,test).ne.0) then
             print *,'locNY or locNZ does not divide evenly'
             print *,locNY,locNZ,numprocs
             stop
          endif

!--- input arrays
          ALLOCATE(XXA(0:locNX+2),XXB(-1:locNX+2),DXA(locNX+2))
          ALLOCATE(XYA(0:locNY+2),XYB(-1:locNY+2),DYA(locNY+2))
          ALLOCATE(XZA(0:locNZ+2),XZB(-1:locNZ+2),DZA(locNZ+2))
          ALLOCATE(V(0:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(G(-1:locNX+2,0:locNY+2,-1:locNZ+2))
          ALLOCATE(H(-1:locNX+2,-1:locNY+2,0:locNZ+2))
          ALLOCATE(RH(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(T(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(ER(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(PG(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(CS(0:locNX+1,0:locNY+1,0:locNZ+1))
          ALLOCATE(XKAPR(1:locNX,1:locNY,1:locNZ))
          ALLOCATE(FLIMX(1:locNX,1:locNY,1:locNZ))
          ALLOCATE(EdepFor(1:locNX,1:locNY,1:locNZ))
          ALLOCATE(EdepAdv(1:locNX,1:locNY,1:locNZ))
          ALLOCATE(EdepRad(1:locNX,1:locNY,1:locNZ))
          ALLOCATE(upperbnd(locNY,locNZ))
          ALLOCATE(xkapA(1:locNX,1:locNY,1:locNZ))
          ALLOCATE(stellarinput(1:locNX,1:locNY,1:locNZ))
          ALLOCATE(delta_ER(1:locNX,1:locNY,1:locNZ))
          ALLOCATE(xkapP(1:locNX,1:locNY,1:locNZ))
!--- IDL combination arrays
          ALLOCATE(global_XXA(global_NX),global_XXB(global_NX))
          ALLOCATE(global_DXA(global_NX))
          ALLOCATE(global_XYA(global_NY),global_XYB(global_NY))
          ALLOCATE(global_DYA(global_NY))
          ALLOCATE(global_XZA(global_NZ),global_XZB(global_NZ))
          ALLOCATE(global_DZA(global_NZ))
          ALLOCATE(global_V(global_NX,global_NY,global_NZ))
          ALLOCATE(global_G(global_NX,global_NY,global_NZ))
          ALLOCATE(global_H(global_NX,global_NY,global_NZ))
          ALLOCATE(global_RH(global_NX,global_NY,global_NZ))
          ALLOCATE(global_T(global_NX,global_NY,global_NZ))
          ALLOCATE(global_ER(global_NX,global_NY,global_NZ))
          ALLOCATE(global_PG(global_NX,global_NY,global_NZ))
          ALLOCATE(global_CS(global_NX,global_NY,global_NZ))
          ALLOCATE(global_upperbnd(global_NY,global_NZ))
          ALLOCATE(global_xkapR(global_NX,global_NY,global_NZ))
          ALLOCATE(global_flimx(global_NX,global_NY,global_NZ))
          ALLOCATE(global_EdepFor(global_NX,global_NY,global_NZ))
          ALLOCATE(global_EdepAdv(global_NX,global_NY,global_NZ))
          ALLOCATE(global_EdepRad(global_NX,global_NY,global_NZ))
          ALLOCATE(global_xkapA(global_NX,global_NY,global_NZ))
          ALLOCATE(global_stellarinput(global_NX,global_NY,global_NZ))
          ALLOCATE(global_delta_ER(global_NX,global_NY,global_NZ))
          ALLOCATE(global_xkapP(global_NX,global_NY,global_NZ))
!--- IFRIT ARRAYS
          ALLOCATE(var1(global_NX,global_NY,global_NZ))
          ALLOCATE(var2(global_NX,global_NY,global_NZ))
          ALLOCATE(var3(global_NX,global_NY,global_NZ))
          ALLOCATE(threeDUPPERBND(global_NX,global_NY,global_NZ))
!--- CARTESIAN MAP ARRAYS
          cartN = 100
          ALLOCATE(cartvar1(cartN,cartN,cartN))
          ALLOCATE(cartvar2(cartN,cartN,cartN))
          ALLOCATE(cartvar3(cartN,cartN,cartN))
          ALLOCATE(cart_T(cartN,cartN,cartN))
          ALLOCATE(cart_RH(cartN,cartN,cartN))
          ALLOCATE(cart_VX(cartN,cartN,cartN))
          ALLOCATE(cart_VY(cartN,cartN,cartN))
          ALLOCATE(cart_VZ(cartN,cartN,cartN))
          ALLOCATE(cart_EdepFor(cartN,cartN,cartN))
          ALLOCATE(cart_EdepAdv(cartN,cartN,cartN))
          ALLOCATE(cart_EdepRad(cartN,cartN,cartN))
      END SUBROUTINE allocatearrays
      END MODULE var_alloc
!---------------------
!
!--- BEGIN PROGRAM
      PROGRAM ppcombine
      USE var_init
      USE var_alloc
      IMPLICIT NONE
      integer :: i,j,k,l,proc,file,startfile,nfiles,fileincrement
      character fnum*3
      character proc_num*3
      character filename*12
      character endname*9
      character IDLoutfile*8
      character IDLoutfile2*13
      character IFRIToutfile*12
      logical :: readrestart,IDLfiles,IFRITfiles
      logical :: EXCLUDEREGION,mapspheretocart,mapcyltocart
!--------------------------------------------------------
!- qqqq
      nfiles=2
      startfile=0
      fileincrement = 1

!-- flag for new write out formats
c      tmpflag = 0  !- oldversions
      tmpflag = 1  !- includes heating+delta_ER

c-- flag to use RSRT_files
c      restart = .true. 
      restart = .false.

c-- flag to read AKTTD
c      readrestart = .true.
      readrestart = .false.

      IDLfiles = .true.
c      IDLfiles = .false.

c      IFRITfiles = .true.
       IFRITfiles = .false.

c      EXCLUDEREGION = .true.
      EXCLUDEREGION = .false.

      mapspheretocart = .true.
c      mapspheretocart = .false.

c      mapcyltocart = .true.
      mapcyltocart = .false.

!--------------------------------------------------------
      call allocatearrays(startfile)

      DO file=startfile,nfiles*fileincrement-fileincrement+startfile,
     &     fileincrement
         if (file.gt.999) then
            print *,'BAD FILE NUM'
            stop
         endif
         write(fnum,'(i3.3)') file
         DO proc=1,numprocs
!---  Read in SNAP files
            write(proc_num,'(i3.3)') proc-1
            filename = 'SNAP_'//TRIM(ADJUSTL(proc_num))//'.'//fnum
            if(restart.and.file.eq.startfile) then
               filename = 'RSRT_'//TRIM(ADJUSTL(proc_num))//'.'//fnum
            endif
            OPEN (12,FILE=filename,STATUS='OLD',FORM='UNFORMATTED')
            rewind(12)
            CALL READFILE
            CLOSE(12)
!---  Fill Global Arrays
            CALL FILLGLOBAL(EXCLUDEREGION)
         ENDDO



!---  Print the GLOB.xxx files
         if(IDLfiles) then
c-use this to offset file: write(fnum,'(i3.3)') file-1
            if(restart.and.file.eq.startfile) then
               idloutfile2 = 'GLOB_RSRT.'//fnum
               print *,'output: ',IDLoutfile2
               OPEN(12,FILE=IDLoutfile2,STATUS='UNKNOWN',
     %              FORM='UNFORMATTED')
            else
               IDLoutfile = 'GLOB.'//fnum
               print *,'output: ',IDLoutfile
               OPEN(12,FILE=IDLoutfile,STATUS='UNKNOWN',
     %              FORM='UNFORMATTED')
            endif
            rewind(12)
            CALL PRINTGLOBAL
            CLOSE(12)
         endif

         
!---  Print the IFRITxxx.bin. files
         if(IFRITfiles) then

            if (mapspheretocart) then
               CALL SPHERETOCARTESIAN(EXCLUDEREGION)
            endif



            IFRIToutfile = 'IFRIT'//fnum//'.bin'
            print *,'output: ',IFRIToutfile
            OPEN(1,FILE=IFRIToutfile,STATUS='UNKNOWN',
     %           FORM='UNFORMATTED')
            rewind(1)

            CALL CALCUALTESSTUFF

            if(mapspheretocart) then
               cartvar1 = cart_T
               cartvar2 = cart_RH
               cartvar3 = cart_VX

               WRITE(1) cartN,cartN,cartN
               WRITE(1) (((cartvar1(I,J,K),I=1,cartN),J=1,cartN),
     %              K=1,cartN)
               WRITE(1) (((cartvar2(I,J,K),I=1,cartN),J=1,cartN),
     %              K=1,cartN)
               WRITE(1) (((cartvar3(I,J,K),I=1,cartN),J=1,cartN),
     %              K=1,cartN)
            else
               var1 = global_T
               var2 = global_RH
               var3 = global_H
               WRITE(1) global_NX,global_NY,global_NZ
               WRITE(1) (((var1(I,J,K),I=1,global_NX),J=1,global_NY),
     %              K=1,global_NZ)
               WRITE(1) (((var2(I,J,K),I=1,global_NX),J=1,global_NY),
     %              K=1,global_NZ)
               WRITE(1) (((var3(I,J,K),I=1,global_NX),J=1,global_NY),
     %              K=1,global_NZ)
            endif
            CLOSE(1) 
         endif
      ENDDO

      IF(readrestart) then
         DO proc=1,numprocs
            write(proc_num,'(i3.3)') proc-1
            endname = 'AKTTD_'//TRIM(ADJUSTL(proc_num))
            OPEN (12,FILE=endname,STATUS='OLD',FORM='UNFORMATTED')
            rewind(12)
            CALL READFILE
            close(12)
            CALL FILLGLOBAL(.false.)
         ENDDO
         OPEN (12,FILE='AKTTD',STATUS='UNKNOWN',FORM='UNFORMATTED')
         rewind(12)
         print *,'output: AKTTD'
         CALL PRINTGLOBAL
         CLOSE(12)
      ENDIF
      END

      SUBROUTINE CALCUALTESSTUFF
      USE var_init
      IMPLICIT NONE
      integer :: i,j,k
!-- transfer the upperbound into a 3d array
      DO K=1,global_NZ            
         DO J=1,global_NY
            DO I=1,global_NX-1
               if(global_upperbnd(j,k).eq.i) then
                  threeDUPPERBND(i,j,k) = 1.d0
               else
                  threeDUPPERBND(i,j,k) = 0.d0
               endif
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CALCUALTESSTUFF


      SUBROUTINE FILLGLOBAL(EXCLUDEREGION)
      USE var_init
      IMPLICIT NONE
      integer :: i,j,k,globali,globalj,globalk
      logical :: EXCLUDEREGION
      if(proc_coords(2).eq.0.and.proc_coords(3).eq.0) then
         DO i=1,locNX
            globali = i+locNX*proc_coords(1)
            global_XXA(globali) = XXA(I)
            global_XXB(globali) = XXB(I)
            global_DXA(globali) = DXA(I)
         enddo
      endif
      IF(proc_coords(1).eq.0.and.proc_coords(3).eq.0) then
         DO j=1,locNY
            globalj = j+locNY*proc_coords(2)
            global_XYA(globalj) = XYA(J)
            global_XYB(globalj) = XYB(J)
            global_DYA(globalj) = DYA(J)
         ENDDO
      ENDIF
      IF(proc_coords(1).eq.0.and.proc_coords(2).eq.0) then
         DO k=1,locNZ
            globalk = k+locNZ*proc_coords(3)
            global_XZA(globalk) = XZA(k)
            global_XZB(globalk) = XZB(k)
            global_DZA(globalk) = DZA(k)
         ENDDO
      ENDIF
      DO k=1,locNZ
         globalk = k+locNZ*proc_coords(3)
         DO j=1,locNY
            globalj = j+locNY*proc_coords(2)
            DO i=1,locNX
               globali = i+locNX*proc_coords(1)
               global_V(globali,globalj,globalk) = V(I,J,K)
               global_G(globali,globalj,globalk) = G(I,J,K)
               global_H(globali,globalj,globalk) = H(I,J,K)
               global_RH(globali,globalj,globalk) = RH(I,J,K)
               global_T(globali,globalj,globalk) = T(I,J,K)
               if(tmpflag.eq.1) then
                  global_ER(globali,globalj,globalk) = ER(I,J,K)
               else
                  global_ER(globali,globalj,globalk) = 0.d0
               endif
               global_PG(globali,globalj,globalk) = PG(I,J,K)
               global_CS(globali,globalj,globalk) = CS(I,J,K)
               global_xkapR(globali,globalj,globalk) = xkapR(I,J,K)
               global_flimx(globali,globalj,globalk) = flimx(I,J,K)
               global_EdepFor(globali,globalj,globalk) = EdepFor(I,J,K)
               global_EdepAdv(globali,globalj,globalk) = EdepAdv(I,J,K)
               global_EdepRad(globali,globalj,globalk) = EdepRad(I,J,K)
               global_xkapA(globali,globalj,globalk) = xkapA(I,J,K)
               if(tmpflag.gt.0) then
                  global_stellarinput(globali,globalj,globalk) = 
     $                 stellarinput(I,J,K)
                  global_delta_ER(globali,globalj,globalk) = 
     $                 delta_ER(I,J,K)
                  global_xkapP(globali,globalj,globalk) = xkapP(I,J,K)
               else
                  global_stellarinput(globali,globalj,globalk) = 0.d0
                  global_delta_ER(globali,globalj,globalk) = 0.d0
                  global_xkapP(globali,globalj,globalk) = 0.d0
               endif
            ENDDO
            global_upperbnd(globalj,globalk)=upperbnd(j,k)
         ENDDO
      ENDDO
!-- Exclude Outer Region
      IF(EXCLUDEREGION) then
         DO k=1,locNZ
            DO j=1,locNY
               if(upperbnd(j,k).lt.locNX+1) then
                  globalj = j+locNY*proc_coords(2)
                  globalk = k+locNZ*proc_coords(3)
                  DO i=upperbnd(j,k),locNX
                     globali = i+locNX*proc_coords(1)
                     global_T(globali,globalj,globalk) = 90.d0
                     global_xkapR(globali,globalj,globalk) = 0.d0
                     global_xkapA(globali,globalj,globalk) = 0.d0
                     global_xkapP(globali,globalj,globalk) = 0.d0
                  ENDDO
               endif
            ENDDO
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE FILLGLOBAL

!-----------------

      SUBROUTINE SPHERETOCARTESIAN(EXCLUDEREGION)
      USE var_init
      IMPLICIT NONE
      integer :: i,j,k,center,xc,yc,zc
      integer :: rindex,phiindex,zindex
      integer :: rinner
      double precision :: rmin,rmax,rc,S
      double precision :: rvalue,phivalue,zvalue
      logical :: EXCLUDEREGION
      center = cartN/2
      rinner = 15
      rmin = global_XXB(1)
      rmax = global_XXB(global_NX)
!--   initilize everything
      cart_T = 90.d0
      cart_VX = 0.d0
      cart_VY = 0.d0
      cart_VZ = 0.d0
!-- serach through the region, matching to data
      do k=1,cartN
         zc = k-center
         do j=1,cartN
            yc = j-center
            do i=1,cartN
               xc = i-center
               rc = (xc*xc + yc*yc + zc*zc)**0.5
               IF (rc.LE.1.d0*center.AND.rc.GE.1.d0*rinner) THEN
                  rvalue = rmin+(rmax-rmin)*
     %                 ((rc-rinner)/(center-rinner))
                  call locateRvalue(rvalue,rindex)
!-- see http://www.math.montana.edu/frankw/ccp/multiworld/multipleIVP/
!---- spherical/body.htm#skip2
                  S = (xc*xc + yc*yc)**0.5
                  if(xc.gt.0) phivalue = asin(yc/S)
                  if(xc.lt.0) phivalue = PI - asin(yc/S)
                  if(phivalue.lt.0) phivalue=2.d0*PI+phivalue
                  call locatePHIvalue(phivalue,phiindex)
                  zvalue = acos(zc/rc) - PI/2.d0
                  if(zvalue.le.global_XZB(global_NZ).and.
     %                 zvalue.ge.global_XZB(1)) then
                     call locateZvalue(zvalue,zindex)
                     cart_T(i,j,k) = global_T(rindex,phiindex,zindex) 
                     cart_RH(i,j,k) = global_RH(rindex,phiindex,zindex)
                     cart_VX(i,j,k) = global_V(rindex,phiindex,zindex)
                     if(ncosys.eq.0) then
                        cart_VY(i,j,k) = global_G(rindex,phiindex,
     %                       zindex)
                        cart_VZ(i,j,k) = global_H(rindex,phiindex,
     %                       zindex)
                     elseif(ncosys.eq.2) then
                        cart_VY(i,j,k) = global_G(rindex,phiindex,
     %                       zindex) *
     %                       global_XXB(rindex)*cos(global_XZB(zindex))
                        cart_VZ(i,j,k) = global_H(rindex,phiindex,
     %                       zindex)*
     %                       global_XXB(rindex)
                     cart_EdepFor(i,j,k) = 
     %                    global_EdepFor(rindex,phiindex,zindex) 
                     cart_EdepAdv(i,j,k) = 
     %                    global_EdepAdv(rindex,phiindex,zindex) 
                     cart_EdepRad(i,j,k) = 
     %                    global_EdepRad(rindex,phiindex,zindex) 
                     else
                        print *,'unknown ncosys for converting vel'
                        stop
                     endif
                  endif
c                  if(EXCLUDEREGION) then
c                     if(rindex.gt.UPPERBND(phiindex,zindex)) then
c                        cart_T(i,j,k) = 90.d0
c                        cart_RH(i,j,k) = 10.**(-10.d0)
c                     endif
c                  endif
!-- when doing velocities, convert to actual vel
c                Vy = Vy*RADIUS*cos(z)
c                Vz = Vz*radius
               ENDIF
            enddo
         enddo
      enddo

      RETURN
      END SUBROUTINE SPHERETOCARTESIAN
      
      SUBROUTINE locateRvalue(rvalue,rindex)
      USE var_init
      IMPLICIT NONE
      integer ::i,rindex
      double precision :: rvalue,delta1,delta2
      DO I=2,global_NX
         if(global_XXB(I).gt.rvalue) then
            delta1 = abs(rvalue-global_XXB(I))
            delta2 = abs(rvalue-global_XXB(I-1))
            if(delta1.gt.delta2) rindex = I-1
            if(delta1.le.delta2) rindex = I
            goto 10
         endif
      ENDDO
      rindex = global_NX
 10   CONTINUE
      RETURN
      END SUBROUTINE locateRvalue

      SUBROUTINE locatePHIvalue(phivalue,phiindex)
      USE var_init
      IMPLICIT NONE
      integer ::j,phiindex
      double precision :: phivalue,delta1,delta2
      DO J=2,global_NY
         if(global_XYB(J).gt.phivalue) then
            delta1 = abs(phivalue-global_XYB(J))
            delta2 = abs(phivalue-global_XYB(J-1))
            if(delta1.gt.delta2) phiindex = J-1
            if(delta1.le.delta2) phiindex = J
            goto 11
         endif
      ENDDO
      phiindex = global_NY
c      print *,'global_NY was exceeded',phivalue,global_XYB(global_NY)
 11   CONTINUE
      RETURN
      END SUBROUTINE locatePHIvalue

      SUBROUTINE locateZvalue(zvalue,zindex)
      USE var_init
      IMPLICIT NONE
      integer ::k,zindex
      double precision :: zvalue,delta1,delta2
      DO k=2,global_NZ
         if(global_XZB(K).gt.zvalue) then
            delta1 = abs(zvalue-global_XZB(K))
            delta2 = abs(zvalue-global_XZB(K-1))
            if(delta1.gt.delta2) zindex = K-1
            if(delta1.le.delta2) zindex = K
            goto 12
         endif
      ENDDO
      zindex = global_NZ
 12   CONTINUE
      RETURN
      END SUBROUTINE locateZvalue

      SUBROUTINE PRINTGLOBAL
      USE var_init
      IMPLICIT NONE
      integer :: i,j,K
      WRITE (12) global_NX
      WRITE (12) global_NY
      WRITE (12) global_NZ
      WRITE (12) locNX,locNY,locNZ,numprocs
      WRITE (12) ZEIT,DELT,NUM_ITER
      WRITE (12) (global_XXB(I),I=1,global_NX)
      WRITE (12) (global_XYB(J),J=1,global_NY)
      WRITE (12) (global_XZB(K),K=1,global_NZ)
      WRITE (12) (((global_V(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_G(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_H(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_RH(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_T(I,J,K) ,I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_ER(I,J,K) ,I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_PG(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_CS(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) ((global_upperbnd(J,K),J=1,global_NY),K=1,global_NZ)
      WRITE (12) TinitBottom,Tinittop,rhinittop,rhinitbottom,
     %     hydrotime,Maccreated
      WRITE (12) NCOSYS
      WRITE (12) (((global_xkapR(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_flimx(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_EdepFor(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_EdepAdv(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_EdepRad(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      WRITE (12) (((global_xkapA(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)

      if(tmpflag.gt.0) then
         WRITE (12) (((global_stellarinput(I,J,K),I=1,global_NX),
     %        J=1,global_NY),K=1,global_NZ)
         WRITE (12) (((global_delta_ER(I,J,K),I=1,global_NX),
     %        J=1,global_NY),K=1,global_NZ)
         WRITE (12) (((global_xkapP(I,J,K),I=1,global_NX),
     %        J=1,global_NY),K=1,global_NZ)
      endif

      RETURN
      END SUBROUTINE PRINTGLOBAL

      SUBROUTINE READFILE
      USE var_init
      IMPLICIT NONE
      integer :: i,j,k
      READ (12) global_NX,global_NY,global_NZ
      READ (12) locNX,locNY,locNZ,numprocs
      READ (12) myid,proc_coords(1),proc_coords(2),proc_coords(3)
      READ (12) proc_dims(1),proc_dims(2),proc_dims(3)
      READ (12) ZEIT,DELT,NUM_ITER
      READ (12) (XXA(I),I=0,locNX+2)
      READ (12) (XXB(I),I=-1,locNX+2)
      READ (12) (DXA(I),I=1,locNX+2)
      READ (12) (XYA(J),J=0,locNY+2)
      READ (12) (XYB(J),J=-1,locNY+2)
      READ (12) (DYA(J),J=1,locNY+2)
      READ (12) (XZA(K),K=0,locNZ+2)
      READ (12) (XZB(K),K=-1,locNZ+2)
      READ (12) (DZA(K),K=1,locNZ+2)
      READ (12) (((V(I,J,K),I=0,locNX+2),J=-1,locNY+2),K=-1,locNZ+2)
      READ (12) (((G(I,J,k),I=-1,locNX+2),J=0,locNY+2),K=-1,locNZ+2)
      READ (12) (((H(I,J,k),I=-1,locNX+2),J=-1,locNY+2),K=0,locNZ+2)
      READ (12) (((RH(I,J,K),I=-1,locNX+2),J=-1,locNY+2),K=-1,locNZ+2)
      READ (12) (((T(I,J,K),I=-1,locNX+2),J=-1,locNY+2),K=-1,locNZ+2)
!--   new stuff
      if(tmpflag.eq.1) then
         READ (12) (((ER(I,J,K),I=-1,locNX+2),J=-1,locNY+2),
     $        K=-1,locNZ+2)
      endif
      READ (12) (((PG(I,J,K),I=-1,locNX+2),J=-1,locNY+2),K=-1,locNZ+2)
      READ (12) (((CS(I,J,K),I=0,locNX+1),J=0,locNY+1),K=0,locNZ+1)
      READ (12) ((UPPERBND(J,K),J=1,locNY),K=1,locNZ)
      READ (12) TinitBottom,Tinittop,rhinittop,rhinitbottom,
     %     hydrotime,Maccreated
!--- THIS STUFF ISN'T READ BACK INTO A RESTART      
      READ (12) NCOSYS
      READ (12) (((XKAPR(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
      READ (12) (((flimx(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
      READ (12) (((EdepFor(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
      READ (12) (((EdepAdv(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
      READ (12) (((EdepRad(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
      READ (12) (((xkapA(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
!--   new stuff
      if(tmpflag.eq.1) then
         READ (12) (((stellarinput(I,J,K),I=1,locNX),J=1,locNY),
     $        K=1,locNZ)
         READ (12) (((delta_ER(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
C         READ (12) (((xkapP(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
         xkapP = 0.D0
      endif


      RETURN
      END

