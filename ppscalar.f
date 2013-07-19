      SUBROUTINE SCALARINIT
**  Initilization for the scalars or read-in routine for a restart
      USE input_init
      USE scalar_var_init
      USE mpi_var_init
      USE fluid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: stat(MPI_STATUS_SIZE)
      integer :: n,i,j,k,loci,locj,lock
      integer :: dummy_nscalar
      double precision global_SCALAR(NSCALAR,global_NX,global_NY,
     %     global_NZ)
!--- INFORMATION FOR A NEW RUN (ALWAYS USE THIS FOR ISCAL=1:2 test cases)
      IF (MOD(IREA,10).EQ.0.or.iscal.eq.1.or.iscal.eq.2) THEN
!--density test case
         if(iscal.eq.1) then
            DO K=-1,locNZ+2
               DO J=-1,locNY+2
                  DO I=-1,locNX+2
                     SCALAR(1,i,j,k) = RH(i,j,k)
                  ENDDO
               ENDDO
            ENDDO
!--diffusion testcase
         elseif(iscal.eq.2) then
!--global gaussian
            do k=-1,locNZ+2
               lock = k+locNZ*proc_coords(3)
               if(lock.gt.global_NZ+2.or.lock.lt.-1) then
                  print *,'ERROR IN SCALARINIT-k'
                  call clean_stop
               endif
               do j=-1,locNY+2
                  locj = j+locNY*proc_coords(2)
                  if(locj.gt.global_NY+2.or.locj.lt.-1) then
                     print *,'ERROR IN SCALARINIT-j'
                     call clean_stop
                  endif
                  do i=-1,locNX+2
                     loci = i+locNX*proc_coords(1)
                     if(loci.gt.global_NX+2.or.loci.lt.-1) then
                        print *,'ERROR IN SCALARINIT-i'
                        call clean_stop
                     endif
                     SCALAR(1,i,j,k) = 1000.0*exp(-(
     %                    (((loci-INT(global_NX/2))**2.d0)/45.d0)))
                  enddo
               enddo
            enddo
         elseif(iscal.eq.3) then
            if(icstyp.eq.1.or.icstyp.eq.2) CALL SOUND
            CALL PRESSG
            CALL CO_SCALAR_init         
         endif
!--- READ IN INFORMATION FOR A RESTART (never for ISCAL=1:2.. test cases)
      ELSEIF (mod(IREA,10).EQ.1) THEN
         OPEN (12,FILE='SCALAR',STATUS='OLD',FORM='UNFORMATTED')
         if(myid.eq.0) print *,'reading SCALAR'
         rewind(12)
         READ (12) dummy_nscalar
!--   check that nscalar is the same
         if(dummy_nscalar.ne.nscalar) then
            print *,'nscalar incorrect',dummy_nscalar,nscalar
            call CLEAN_STOP
         endif
         do n=1,nscalar
            READ (12) (((global_SCALAR(N,I,J,K),I=1,global_NX),
     %           J=1,global_NY),K=1,global_NZ)
         enddo
         CLOSE(12)
!---  Initilize local copies of the variables
         do k=1,locNZ
            lock = k+locNZ*proc_coords(3)
            if(lock.gt.global_NZ.or.lock.lt.1) then
               print *,'ERROR IN INIT-k:initscalar'
               stop
            endif
            do j=1,locNY
               locj = j+locNY*proc_coords(2)
               if(locj.gt.global_NY.or.locj.lt.1) then
                  print *,'ERROR IN INIT-j:initscalar'
                  stop
               endif
               do i=1,locNX
                  loci = i+locNX*proc_coords(1)
                  if(loci.gt.global_NX.or.loci.lt.1) then
                     print *,'ERROR IN INIT-i:initscalar'
                     stop
                  endif
                  do n=1,nscalar
                     SCALAR(n,i,j,k)=global_SCALAR(n,loci,locj,lock)
                  enddo
               enddo
            enddo
         enddo
      ELSE
         print *,'irea',irea,' is undefined:SCALARINIT'
         CALL CLEAN_STOP
      END IF
      CALL SCALAR_BOUNDS
      RETURN
      END SUBROUTINE SCALARINIT
      
      SUBROUTINE SCALAR_BOUNDS
!--routine for the boundary conditions on scalars
      USE scalar_var_init
      USE input_init
      USE fluid_var_init
      USE grid_var_init
      USE mpi_grid_init
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: stat(MPI_STATUS_SIZE)
      integer :: i,j,k,n
!---  CALL SHIFTING ROUTINE
      do n=1,nscalar
         pass_sc(:,:,:) = scalar(n,:,:,:) 
         if(modtyp.eq.2) then
            CALL PLANETSHIFT(9,n)
!--   now correct the polar shifts if necessary
            if(poles) then
               if(proc_coords(3).eq.proc_dims(3)-1) then
                  CALL NORTH_POLE_SHIFT(9,n)
               endif
               if(proc_coords(3).eq.0) then
                  CALL SOUTH_POLE_SHIFT(9,n)
               endif
            endif            
         else
            CALL NEWSHIFTVAR3D(9)
         endif
         scalar(n,:,:,:) = pass_sc(:,:,:)
      enddo
!---  IMPLEMENT BOUNDARY CONDITIONS
      SELECT CASE(iscal)
      CASE (1)
!  1=density test case
         do n=1,nscalar
!--XMIN
            if(MPIleft.eq.MPI_PROC_NULL) then
               DO K=-1,locNZ+2
                  DO J=-1,locNY+2
                     SCALAR(n,1,J,K) = RHinitbottom
                     SCALAR(n,0,J,K) = SCALAR(n,1,J,K)-D_RH*
     $                    (xxb(1)-xxb(0))
                     SCALAR(n,-1,J,K)= SCALAR(n,1,J,K)-D_RH*
     $                    (xxb(1)-xxb(-1))
                  ENDDO
               ENDDO
            endif
!--   XMAX
            if(MPIright.eq.MPI_PROC_NULL) then
               DO K=-1,locNZ+2
                  DO J=-1,locNY+2
                     DO I=upperbnd(j,k)+1,locNX+2
                        SCALAR(n,I,J,K) = 1.0e-10
                     ENDDO
                     IF(upperbnd(j,k).eq.locNX) THEN
                        SCALAR(n,locNX,j,k) = 1.0e-10
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
!--   ZMIN
            if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
               DO I=-1,locNX+2
                  DO J=-1,locNY+2
                     SCALAR(n,I,j,1)=SCALAR(n,I,j,2)
                     SCALAR(n,I,j,0)=SCALAR(n,I,j,1)
                     SCALAR(n,I,j,-1)=SCALAR(n,I,j,1)
                  ENDDO
               ENDDO
            endif
!--   ZMAX
            if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL) then
               DO I=-1,locNX+2
                  DO J=-1,locNY+2
                     SCALAR(n,I,j,locNZ)=SCALAR(n,I,j,locNZ-1)
                     SCALAR(n,I,j,locNZ+1)=SCALAR(n,I,j,locNZ)
                     SCALAR(n,I,j,locNZ+2)=SCALAR(n,I,j,locNZ)
                  ENDDO
               ENDDO
            endif
         enddo
      CASE (2:3)
!  2=diffusion test case
!  3=Xco, zero derivatives at boundaries
         do n=1,nscalar
!--XMIN
            if(MPIleft.eq.MPI_PROC_NULL) then
               DO K=-1,locNZ+2
                  DO J=-1,locNY+2
                     SCALAR(n,0,J,K) =SCALAR(n,1,J,K)
                     SCALAR(n,-1,J,K)=SCALAR(n,1,J,K)
                  ENDDO
               ENDDO
            endif
!--   XMAX
            if(MPIright.eq.MPI_PROC_NULL) then
               DO K=-1,locNZ+2
                  DO J=-1,locNY+2
                     DO I=upperbnd(j,k),locNX+2
                        SCALAR(n,I,J,K) = SCALAR(n,upperbnd(j,k)-1,J,K)
                     ENDDO
                  ENDDO
               ENDDO
            endif
!--   ZMIN
            if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
               DO I=-1,locNX+2
                  DO J=-1,locNY+2
                     SCALAR(n,I,j,1)=SCALAR(n,I,j,2)
                     SCALAR(n,I,j,0)=SCALAR(n,I,j,1)
                     SCALAR(n,I,j,-1)=SCALAR(n,I,j,1)
                  ENDDO
               ENDDO
            endif
!--   ZMAX
            if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL) then
               DO I=-1,locNX+2
                  DO J=-1,locNY+2
                     SCALAR(n,I,j,locNZ)=SCALAR(n,I,j,locNZ-1)
                     SCALAR(n,I,j,locNZ+1)=SCALAR(n,I,j,locNZ)
                     SCALAR(n,I,j,locNZ+2)=SCALAR(n,I,j,locNZ)
                  ENDDO
               ENDDO
            endif
         enddo
      CASE (4:)
         print *,'unknown ISCAL=',ISCAL
         stop
      CASE DEFAULT
         print *,'also unknown ISCAL=',ISCAL
         stop
      END SELECT
      RETURN
      END SUBROUTINE SCALAR_BOUNDS

      SUBROUTINE ADVECTSCALAR3D(scalar_num)
**  3D Advection Terms for scalars
      USE input_init
      USE fluid_var_init
      USE grid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE scalar_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k,scalar_num
      double precision :: FLSX(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLSY(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLSZ(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: S(NTOP),GRAD(0:NTOP),GRAD1(0:NTOP)
      double precision :: mu_co
      integer :: IND(NTOP)
      double precision,parameter ::half=0.5d0,one=1.0d0,two=2.0d0
      double precision,parameter :: xdmin=1.0d-50
      integer :: stat(MPI_STATUS_SIZE)
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: SCAL
!-- point to scalar
c      SCAL => scalar(scalar_num,:,:,:)
      SCAL => pass_sc

!-If using iscal=3 (XCO), change scal from X to rh_co
      if(iscal.eq.3) then 
!-mean molecular weight of Carbon Monoxide (CO)
c         mu_co = 2.3
c         SCAL=SCAL*RH*(mu_co/mu_gas)
         SCAL=SCAL*RH
      endif
!-call bounds to get the velocities correct
      CALL BOUNDS

**  X-Direction
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX+1
               IND(I)=I
               S(I)=V(I,J,K)*DELT
               IF (V(I,J,K).GT..0) IND(I)=I-1
            ENDDO
            DO I=0,locNX+1
               GRAD1(I)=(SCAL(I,J,K)-SCAL(I-1,J,K))*
     &              (SCAL(I+1,J,K)-SCAL(I,J,K))
               GRAD(I)=0.d0
               IF (GRAD1(I).GT.0.d0) GRAD(I)=two*GRAD1(I)/
     &              ((SCAL(I,J,K)-SCAL(I-1,J,K))*DXB(I+1)+
     &              (SCAL(I+1,J,K)-SCAL(I,J,K))*DXB(I))
            ENDDO               
            DO I=1,locNX+1
               FLSX(I,J,K)=(surxa(i)-ccxa(i)*s(i)) * S(I)
     &              * (SCAL(IND(I),J,K)
     &              +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
            ENDDO
         ENDDO
      ENDDO
      
**  Y-Direction
      DO K=1,locNZ
         DO I=1,locNX
            DO J=1,locNY+1
               S(J)=G(I,J,K)*DELT 
               IND(J)=J
               IF (G(I,J,K).GT.0.d0) IND(J)=J-1
            ENDDO
            DO J=0,locNY+1
               GRAD1(J)=(SCAL(I,J,K)-SCAL(I,J-1,K))*
     &              (SCAL(I,J+1,K)-SCAL(I,J,K))
               GRAD(J)=0.d0
               IF (GRAD1(J).GT..0) GRAD(J)=two*GRAD1(J)/
     &              ((SCAL(I,J,K)-SCAL(I,J-1,K))*DYB(J+1)+
     &              (SCAL(I,J+1,K)-SCAL(I,J,K))*DYB(J))
            ENDDO
            DO J=1,locNY+1
               FLSY(I,J,K)=S(J) * (surya(j)+ccya(j)*S(J)) *
     &              (SCAL(I,IND(J),K)
     &              +(XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)))
            ENDDO
         ENDDO
      ENDDO
      
**  Z-Direction
      DO J=1,locNY
         DO I=1,locNX
            DO K=1,locNZ+1
               IND(K)=K
               S(K)=H(I,J,K)*DELT
               IF (H(I,J,K).GT..0) IND(K)=K-1
            ENDDO              
            DO K=0,locNZ+1
               GRAD1(K)=(SCAL(I,J,K)-SCAL(I,J,K-1))*
     &              (SCAL(I,J,K+1)-SCAL(I,J,K))
               GRAD(K)=0.d0
               IF (GRAD1(K).GT.0.d0) GRAD(K)=two*GRAD1(K)/
     &              ((SCAL(I,J,K)-SCAL(I,J,K-1))*DZB(K+1)+
     &              (SCAL(I,J,K+1)-SCAL(I,J,K))*DZB(K))
            ENDDO
            DO K=1,locNZ+1
               FLSZ(I,J,K)=(surza(k) + ccza(k)*s(k)) * S(k)
     &              * (SCAL(I,J,ind(k))
     &              +(XZA(k)-XZB(IND(k))-half*S(k))*GRAD(IND(k)))
            ENDDO
         ENDDO
      ENDDO

!-MPI FOR FLSX,FLSY,FLSZ
!--- Pass FLSX(*,locNY,*)->FLSX(*,0,*)
!----- Y-COORD, + DIRECTION
      call MPI_SENDRECV(FLSX(1:locNX+1,locNY,1:locNZ),
     %     (locNX+1)*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     FLSX(1:locNX+1,0,1:locNZ),(locNX+1)*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIlower,1,COMM_CART,stat,ierr)
!--- Pass FLSX(*,*,locNZ)->FLSX(*,*,0)
!----- Z-COORD, + DIRECTION
      call MPI_SENDRECV(FLSX(1:locNX+1,1:locNY,locNZ),
     %     (locNX+1)*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     FLSX(1:locNX+1,1:locNY,0),(locNX+1)*locNY,
     %     MPI_DOUBLE_PRECISION,MPIbelow,1,COMM_CART,stat,ierr)
      if(poles.and.proc_coords(3).eq.0) then
!--  redo South pole, flmx defined on xzb, no sign flip
         call MPI_SENDRECV(FLSX(1:locNX+1,0:locNY,1),
     %        (locNX+1)*(1+locNY),MPI_DOUBLE_PRECISION,MPIbelow_pole,1,
     %        FLSX(1:locNX+1,0:locNY,0),(locNX+1)*(1+locNY),
     %        MPI_DOUBLE_PRECISION,MPIbelow_pole,1,COMM_CART,stat,ierr)
      endif

!--- Pass FLSY(locNX,*,*)->FLSY(0,*,*)
!----- X-COORD, + DIRECTION
      call MPI_SENDRECV(FLSY(locNX,1:locNY+1,1:locNZ),
     %     (locNY+1)*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %     FLSY(0,1:locNY+1,1:locNZ),(locNY+1)*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
!--- Pass FLSY(*,*,locNZ)->FLSY(*,*,0)
!----- Z-COORD, + DIRECTION
      call MPI_SENDRECV(FLSY(1:locNX,1:locNY+1,locNZ),
     %     locNX*(locNY+1),MPI_DOUBLE_PRECISION,MPIabove,1,
     %     FLSY(1:locNX,1:locNY+1,0),locNX*(locNY+1),
     %     MPI_DOUBLE_PRECISION,MPIbelow,1,COMM_CART,stat,ierr)
!--  redo South pole, flmy defined on xzb, use sign flip
      if(poles.and.proc_coords(3).eq.0) then
         call MPI_SENDRECV(-FLSY(1:locNX,1:locNY+1,1),
     %        locNX*(locNY+1),MPI_DOUBLE_PRECISION,MPIbelow_pole,1,
     %        FLSY(1:locNX,1:locNY+1,0),locNX*(locNY+1),
     %        MPI_DOUBLE_PRECISION,MPIbelow_pole,1,COMM_CART,stat,ierr)
      endif
!--- Pass FLSZ(locNX,*,*)->FLSZ(0,*,*)
!----- X-COORD, + DIRECTION
      call MPI_SENDRECV(FLSZ(locNX,1:locNY,1:locNZ+1),
     %     locNY*(locNZ+1),MPI_DOUBLE_PRECISION,MPIright,1,
     %     FLSZ(0,1:locNY,1:locNZ+1),locNY*(locNZ+1),
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
!--- Pass FLSZ(*,locNY,*)->FLSZ(*,0,*)
!----- Y-COORD, + DIRECTION
      call MPI_SENDRECV(FLSZ(1:locNX,locNY,1:locNZ+1),
     %     locNX*(locNZ+1),MPI_DOUBLE_PRECISION,MPIupper,1,
     %     FLSZ(1:locNX,0,1:locNZ+1),locNX*(locNZ+1),
     %     MPI_DOUBLE_PRECISION,MPIlower,1,COMM_CART,stat,ierr)

!--- Pass FLSZ(*,*,locNZ-1)->-FLSZ(*,*,locNZ+1)
!----- Z-COORD, + DIRECTION (OVER NORTHPOLE)
      if(poles.and.proc_coords(3).eq.proc_dims(3)-1) then
         call MPI_SENDRECV(-FLSZ(0:locNX+1,0:locNY+1,locNZ-1),
     %        (locNX+2)*(locNY+2),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,
     %        FLSZ(0:locNX+1,0:locNY+1,locNZ+1),
     %        (locNX+2)*(locNY+2),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,COMM_CART,stat,ierr)
      endif

******************************************************************
* * CALCULATE NEW VALUES OF SCALAR
****************************************************************** for
!     XCO, advection was done in rh_co, so transfer back to XCO while
!     using calculated fluxes
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               SCAL(I,J,K)=SCAL(I,J,K)
     &              -(FLSX(I+1,J,K)-FLSX(I,J,K))/volxb(i)
     &              -(FLSY(I,J+1,K)-FLSY(I,J,K))/volyb(j)
     &              -(FLSZ(I,J,K+1)-FLSZ(I,J,K))/volzb(k)
            ENDDO
         ENDDO
      ENDDO
      if(iscal.eq.3) then 
         SCAL=SCAL/RH
      endif
!--- RELEASE THE POINTER
      NULLIFY(SCAL)
      RETURN
      END SUBROUTINE ADVECTSCALAR3D


      SUBROUTINE SCALAR_ADJUST_XCO(scalar_num)
**  3D Advection Terms for scalars
      USE input_init
      USE fluid_var_init
      USE scalar_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k,scalar_num
      double precision :: mu_co
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: SCAL
!-- point to scalar
c      SCAL => scalar(scalar_num,:,:,:)
      SCAL => pass_sc
!-mean molecular weight of Carbon Monoxide (CO)
      mu_co = mu_gas
!  Subtract off the continuity equation for the total density because
!  SCALAR is defined as a mass fraction, not a density
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               SCAL(I,J,K) = SCAL(I,J,K)*
     %              (2.d0 - RH(i,j,k)/rh_pre_adv(i,j,k))
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE SCALAR_ADJUST_XCO


