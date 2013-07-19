      SUBROUTINE TESTADVECT
**  2D Advection Terms  (Mass X variable, includes Coriolis forces (for Vphi)
      USE input_init
      USE fluid_var_init
      USE grid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: ED(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: VD(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: FLX(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLY(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLMX(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLMY(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: S(NTOP),GRAD(0:NTOP),GRAD1(0:NTOP)
      integer :: IND(NTOP)
      integer :: irht
      double precision,parameter ::half=0.5d0,one=1.0d0,two=2.0d0
      double precision,parameter :: xdmin=1.0d-50
      integer :: stat(MPI_STATUS_SIZE)

!--- SET THERMAL ENERGY TO BE ADVECTED: 0=T: 1=T*CV*RH
      irht = 1
**  General Boundary Conditions
      CALL BOUNDS
******************************************
** CALCULATE CHANGE IN DENSITY: MASSFLOW
******************************************

      if (isoden.eq.0) then
**  X-Direction
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX+1
                  IND(I)=I
                  S(I)=V(I,J,K)*DELT
                  IF (V(I,J,K).GT.0.d0) IND(I)=I-1
               ENDDO
            
               DO I=0,locNX+1
                  GRAD1(I)=(RH(I,J,K)-RH(I-1,J,K))*
     &                 (RH(I+1,J,K)-RH(I,J,K))
                  GRAD(I)=0.d0
                  IF (GRAD1(I).GT.0.d0) GRAD(I)=two*GRAD1(I)/
     &                 ((RH(I,J,K)-RH(I-1,J,K))*DXB(I+1)+
     &                 (RH(I+1,J,K)-RH(I,J,K))*DXB(I))
               ENDDO
               
               DO I=1,locNX+1
                  FLMX(I,J,K)=(surxa(i)-ccxa(i)*s(i)) * S(I)
     &                 * (RH(IND(I),J,K)
     &                 +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
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
                  GRAD1(J)=(RH(I,J,K)-RH(I,J-1,K))*
     &                 (RH(I,J+1,K)-RH(I,J,K))
                  GRAD(J)=0.d0
                  IF (GRAD1(J).GT.0.d0) GRAD(J)=two*GRAD1(J)/
     &                 ((RH(I,J,K)-RH(I,J-1,K))*DYB(J+1)+
     &                 (RH(I,J+1,K)-RH(I,J,K))*DYB(J))
               ENDDO
               DO J=1,locNY+1
                  FLMY(I,J,K)=S(J) * (surya(j)+ccya(j)*S(J)) *
     &                 (RH(I,IND(J),K)
     &                 +(XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)))
               ENDDO
            ENDDO
         ENDDO

      else if (isoden.eq.1) then
         FLMX = 0.d0
         FLMY = 0.d0
      end if

!-MPI FOR FLMX,FLMY
!--- Pass FLMX(*,locNY,*)->FLMX(*,0,*)
!---- Y-COORD, + DIRECTION
      call MPI_SENDRECV(FLMX(1:locNX+1,locNY,1:locNZ),
     %     (locNX+1)*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     FLMX(1:locNX+1,0,1:locNZ),(locNX+1)*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIlower,1,COMM_CART,stat,ierr)
!--- Pass FLMY(locNX,*)->FLMY(0,*)
!---- X-COORD, + DIRECTION
      call MPI_SENDRECV(FLMY(locNX,1:locNY+1,1:locNZ),
     %     (locNY+1)*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %     FLMY(0,1:locNY+1,1:locNZ),(locNY+1)*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
c   moved down here from about after 'icell' modification
      CALL QDENS


******************************************
**  CALCULATE CHANGE IN X-VELOCITY
******************************************
**  X-Direction
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               S(I)=(V(I,J,K)+V(I+1,J,K))*DELT/2.
               IND(I)=I+1
               IF (S(I).GT.0.d0) IND(I)=I
            ENDDO
            DO I=1,locNX+1
               GRAD1(I)=(V(I,J,K)-V(I-1,J,K))*
     &              (V(I+1,J,K)-V(I,J,K))
               GRAD(I)=.0
               IF (GRAD1(I).GT.0.d0) GRAD(I)=two*GRAD1(I)/
     &              ((V(I+1,J,K)-V(I,J,K))*DXA(I)+
     &              (V(I,J,K)-V(I-1,J,K))*DXA(I+1))
            ENDDO
            DO I=1,locNX
               FLX(I,J,K)=half*(flmx(i,j,K)+flmx(i+1,j,K)) *
     &              ( V(IND(I),j,K)+
     &              (XXB(I)-XXA(IND(I))-half*S(I))*GRAD(IND(I)))
            ENDDO
         ENDDO
      ENDDO
**  Y-Direction
      if(MPIleft.eq.MPI_PROC_NULL) then
         DO K=1,locNZ
            DO I=2,locNX
               DO J=1,locNY+1
                  S(J)=(DDX1(I)*G(I,J,K)+DDX0(I)*G(I-1,J,K))*DELT
                  IND(J)=J
                  IF (S(J).GT.0.d0) IND(J)=J-1
               ENDDO
               DO J=0,locNY+1
                  GRAD1(J)=(V(I,J,K)-V(I,J-1,K))*
     &                 (V(I,J+1,K)-V(I,J,K))
                  GRAD(J)=.0
                  IF (GRAD1(J).GT.0.d0) GRAD(J)=two*GRAD1(J)/
     &                 ((V(I,J+1,K)-V(I,J,K))*DYB(J)+
     &                 (V(I,J,K)-V(I,J-1,K))*DYB(J+1))
               ENDDO
               DO J=1,locNY+1
                  FLY(I,J,k)=(ddx1(i)*FLMY(i,j,K)+DDX0(i)*FLMY(i-1,j,K))
     &                 *( V(I,IND(J),K)+
     &                 (XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)) )
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO I=1,locNX
               DO J=1,locNY+1
                  S(J)=(DDX1(I)*G(I,J,K)+DDX0(I)*G(I-1,J,K))*DELT
                  IND(J)=J
                  IF (S(J).GT..0) IND(J)=J-1
               ENDDO
               DO J=0,locNY+1
                  GRAD1(J)=(V(I,J,K)-V(I,J-1,K))*
     &                 (V(I,J+1,K)-V(I,J,K))
                  GRAD(J)=.0
                  IF (GRAD1(J).GT..0) GRAD(J)=two*GRAD1(J)/
     &                 ((V(I,J+1,K)-V(I,J,K))*DYB(J)+
     &                 (V(I,J,K)-V(I,J-1,K))*DYB(J+1))
               ENDDO
               DO J=1,locNY+1
                  FLY(I,J,k)=(ddx1(i)*FLMY(i,j,K)+DDX0(i)*FLMY(i-1,j,K))
     &                 *( V(I,IND(J),K)+
     &                 (XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)) )
               ENDDO
            ENDDO
         ENDDO
      endif
!--- Pass FLX(locNX,*)->FLX(0,*) for below calculation of VD from i=1,locNX
!---- X-COORD, + DIRECTION
!FLX(0:locNX+1,0:locNY+1,0:locNZ+1)

      call MPI_SENDRECV(FLX(locNX,1:locNY,1:locNZ),
     %     locNY*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %     FLX(0,1:locNY,1:locNZ),locNY*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)

******************************************      
**  Calculate New X-Momentum density based on calculated fluxes
******************************************
      if(MPIleft.eq.MPI_PROC_NULL) then
         DO K=1,locNZ
            DO J=1,locNY
               DO I=2,locNX
                  VD(I,J,K)=V(I,J,K)*rhqx(i,j,K)
     &                 - (FLX(I,J,K)-FLX(I-1,J,K))/volxa(i)
     &                 - (FLY(I,J+1,K)-FLY(I,J,K))/volyb(j)
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  VD(I,J,K)=V(I,J,K)*rhqx(i,j,K)
     &                 - (FLX(I,J,K)-FLX(I-1,J,K))/volxa(i)
     &                 - (FLY(I,J+1,K)-FLY(I,J,K))/volyb(j)
               ENDDO
            ENDDO
         ENDDO
      endif

******************************************
**  CALCULATE CHANGE IN Y-MOMENTUM DENSITY (eg. MERIDIONAL)
*     note that the order of y and z directions are switched from original
******************************************
      DO K=-1,locNZ+2
         DO J=0,locNY+2
            DO I=-1,locNX+2
               GD(I,J,K)= (G(I,J,K)+OMROT) * geoxg(i)**2 * geozg(k)**2
            ENDDO
         ENDDO
      ENDDO
      
**  X-Direction
      if(MPIlower.eq.MPI_PROC_NULL) then
c         print *,'INHERE'
c         stop
         DO K=1,locNZ
            DO J=2,locNY
               DO I=1,locNX+1
                  S(I)=(DDY1(J)*V(I,J,k)+DDY0(J)*V(I,J-1,K))*DELT
                  IND(I)=I
                  IF (S(I).GT.0.d0) IND(I)=I-1
               ENDDO
               DO I=0,locNX+1
                  GRAD1(I)=(GD(I,J,K)-GD(I-1,J,K))*
     &                 (GD(I+1,J,K)-GD(I,J,K))
                  GRAD(I)=0.d0
                  IF (GRAD1(I).GT.0.d0) GRAD(I)=two*GRAD1(I)/
     &                 ((GD(I+1,J,K)-GD(I,J,K))*DXB(I)+
     &                 (GD(I,J,K)-GD(I-1,J,K))*DXB(I+1))
               ENDDO
               DO I=1,locNX+1
                  FLX(I,J,K)= (ddy1(j)*FLMX(i,j,K)+ddy0(j)*
     &                 FLMX(i,j-1,K))*(GD(IND(I),J,K)
     &                 +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX+1
                  S(I)=(DDY1(J)*V(I,J,k)+DDY0(J)*V(I,J-1,K))*DELT
                  IND(I)=I
                  IF (S(I).GT.0.d0) IND(I)=I-1
               ENDDO
               DO I=0,locNX+1
                  GRAD1(I)=(GD(I,J,K)-GD(I-1,J,K))*
     &                 (GD(I+1,J,K)-GD(I,J,K))
                  GRAD(I)=0.d0
                  IF (GRAD1(I).GT.0.d0) GRAD(I)=two*GRAD1(I)/
     &                 ((GD(I+1,J,K)-GD(I,J,K))*DXB(I)+
     &                 (GD(I,J,K)-GD(I-1,J,K))*DXB(I+1))
               ENDDO
               DO I=1,locNX+1
                  FLX(I,J,K)= (ddy1(j)*FLMX(i,j,K)+ddy0(j)*
     &                 FLMX(i,j-1,K))*(GD(IND(I),J,K)
     &                 +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
               ENDDO               
            ENDDO
         ENDDO
      endif
**  Y-Direction
      DO K=1,locNZ
         DO I=1,locNX
            DO J=1,locNY
               S(J)=(DCY0(J+1)*G(I,J,K)+DCY1(J+1)*G(I,J+1,K))*DELT
               IND(J)=J+1
               IF (S(J).GT.0.d0) IND(J)=J
            ENDDO
            
            DO J=1,locNY+1
               GRAD1(J)=(GD(I,J,K)-GD(I,J-1,K))*
     &              (GD(I,J+1,K)-GD(I,J,K))
               GRAD(J)=0.d0
               IF (GRAD1(J).GT..0) GRAD(J)=two*GRAD1(J)/
     &              ((GD(I,J+1,K)-GD(I,J,K))*DYA(J)+
     &              (GD(I,J,K)-GD(I,J-1,K))*DYA(J+1))
            ENDDO
            DO J=1,locNY
               FLY(I,J,K)= (DCY0(j+1)*flmy(i,j,K)+DCY1(j+1)*
     %              flmy(i,j+1,K))
     1              * (GD(I,IND(J),K) 
     &              +(XYB(J)-XYA(IND(J))-half*S(J))*GRAD(IND(J)))
            ENDDO
         ENDDO
      ENDDO

!--- PASS FLY(*,locNY,*)->FLY(*,0,*) for calculation of GD
!---- Y-COORD, + DIRECTION
      call MPI_SENDRECV(FLY(1:locNX,locNY,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     FLY(1:locNX,0,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     COMM_CART,stat,ierr)

******************************************
**  Calculate New Y-Momentum density based on calculated fluxes
******************************************
      if(MPIlower.eq.MPI_PROC_NULL) then
c         print *,'INHERE'
c         stop
         DO K=1,locNZ
            DO J=2,locNY
               DO I=1,locNX
                  GD(I,J,K)=GD(I,J,K)*rhqy(i,j,K) 
     &                 -(FLX(I+1,J,K)-FLX(I,J,K))/volxb(i)
     &                 -(FLY(I,J,K)-FLY(I,J-1,K))/volya(j)
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  GD(I,J,K)=GD(I,J,K)*rhqy(i,j,K) 
     &                 -(FLX(I+1,J,K)-FLX(I,J,K))/volxb(i)
     &                 -(FLY(I,J,K)-FLY(I,J-1,K))/volya(j)
               ENDDO
            ENDDO
         ENDDO
      endif

******************************************************************
**  CALCULATE NEW VALUES OF PHYSICAL VARIABLES USED IN ONTHER ROUTINES
**   based on previous calculations of ED,VD, and GD 
**   However, this is the first calculation of RH
******************************************************************
**  DENSITY (has to come before the other variables)
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               RH(I,J,K)=RH(I,J,K)
     &              -(FLMX(I+1,J,K)-FLMX(I,J,K))/volxb(i)
     &              -(FLMY(I,J+1,K)-FLMY(I,J,K))/volyb(j)
            ENDDO
         ENDDO
      ENDDO

**  Boundary conditions (especially for densities, or velocities, if ietot=1)
      CALL BOUNDS
**  Density at Cell-Faces (and Boundary Conditions for density)
      CALL QDENS

**  Velocity X,Y,Z
      if(MPIleft.eq.MPI_PROC_NULL) then
         DO K=1,locNZ
            DO J=1,locNY
               DO I=2,locNX
                  V(I,J,K)=VD(I,J,K)/RHQX(I,J,K)
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  V(I,J,K)=VD(I,J,K)/RHQX(I,J,K)
               ENDDO
            ENDDO
         ENDDO
      endif

      if(MPIlower.eq.MPI_PROC_NULL) then
c         print *,'INHERE'
c         stop
         DO K=1,locNZ
            DO J=2,locNY
               DO I=1,locNX
                  G(I,J,K)=GD(I,J,K)/(RHQY(i,j,K)*geoxg(i)**2*
     %                 geozg(k)**2) - OMROT
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  G(I,J,K)=GD(I,J,K)/(RHQY(i,j,K)*geoxg(i)**2*
     %                 geozg(k)**2) - OMROT
               ENDDO
            ENDDO
         ENDDO
      endif

      RETURN
      END SUBROUTINE TESTADVECT

