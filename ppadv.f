      SUBROUTINE ADVECT2D
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
               
                  IF (V(I,J,K).GT..0) IND(I)=I-1
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
                  IF (GRAD1(J).GT..0) GRAD(J)=two*GRAD1(J)/
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
         stop 'this is dumb.  set iadv=0'
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
**  CALCULATE CHANGE IN INTERNAL ENERGY - IE THERMAL ENERGY DENSITY
******************************************
*** set isoth=0 to include energy advection
*** set isoth=1 to exclude energy advection
*** set isoth=3 to include the advective flow of thermal energy ...
***             but to exclude the work (P.dV) compressional heating in ttfor
      if (isoth.eq.0 .or. isoth.eq.3) then
**  Thermal energy density
         if (irht.eq.0) then
            ED=T
         else
            ED=T*RH*CV
         endif

**  X-Direction
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX+1
                  IND(I)=I
                  S(I)=V(I,J,K)*DELT
                  IF (V(I,J,K).GT..0) IND(I)=I-1
               ENDDO
               
               DO I=0,locNX+1
                  GRAD1(I)=(ED(I,J,K)-ED(I-1,J,K))*
     1                 (ED(I+1,J,K)-ED(I,J,K))
                  GRAD(I)=.0
                  IF (GRAD1(I).GT..0) GRAD(I)=two*GRAD1(I)/
     1                 ((ED(I,J,K)-ED(I-1,J,K))*DXB(I+1)+
     2                 (ED(I+1,J,K)-ED(I,J,K))*DXB(I))
               ENDDO            
c     these give the same results (if ed contains rh) (label T*RH)
               if (irht.eq.0) then
                  DO I=1,locNX+1
                     FLX(I,J,K)= flmx(i,j,k)*
     &                    (ED(IND(I),J,K)
     &                    +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
                  ENDDO
               else
                  DO I=1,locNX+1
                     FLX(I,J,K)= (surxa(i)-ccxa(i)*s(i)) * S(i) *
     &                    (ED(IND(I),J,K)
     &                    +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
                  ENDDO
               endif
            ENDDO
         ENDDO
         
**  Y-Direction
         DO K=1,locNZ
            DO I=1,locNX
               DO J=1,locNY+1
                  IND(J)=J
                  S(J)=G(I,J,K)*DELT
                  IF (G(I,J,K).GT..0) IND(J)=J-1
               ENDDO
               DO J=0,locNY+1
                  GRAD1(J)=(ED(I,J,K)-ED(I,J-1,K))*
     1                 (ED(I,J+1,K)-ED(I,J,K))
                  GRAD(J)=.0
                  IF (GRAD1(J).GT..0) GRAD(J)=two*GRAD1(J)/
     1                 ((ED(I,J,K)-ED(I,J-1,K))*DYB(J+1)+
     2                 (ED(I,J+1,K)-ED(I,J,K))*DYB(J))
               ENDDO
c     these give the same results (if ed contains rh) (label T*RH)
               if (irht.eq.0) then
                  DO J=1,locNY+1
                     FLY(I,J,K)= flmy(i,j,K)*
     &                    (ED(I,IND(J),K)
     &                    +(XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)))
                  ENDDO
               else
                  DO J=1,locNY+1
                     FLY(I,J,K)= (surya(j)-ccya(j)*s(j)) * S(j) *
     &                    (ED(I,IND(J),K)
     &                    +(XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)))
                  ENDDO
               endif 
            ENDDO
         ENDDO
******************************************
**  Calculate New Temperature (Density) based on the calculated fluxes
******************************************
         if (irht.eq.0) then
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     ED(I,J,K)=ED(I,J,K)*rh(i,j,k)
     &                    -(FLX(I+1,J,K)-FLX(I,J,K))/volxb(i)
     &                    -(FLY(I,J+1,K)-FLY(I,J,K))/volyb(j)
                  ENDDO
               ENDDO
            ENDDO   
         else
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     ED(I,J,K)=ED(I,J,K)
     &                    -(FLX(I+1,J,K)-FLX(I,J,K))/volxb(i)
     &                    -(FLY(I,J+1,K)-FLY(I,J,K))/volyb(j)
                  ENDDO
               ENDDO
            ENDDO
         endif
ci end of isoth if statement         
      ENDIF

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

**  Temperature
      if (isoth.eq.0 .or. isoth.eq.3) then
         if (irht.eq.0) then
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     T(i,j,k)=ED(i,j,k)/RH(i,j,k)
                  ENDDO
               ENDDO
            ENDDO
         else
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     T(i,j,k)=ED(i,j,k)/(RH(i,j,k)*CV(i,j,k))
                  ENDDO
               ENDDO
            ENDDO
         endif            
      endif
      RETURN
      END SUBROUTINE ADVECT2D

      SUBROUTINE ADVECT3D
**  3D Advection Terms  (Mass X variable, includes Coriolis forces (for Vphi)
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
      double precision :: HD(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: FLX(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLY(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLZ(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLMX(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLMY(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: FLMZ(0:locNX+1,0:locNY+1,0:locNZ+1)
      double precision :: S(NTOP),GRAD(0:NTOP),GRAD1(0:NTOP)
      double precision :: S_POLE,HD_NEG1(locNX,locNY)
      integer :: IND(NTOP),IND_POLE
      integer :: irht
      double precision,parameter ::half=0.5d0,one=1.0d0,two=2.0d0
      double precision,parameter :: xdmin=1.0d-50
      integer :: stat(MPI_STATUS_SIZE)
!--- SET THERMAL ENERGY TO BE ADVECTED: T*CV*RH
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
                  IF (V(I,J,K).GT..0) IND(I)=I-1
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
                  IF (GRAD1(J).GT..0) GRAD(J)=two*GRAD1(J)/
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

**  Z-Direction
         DO J=1,locNY
            DO I=1,locNX
               DO K=1,locNZ+1
                  IND(K)=K
                  S(K)=H(I,J,K)*DELT
                  IF (H(I,J,K).GT..0) IND(K)=K-1
               ENDDO              
               DO K=0,locNZ+1
                  GRAD1(K)=(RH(I,J,K)-RH(I,J,K-1))*
     &                 (RH(I,J,K+1)-RH(I,J,K))
                  GRAD(K)=0.d0
                  IF (GRAD1(K).GT.0.d0) GRAD(K)=two*GRAD1(K)/
     &                 ((RH(I,J,K)-RH(I,J,K-1))*DZB(K+1)+
     &                 (RH(I,J,K+1)-RH(I,J,K))*DZB(K))
               ENDDO
               DO K=1,locNZ+1
                  if(ncosys.eq.0.or.ncosys.eq.2) then
                     FLMZ(I,J,K)=(surza(k) + ccza(k)*s(k)) * S(k)
     &                    * (RH(I,J,ind(k))
     &                    +(XZA(k)-XZB(IND(k))-half*S(k))*GRAD(IND(k)))
                   else
c     Original version; This is a trig idenity, then a taylor expansion
c     utilizing the small value H*dt=S. The original appears to have a
c     sign error.. CHECK THIS FOR NON_SPHERICAL COORDINATES..
c     Note for ncosys=0 ccza=0, so this switch doesn't matter (this is true for ncosys=1)
                     print *,'Check the sign of surza(k)-ccza(k)'
                     print *,'STOP: advect3d'
                     call clean_stop
                     FLMZ(I,J,K)=(surza(k)-ccza(k)*s(k)) * S(k)
     &                    * (RH(I,J,ind(k))
     &                    +(XZA(k)-XZB(IND(k))-half*S(k))*GRAD(IND(k)))
                  endif 
               ENDDO
            ENDDO
         ENDDO

      else if (isoden.eq.1) then
         stop 'this is dumb.  set iadv=0'
      end if


!-MPI FOR FLMX,FLMY,FLMZ
!--- Pass FLMX(*,locNY,*)->FLMX(*,0,*)
!----- Y-COORD, + DIRECTION
      call MPI_SENDRECV(FLMX(1:locNX+1,locNY,1:locNZ),
     %     (locNX+1)*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     FLMX(1:locNX+1,0,1:locNZ),(locNX+1)*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIlower,1,COMM_CART,stat,ierr)
!--- Pass FLMX(*,*,locNZ)->FLMX(*,*,0)
!----- Z-COORD, + DIRECTION
      call MPI_SENDRECV(FLMX(1:locNX+1,1:locNY,locNZ),
     %     (locNX+1)*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     FLMX(1:locNX+1,1:locNY,0),(locNX+1)*locNY,
     %     MPI_DOUBLE_PRECISION,MPIbelow,1,COMM_CART,stat,ierr)
      if(poles.and.proc_coords(3).eq.0) then
!--  redo South pole, flmx defined on xzb, no sign flip
         call MPI_SENDRECV(FLMX(1:locNX+1,0:locNY,1),
     %        (locNX+1)*(1+locNY),MPI_DOUBLE_PRECISION,MPIbelow_pole,1,
     %        FLMX(1:locNX+1,0:locNY,0),(locNX+1)*(1+locNY),
     %        MPI_DOUBLE_PRECISION,MPIbelow_pole,1,COMM_CART,stat,ierr)
      endif

!--- Pass FLMY(locNX,*,*)->FLMY(0,*,*)
!----- X-COORD, + DIRECTION
      call MPI_SENDRECV(FLMY(locNX,1:locNY+1,1:locNZ),
     %     (locNY+1)*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %     FLMY(0,1:locNY+1,1:locNZ),(locNY+1)*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
!--- Pass FLMY(*,*,locNZ)->FLMY(*,*,0)
!----- Z-COORD, + DIRECTION
      call MPI_SENDRECV(FLMY(1:locNX,1:locNY+1,locNZ),
     %     locNX*(locNY+1),MPI_DOUBLE_PRECISION,MPIabove,1,
     %     FLMY(1:locNX,1:locNY+1,0),locNX*(locNY+1),
     %     MPI_DOUBLE_PRECISION,MPIbelow,1,COMM_CART,stat,ierr)
!--  redo South pole, flmy defined on xzb, use sign flip
      if(poles.and.proc_coords(3).eq.0) then
         call MPI_SENDRECV(-FLMY(1:locNX,1:locNY+1,1),
     %        locNX*(locNY+1),MPI_DOUBLE_PRECISION,MPIbelow_pole,1,
     %        FLMY(1:locNX,1:locNY+1,0),locNX*(locNY+1),
     %        MPI_DOUBLE_PRECISION,MPIbelow_pole,1,COMM_CART,stat,ierr)
      endif
c!added 4/20/2011
c!--  redo North pole, flmy defined on xzb, use sign flip
c      if(poles.and.proc_coords(3).eq.proc_dims(3)-1) then
c         call MPI_SENDRECV(-FLMY(1:locNX,1:locNY+1,locNZ-1),
c     %        locNX*(locNY+1),MPI_DOUBLE_PRECISION,MPIabove_pole,1,
c     %        FLMY(1:locNX,1:locNY+1,locNZ),locNX*(locNY+1),
c     %        MPI_DOUBLE_PRECISION,MPIabove_pole,1,COMM_CART,stat,ierr)
c      endif
!--- Pass FLMZ(locNX,*,*)->FLMZ(0,*,*)
!----- X-COORD, + DIRECTION
      call MPI_SENDRECV(FLMZ(locNX,1:locNY,1:locNZ+1),
     %     locNY*(locNZ+1),MPI_DOUBLE_PRECISION,MPIright,1,
     %     FLMZ(0,1:locNY,1:locNZ+1),locNY*(locNZ+1),
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
!--- Pass FLMZ(*,locNY,*)->FLMZ(*,0,*)
!----- Y-COORD, + DIRECTION
      call MPI_SENDRECV(FLMZ(1:locNX,locNY,1:locNZ+1),
     %     locNX*(locNZ+1),MPI_DOUBLE_PRECISION,MPIupper,1,
     %     FLMZ(1:locNX,0,1:locNZ+1),locNX*(locNZ+1),
     %     MPI_DOUBLE_PRECISION,MPIlower,1,COMM_CART,stat,ierr)

!--- Pass FLMZ(*,*,locNZ-1)->-FLMZ(*,*,locNZ+1)
!----- Z-COORD, + DIRECTION (OVER NORTHPOLE)
      if(poles.and.proc_coords(3).eq.proc_dims(3)-1) then
         call MPI_SENDRECV(-FLMZ(0:locNX+1,0:locNY+1,locNZ-1),
     %        (locNX+2)*(locNY+2),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,
     %        FLMZ(0:locNX+1,0:locNY+1,locNZ+1),
     %        (locNX+2)*(locNY+2),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,COMM_CART,stat,ierr)
      endif



c      DO K=0,locNZ+1
c         DO J=0,locNY+1
c            DO I=0,locNX+1
c               test_vadv(1,I,J,K)=FLMX(i,j,k)
c               test_vadv(2,I,J,K)=FLMY(i,j,k)
c               test_vadv(3,I,J,K)=FLMZ(i,j,k)
c            enddo
c         enddo
c      enddo


c   moved down here from about after 'icell' modification
      CALL QDENS

******************************************      
**  CALCULATE CHANGE IN INTERNAL ENERGY - IE THERMAL ENERGY DENSITY
******************************************

*** set isoth=3 to include the advective flow of thermal energy
***  but to exclude the work (P.dV) compressional heating in ttfor
      if (isoth.eq.0 .or. isoth.eq.3) then
         if (ietot.eq.0) then
**  Thermal energy density
            if (irht.eq.0) then
               ED=T
            else
               ED=T*RH*CV
            endif
         endif

**  X-Direction
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX+1
                  IND(I)=I
                  S(I)=V(I,J,K)*DELT
                  IF (V(I,J,K).GT..0) IND(I)=I-1
               ENDDO
               
               DO I=0,locNX+1
                  GRAD1(I)=(ED(I,J,K)-ED(I-1,J,K))*
     1                 (ED(I+1,J,K)-ED(I,J,K))
                  GRAD(I)=.0
                  IF (GRAD1(I).GT..0) GRAD(I)=two*GRAD1(I)/
     1                 ((ED(I,J,K)-ED(I-1,J,K))*DXB(I+1)+
     2                 (ED(I+1,J,K)-ED(I,J,K))*DXB(I))
               ENDDO
c     these give the same results (if ed contains rh) (label T*RH)
               if (irht.eq.0) then
                  DO I=1,locNX+1
                     FLX(I,J,K)= flmx(i,j,k)*
     &                    (ED(IND(I),J,K)
     &                    +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
                  ENDDO
               else
                  DO I=1,locNX+1
                     FLX(I,J,K)= (surxa(i)-ccxa(i)*s(i)) * S(i) *
     &                    (ED(IND(I),J,K)
     &                    +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
                  ENDDO
               endif
            ENDDO
         ENDDO
         
**  Y-Direction
         DO K=1,locNZ
            DO I=1,locNX
               DO J=1,locNY+1
                  IND(J)=J
                  S(J)=G(I,J,K)*DELT
                  IF (G(I,J,K).GT..0) IND(J)=J-1
               ENDDO
               DO J=0,locNY+1
                  GRAD1(J)=(ED(I,J,K)-ED(I,J-1,K))*
     1                 (ED(I,J+1,K)-ED(I,J,K))
                  GRAD(J)=.0
                  IF (GRAD1(J).GT..0) GRAD(J)=two*GRAD1(J)/
     1                 ((ED(I,J,K)-ED(I,J-1,K))*DYB(J+1)+
     2                 (ED(I,J+1,K)-ED(I,J,K))*DYB(J))
               ENDDO
c     these give the same results (if ed contains rh) (label T*RH)
               if (irht.eq.0) then
                  DO J=1,locNY+1
                     FLY(I,J,K)= flmy(i,j,K)*
     &                    (ED(I,IND(J),K)
     &                    +(XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)))
                  ENDDO
               else
                  DO J=1,locNY+1
                     FLY(I,J,K)= (surya(j)-ccya(j)*s(j)) * S(j) *
     &                    (ED(I,IND(J),K)
     &                    +(XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)))
                  ENDDO
               endif 
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
                  GRAD1(K)=(ED(I,J,K)-ED(I,J,K-1))*
     1                 (ED(I,J,K+1)-ED(I,J,K))
                  GRAD(K)=.0
                  IF (GRAD1(K).GT..0) GRAD(K)=two*GRAD1(K)/
     1                 ((ED(I,J,K)-ED(I,J,K-1))*DZB(K+1)+
     2                 (ED(I,J,K+1)-ED(I,J,K))*DZB(K))
               ENDDO
c     these give the same results (if ed contains rh) (label T*RH)
               if (irht.eq.0) then 
                  DO K=1,locNZ+1
                     FLZ(I,J,K)= flmz(i,j,k)*
     &                    (ED(I,J,ind(k))
     &                    +(XZA(k)-XZB(IND(k))-half*S(k))*GRAD(IND(k)))
                  ENDDO
               else
                  DO K=1,locNZ+1
                     FLZ(I,J,K)= (surza(k)-ccza(k)*s(k)) * S(k) *
     &                    (ED(I,J,ind(k))
     &                    +(XZA(k)-XZB(IND(k))-half*S(k))*GRAD(IND(k)))
                  ENDDO
               endif               
            ENDDO
         ENDDO


!-Save the thermal advection flux
         DO K=0,locNZ+1
            DO J=0,locNY+1
               DO I=0,locNX+1
                  conv_flux(1,I,J,K)=FLX(i,j,k)
                  conv_flux(2,I,J,K)=FLY(i,j,k)
                  conv_flux(3,I,J,K)=FLZ(i,j,k)
               enddo
            enddo
         enddo

******************************************
**  Calculate New Temperature (Density) based on the calculated fluxes
******************************************
         if (irht.eq.0) then
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     ED(I,J,K)=ED(I,J,K)*rh(i,j,k)
     &                    -(FLX(I+1,J,K)-FLX(I,J,K))/volxb(i)
     &                    -(FLY(I,J+1,K)-FLY(I,J,K))/volyb(j)
     &                    -(FLZ(I,J,K+1)-FLZ(I,J,K))/volzb(k)
                  ENDDO
               ENDDO
            ENDDO   
         else
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     ED(I,J,K)=ED(I,J,K)
     &                    -(FLX(I+1,J,K)-FLX(I,J,K))/volxb(i)
     &                    -(FLY(I,J+1,K)-FLY(I,J,K))/volyb(j)
     &                    -(FLZ(I,J,K+1)-FLZ(I,J,K))/volzb(k)
                  ENDDO
               ENDDO
            ENDDO
         endif
ci end of isoth if statement         
      ENDIF

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

**  Z-Direction
      if(MPIleft.eq.MPI_PROC_NULL) then
         DO J=1,locNY
            DO I=2,locNX
               DO K=1,locNZ+1
                  S(K)=(DDX1(i)*H(I,J,K)+DDX0(i)*H(I-1,J,K))*DELT
                  IND(K)=K
                  IF (S(K).GT..0) IND(K)=K-1
               ENDDO
               DO K=0,locNZ+1
                  GRAD1(K)=(V(I,J,K)-V(I,J,K-1))*
     1                 (V(I,J,K+1)-V(I,J,K))
                  GRAD(K)=.0
                  IF (GRAD1(K).GT..0) GRAD(K)=two*GRAD1(K)/
     1                 ((V(I,J,K+1)-V(I,J,K))*DZB(K)+
     2                 (V(I,J,K)-V(I,J,K-1))*DZB(K+1))
               ENDDO
               DO K=1,locNZ+1
                  FLZ(I,J,K)=(ddx1(i)*flmz(i,j,k)+ddx0(i)*flmz(i-1,j,k))
     &                 *( V(I,J,ind(k))+
     1                 (XZA(K)-XZB(IND(K))-half*S(K))*GRAD(IND(K)) )
               ENDDO
            ENDDO
         ENDDO
      else
         DO J=1,locNY
            DO I=1,locNX
               DO K=1,locNZ+1
                  S(K)=(DDX1(i)*H(I,J,K)+DDX0(i)*H(I-1,J,K))*DELT
                  IND(K)=K
                  IF (S(K).GT..0) IND(K)=K-1
               ENDDO
               DO K=0,locNZ+1
                  GRAD1(K)=(V(I,J,K)-V(I,J,K-1))*
     1                 (V(I,J,K+1)-V(I,J,K))
                  GRAD(K)=.0
                  IF (GRAD1(K).GT..0) GRAD(K)=two*GRAD1(K)/
     1                 ((V(I,J,K+1)-V(I,J,K))*DZB(K)+
     2                 (V(I,J,K)-V(I,J,K-1))*DZB(K+1))
               ENDDO
               DO K=1,locNZ+1
                  FLZ(I,J,K)=(ddx1(i)*flmz(i,j,k)+ddx0(i)*flmz(i-1,j,k))
     &                 *( V(I,J,ind(k))+
     1                 (XZA(K)-XZB(IND(K))-half*S(K))*GRAD(IND(K)) )
               ENDDO
            ENDDO
         ENDDO
      endif

!--- Pass FLX(locNX,*,*)->FLX(0,*,*) for below calculation of VD from i=1,locNX
!---- X-COORD, + DIRECTION
      call MPI_SENDRECV(FLX(locNX,1:locNY,1:locNZ),
     %     locNY*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %     FLX(0,1:locNY,1:locNZ),locNY*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)

c      DO K=0,locNZ+1
c         DO J=0,locNY+1
c            DO I=0,locNX+1
c               test_vadv(1,I,J,K)=FLX(i,j,k)
c               test_vadv(2,I,J,K)=FLY(i,j,k)
c               test_vadv(3,I,J,K)=FLZ(i,j,k)
c            enddo
c         enddo
c      enddo

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
     &                 - (FLZ(I,J,K+1)-FLZ(I,J,K))/volzb(k)

c                  if(myid.eq.19.and.j.eq.5.and.k.eq.3.and.I.eq.93) then
c                     write(*,'(A,I8,I8,10(1x,e11.5))') 'VD:',
c     $                    iii,i,VD(I,J,K),V(I,J,K),
c     $                    VD(I,J,K)/RHQX(I,J,K),
c     &                    rhqx(i,j,K),FLX(I,J,K),FLX(I-1,J,K),
c     &                    FLY(I,J+1,K),FLY(I,J,K),
c     &                    FLZ(I,J,K+1),FLZ(I,J,K)
c                  endif


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
     &                 - (FLZ(I,J,K+1)-FLZ(I,J,K))/volzb(k)
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

!--pass info for GD
      if(modtyp.eq.2) CALL PLANETSHIFT(11,-1)
!--pass info for GD across the poles
      if(poles.and.proc_coords(3).eq.proc_dims(3)-1) then
         CALL NORTH_POLE_SHIFT(11,-1)
      endif
      if(poles.and.proc_coords(3).eq.0) then
         CALL SOUTH_POLE_SHIFT(11,-1)
      endif

**  X-Direction
      if(MPIlower.eq.MPI_PROC_NULL) then
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


c               if(myid.eq.12.and.k.eq.1.and.i.eq.150) then
c                  write(*,'(A,I8,I8,10(1x,e12.3))') 'FLY12=',III,ind(j),
c     $                 FLY(I,J,K),flmy(i,j,K),
c     %                 flmy(i,j+1,K),
c     1                 GD(I,IND(J),K), 
c     &                 XYB(J),XYA(IND(J)),S(J),GRAD(IND(J)),
c     $                 G(I,J,K),G(I,J+1,K)
c               endif
            ENDDO
         ENDDO
      ENDDO

**  Z-Direction
      if(MPIlower.eq.MPI_PROC_NULL) then
         DO I=1,locNX
            DO J=2,locNY
               DO K=1,locNZ+1
                  S(K)=(DDY1(J)*H(I,J,K)+DDY0(J)*H(I,J-1,K))*DELT
                  IND(K)=K
                  IF (S(k).GT..0) IND(k)=k-1
               ENDDO
               DO K=0,locNZ+1
                  GRAD1(k)=(GD(I,J,K)-GD(I,J,K-1))*
     1                 (GD(I,J,K+1)-GD(I,J,K))
                  GRAD(K)=.0
                  IF (GRAD1(K).GT..0) GRAD(K)=two*GRAD1(K)/
     1                 ((GD(I,J,K+1)-GD(I,J,K))*DZB(K)+
     2                 (GD(I,J,K)-GD(I,J,K-1))*DZB(K+1))
               ENDDO
               DO K=1,locNZ+1
                  FLZ(I,J,K)= (ddy1(j)*FLMZ(i,j,k)+ddy0(j)*
     &                 FLMZ(i,j-1,k))*(GD(I,J,IND(k))
     &                 +(XZA(K)-XZB(IND(K))-half*S(K))*GRAD(IND(K)))
               ENDDO
            ENDDO
         ENDDO
      else
         DO I=1,locNX
            DO J=1,locNY
               DO K=1,locNZ+1
                  S(K)=(DDY1(J)*H(I,J,K)+DDY0(J)*H(I,J-1,K))*DELT
                  IND(K)=K
                  IF (S(k).GT..0) IND(k)=k-1
               ENDDO
               DO K=0,locNZ+1
                  GRAD1(k)=(GD(I,J,K)-GD(I,J,K-1))*
     1                 (GD(I,J,K+1)-GD(I,J,K))
                  GRAD(K)=.0
                  IF (GRAD1(K).GT..0) GRAD(K)=two*GRAD1(K)/
     1                 ((GD(I,J,K+1)-GD(I,J,K))*DZB(K)+
     2                 (GD(I,J,K)-GD(I,J,K-1))*DZB(K+1))
               ENDDO
               DO K=1,locNZ+1
                  FLZ(I,J,K)= (ddy1(j)*FLMZ(i,j,k)+ddy0(j)*
     &                 FLMZ(i,j-1,k))*(GD(I,J,IND(k))
     &                 +(XZA(K)-XZB(IND(K))-half*S(K))*GRAD(IND(K)))
               ENDDO
            ENDDO
         ENDDO
      endif

!--- PASS FLY(*,locNY,*)->FLY(*,0,*) for calculation of GD
!---- Y-COORD, + DIRECTION
      call MPI_SENDRECV(FLY(1:locNX,locNY,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     FLY(1:locNX,0,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     COMM_CART,stat,ierr)

c      DO K=0,locNZ+1
c         DO J=0,locNY+1
c            DO I=0,locNX+1
c               test_vadv(1,I,J,K)=FLX(i,j,k)
c               test_vadv(2,I,J,K)=FLY(i,j,k)
c               test_vadv(3,I,J,K)=FLZ(i,j,k)
c            enddo
c         enddo
c      enddo


******************************************
**  Calculate New Y-Momentum density based on calculated fluxes
******************************************
      if(MPIlower.eq.MPI_PROC_NULL) then
         DO K=1,locNZ
            DO J=2,locNY
               DO I=1,locNX
                  GD(I,J,K)=GD(I,J,K)*rhqy(i,j,K) 
     &                 -(FLX(I+1,J,K)-FLX(I,J,K))/volxb(i)
     &                 -(FLY(I,J,K)-FLY(I,J-1,K))/volya(j)
     &                 -(FLZ(I,J,K+1)-FLZ(I,J,K))/volzb(k)
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
     &                 -(FLZ(I,J,K+1)-FLZ(I,J,K))/volzb(k)
               ENDDO
            ENDDO
         ENDDO
      endif
      
******************************************
**  CALCULATE CHANGE IN Z-MOMENTUM DENSITY (eg. MERIDIONAL)
******************************************
      DO K=0,locNZ+2
         DO J=-1,locNY+2
            DO I=-1,locNX+2
               HD(I,J,K)=H(I,J,K) *geoxh(i)**2
            ENDDO
         ENDDO
      ENDDO
      
**  X-Direction
      if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
         DO K=2,locNZ
            DO J=1,locNY            
               DO I=1,locNX+1
                  S(I)=(DDZ1(K)*V(I,J,K)+DDZ0(K)*V(I,J,K-1))*DELT
                  IND(I)=I
                  IF (S(I).GT..0) IND(I)=I-1
               ENDDO
               DO I=0,locNX+1
                  GRAD1(I)=(HD(I,J,K)-HD(I-1,J,K))*
     1                 (HD(I+1,J,K)-HD(I,J,K)) 
                  GRAD(I)=.0
                  IF (GRAD1(I).GT..0) GRAD(I)=two*GRAD1(I)/
     1                 ((HD(I+1,J,K)-HD(I,J,K))*DXB(I)+
     2                 (HD(I,J,K)-HD(I-1,J,K))*DXB(I+1))
               ENDDO
               DO I=1,locNX+1
                  FLX(I,J,K)= (ddz1(k)*FLMX(i,j,k)+ddz0(k)*
     &                 FLMX(i,j,k-1))*(HD(IND(I),J,K)
     &                 +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO J=1,locNY            
               DO I=1,locNX+1
                  S(I)=(DDZ1(K)*V(I,J,K)+DDZ0(K)*V(I,J,K-1))*DELT
                  IND(I)=I
                  IF (S(I).GT..0) IND(I)=I-1
               ENDDO
               DO I=0,locNX+1
                  GRAD1(I)=(HD(I,J,K)-HD(I-1,J,K))*
     1                 (HD(I+1,J,K)-HD(I,J,K)) 
                  GRAD(I)=.0
                  IF (GRAD1(I).GT..0) GRAD(I)=two*GRAD1(I)/
     1                 ((HD(I+1,J,K)-HD(I,J,K))*DXB(I)+
     2                 (HD(I,J,K)-HD(I-1,J,K))*DXB(I+1))
               ENDDO
               DO I=1,locNX+1
                  FLX(I,J,K)= (ddz1(k)*FLMX(i,j,k)+ddz0(k)*
     &                 FLMX(i,j,k-1))*(HD(IND(I),J,K)
     &                 +(XXA(I)-XXB(IND(I))-half*S(I))*GRAD(IND(I)))
               ENDDO
            ENDDO
         ENDDO
      endif

**  Y-Direction
      if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
         DO K=2,locNZ
            DO I=1,locNX
               DO J=1,locNY+1
                  S(J)=(DDZ1(K)*G(I,J,K)+DDZ0(K)*G(I,J,K-1))*DELT
                  IND(J)=J
                  IF (S(J).GT..0) IND(J)=J-1
               ENDDO
               DO J=0,locNY+1
                  GRAD1(J)=(HD(I,J,K)-HD(I,J-1,K))*
     1                 (HD(I,J+1,K)-HD(I,J,K))
                  GRAD(J)=.0
                  IF (GRAD1(J).GT..0) GRAD(J)=two*GRAD1(J)/
     1                 ((HD(I,J+1,K)-HD(I,J,k))*DYB(J)+
     2                 (HD(I,J,K)-HD(I,J-1,K))*DYB(J+1))
               ENDDO            
               DO J=1,locNY+1
                  FLY(I,J,K)=(ddz1(k)*flmy(i,j,k)+ddz0(k)*flmy(i,j,k-1))
     &                 *( HD(I,IND(J),k)+
     1                 (XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)) )
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO I=1,locNX
               DO J=1,locNY+1
                  S(J)=(DDZ1(K)*G(I,J,K)+DDZ0(K)*G(I,J,K-1))*DELT
                  IND(J)=J
                  IF (S(J).GT..0) IND(J)=J-1
               ENDDO
               DO J=0,locNY+1
                  GRAD1(J)=(HD(I,J,K)-HD(I,J-1,K))*
     1                 (HD(I,J+1,K)-HD(I,J,K))
                  GRAD(J)=.0
                  IF (GRAD1(J).GT..0) GRAD(J)=two*GRAD1(J)/
     1                 ((HD(I,J+1,K)-HD(I,J,k))*DYB(J)+
     2                 (HD(I,J,K)-HD(I,J-1,K))*DYB(J+1))
               ENDDO            
               DO J=1,locNY+1
                  FLY(I,J,K)=(ddz1(k)*flmy(i,j,k)+ddz0(k)*flmy(i,j,k-1))
     &                 *( HD(I,IND(J),k)+
     1                 (XYA(J)-XYB(IND(J))-half*S(J))*GRAD(IND(J)) )
!     At the northern/southern most grid, force FLY=0. If this was really
!     +/-(pi/2) then G=0, so both there is no advection of H in by G
                  if(poles.and.proc_coords(3).eq.proc_dims(3)-1.and.
     $                 k.eq.locNZ) then
                     FLY(i,j,k) = 0.d0
                  endif
                  if(poles.and.proc_coords(3).eq.0.and.k.eq.1) then
                     FLY(i,j,k)= 0.d0
                  endif
               ENDDO
            ENDDO
         ENDDO
      endif

**  Z-Direction
      DO J=1,locNY
         DO I=1,locNX
            DO K=1,locNZ
               S(K)=half*(H(I,J,K)+H(I,J,K+1))*DELT
               IND(K)=K+1
               IF (S(K).GT..0) IND(K)=K
            ENDDO
            DO K=1,locNZ+1
               GRAD1(K)=(HD(I,J,K)-HD(I,J,K-1))*
     1              (HD(I,J,K+1)-HD(I,J,K))
               GRAD(K)=.0
               IF (GRAD1(K).GT..0) GRAD(K)=two*GRAD1(K)/
     1              ((HD(I,J,K+1)-HD(I,J,K))*DZA(K)+
     2              (HD(I,J,K)-HD(I,J,K-1))*DZA(K+1))
            ENDDO
            DO K=1,locNZ
               FLZ(I,J,K)= half*(flmz(i,j,k)+flmz(i,j,k+1)) *
     1              (HD(I,J,ind(k))
     &              +(XZB(K)-XZA(IND(K))-half*S(K))*GRAD(IND(K)))
            ENDDO
         ENDDO
      ENDDO


!--- PASS FLZ(*,*,locNZ)->FLZ(*,*,0) for calculation of HD
!---- Z-COORD, + DIRECTION
      call MPI_SENDRECV(FLZ(0:locNX+1,0:locNY+1,locNZ),
     %     (locNX+2)*(locNY+2),MPI_DOUBLE_PRECISION,MPIabove,1,
     %     FLZ(0:locNX+1,0:locNY+1,0),
     %     (locNX+2)*(locNY+2),MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     COMM_CART,stat,ierr)
!--  redo South Pole. For H, FLZ defined on xzb. NO sign flip      
      if(poles.and.proc_coords(3).eq.0) then
         call MPI_SENDRECV(FLZ(0:locNX+1,0:locNY+1,1),
     %        (locNX+2)*(locNY+2),MPI_DOUBLE_PRECISION,MPIbelow_pole,1,
     %        FLZ(0:locNX+1,0:locNY+1,0),
     %        (locNX+1)*(locNY+1),MPI_DOUBLE_PRECISION,MPIbelow_pole,1,
     %        COMM_CART,stat,ierr)
      endif

******************************************
**  Calculate New Z-Momentum density based on calculated fluxes
******************************************
      if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
         DO K=2,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  HD(I,J,K)=HD(I,J,K)*rhqz(i,j,k) 
     &                 -(FLX(I+1,J,K)-FLX(I,J,K))/volxb(i)
     &                 -(FLY(I,J+1,K)-FLY(I,J,K))/volyb(j)
     &                 -(FLZ(I,J,K)-FLZ(I,J,K-1))/volza(k)
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  HD(I,J,K)=HD(I,J,K)*rhqz(i,j,k) 
     &                 -(FLX(I+1,J,K)-FLX(I,J,K))/volxb(i)
     &                 -(FLY(I,J+1,K)-FLY(I,J,K))/volyb(j)
     &                 -(FLZ(I,J,K)-FLZ(I,J,K-1))/volza(k)
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
c               if(myid.eq.19.and.j.eq.6.and.k.eq.6.and.i.eq.100) then
c                  write(*,'(A,I8,6(1x,e12.6))') ,'IN ADV',III,
c     $                 RH(I,J,K),
c     &                 -(FLMX(I+1,J,K)-FLMX(I,J,K))/volxb(i),
c     &                 FLMX(I+1,J,K),FLMX(I,J,K),
c     &                 -(FLMY(I,J+1,K)-FLMY(I,J,K))/volyb(j),
c     &                 -(FLMZ(I,J,K+1)-FLMZ(I,J,K))/volzb(k)
c               endif
               RH(I,J,K)=RH(I,J,K)
     &              -(FLMX(I+1,J,K)-FLMX(I,J,K))/volxb(i)
     &              -(FLMY(I,J+1,K)-FLMY(I,J,K))/volyb(j)
     &              -(FLMZ(I,J,K+1)-FLMZ(I,J,K))/volzb(k)
            ENDDO
         ENDDO
      ENDDO 

**  Boundary conditions (especially for densities, or velocities, if ietot=1)
      CALL BOUNDS
**  Density at Cell-Faces (and Boundary Conditions for density)
      CALL QDENS

**  Final Velocity (X)
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
**  Final Velocity (Y)
      if(MPIlower.eq.MPI_PROC_NULL) then
         DO J=2,locNY
            DO I=1,locNX
               DO K=1,locNZ
                  G(I,J,K)=GD(I,J,K)/(RHQY(i,j,K)*geoxg(i)**2*
c     %                 geozg(k)**2)
     %                 geozg(k)**2) - OMROT
               ENDDO
            ENDDO
         ENDDO
      else
         DO I=1,locNX
            DO J=1,locNY
               DO K=1,locNZ
                  G(I,J,K)=GD(I,J,K)/(RHQY(i,j,K)*geoxg(i)**2*
c     %                 geozg(k)**2)
     %                 geozg(k)**2) - OMROT
               ENDDO
            ENDDO
         ENDDO
      endif
**  Final Velocity (Z)
      if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
         DO K=2,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  H(I,J,K)=HD(I,J,K)/(RHQZ(I,J,K)*geoxh(i)**2)
               ENDDO
            ENDDO
         ENDDO
      else
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  H(I,J,K)=HD(I,J,K)/(RHQZ(I,J,K)*geoxh(i)**2)
               ENDDO
            ENDDO
         ENDDO
      endif
**  Temperature
      if (isoth.eq.0 .or. isoth.eq.3) then
         if (ietot.eq.0) then
            if (irht.eq.0) then
               DO K=1,locNZ
                  DO J=1,locNY
                     DO I=1,locNX
                        T(i,j,k)=ED(i,j,k)/RH(i,j,k)
                     ENDDO
                  ENDDO
               ENDDO
            else
               DO K=1,locNZ
                  DO J=1,locNY
                     DO I=1,locNX
                        T(i,j,k)=ED(i,j,k)/(RH(i,j,k)*CV(i,j,k))
                     ENDDO
                  ENDDO
               ENDDO
            endif            
         else if (ietot.eq.2.or.ietot.eq.1) then
            print *,'removed from advction'
            stop
         endif
      endif
      RETURN
      END SUBROUTINE ADVECT3D
