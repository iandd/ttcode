
      SUBROUTINE VISC
!-- VISCOSITY ROUTINE
      USE input_init
      USE force_var_init
      USE fluid_var_init
      USE grid_var_init
      USE mpi_var_init
      IMPLICIT NONE
      INTEGER I,J,K
      DOUBLE PRECISION :: DIVV(0:locNX+1,0:locNY+1,0:locNZ+1)
      DOUBLE PRECISION :: TXYAV(0:locNX+1,0:locNY+1,0:locNZ+1)
      DOUBLE PRECISION :: TYZAV(0:locNX+1,0:locNY+1,0:locNZ+1)
      DOUBLE PRECISION :: TXZAV(0:locNX+1,0:locNY+1,0:locNZ+1)
      DOUBLE PRECISION :: TXX(0:locNX+1,0:locNY+1,0:locNZ+1)
      DOUBLE PRECISION :: TXY(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      DOUBLE PRECISION :: TYY(0:locNX+1,0:locNY+1,0:locNZ+1)
      DOUBLE PRECISION :: TYZ(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      DOUBLE PRECISION :: TZZ(0:locNX+1,0:locNY+1,0:locNZ+1)
      DOUBLE PRECISION :: TXZ(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      DOUBLE PRECISION :: alprat
      character visc_file*8
      character proc_num*3
      SELECT CASE(IVISTYP)
      CASE (0)
         DIVTX = 0.d0
         DIVTY = 0.d0
         DIVTZ = 0.d0
      CASE (1)
         alprat=1.d0
*  Cartesian
*    TXX can be defined (0:NX+1,-1:NY+2,-1:NZ+2) (not dimensioned though)
*    DIVV is defined (0:NX+1,0:NY+1,0:NZ+1)
         DO K=0,locNZ+1
            DO J=0,locNY+1
               DO I=0,locNX+1
                  TXX(I,J,K) = 2.d0*ETA(I,J,K)*
     &                 (V(I+1,J,K)-V(I,J,K))/(XXA(I+1)-XXA(I))
                  DIVV(I,J,K) = 2.d0/3.d0 * ETA(I,J,K)*(
     &               (surxa(i+1)*V(I+1,J,k)-surxa(i)*V(I,J,k))/volxb(i)+
     &               (surya(j+1)*G(I,J+1,k)-surya(j)*G(I,J,k))/volyb(j)+
     &               (surza(k+1)*H(I,J,k+1)-surza(k)*H(I,J,k))/volzb(k))
               ENDDO
            ENDDO
         ENDDO
      
c--------------------------
c CARTESIAN COORDINATES
c--------------------------
         if (ncosys.eq.0) then
            DO K=0,locNZ+1
               DO J=0,locNY+1
                  DO I=0,locNX+1
                     TYY(I,J,K) = 2.d0*ETA(I,J,K)*
     &                    (G(I,J+1,K)-G(I,J,K))/(XYA(J+1)-XYA(J))
                     TZZ(I,J,K) =  2.d0*ETA(I,J,K)*
     &                    (H(I,J,K+1)-H(I,J,K))/(XZA(K+1)-XZA(K))
                  ENDDO
               ENDDO
            ENDDO
c--------------------------
c CYLINDRICAL COORDINATES
c--------------------------
         else if (ncosys.eq.1) then
            DO K=0,locNZ+1
               DO J=0,locNY+1
                  DO I=0,locNX+1
                     TYY(I,J,K) = 2.d0*ETA(I,J,K)*(
     &                    (G(I,J+1,K)-G(I,J,K))/(XYA(J+1)-XYA(J))+
     &                    (V(I+1,J,K)+V(I,J,K))/2.d0/XXB(I))
                     TZZ(I,J,K) =  2.d0*ETA(I,J,K)*
     &                    (H(I,J,K+1)-H(I,J,K))/(XZA(K+1)-XZA(K))
                  ENDDO
               ENDDO
            ENDDO
c--------------------------
c SPHERICAL COORDINATES
c--------------------------
         else if (ncosys.eq.2) then
            DO K=0,locNZ+1
               DO J=0,locNY+1
                  DO I=0,locNX+1
                     TYY(I,J,K) = 2.d0*ETA(I,J,K)*(
     &                    (G(I,J+1,K)-G(I,J,K))/(XYA(J+1)-XYA(J))+
     &                    (V(I+1,J,K)+V(I,J,K))/2.d0/XXB(I)+
     &                    (H(I,J,K+1)+H(I,J,K))/2.d0*(-TAN(XZB(K))))
                     TZZ(I,J,K) =  2.d0*ETA(I,J,K)*(
     &                    (H(I,J,K+1)-H(I,J,K))/(XZA(K+1)-XZA(K))+
     &                    (V(I+1,J,K)+V(I,J,K))/2.d0/XXB(I))
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

c--------------------------
c ALL COORDINATES
c--------------------------
*    TXY can be defined (0:NX+2,0:NY+2,-1:NZ+2)
         DO K=0,locNZ+2
            DO J=0,locNY+2
               DO I=0,locNX+2
                  TXY(I,J,K) = 0.5 * (
     &                 (V(I,J,K)-V(I,J-1,K))/(XYB(J)-XYB(J-1))
     &                 /geoxga(i)/geozga(k)+
     &                 (G(I,J,K)-G(I-1,J,K))/(XXB(I)-XXB(I-1))
     &                 *geoxga(i)*geozga(k)/alprat)
                  TXZ(I,J,K) = 0.5 * (
     &                 (V(I,J,K)-V(I,J,K-1))/(XZB(K)-XZB(K-1))
     &                 /geoxha(i)+
     &                 (H(I,J,K)-H(I-1,J,K))/(XXB(I)-XXB(I-1))
     &                 *geoxha(i))
                  if (geoxha(i).lt.1.e-10) then
                     TXZ(I,J,K) = 0.d0
                  endif
                  TYZ(I,J,K) = 0.5 * (
     &                 (H(I,J,K)-H(I,J-1,K))/(XYB(J)-XYB(J-1))
     &                 /geoxg(i)/geozga(k)*geoxh(i)+
     &                 (G(I,J,K)-G(I,J,K-1))/(XZB(K)-XZB(K-1))
     &                 *geoxg(i)*geozga(k)/geoxh(i)
     &                 /alprat)
               ENDDO
            ENDDO
         ENDDO
         
c         if (ivis.eq.5) then
c* no shear viscosity
c            TXY = 0.d0
c            TYZ = 0.d0
c         endif
      
*   TXYAV can be defined (0:NX+1,0:NY+1,-1:NZ+2)
         DO K=0,locNZ+1
            DO J=0,locNY+1
               DO I=0,locNX+1
                  TXYAV(I,J,K)= ETA(I,J,K)*(
     &                 TXY(I,J,K)+TXY(I+1,J,K)+TXY(I,J+1,K)+
     &                 TXY(I+1,J+1,K))/2.d0
                  TYZAV(I,J,K)= ETA(I,J,K)*(
     &                 TYZ(I,J,K)+TYZ(I,J+1,K)+TYZ(I,J,K+1)+
     &                 TYZ(I,J+1,K+1))/2.d0
                  TXZAV(I,J,K)= ETA(I,J,K)*(
     &                 TXZ(I,J,K)+TXZ(I+1,J,K)+TXZ(I,J,K+1)+
     &                 TXZ(I+1,J,K+1))/2.d0
               ENDDO
            ENDDO
         ENDDO
         
         DO K=0,locNZ+2
            DO J=0,locNY+2
               DO I=0,locNX+2
                  TXY(I,J,K)=TXY(I,J,K)/2.d0*(
     &              ETA(I,J,K)+ETA(I-1,J,K)+ETA(I,J-1,K)+ETA(I-1,J-1,K))
                  TXZ(I,J,K)=TXZ(I,J,K)/2.d0*(
     &              ETA(I,J,K)+ETA(I-1,J,K)+ETA(I,J,K-1)+ETA(I-1,J,K-1))
                  TYZ(I,J,K)=TYZ(I,J,K)/2.d0*(
     &              ETA(I,J,K)+ETA(I,J,K-1)+ETA(I,J-1,K)+ETA(I,J-1,K-1))
               ENDDO
            ENDDO
         ENDDO
         
*     DISFN can be defined (0:NX+1,0:NY+1,0:NZ+1)
         if(IDISS.ne.0) then
            DO K=0,locNZ+1
               DO J=0,locNY+1
                  DO I=0,locNX+1
                     DISFN(I,J,K)= (3.d0*DIVV(I,J,K)*DIVV(I,J,K)
     &                    +TXX(I,J,K)*TXX(I,J,K)+2.d0*TXYAV(I,J,K)*
     &                    TXYAV(I,J,K)*alprat
     &                    +TYY(I,J,K)*TYY(I,J,K)+2.d0*TYZAV(I,J,K)*
     &                    TYZAV(I,J,K)*alprat
     &                    +TZZ(I,J,K)*TZZ(I,J,K)+2.d0*TXZAV(I,J,K)*
     &                    TXZAV(I,J,K))/2.d0/ETA(I,J,K)              
                  ENDDO
               ENDDO
            ENDDO
         endif

c--------------------------
c CARTESIAN COORDINATES
c--------------------------      
*    DIVTX can be defined (1:NX+1,0:NY+1,0:NZ+1)
* All quantities (divv, txx, etc) are defined in grid centers. Averages
*   are used to find it at grid edges
* Note: negative sign for divv is incorperated by fliping the difference
         if (ncosys.eq.0) then
            DO K=0,locNZ+1
               DO J=0,locNY+1
                  DO I=1,locNX+1
                     DIVTX(I,J,K)=
     &          (surxb(i-1)*DIVV(I-1,J,K)-surxb(i)*DIVV(I,J,K))/volxa(i)
     &          + (surxb(i)*TXX(I,J,K)-surxb(i-1)*TXX(I-1,J,K))/volxa(i)
     &          + (surya(j+1)*TXY(I,J+1,K)-surya(j)*TXY(I,J,K))/volyb(j)
     &                    /geoxga(i)/geozg(k)
     &          + (surza(k+1)*TXZ(I,J,K+1)-surza(k)*TXZ(I,J,K))/volzb(k)
     &                    /geoxha(i)
                  ENDDO
               ENDDO
            ENDDO
            DO K=0,locNZ+1
               DO J=1,locNY+1
                  DO I=0,locNX+1
                     DIVTY(I,J,K)=
     &          (suryb(j-1)*DIVV(I,J-1,K)-suryb(j)*DIVV(I,J,K))/volya(j)
     &                    /geoxg(i)/geozg(k)
     &          + (surxa(i+1)*TXY(I+1,J,K)-surxa(i)*TXY(I,J,K))/volxb(i)
     &          + (suryb(j)*TYY(I,J,K)-suryb(j-1)*TYY(I,J-1,K))/volya(j)
     &                    /geoxg(i)/geozg(k)
     &          + (surza(k+1)*TYZ(I,J,K+1)-surza(k)*TYZ(I,J,K))/volzb(k) 
     &                    /geoxh(i)
                  ENDDO
               ENDDO
            ENDDO
            DO K=1,locNZ+1
               DO J=0,locNY+1
                  DO I=0,locNX+1
                     DIVTZ(I,J,K)=
     &          (surzb(k-1)*DIVV(I,J,K-1)-surzb(k)*DIVV(I,J,K))/volza(k)
     &                    /geoxh(i)
     &          + (surxa(i+1)*TXZ(I+1,J,K)-surxa(i)*TXZ(I,J,K))/volxb(i)
     &          + (surya(j+1)*TYZ(I,J+1,K)-surya(j)*TYZ(I,J,K))/volyb(j)
     &                    /geoxg(i)/geozga(k)
     &          + (surzb(k)*TZZ(I,J,K)-surzb(k-1)*TZZ(I,J,K-1))/volza(k)
     &                    /geoxh(i)
                  ENDDO
               ENDDO
            ENDDO
 
c--------------------------
c CYLINDRICAL COORDINATES
c--------------------------
         else if (ncosys.eq.1) then
            DO K=0,locNZ+1
               DO J=0,locNY+1
                  DO I=1,locNX+1
                     DIVTX(I,J,K)=
     &          (surxb(i-1)*DIVV(I-1,J,K)-surxb(i)*DIVV(I,J,K))/volxa(i)
     &          + (surxb(i)*TXX(I,J,K)-surxb(i-1)*TXX(I-1,J,K))/volxa(i)
     &          + (surya(j+1)*TXY(I,J+1,K)-surya(j)*TXY(I,J,K))/volyb(j)
     &                    /geoxga(i)/geozg(k)
     &          + (surza(k+1)*TXZ(I,J,K+1)-surza(k)*TXZ(I,J,K))/volzb(k)
     &                    /geoxha(i)
     &          -((TYY(I,J,K)+TYY(I-1,J,K))-(DIVV(I,J,K)+DIVV(I-1,J,K)))
     &                    /2.d0/XXA(I)
                  ENDDO
               ENDDO
            ENDDO
            DO K=0,locNZ+1
               DO J=1,locNY+1
                  DO I=0,locNX+1
                     DIVTY(I,J,K)=
     &          (suryb(j-1)*DIVV(I,J-1,K)-suryb(j)*DIVV(I,J,K))/volya(j)
     &                    /geoxg(i)/geozg(k)
     &          + (surxa(i+1)*TXY(I+1,J,K)-surxa(i)*TXY(I,J,K))/volxb(i)
     &          + (suryb(j)*TYY(I,J,K)-suryb(j-1)*TYY(I,J-1,K))/volya(j)
     &                    /geoxg(i)/geozg(k)
     &          + (surza(k+1)*TYZ(I,J,K+1)-surza(k)*TYZ(I,J,K))/volzb(k) 
     &                    /geoxh(i)
     &                    + (TXYAV(I,J,K)+TXYAV(I,J-1,K))/2.d0/XXB(I)
                  ENDDO
               ENDDO
            ENDDO
            DO K=1,locNZ+1
               DO J=0,locNY+1
                  DO I=0,locNX+1
                     DIVTZ(I,J,K)=
     &          (surzb(k-1)*DIVV(I,J,K-1)-surzb(k)*DIVV(I,J,K))/volza(k)
     &                    /geoxh(i)
     &          + (surxa(i+1)*TXZ(I+1,J,K)-surxa(i)*TXZ(I,J,K))/volxb(i)
     &          + (surya(j+1)*TYZ(I,J+1,K)-surya(j)*TYZ(I,J,K))/volyb(j)
     &                    /geoxg(i)/geozga(k)
     &          + (surzb(k)*TZZ(I,J,K)-surzb(k-1)*TZZ(I,J,K-1))/volza(k)
     &                    /geoxh(i)
                  ENDDO
               ENDDO
            ENDDO
            
c--------------------------
c SPHERICAL COORDINATES
c--------------------------
         else if (ncosys.eq.2) then
            DO K=0,locNZ+1
               DO J=0,locNY+1
                  DO I=1,locNX+1
                     DIVTX(I,J,K)=
     &          (surxb(i-1)*DIVV(I-1,J,K)-surxb(i)*DIVV(I,J,K))/volxa(i)
     &          + (surxb(i)*TXX(I,J,K)-surxb(i-1)*TXX(I-1,J,K))/volxa(i)
     &          + (surya(j+1)*TXY(I,J+1,K)-surya(j)*TXY(I,J,K))/volyb(j)
     &                    /geoxga(i)/geozg(k)
     &          + (surza(k+1)*TXZ(I,J,K+1)-surza(k)*TXZ(I,J,K))/volzb(k)
     &                    /geoxha(i)- (TYY(I,J,K)+TYY(I-1,J,K) +
     &                    TZZ(I,J,K)+TZZ(I-1,J,K) -
     &                    (DIVV(I,J,K)+DIVV(I-1,J,K))*2.d0)/2.d0/XXA(I)
                  ENDDO
               ENDDO
            ENDDO
            DO K=0,locNZ+1
               DO J=1,locNY+1
                  DO I=0,locNX+1
                     DIVTY(I,J,K)=
     &          (suryb(j-1)*DIVV(I,J-1,K)-suryb(j)*DIVV(I,J,K))/volya(j)
     &                    /geoxg(i)/geozg(k)
     &          + (surxa(i+1)*TXY(I+1,J,K)-surxa(i)*TXY(I,J,K))/volxb(i)
     &          + (suryb(j)*TYY(I,J,K)-suryb(j-1)*TYY(I,J-1,K))/volya(j)
     &                    /geoxg(i)/geozg(k)
     &          + (surza(k+1)*TYZ(I,J,K+1)-surza(k)*TYZ(I,J,K))/volzb(k) 
     &                    /geoxh(i)
     &                 + ((TYZAV(I,J,K)+TYZAV(I,J-1,K))*(-TAN(XZB(K))) + 
     &                    (TXYAV(I,J,K)+TXYAV(I,J-1,K)))/2.d0/XXB(I)
                  ENDDO
               ENDDO
            ENDDO
            DO K=1,locNZ+1
               DO J=0,locNY+1
                  DO I=0,locNX+1
                     DIVTZ(I,J,K)=
     &          (surzb(k-1)*DIVV(I,J,K-1)-surzb(k)*DIVV(I,J,K))/volza(k)
     &                    /geoxh(i)
     &          + (surxa(i+1)*TXZ(I+1,J,K)-surxa(i)*TXZ(I,J,K))/volxb(i)
     &          + (surya(j+1)*TYZ(I,J+1,K)-surya(j)*TYZ(I,J,K))/volyb(j)
     &                    /geoxg(i)/geozga(k)
     &          + (surzb(k)*TZZ(I,J,K)-surzb(k-1)*TZZ(I,J,K-1))/volza(k)
     &                    /geoxh(i)
     &                    + ((TYY(I,J,K)+TYY(I,J,K-1)-
     &                    (DIVV(I,J,K)+DIVV(I,J,K-1)))*TAN(XZA(K))
     &                    +TXZAV(I,J,K)+TXZAV(I,J,K-1))/2.d0/XXB(I)

!     At the northern/southern most grid in the pole simulations, force
!     divtz=0. This is a fudge... I'm not sure how to correctly write
!     this term
                     if(poles.and.proc_coords(3).eq.proc_dims(3)-1.and.
     $                    k.eq.locNZ) then
                        divtz(i,j,k)= 0.d0
                     endif
                     if(poles.and.proc_coords(3).eq.0.and.k.eq.1) then
                        divtz(i,j,k)= 0.d0
                     endif                     
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

!     For the convecting box save the viscous flux. Because this is a
!     flux for velocity it is defined at grid centers. Recall that
!     T_{i,j} is symmetric, so I use TXY (instead of TYX), TXZ (instead
!     of TZX), and TYZ (instead of TZY)
         if(MODTYP.eq.1) then
            DO K=1,locNZ
               DO J=1,locNY+1
                  DO I=1,locNX
                     visc_flux(1,i,j,k)=
     $                    V(i,j,k)*(TXX(i,j,k)-divv(i,j,k)) + 
     $                    G(i,j,k)*TXY(i,j,k) + 
     $                    H(i,j,k)*TXZ(i,j,k)
                     visc_flux(2,i,j,k)=
     $                    V(i,j,k)*TXY(i,j,k) + 
     $                    G(i,j,k)*(TYY(i,j,k)-DIVV(i,j,k)) + 
     $                    H(i,j,k)*TYZ(i,j,k)
                     visc_flux(3,i,j,k)=
     $                    V(i,j,k)*TXZ(i,j,k) + 
     $                    G(i,j,k)*TYZ(i,j,k) + 
     $                    H(i,j,k)*(TZZ(i,j,k)-DIVV(i,j,k)) 
                  enddo
               enddo
            enddo
         endif
      CASE (2:)
         print *,'unknown IVISTYP=',IVISTYP
         stop
      CASE DEFAULT
         print *,'also not defined IVISTYP',IVISTYP
      END SELECT
      RETURN
      END SUBROUTINE VISC

      SUBROUTINE ETACALC
!-- eta is the dynamic viscosity, related to the kinematic 
!--  viscosity (mu) via eta = mu*rh
! ieta = 0:  constant kinematic viscosity
! ieta = 1:  alpha viscosity, eta=rh*(alpha*cs^2/om) = rh*(alpha*cs*H)
! ieta = 2:  ramp up kinematic viscosity in time
! ieta = 3:  increase kinematic viscosity at the bottom boundary
! ieta = 4:  use a molecular viscosity. visc=cs/(n*sigma)
! ieta = 5:  constant kinematic viscosity, maximum Pr_num (for difrx)
!----------
      USE input_init
      USE force_var_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      integer :: I,J,K,STEPS,ivisc_bndry
      double precision :: om,XNUE_RAMP,XNUE_VAL,XNUE_START
      double precision :: TIME_XNU_START,TIME_XNU_END,mslope
      double precision :: xnue_extra
**  Calculation of the DYNAMIC Viscosity
      SELECT CASE(IETA)
      CASE (0) !-- constant kinematic viscosity
         ETA=RH*XNUE
      CASE (1) !--  Calculation of the ALPHA Viscosity
         DO J=-1,locNY+2
            DO I=-1,locNX+2
               om = (GRAV*MSOL/(xxb(i)**3.d0))**0.5
               ETA(I,J,K)=alpha*RH(I,J,K)*(CS(I,J,K)**2.d0)/OM
            ENDDO
         ENDDO
      CASE (2) !-- Ramp of kinematic visc from XNUE_START to XNUE
c         STEPS = 100000
c         xnue_ramp = (1.d0*III/(1.d0*STEPS))
c         if(xnue_ramp.ge.1.0) xnue_ramp = 1.d0
c-linear ramp
c         XNUE_START = 10.d0**(12.d0)
c         xnue_val = XNUE_START-(XNUE_START-XNUE)*xnue_ramp
c-linear in log ramp over STEPS
c         XNUE_START = 10.d0
c         xnue_val=10.d0**(XNUE_START-(XNUE_START-LOG10(XNUE))*xnue_ramp)
c         IF(((III-1)/(IREA/10))*(IREA/10).eq.(III-1).and.myid.eq.0) THEN
c            write(*,'(A,3(1x,e12.6),1x,I8,1x,I8)') 
c     $           'viscosity ramp %=',xnue_ramp,xnue_val,xnue_start,
c     $           steps,III
c         ENDIF

!-ramp that is linear in log, but uses a fixed time rather then step number
         XNUE_START = 9.d0 !(defined in log-space)
         TIME_XNU_START = 282.5*DAY
         TIME_XNU_END   = TIME_XNU_START + 10.d0*DAY
         if(hydrotime.lt.time_xnu_start) then
            xnue_val = 10.d0**(XNUE_START)
         elseif (hydrotime.ge.time_xnu_start.and.
     %           hydrotime.lt.time_xnu_end) then
            mslope = (LOG10(xnue)-XNUE_START)/
     $           (TIME_XNU_END-TIME_XNU_START)
            xnue_val=10.d0**(mslope*(hydrotime-TIME_XNU_START)+
     $           XNUE_START)
         elseif(hydrotime.ge.TIME_XNU_END) then
            xnue_val = XNUE
         endif
         IF(((III-1)/(IREA/10))*(IREA/10).eq.(III-1).and.myid.eq.0) THEN
            write(*,'(A,4(1x,e12.6))') 
     $           'viscosity ramp =',hydrotime/DAY,xnue_val,
     $           xnue_start,xnue
         ENDIF
         ETA=RH*XNUE_VAL
      CASE (3) ! ramp up kinematic viscosity at the bottom boundary
         ETA=RH*XNUE
         xnue_extra = 4.d0
         ivisc_bndry= 15
         mslope = xnue_extra/(1.d0-1.d0*ivisc_bndry)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               DO I=1,ivisc_bndry
                  ETA(i,j,k)=RH(i,j,k)*
     $                 10.d0**(mslope*(i-ivisc_bndry)+LOG10(XNUE))
               enddo
               DO I=-1,0
                  ETA(i,j,k)=RH(i,j,k)*10.d0**(xnue_extra+LOG10(XNUE))
               enddo
            enddo
         enddo
         IF(III.eq.1.and.myid.eq.0) then
            write(*,'(A,2(1x,e12.6))') 
     $           'viscosity profile (ieta=3):',xnue,xnue_extra
            DO I=-1,locNX+2
               write(*,'(A,1x,I8,1x,e12.6)') 
     $              'log10(\nu)=',i,log10(eta(i,1,1)/rh(i,1,1))
            ENDDO
         ENDIF
      CASE (4) !-molecular viscosity only = cs*(1/(n\sigma) where sigma=1Ang^2
         CALL SOUND
c         ETA = (10.d0**5.d0)*CS*mu_gas*(1.67*10.d0**(-8.d0))/RH
         ETA = CS*mu_gas*(1.67*10.d0**(-8.d0))/RH
         ETA = min(ETA,XNUE)
         ETA = ETA*RH
      CASE (5) !-- constant kinematic viscosity, and maximum limit on Pr
               !-- (important for difrx/y/z calculation)
         ETA=RH*XNUE
      CASE (6:)
         print *,'unknown IETA=',IETA
         stop
      CASE DEFAULT
         print *,'also not defined IETA',IETA
      END SELECT
      RETURN
      END SUBROUTINE ETACALC



      SUBROUTINE RAYLEIGHDRAG
      USE input_init
      USE fluid_var_init
      USE newton_heat_init
      IMPLICIT NONE
      INTEGER :: I,J,K
      double precision :: kv(locNX,locNY,locNZ),CALCULATE_HELD_KV
      if(modtyp.eq.2.and.modver.eq.6) then
         kV = CALCULATE_HELD_KV
      else
         print *,'RAYLEIGHDRAG IS NOT AVAILABLE FOR THIS MODTYP/VER'
         call clean_stop
      endif
      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX
               V(i,j,k)=V(i,j,k) - DELT*kV(i,j,k)*V(i,j,k)
               G(i,j,k)=V(i,j,k) - DELT*kV(i,j,k)*G(i,j,k)
               H(i,j,k)=V(i,j,k) - DELT*kV(i,j,k)*H(i,j,k)
            enddo
         enddo
      enddo

      RETURN
      END SUBROUTINE RAYLEIGHDRAG
