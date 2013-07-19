
      SUBROUTINE FORGAS2D
**  Calculates 2D Force terms (centrifugal, graviation, pressure, viscosity)
ci  Also calculates the Pdiv(u) contribution for the energy equation
      USE input_init
      USE force_var_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer i,j,k
      double precision :: TOLD(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: VH(locNX,locNY,locNZ)
      double precision :: GH(locNX,locNY,locNZ)
**  Boundary Conditions
      CALL BOUNDS
**  Densities at Cell Boundaries - Only need to call if advection is not on
      if(IADV.EQ.0) CALL QDENS
**  Equation of state
      CALL PRESSG
**  Viscosity
c      IF (IVIS.GT.0) THEN
c         CALL ETACALC
c         CALL VISC2D
c      ENDIF
**  update Gravity force (only necessary if adding mass, eccentric or
**    non-syncronized roche potential, or ramping up gravity)
**  The other cases are set during initialization
      if(IGRAV.eq.3.or.(IGRAV.eq.5.and.(f_incident.eq.6.or.
     $     f_incident.eq.7)).or.IGRAV.eq.6) CALL FGRAV
**  Centrifugal and Coriolis Terms
      CALL CENTRIFCORIOLIS

**  Old Temperature
      TOLD=T

**  Temperature at t(n+1/2), by p div(v) Work (and dissipative energy)
**  (also in the case of total Energy)
      IF (ISOTH.EQ.0) THEN
         DO K=0,locNZ+1
            DO J=0,locNY+1
               DO I=0,locNX+1
                  T(I,J,K)=T(I,J,K) - DELT/2.d0
     &                 * PG(I,J,K) / RH(I,J,K) / CV(I,J,K)
     &             * ((surxa(i+1)*V(I+1,J,K)-surxa(i)*V(I,J,K))/volxb(i)
     &             + (surya(j+1)*G(I,J+1,K)-surya(j)*G(I,J,K))/volyb(j)
     &             + (surza(k+1)*H(I,J,k+1)-surza(k)*H(I,J,k))/volzb(k))
               ENDDO
            ENDDO
 3       ENDDO
**  Pressure at t(n+1/2) (Approximation, i.e. don't call PRESSG)
         CALL BOUNDS
         CALL PRESSG
      ENDIF 

      
**  Change in Momentum
**  X-Direction
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               VH(I,J,K)=V(I,J,K) + DELT * (
     &              ( -(PG(I,J,K)-PG(I-1,J,K)) /DXB(I)
     &              + DIVTX(I,J,K) )/RHQX(I,J,K)
     &              + centcorx(i,j,K) + gravx(i,j,K) )
            ENDDO
         ENDDO
      ENDDO

**  Y-Direction
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               GH(I,J,K)=G(I,J,K) + DELT * ((
     &              ( -(PG(I,J,K)-PG(I,J-1,K))
     &              /DYB(J)/geoxg(i)/geozg(k)
     &              + DIVTY(I,J,K) )/RHQY(I,J,K)
     &              + centcory(i,j,K)+ gravy(i,j,K))
     &              /geoxg(i)/geozg(k))
            ENDDO
         ENDDO
      ENDDO

**  reSet T=T(old) 
      T=TOLD
**  Calculate new temperature at T(N+1) using Velocities at t(n+1/2)
**  Boundary Conditions (for velocities) (T is reset, RH is unchanged)
**   don't have to call BOURAD here.  only velocities have changed. 
      IF (ISOTH.eq.0) THEN
!-----  Velocities at t(n+1/2)
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  V(I,J,K)=(V(I,J,K)+VH(I,J,K))/2.d0
                  G(I,J,K)=(G(I,J,K)+GH(I,J,K))/2.
               ENDDO
            ENDDO
         ENDDO
         CALL BOUNDS         
         if (ietot.eq.1) then
            print *,'this was changed'
            stop
         elseif (ietot.eq.0) then
**  p div(v) Work (Plus the dissipation function)
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     T(I,J,K)=T(I,J,K) - delt
     &                    * PG(I,J,K) / RH(I,J,K) / CV(I,J,K) *
     &             ( (surxa(i+1)*V(I+1,J,K)-surxa(i)*V(I,J,K))/volxb(i)
     &             + (surya(j+1)*G(I,J+1,K)-surya(j)*G(I,J,K))/volyb(j)
     &             + (surza(k+1)*H(I,J,k+1)-surza(k)*H(I,J,k))/volzb(k))
                  ENDDO
               ENDDO
            ENDDO
            if (idiss.ne.0) then
               print *,'diss N/A: 2d force'
               stop
            endif    
         else if (ietot.ge.2) then
*     dissipation and the v dot delP term (opposite of P div v)
            stop 'ietot N/A in force5C'
         endif
      ENDIF

**  New Velocities at t(n+1)
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               V(i,J,K) = VH(i,J,K)
               G(i,J,K) = GH(i,J,K)
            ENDDO
         ENDDO
      ENDDO
      CALL BOUNDS
      RETURN
      END SUBROUTINE FORGAS2D


      SUBROUTINE FORGAS3D
**  Calculates 3D Force terms (centrifugal, graviation, pressure, viscosity)
ci  Also calculates the Pdiv(u) contribution for the energy equation
      USE input_init
      USE force_var_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer i,j,k
      integer :: vdamp_steps
      double precision :: TOLD(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: VH(locNX,locNY,locNZ)
      double precision :: GH(locNX,locNY,locNZ)
      double precision :: HH(locNX,locNY,locNZ)
      double precision :: vdamp_ramp,depth_vrdamp
      double precision :: test_T(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: test_V(0:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: test_G(-1:locNX+2,0:locNY+2,-1:locNZ+2)
      double precision :: test_H(-1:locNX+2,-1:locNY+2,0:locNZ+2)
**  Boundary Conditions
      CALL BOUNDS
**  Densities at Cell Boundaries - Only need to call if advection is not on
      if(IADV.EQ.0) CALL QDENS
**  Equation of state
      CALL PRESSG
**  Viscosity
      IF (IVIS.GT.0) THEN
         CALL ETACALC
         CALL VISC
      ENDIF
**  update Gravity force (only necessary if adding mass, eccentric or
**    non-syncronized roche potential, or ramping up gravity)
**  The other cases are set during initialization
      if(IGRAV.eq.3.or.(IGRAV.eq.5.and.(f_incident.eq.6.or.
     $     f_incident.eq.7)).or.IGRAV.eq.6) CALL FGRAV
**  Centrifugal and Coriolis Terms + 2 advection terms
      CALL CENTRIFCORIOLIS

**  Old Temperature
      TOLD=T

**  Temperature at t(n+1/2), by p div(v) Work (and dissipative energy)
**  (also in the case of total Energy)

      IF (ISOTH.EQ.0) THEN
         DO K=0,locNZ+1
            DO J=0,locNY+1
               DO I=0,locNX+1
                  T(I,J,K)=T(I,J,K) - DELT/2.d0
     &                * PG(I,J,K) / RH(I,J,K) / CV(I,J,K)
     &             * ((surxa(i+1)*V(I+1,J,K)-surxa(i)*V(I,J,K))/volxb(i)
     &             + (surya(j+1)*G(I,J+1,K)-surya(j)*G(I,J,K))/volyb(j)
     &             + (surza(k+1)*H(I,J,k+1)-surza(k)*H(I,J,k))/volzb(k))
               ENDDO
            ENDDO
         ENDDO
         if (idiss.ne.0.and.ivis.ne.0) then
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     T(I,J,K)=T(I,J,K) + DELT/2.d0
     &                    * DISFN(I,J,K) / CV(I,J,K) / RH(I,J,K)
                  ENDDO
               ENDDO
            ENDDO
         endif
**  Pressure at t(n+1/2) (Approximation, i.e. don't call PRESSG)
         CALL BOUNDS
         CALL PRESSG
      ENDIF 

**  Change in Momentum
**  X-Direction
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               VH(I,J,K)=V(I,J,K) + DELT * (
     &              ( -(PG(I,J,K)-PG(I-1,J,K)) /DXB(I)
     &              + DIVTX(I,J,K) )/RHQX(I,J,K)
     &              + centcorx(i,j,K) + gravx(i,j,K) )

c               if(myid.eq.19.and.III.le.10.and.j.eq.5.and.k.eq.3) then
c               IF ((IREA/10).ne.0 .and. 
c     &              ((III-1)/(IREA/10))*(IREA/10).eq.(III-1)) THEN
c                  if(myid.eq.19.and.j.eq.5.and.k.eq.3) then
c                  if(myid.eq.19.and.j.eq.5.and.k.eq.3.and.
c     $              (I.gt.0.and.I.lt.7)) then
c                     write(*,'(A,I8,I8,6(1x,e11.5))') 'VH:',
c     $                    iii,i,VH(I,J,K),V(I,J,K),
c     &                    -(PG(I,J,K)-PG(I-1,J,K))/DXB(I)/RHQX(I,J,K),
cc     &                    PG(I,J,K),PG(I-1,J,K),
c     $                    RHQX(I,J,K),gravx(i,j,K),
c     &                    (-(PG(I,J,K)-PG(I-1,J,K)) /DXB(I)
c     &                    )/RHQX(I,J,K)+ gravx(i,j,K)
c                  endif
               

            ENDDO
         ENDDO
      ENDDO

**  Y-Direction
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               GH(I,J,K)=G(I,J,K) + DELT * ((
     &              ( -(PG(I,J,K)-PG(I,J-1,K))
     &              /DYB(J)/geoxg(i)/geozg(k)
     &              + DIVTY(I,J,K) )/RHQY(I,J,K)
     &              + centcory(i,j,K)+ gravy(i,j,K))
     &              /geoxg(i)/geozg(k))
c               if(myid.eq.19.and.III.le.10.and.i.eq.5.and.k.eq.1) then
c                  write(*,'(A,I8,I8,8(1x,e15.9))') 'GH:',
c     $                 iii,j,GH(I,J,K),G(I,J,K),
c     &                 -(PG(I,J,K)-PG(I,J-1,K))
c     &                 /DYB(J)/geoxg(i)/geozg(k)/RHQY(I,J,K),
c     &                 DIVTY(I,J,K)/RHQY(I,J,K),RHQY(I,J,K),PG(i,j,k),
c     &                 RH(I,J,K),T(i,j,k)
c               endif
            ENDDO
         ENDDO
      ENDDO

**  Z-Direction
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               HH(I,J,K)=H(I,J,K) + DELT * ((
     &              (-(PG(I,J,k)-PG(I,J,K-1)) /DZB(K)/geoxh(i)
     &              + DIVTZ(I,J,K) )/RHQZ(I,J,k)
     &              + centcorz(i,j,k) + gravz(i,j,k) )/geoxh(i))
            ENDDO
         ENDDO
      ENDDO


!--- DAMP VELOCITIES
      IF(VELDAMP.GT.0.d0) 
     $     CALL VELOCITYDAMPING(III,hydrotime,locNX,locNY,locNZ,
     $     VH,GH,HH)


**  reSet T=T(old) 
      T=TOLD
**  Calculate new temperature at T(N+1) using Velocities at t(n+1/2)
**  Boundary Conditions (for velocities) (T is reset, RH is unchanged)
**   don't have to call BOURAD here.  only velocities have changed. 
      if (isoth.eq.0) then
!-----  Velocities at t(n+1/2)
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  V(I,J,K)=(V(I,J,K)+VH(I,J,K))/2.d0
                  G(I,J,K)=(G(I,J,K)+GH(I,J,K))/2.d0
                  H(I,J,K)=(H(I,J,K)+HH(I,J,K))/2.d0
               ENDDO
            ENDDO
         ENDDO

!---- LIMIT ON VELOCITIES FOR INTERMEDIATE STEP
         if (abs(IVEL).gt.0.d0) then
            if(ncosys.eq.0) then
               CALL VARLIMIT (IVEL,V,-VELMAX,VELMAX)
               CALL VARLIMIT (IVEL,G,-VELMAX,VELMAX)
               CALL VARLIMIT (IVEL,H,-VELMAX,VELMAX)
            elseif(ncosys.eq.2) then
               CALL SPHERICALVELLIMIT
            else
               print *,'ncosys not defined:pp.f'
               stop
            endif 
         endif
         CALL BOUNDS         
         if (ietot.eq.1) then
            print *,'this was changed'
            stop
         elseif (ietot.eq.0) then
**  p div(v) Work (Plus the dissipation function)
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     T(I,J,K)=T(I,J,K) - delt
     &              * PG(I,J,K) /RH(I,J,K) / CV(I,J,K) *
     &             ( (surxa(i+1)*V(I+1,J,K)-surxa(i)*V(I,J,K))/volxb(i)
     &             + (surya(j+1)*G(I,J+1,K)-surya(j)*G(I,J,K))/volyb(j)
     &             + (surza(k+1)*H(I,J,k+1)-surza(k)*H(I,J,k))/volzb(k))
                   ENDDO
               ENDDO
            ENDDO
            if (idiss.ne.0.and.ivis.ne.0) then
               DO K=1,locNZ
                  DO J=1,locNY
                     DO I=1,locNX
                        T(I,J,K)=T(I,J,K) + DELT
     &                       * DISFN(I,J,K) / CV(I,J,K) / RH(I,J,K)
                     enddo
                  enddo
               enddo
            endif
         else if (ietot.ge.2) then
*     dissipation and the v dot delP term (opposite of P div v)
            stop 'ietot N/A in force'
         endif
      endif
**  New Velocities at t(n+1)
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               V(i,J,K) = VH(i,J,K)
               G(i,J,K) = GH(i,J,K)
               H(i,J,K) = HH(i,J,K)
            ENDDO
         ENDDO
      ENDDO

!---- LIMIT ON VELOCITIES FOR INTERMEDIATE STEP
      if (abs(IVEL).gt.0.d0) then
         if(ncosys.eq.0) then
            CALL VARLIMIT (IVEL,V,-VELMAX,VELMAX)
            CALL VARLIMIT (IVEL,G,-VELMAX,VELMAX)
            CALL VARLIMIT (IVEL,H,-VELMAX,VELMAX)
         elseif(ncosys.eq.2) then
            CALL SPHERICALVELLIMIT
         else
            print *,'ncosys not defined:pp.f'
            stop
         endif 
      endif

!---- Artifical Viscosity
      if(IARTVIS.gt.0) then
         CALL BOUNDS
         SELECT CASE(IARTVIS)
         CASE (1)
c            test_T = T
c            test_V = V
c            test_G = G
c            test_H = H
            CALL NR_ARTIFICIAL_VISC
c            ARTVIS_T = T - test_T
c            ARTVIS_V = V - test_V
c            ARTVIS_G = G - test_G
c            ARTVIS_H = H - test_H
         CASE (2:)
            print *,'unknown IARTVIS=',IARTVIS
            stop
         CASE DEFAULT
            print *,'also not defined IARTVIS=',IARTVIS
         END SELECT
      ENDIF
      CALL BOUNDS
      RETURN
      END SUBROUTINE FORGAS3D


      SUBROUTINE CENTRIFCORIOLIS
**  Centrifugal and Coriolis Forces and also one advection term for x and z-directions
**  The y-direction Coriolis Force is included in ADVECT
      USE force_var_init
      USE fluid_var_init
      USE grid_var_init
      USE input_init
      USE mpi_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: GQ,H1,H2,HQ
      double precision :: G1,G2,GQkep2
      double precision :: o2r2
      if (ncosys.eq.0) then
         centcorx = 0.d0
         centcory = 0.d0
         centcorz = 0.d0
      elseif(ncosys.eq.1) then
**  Cylindrical (r-z)
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
** (V_phi^2)/r + OMROT^2*r + 2OMROT*V_phi (Centrifugal + corriolis)
** the Centrifugal + corriolis are included with the OMROT in G1 and G2
                  G1=DCY0(J+1)*G(I-1,J,K)+DCY1(J+1)*G(I-1,J+1,K)+omrot
                  G2=DCY0(J+1)*G(I,J,K)  +DCY1(J+1)*G(I,J+1,K)+omrot
                  GQkep2=DDX0kep(I)*G1 + DDX1kep(I)*G2
                  centcorx(i,j,k) = (GQkep2/XXA(I))**2.d0
               ENDDO
            ENDDO
         ENDDO
         centcory = 0.d0
         centcorz = 0.d0
      elseif (ncosys.eq.2) then
         if(igrav.eq.5.or.igrav.eq.6) then
            call ROCHE_CENTRIFCORIOLIS
         else
!--  Centrifugal + Corriolis + 2 advection terms
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
*     *  X-Direction: GQkep2 is G+omega centered where V is defined.. with an
*     *  extra r**1.5. HQ is H centered where V is defined
                     G1=DCY0(J+1)*G(I-1,J,K)+DCY1(J+1)*G(I-1,J+1,K)+
     $                    omrot
                     G2=DCY0(J+1)*G(I,J,K)  +DCY1(J+1)*G(I,J+1,K)+omrot
                     GQkep2=DDX0kep(I)*G1 + DDX1kep(I)*G2
                     H1=DDZ0(K+1)*H(I-1,J,K)+DDZ1(K+1)*H(I-1,J,K+1)
                     H2=DDZ0(K+1)*H(I  ,J,K)+DDZ1(K+1)*H(I  ,J,K+1)
                     HQ=DDX0(I)*H1 + DDX1(I)*H2
                     centcorx(i,j,k) = (GQkep2/XXA(I)*geozg(k))**2 
     1                    + XXA(I)*HQ**2
c                     centcorx(i,j,k) = 0.d0
*     *  Y-Direction
                     centcory(i,j,k)= 0.d0
*     *  Z-Direction
                     G1=DCY0(J+1)*G(I,J,K-1)+DCY1(J+1)*G(I,J+1,K-1)
                     G2=DCY0(J+1)*G(I,J,K)  +DCY1(J+1)*G(I,J+1,K)  
                     GQ=DDZ1(K)*G2 +  DDZ0(K)*G1 + omrot
                     centcorz(i,j,k)= - XXB(I) * GQ*GQ 
     &                    *SIN(XZA(K))*COS(XZA(K))
c                     centcorz(i,j,k)= 0.d0
!     At the northern/southern most grid, force centcorz=0. If this was
!     really +/-(pi/2) then the centrifugal force=0, and G=0, so both the
!     advection term and the corriolis term would also be zero
                     if(poles.and.proc_coords(3).eq.proc_dims(3)-1.and.
     $                    k.eq.locNZ) then
                        centcorz(i,j,k)= 0.d0
                     endif
                     if(poles.and.proc_coords(3).eq.0.and.k.eq.1) then
                        centcorz(i,j,k)= 0.d0
                     endif                     

!-try replacing the v=0 at the poles with this only
c                     if(poles.and.proc_coords(3).eq.0.and.k.eq.1) then
c                        centcorx(i,j,k) = 0.d0
c                     endif         
c                     if(poles.and.proc_coords(3).eq.proc_dims(3)-1.and.
c     $                    (k.eq.locNZ.or.k.eq.locNZ-1)) then
c                        centcorx(i,j,k) = 0.d0
c                     endif
                  enddo
               enddo
            enddo
         endif
      else
         print *,'CENTRIF N/A for ncosys=',ncosys
         stop
      endif
      RETURN
      END SUBROUTINE CENTRIFCORIOLIS


      SUBROUTINE ROCHE_CENTRIFCORIOLIS
** Routine to replace normal CENTRIFCORIOLIS function for roche_gravity
!   case. As with the other, this routine contains both Centrifugal and
!   Coriolis Forces and also one advection term for x and z-directions
!   The y-direction Coriolis Force is still included in ADVECT. Note
!   that the form of the Coriolis Force is unchanged with a simple
!   translation of omega from the center of the planet to the center of
!   mass
!----
      USE force_var_init
      USE fluid_var_init
      USE grid_var_init
      USE input_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,j,k
      double precision :: rcom,beta,qratio,u_sq
      double precision :: uphi_bar,utheta_bar
      double precision :: phi_bar,theta_bar
      double precision :: drcom_dr,dbeta_dr,drcom_dphi,dbeta_dphi
      double precision :: drcom_dtheta,dbeta_dtheta
      double precision :: CALC_ROCHE_RCOM,CALC_ROCHE_BETA
      qratio = mstar2/(mstar1_ramp*mstar1) !--mplanet/mstar
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
!-- THE R-TERM ---------- ---------- ----------
               rcom = CALC_ROCHE_RCOM(xxa(i),xyb(j),xzb(k))
               beta = CALC_ROCHE_BETA(xxa(i),xzb(k),rcom)
!-- phi-velcocity at (ia,jb,kb)
               uphi_bar = (0.25*cos(xzb(k)))*( xxb(i-1)*(G(i-1,j,k)+
     $              G(i-1,j+1,k)) + xxb(i)*(G(i,j,k)+G(i,j+1,k)))
!-- theta-velcocity at (ia,jb,kb)
               utheta_bar = 0.25*( xxb(i-1)*(H(i-1,j,k)+H(i-1,j,k+1)) + 
     $              xxb(i)*(H(i,j,k)+H(i,j,k+1)) )
!-- derivative: d(R_com)/dr
               drcom_dr = ( xxa(i) - radius_orb*(cos(xyb(j)))*
     $              (cos(xzb(k)))/(1.d0+qratio) )/rcom
!-- derivative: d(beta)/dr
               u_sq = (xxa(i)*(sin(xzb(k)))/rcom)**2
               dbeta_dr = -( (sin(xzb(k)))/((1.d0-u_sq)**(0.5)))*
     $              ( (1.d0/rcom) - (xxa(i)*drcom_dr/(rcom*rcom)) )
!-- put it all together:
!--------
!-   = (uphi^2/r) + (utheta^2/r) + 2*omega*uphi*cos(theta) + centrifugal
!      -advection-   -advection-      -coriolis-            -centrifugal-
!--------
               centcorx(i,j,k) = ((uphi_bar**2)/xxa(i)) +
     $              ((utheta_bar**2)/xxa(i)) + 
     $              2*OMROT*uphi_bar*cos(xzb(k)) + 
     $              (OMROT**2)*( rcom*drcom_dr*((sin(beta))**2) + 
     $              (rcom**2)*(sin(beta))*(cos(beta))*dbeta_dr )


               test_roche(1,i,j,k) = gravx(i,j,k) +
     $              (OMROT**2)*( rcom*drcom_dr*((sin(beta))**2) + 
     $              (rcom**2)*(sin(beta))*(cos(beta))*dbeta_dr )


!--THE PHI-TERM (note coriolis portion is in advection routine) ----------
               rcom = CALC_ROCHE_RCOM(xxb(i),xya(j),xzb(k))
               beta = CALC_ROCHE_BETA(xxb(i),xzb(k),rcom)
!-- derivative: d(R_com)/dphi
               drcom_dphi = radius_orb*xxb(i)*(sin(xya(j)))*
     $              (cos(xzb(k)))/(rcom*(1.d0+qratio))
!-- derivative: d(beta)/dphi
               u_sq = (xxb(i)*(sin(xzb(k)))/rcom)**2
               dbeta_dphi = (1.d0/((1.d0-u_sq)**(0.5)))*
     $              (xxb(i)*(sin(xzb(k)))/(rcom**2))*drcom_dphi
!-- put it all together:
!--------
! = just centrifugal
!--------
               centcory(i,j,k) = (OMROT**2)*( rcom*drcom_dphi*
     $              ((sin(beta))**2) + 
     $              (rcom**2)*(sin(beta))*(cos(beta))*dbeta_dphi )/
     $              (xxb(i)*cos(xzb(k)))
               test_roche(2,i,j,k) = gravy(i,j,k) +
     $              (OMROT**2)*( rcom*drcom_dphi*
     $              ((sin(beta))**2) + 
     $              (rcom**2)*(sin(beta))*(cos(beta))*dbeta_dphi )/
     $              (xxb(i)*cos(xzb(k)))
               
!--THE THETA-TERM ---------- ---------- ----------
               rcom = CALC_ROCHE_RCOM(xxb(i),xyb(j),xza(k))
               beta = CALC_ROCHE_BETA(xxb(i),xza(k),rcom)
!-- phi-velcocity at (ib,jb,ka)
               uphi_bar = (0.25*xxb(i))*( (cos(xzb(k-1)))*(G(i,j,k-1)+
     $              G(i,j+1,k-1)) + (cos(xzb(k)))*(G(i,j,k)+G(i,j+1,k)))
!-- derivative: d(R_com)/dtheta
               drcom_dtheta = radius_orb*xxb(i)*(cos(xyb(j)))*
     $              (sin(xza(k)))/(rcom*(1.d0+qratio))
!-- derivative: d(beta)/dtheta
               u_sq = (xxb(i)*(sin(xza(k)))/rcom)**2
               dbeta_dtheta = -(xxb(i)/((1.d0-u_sq)**(0.5)))*
     $              ( ((cos(xza(k)))/rcom) - 
     $              ((sin(xza(k)))*drcom_dtheta/(rcom*rcom)) )
!-- put it all together:
!--------
! = -(uphi^2/r)tan(theta) - 2*omega*uphi*sin(theta) - centrifugal
!      -advection-              -coriolis-           -centrifugal-
!--------
               centcorz(i,j,k) = -((uphi_bar**2)*(tan(xza(k)))/xxb(i))
     $              -(2.d0*OMROT*uphi_bar*sin(xza(k))) + 
     $              ((OMROT**2)*( rcom*drcom_dtheta*((sin(beta))**2) +
     $              (rcom**2)*(sin(beta))*(cos(beta))*
     $              dbeta_dtheta )/xxb(i))

               test_roche(3,i,j,k) = gravz(i,j,k) +
     $              ((OMROT**2)*( rcom*drcom_dtheta*((sin(beta))**2) +
     $              (rcom**2)*(sin(beta))*(cos(beta))*
     $              dbeta_dtheta )/xxb(i))

            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE ROCHE_CENTRIFCORIOLIS


      FUNCTION CALC_ROCHE_RCOM(r_internal,phi_internal,theta_internal)
!-distance to the center of mass
      USE fluid_var_init
      USE force_var_init
      USE input_init
      IMPLICIT NONE
      double precision :: CALC_ROCHE_RCOM,qratio
      double precision :: r_internal,phi_internal,theta_internal
      qratio = mstar2/(mstar1_ramp*mstar1) !--mplanet/mstar
      CALC_ROCHE_RCOM = (  ((radius_orb/(1.d0+qratio))**2) + 
     $     (r_internal**2) - 
     $     2.d0*r_internal*radius_orb*(cos(theta_internal))*
     $     (cos(phi_internal))/(1.d0+qratio)  )**(0.5)
      RETURN
      END FUNCTION CALC_ROCHE_RCOM

      FUNCTION CALC_ROCHE_BETA(r_internal,theta_internal,rcom_internal)
!-angle between r_com and omega
      USE global_constants
      IMPLICIT NONE
      double precision :: CALC_ROCHE_BETA
      double precision :: r_internal,theta_internal
      double precision :: rcom_internal
      CALC_ROCHE_BETA = (PI/2.d0) - 
     $     asin(r_internal*(sin(theta_internal))/rcom_internal)
      RETURN
      END FUNCTION CALC_ROCHE_BETA
      


      SUBROUTINE FGRAV
**  Gravitational Force
!
!     IGRAV = 0: All grav=0
!     IGRAV = 1: Constant Radial Gravity
!     IGRAV = 2: CENTRAL STAR (MSOL*MSTAR1) + PLANET (MSOL*MSTAR2)
!     IGRAV = 3: CENTRAL STAR (MSOL*MSTAR1) + PLANET (MSOL*MSTAR2) + MASS ACCRETED
!     IGRAV = 4: POINT MASS POTENTIAL (1/r^2) OF MPLANET = MSOL*MSTAR2
!     IGRAV = 5: ROCHE POTENTIAL (BOTH CIRCULAR AND ECCENTRIC ORBITS)
!     IGRAV = 6: RAMP UP M1 OF ROCHE POTENTIAL
!     IGRAV = 7: UNIFORM GRAVITY RAMP, SET PRARAMETERS HERE
!     IGRAV = 8: 
!
      USE force_var_init
      USE fluid_var_init
      USE input_init
      USE global_constants
      USE grid_var_init
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,j,k
      integer :: grav_STEPS
      double precision :: rplanet,phiplanet,Mplanet,d
      double precision :: r1,GM1,GM2
      double precision :: x_dist,z_dist,dist
      double precision :: init_gravx,final_gravx,grav_ramp
      double precision :: SETPLANETXGRAV
      double precision :: mstar1_init,m1ramp_time_start,m1ramp_time_end
      double precision :: mslope
!--- NO GRAVITY
      SELECT CASE(igrav)
      CASE (0)
         gravx=0.d0
         gravy=0.d0
         gravz=0.d0
!--- CONSTANT GRAVITY
      CASE (1)
         gravx = SETPLANETXGRAV(planet_num)
         gravy=0.d0
         gravz=0.d0
!--- 2 BODY GRAVITY: CENTRAL STAR (MSTAR1) + PLANET (MSTAR2)
!--- 3 BODY GRAVITY: CENTRAL STAR + PLANET + MASS ACCRETED
      CASE (2:3)
         rplanet = rp*AU
         phiplanet = phip*PI
         Mplanet = MSOL*MSTAR2
         if(locNZ.ne.1) then
            print *,'THIS GRAV IS ONLY FOR 2D'
            stop
         endif
**    Normal One-Mass Potential
         DO i=-1,locNX+2
            gravx(i,1,1)=-GRAV*MSOL*MSTAR1/(xxa(i)**2)
         ENDDO
**    Normal One-Mass Potential
         DO k=-1,locNZ+2
            DO j=-1,locNY+2
               DO i=-1,locNX+2
                  gravx(i,j,k)=gravx(i,1,1)
               ENDDO
            ENDDO
         ENDDO
         gravy = 0.d0
         gravz = 0.d0
** GRAVITATIONAL FORCE FROM SECOND BODY         
         if (MSTAR2.gt.0.0) then
            if(myid.eq.0.and.num_iter.eq.0) then
               write(*,'(A,e12.3)') 'a_planet (AU)=',rplanet/AU
               write(*,'(A,e12.3)') 'phi(planet) =',phiplanet
            endif
            CALL JUPGRAV2D(rplanet,phiplanet)
         endif
      CASE (4)
!--- SINGLE POINT MASS POTENTIAL OF  MPLANET = MSOL*MSTAR2
         MPLANET = MSOL*MSTAR2
         if(myid.eq.0) then
            print *,'igravity=4'
            write(*,'(A,e12.6)') '   mplanet=',mplanet
         endif
         if(ncosys.eq.1) then
            DO K=-1,locNZ+2
               DO J=-1,locNY+2
                  DO I=-1,locNX+2
!-----  define distance with xxa and xzb
                     if(i.eq.-1) then !--xxa(-1) is not defined
                        x_dist = xxa(0)-dxa(1)
                     else
                        x_dist = xxa(i)
                     endif
                     z_dist = xzb(k)
                     dist = ( (x_dist**(2.d0))+(z_dist**(2.d0)) )**(0.5)
                     gravx(i,1,k)=-(GRAV*MPLANET/(dist**2.d0))*
     %                    (x_dist/dist)
!-----  now define distance with xxb and xza
                     x_dist = xxb(i)
                     if(k.eq.-1) then !--xza(-1) is not defined
                        z_dist = xza(0)-dza(1)
                     else
                        z_dist = xza(k)
                     endif
                     dist = ( (x_dist**(2.d0))+(z_dist**(2.d0)) )**(0.5)
                     gravz(i,1,k)=-(GRAV*MPLANET/(dist**2.d0))*
     %                    (z_dist/dist)
                  ENDDO
               ENDDO
            ENDDO
            gravy = 0.d0
         elseif(ncosys.eq.2) then
            gravx(-1,1,1)=-GRAV*MPLANET/((xxa(0)-dxa(1))**2) !-xxa(-1) and dxa(0) not defined
            DO i=0,locNX+2
               gravx(i,1,1)=-GRAV*MPLANET/(xxa(i)**2)
            ENDDO
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,locNX+2
                     gravx(i,j,k)=gravx(i,1,1)
                  ENDDO
               ENDDO
            ENDDO
            gravy = 0.d0
            gravz = 0.d0
         else
            print *,'fgrav=4, ncosys undefined',ncosys
         endif
      CASE (5:6)
!--- ROCHE POTENTIAL (Centrifugal portion in CENTRIFCORIOLIS)
!-- note radius_orb and phi0_ecc are defined (and updated) elsewhere
         if(ncosys.eq.2) then
            GM2 = GRAV*MSOL*mstar2 !--planet
            if(igrav.eq.6) then ! ramp up mstar1 (from 0.0001 to mstar1)
! Initial m1 mass (in solar units). Cannot be 0 because of qratio. Also
!  should be large enough such that the center of mass of the two body
!  system is at a distance from the planetary center more then twice XMAX
               mstar1_init = 5.e-3
               m1ramp_time_start = 71*DAY !-start of ramp
               m1ramp_time_end = m1ramp_time_start + 25.*DAY !-end of ramp
               if(hydrotime.lt.m1ramp_time_start) then
                  mstar1_ramp = mstar1_init/mstar1
               elseif(hydrotime.ge.m1ramp_time_start.and.
     $                 hydrotime.lt.m1ramp_time_end) then
                  mslope = (LOG10(mstar1/mstar1_init))/
     $                 (m1ramp_time_end-m1ramp_time_start)
                  mstar1_ramp = (mstar1_init/mstar1)*(
     $                 10.d0**(mslope*(hydrotime-m1ramp_time_start)))
               elseif(hydrotime.ge.m1ramp_time_end) then
                  mstar1_ramp = 1.0
               endif
            else
               mstar1_ramp = 1.d0
            endif
            
            IF(((III-1)/(IREA/10))*(IREA/10).eq.(III-1).and.
     $           myid.eq.0) THEN
               write(*,'(A,4(1x,e12.6))') 
     $              ' Primary Mass Ramp =',hydrotime/DAY,mstar1_ramp,
     $              mstar1_init/mstar1,mstar1
            ENDIF
            

            GM1 = GRAV*MSOL*mstar1*mstar1_ramp !--star
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,locNX+2
!-- distance to primary (M1) using xxa,xyb,xzb and secondary (M2) using xxa
                     if(i.eq.-1) then !--xxa(-1) is not defined
                        x_dist = xxa(0)-dxa(1)
                     else
                        x_dist = xxa(i)
                     endif
                     r1 = ( (x_dist**2.d0)+(radius_orb**2.)-
     1                    (2.d0*radius_orb*x_dist*cos(xzb(k))*
     2                    cos(xyb(j)-phi0_ecc)))**(0.5)
                     gravx(i,j,k)=-(GM1/(r1**3.d0))*(x_dist-
     1                    radius_orb*cos(xzb(k))*cos(xyb(j)-phi0_ecc)) - 
     2                    (GM2/(x_dist**2.d0))
!-- distance to primary (M1) using xxb,xya,xzb
                     r1 = ( (xxb(i)**2.d0)+(radius_orb**2.)-
     1                    (2.d0*radius_orb*xxb(i)*cos(xzb(k))*
     2                    cos(xya(j)-phi0_ecc)) )**(0.5)
                     gravy(i,j,k)=-(GM1/(r1**3.d0))*radius_orb*
     1                    sin(xya(j)-phi0_ecc)
!-- distance to primary (M1) using xxb,xyb,xza
                     r1 = ( (xxb(i)**2.d0)+(radius_orb**2.)-
     1                    (2.d0*radius_orb*xxb(i)*cos(xza(k))*
     2                    cos(xyb(j)-phi0_ecc)) )**(0.5)
                     gravz(i,j,k)=-(GM1/(r1**3.d0))*(
     1                    radius_orb*sin(xza(k))*cos(xyb(j)-phi0_ecc) )
                  ENDDO
               ENDDO
            ENDDO
         else
            print *,'ROCHE POTENTIAL NOT AVALIBLE IN ncosys=',ncosys
            stop
         endif
      CASE (7) !--- Gradual Linear change in uniform radial gravity
         init_gravx = -910.d0
         final_gravx = -2209.44 !-hatp7
         grav_STEPS = 100000
         grav_ramp = (1.d0*III/(1.d0*grav_STEPS))
         if(grav_ramp.ge.1.0) grav_ramp = 1.d0
c-linear in log10 space
         gravx=-10.d0**(LOG10(-init_gravx)-(LOG10(-init_gravx)-
     $        LOG10(-final_gravx))*grav_ramp)
c         gravx = init_gravx - (init_gravx-final_gravx)*grav_ramp
         IF(((III-1)/(IREA/10))*(IREA/10).eq.(III-1).and.
     $        myid.eq.0) THEN
            write(*,'(A,2(1x,e12.6),1x,I10)') 
     $           'gravitational ramp %=',grav_ramp,gravx(1,1,1),III
         ENDIF        
         gravy=0.d0
         gravz=0.d0
      CASE (8) !--sudo roche potential
         print *,'unknown IGRAV=',IGRAV
         stop
      CASE (9:)
         print *,'unknown IGRAV=',IGRAV
         stop
      CASE DEFAULT
         print *,'also not defined IGRAV ',IGRAV
      END SELECT
      RETURN
      END SUBROUTINE FGRAV


      SUBROUTINE JUPGRAV2D(rplanet,phiplanet)
**  Calculate the force on the gas due to a protoplanet/star at r=rplanet,phi=phiplanet
      USE input_init
      USE force_var_init
      USE global_constants
      USE mpi_var_init
      USE fluid_var_init
      USE grid_var_init
      IMPLICIT NONE
      double precision :: rplanet,phiplanet,rsoft
      integer :: i,j
      double precision :: GM,Acosb,Acosa,Asina,Rtwo,slowfac,t0
      GM = GRAV*MSOL
      rsoft = psoft*rplanet*(mstar2/3.d0)**(1.d0/3.d0)
      if(myid.eq.0.and.num_iter.eq.0) write(*,'(A,e12.3)')
     %     'Gravity softening param (RH)=',psoft
c      if(myid.eq.0) print *,'deltax=',XXA(2)-XXA(1)
      if(myid.eq.0.and.num_iter.eq.0) write(*,'(A,e12.3)')
     %     ' rhill=',rplanet*(mstar2/3.d0)**(1.d0/3.d0)
      if(myid.eq.0.and.num_iter.eq.0) write(*,'(A,2(e12.3))')
     %     ' rhill/dx, rhill/dy=',
     %     rplanet*(mstar2/3.d0)**(1.d0/3.d0)/(XXA(2)-XXA(1)),
     %     rplanet*(mstar2/3.d0)**(1.d0/3.d0)/(rplanet*(XYA(2)-XYA(1)))
      if(myid.eq.0.and.num_iter.eq.0) 
     %     write(*,'(A,2(e12.3))') ' rsoft/dx, rsoft/dy, =',
     %     rsoft/(XXA(2)-XXA(1)), rsoft/(rplanet*(XYA(2)-XYA(1)))
* 3/14/03 Pawel ramp up
c      if (iartvis.eq.5) then
c      slowfac=0.9
c      t0=5. * 2.*pi/omrot
c      if (zeit.gt.t0) then
c         slowfac=1.
c      else
c         slowfac=sin(zeit/t0*pi/2.)**2
c      endif
c      else
      slowfac=1.d0
c      endif         
      DO J=1,locNY
         Acosb=rplanet*cos(XYB(J)-phiplanet)
         Acosa=rplanet*cos(XYA(J)-phiplanet)
         Asina=rplanet*sin(XYA(J)-phiplanet)
         DO I=1,locNX
            Rtwo=dsqrt(rplanet*rplanet + XXA(I)*XXA(I)
     &           -2.*XXA(I)*Acosb
     &           +rsoft*rsoft)
            
            gravx(I,J,1)=gravx(I,J,1) + slowfac* (
     &           - GM*MSTAR2 /Rtwo**3
     &           *(XXA(I) - Acosb) )

            Rtwo=sqrt(rplanet*rplanet + XXB(I)*XXB(I)
     &           -2.*XXB(I)*Acosa
     &           +rsoft*rsoft)

            IF(IGRAV.EQ.2) THEN !-- PLANET MASS ONLY
               gravy(I,J,1)=gravy(I,J,1)+ slowfac* ( 
     &              - GM*MSTAR2/Rtwo**3
     &              *Asina )
            ELSEIF(IGRAV.EQ.3) THEN ! INCLUDE ACCRETED MASS
               gravy(I,J,1)=gravy(I,J,1)+ slowfac* ( 
     &              - (GM*MSTAR2+GRAV*Maccreated)/Rtwo**3
     &              *Asina )
            ELSE
               print *,'ill defined igrav:jupgrav2d',igrav
               stop
            ENDIF
c            if(myid.eq.0)print *,gravy(i,j,1)
         ENDDO
      ENDDO
c      if (xMode.ne.13.and.abs(xMode).ne.14) then
c*** add additional inertial term if not in center-of-mass frame (13) 
c*** (or for 14, special case)
c         Aconst=slowfac*GM*MSTAR2/(rplanet*rplanet)
c         DO J=1,NY
c            Acos=-Aconst* cos(XYB(J)-phiplanet)
c            Asin= Aconst* sin(XYA(J)-phiplanet)
c            DO I=1,NX
c               gravx(I,J,1)=gravx(I,J,1) + Acos
c               gravy(I,J,1)=gravy(I,J,1) + Asin
c            ENDDO
c         ENDDO 
c      endif
      RETURN
      END SUBROUTINE JUPGRAV2D


      SUBROUTINE JUPGRAV3D(rplanet,phiplanet)
**  Calculate the force on the gas due to a protoplanet/star at r=rplanet,phi=phiplanet
      USE input_init
      USE force_var_init
      USE global_constants
      USE mpi_var_init
      USE fluid_var_init
      USE grid_var_init
      IMPLICIT NONE
      double precision :: rplanet,phiplanet,rsoft
      integer :: i,j
      double precision :: GM,Acosb,Acosa,Asina,Rtwo,slowfac,t0
      GM = GRAV*MSOL
      rsoft = psoft*rplanet*(mstar2/3.d0)**(1.d0/3.d0)
      if(myid.eq.0.and.num_iter.eq.0) write(*,'(A,e12.3)')
     %     'Gravity softening param (RH)=',psoft
      if(myid.eq.0.and.num_iter.eq.0) write(*,'(A,e12.3)')
     %     ' rhill=',rplanet*(mstar2/3.d0)**(1.d0/3.d0)
      if(myid.eq.0.and.num_iter.eq.0) write(*,'(A,2(e12.3))')
     %     ' rhill/dx, rhill/dy=',
     %     rplanet*(mstar2/3.d0)**(1.d0/3.d0)/(XXA(2)-XXA(1)),
     %     rplanet*(mstar2/3.d0)**(1.d0/3.d0)/(rplanet*(XYA(2)-XYA(1)))
      if(myid.eq.0.and.num_iter.eq.0) 
     %     write(*,'(A,2(e12.3))') ' rsoft/dx, rsoft/dy, =',
     %     rsoft/(XXA(2)-XXA(1)), rsoft/(rplanet*(XYA(2)-XYA(1)))
      slowfac=1.d0
      DO J=1,locNY
         Acosb=rplanet*cos(XYB(J)-phiplanet)
         Acosa=rplanet*cos(XYA(J)-phiplanet)
         Asina=rplanet*sin(XYA(J)-phiplanet)
         DO I=1,locNX
            Rtwo=dsqrt(rplanet*rplanet + XXA(I)*XXA(I)
     &           -2.*XXA(I)*Acosb
     &           +rsoft*rsoft)
            
            gravx(I,J,1)=gravx(I,J,1) + slowfac* (
     &           - GM*MSTAR2 /Rtwo**3
     &           *(XXA(I) - Acosb) )

            Rtwo=sqrt(rplanet*rplanet + XXB(I)*XXB(I)
     &           -2.*XXB(I)*Acosa
     &           +rsoft*rsoft)

            IF(IGRAV.EQ.2) THEN !-- PLANET MASS ONLY
               gravy(I,J,1)=gravy(I,J,1)+ slowfac* ( 
     &              - GM*MSTAR2/Rtwo**3
     &              *Asina )
            ELSEIF(IGRAV.EQ.3) THEN ! INCLUDE ACCRETED MASS
               gravy(I,J,1)=gravy(I,J,1)+ slowfac* ( 
     &              - (GM*MSTAR2+GRAV*Maccreated)/Rtwo**3
     &              *Asina )
            ELSE
               print *,'ill defined igrav:jupgrav2d',igrav
               stop
            ENDIF
c            if(myid.eq.0)print *,gravy(i,j,1)
         ENDDO
      ENDDO
c      if (xMode.ne.13.and.abs(xMode).ne.14) then
c*** add additional inertial term if not in center-of-mass frame (13) 
c*** (or for 14, special case)
c         Aconst=slowfac*GM*MSTAR2/(rplanet*rplanet)
c         DO J=1,NY
c            Acos=-Aconst* cos(XYB(J)-phiplanet)
c            Asin= Aconst* sin(XYA(J)-phiplanet)
c            DO I=1,NX
c               gravx(I,J,1)=gravx(I,J,1) + Acos
c               gravy(I,J,1)=gravy(I,J,1) + Asin
c            ENDDO
c         ENDDO 
c      endif
      RETURN
      END SUBROUTINE JUPGRAV3D


      SUBROUTINE NR_ARTIFICIAL_VISC
!-- Standard Neumann-Richtmyer method
!--  use formalisim as in ZEUS code. 
!--  see ZUES3D viscus.src
      USE input_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      USE artvis_var_init
      USE mpi_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      integer i,j,k,zmin_index,zmax_index
      double precision :: diff_vx,diff_vy,diff_vz
      
      integer :: i_av,j_av,k_av,dir_av
      double precision :: dt_av_tstX,dt_av_tstY,dt_av_tstZ,dt_av_all
      dt_artvis = 0.d0

      dt_av_tstX= 0.d0
      dt_av_tstY= 0.d0
      dt_av_tstZ= 0.d0
      dt_av_all = 0.d0


!-set z-range for artifical viscosity (important for dt calculation)
      zmin_index = 1
      zmax_index = locNZ
      if(poles.and.MPIabove.eq.MPI_PROC_NULL) then
         zmin_index = 1
         zmax_index = locNZ-1
      endif
      if(poles.and.MPIbelow.eq.MPI_PROC_NULL) then
         zmin_index = 2
         zmax_index = locNZ
      endif
c**  X-Direction
      DO K=zmin_index,zmax_index
         DO J=1,locNY
            diff_vx = (V(1,J,K)-V(0,J,K))
            if(diff_vx.gt.0.d0) diff_vx = 0.d0
            QviscX(0,j,k)= artvisc_coeff*delt*RH(0,J,K)*
     %           (diff_vx**(2.d0))
            dt_artvis = min(dt_artvis,(diff_vx/dxb(0)))
!----------------------------------------------------------------------
            if((diff_vx/dxb(0)).lt.dt_av_all) then
               i_av = 0
               j_av = j
               k_av = k
               dir_av = 4
               dt_av_all = (diff_vx/dxb(0))
            endif
            dt_av_tstX = min(dt_av_tstX,(diff_vx/dxb(0)))
!----------------------------------------------------------------------
            DO I=1,upperbnd(j,k)-2
               diff_vx = (V(I+1,J,K)-V(I,J,K))
               if(diff_vx.gt.0.d0) diff_vx = 0.d0
               QviscX(i,j,k)= artvisc_coeff*delt*RH(I,J,K)*
     $              (diff_vx**(2.d0))
               V(I,J,K) = V(I,J,K) - (QviscX(i,j,k)-QviscX(i-1,j,k))/
     %              (DXB(I)*RHQX(I,J,k))
               T(i,j,k) = T(i,j,k) - QviscX(i,j,k)*diff_vx/
     %              (DXA(i)*RH(I,J,K)*CV(I,J,K))
               dt_artvis = min(dt_artvis,(diff_vx/dxb(i)))
!----------------------------------------------------------------------
               if((diff_vx/dxb(i)).lt.dt_av_tstX) then
                  i_av = i
                  j_av = j
                  k_av = k
                  dir_av = 1
                  dt_av_all = (diff_vx/dxb(i))
               endif
               dt_av_tstX = min(dt_av_tstX,(diff_vx/dxb(i)))
!----------------------------------------------------------------------
            ENDDO
         ENDDO
      ENDDO         
**   Y-Direction
      DO K=zmin_index,zmax_index
         DO J=1,locNY
            DO I=1,upperbnd(j,k)-2
               if(j.eq.1) then  !-calculate Q(j=0) and timestep constraint
                  diff_vy = (G(I,1,K)-G(i,0,K))
                  if(diff_vy.gt.0.d0) diff_vy = 0.d0
                  QviscY(i,0,k)= artvisc_coeff*delt*RH(i,0,K)*
     %                 (diff_vy**(2.d0))
                  dt_artvis = min(dt_artvis,(diff_vy/dyb(0)))

!----------------------------------------------------------------------
                  if((diff_vy/dyb(0)).lt.dt_av_all) then
                     i_av = i
                     j_av = 0
                     k_av = k
                     dir_av = 5
                     dt_av_all = (diff_vy/dyb(0))
                  endif
                  dt_av_tstY = min(dt_av_tstY,(diff_vy/dyb(0)))
!----------------------------------------------------------------------

               endif
               diff_vy = (G(I,j+1,K)-G(i,j,K))
               if(diff_vy.gt.0.d0) diff_vy = 0.d0
               QviscY(i,j,k)= artvisc_coeff*delt*RH(I,J,K)*
     %              (diff_vy**(2.d0))
               G(I,J,K) = G(I,J,K) - (QviscY(i,j,k)-QviscY(i,j-1,k))/
     %              (DYB(j)*RHQY(I,J,k))
               T(i,j,k) = T(i,j,k) - QviscY(i,j,k)*diff_vy/
     %              (DYA(j)*RH(I,J,K)*CV(I,J,K))
               dt_artvis = min(dt_artvis,(diff_vy/dyb(j)))

!----------------------------------------------------------------------
               if((diff_vy/dyb(j)).lt.dt_av_all) then
                  i_av = i
                  j_av = j
                  k_av = k
                  dir_av = 2
                  dt_av_all = (diff_vy/dyb(j))
               endif
               dt_av_tstY = min(dt_av_tstY,(diff_vy/dyb(j)))
!----------------------------------------------------------------------

            ENDDO
         ENDDO
      ENDDO

      if(poles.and.MPIabove.eq.MPI_PROC_NULL) then
         zmin_index = 1
         zmax_index = locNZ-2
      endif
      if(poles.and.MPIbelow.eq.MPI_PROC_NULL) then
         zmin_index = 3
         zmax_index = locNZ
      endif
**   Z-Direction
      DO K=zmin_index,zmax_index
         DO J=1,locNY
            DO I=1,upperbnd(j,k)-2
               if(k.eq.zmin_index) then  !-calculate Q(k=zmin_index-1) and timestep constraint
                  diff_vz = (H(i,j,zmin_index)-H(i,j,zmin_index-1))
                  if(diff_vz.gt.0.d0) diff_vz = 0.d0
                  QviscZ(i,j,zmin_index-1)= artvisc_coeff*delt*
     %                 RH(i,j,zmin_index-1)*(diff_vz**(2.d0))
                  dt_artvis = min(dt_artvis,(diff_vz/dzb(zmin_index-1)))
!----------------------------------------------------------------------
                  if((diff_vz/dzb(zmin_index-1)).lt.dt_av_all) then
                     i_av = i
                     j_av = j
                     k_av = zmin_index-1
                     dir_av = 6
                     dt_av_all = (diff_vz/dzb(zmin_index-1))
                  endif
                  dt_av_tstZ = min(dt_av_tstZ,
     $                 (diff_vz/dzb(zmin_index-1)))
!----------------------------------------------------------------------
               endif             
               diff_vz = (H(i,j,k+1)-H(i,j,k))
               if(diff_vz.gt.0.d0) diff_vz = 0.d0
               QviscZ(i,j,k)= artvisc_coeff*delt*RH(I,J,K)*
     %              (diff_vz**(2.d0))
               H(I,J,K) = H(I,J,K) - (QviscZ(i,j,k)-QviscZ(i,j,k-1))/
     %              (DZB(k)*RHQZ(I,J,k))
               T(i,j,k) = T(i,j,k) - QviscZ(i,j,k)*diff_vz/
     %              (DZA(k)*RH(I,J,K)*CV(I,J,K))
               dt_artvis = min(dt_artvis,(diff_vz/dzb(k)))
!----------------------------------------------------------------------
               if((diff_vz/dzb(k)).lt.dt_av_all) then
                  i_av = i
                  j_av = j
                  k_av = k
                  dir_av = 3
                  dt_av_all = (diff_vz/dzb(k))
               endif
               dt_av_tstZ = min(dt_av_tstZ,(diff_vz/dzb(k)))
!----------------------------------------------------------------------
            ENDDO
         ENDDO
      ENDDO
!--   final artifical viscosity timestep calculation
      dt_artvis = min(1.d0/(abs(8.d0*artvisc_coeff*(dt_artvis-TINY))),
     %     DTMAX)



      dt_av_tstX = min(1.d0/(abs(8.d0*artvisc_coeff*(dt_av_tstX-TINY))),
     %     DTMAX)
      dt_av_tstY = min(1.d0/(abs(8.d0*artvisc_coeff*(dt_av_tstY-TINY))),
     %     DTMAX)
      dt_av_tstZ = min(1.d0/(abs(8.d0*artvisc_coeff*(dt_av_tstZ-TINY))),
     %     DTMAX)
      dt_av_all = min(1.d0/(abs(8.d0*artvisc_coeff*(dt_av_all-TINY))),
     %     DTMAX)

cc      write(*,'(A,I6,1x,I6,4(1x,e12.3),5(1x,I6))') 'XYZ:',iii,myid,
c     %     dt_av_tstX,dt_av_tstY,dt_av_tstZ,dt_av_all,i_av,j_av,k_av,
c     $     upperbnd(j_av,k_av),dir_av

      RETURN
      END SUBROUTINE NR_ARTIFICIAL_VISC
