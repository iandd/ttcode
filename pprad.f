      SUBROUTINE RADDIF
**     **Implicit** Radiation Transport (SOR)
!       irad=1: standard SOR calculation of div(F)
      USE rad_var_init
      USE fluid_var_init
      USE global_constants
      USE input_init
      IMPLICIT NONE
      integer :: i,j,K
!----- Convert to radiation energy density
      ER=arad*(T**(4.d0))
      if(igastyp.gt.2) then
         print *,'I may need to recalculate CV here:pprad'
         call clean_stop
      endif
**  Matrixelements for Diffusion (A,C only)
      CALL DIFFEL
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               B(I,J,K)=1.d0/delt*CV(I,J,K)*RH(I,J,K)
     &              /(4.d0*arad*T(I,J,K)*T(I,J,K)*T(I,J,K))
     &              -AX(I,J,K)-CX(I,J,K)
     &              -AY(I,J,K)-CY(I,J,K)
     &              -AZ(I,J,K)-CZ(I,J,K)
               RHS(I,J,K)=ER(I,J,K)/delt*CV(I,J,K)*RH(I,J,K)
     &              /(4.d0*arad*T(I,J,K)*T(I,J,K)*T(I,J,K))
            ENDDO
         ENDDO
      ENDDO
**  Boundary Conditions imposed on Matrixelements (also set T,ER in ttbou.f)
      CALL DIFFELBC
**    Solution of Matrixequation (in Array ER(I,J,K) )
      CALL SOLVMAT
!---Convert back into temperature
      T = (ER/arad)**(0.25)
      RETURN
      END SUBROUTINE RADDIF


      SUBROUTINE DIFFEL
**    Sets Matrixelements for Diffusion of Radiation
      USE input_init
      USE rad_var_init
      USE fluid_var_init
      USE grid_var_init
      IMPLICIT NONE
      INTEGER :: I,J,K
**  Calculate diffusion coefficients 
      CALL DIFF
*   Matrixelements (Diffusion of Radiation), 
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               AX(I,J,K) = -difrx(I,J,K)*surxa(i)/volxb(i)
               CX(I,J,K) = -difrx(I+1,J,K)*surxa(i+1)/volxb(i)
               AY(I,J,K) = -difry(I,J,K)*surya(j)/volyb(j)/geoxg(i)
     %              /geozg(k)
               CY(I,J,K) = -difry(I,J+1,K)*surya(j+1)/volyb(j)/geoxg(i)
     %              /geozg(k)
               AZ(I,J,K) = -difrz(I,J,K)*surza(k)/volzb(k)/geoxh(i)
               CZ(I,J,K) = -difrz(I,J,K+1)*surza(k+1)/volzb(k)/geoxh(i)
            ENDDO
         ENDDO
      ENDDO
**  No radiative transfer in z-direction for 2D simulations
      if(locNZ.eq.1) then
         AZ = 0.d0
         CZ = 0.d0
      endif
      RETURN
      END SUBROUTINE DIFFEL

      SUBROUTINE DIFF
**    calculate radiative diffusion coefficients
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      INTEGER :: I,J,k
      double precision :: difr_min
      if(myid.eq.0.and.ieta.eq.5.and.iii.eq.0) then
         write(*,'(A)') 
     $        'Setting a floor for radiative diffusion coefficent'
         write(*,'(A,1(1x,e12.6))') 
     $        ' Utilizing Pr_number=',pr_num
      endif
**  Flux-Limiter (calculates flimx,flimy)
      CALL FLIM
** difrx = (1/DX_b)*lambdax*c/(rho*(kappa+sig)) - defined on iA
** X
      if(irad.eq.5.or.irad.eq.6) then !-only use FLD for horizontal component
         difrx=0.d0
      else
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX+1
                  difrx(I,J,K)= fd*flimx(I,J,K)
     &                 / ( DXB(I)*RHQX(I,J,K)*
     &                 (DDX1(I)*(xkapR(I,J,K)+sig(I,J,K))
     &                 +DDX0(I)*(xkapR(i-1,j,K)+sig(I-1,J,K))) )
                  if(ieta.eq.5) then
                     difr_min = RHQX(I,J,K)*cv(i,j,k)*XNUE
     &                    / ( DXB(I)*4.d0*ARAD*
     &                    ((DDX1(I)*T(I,J,K)+DDX0(I)*T(i-1,j,K))**3.d0)
     &                    *pr_num)
                     difrx(i,j,k) = max(difrx(i,j,k),difr_min)
                  endif
               ENDDO
            ENDDO
         ENDDO
      endif
** Y
** difry = (1/DY*r*cos(theta))*lambday*c/(rho*(kappa+sig)) - defined on jA
      DO K=1,locNZ
         DO J=1,locNY+1
            DO I=1,locNX
               difry(I,J,K)= FD*FLIMY(I,J,K)/ 
     &              (DYB(J)*geoxg(i)*geozg(k)*RHQY(I,J,K)*
     &              (DDY1(J)*(xkapR(i,j,K)+sig(I,J,K))
     &              +DDY0(J)*(xkapR(I,J-1,K)+sig(I,J-1,K))) )
               if(ieta.eq.5) then
                  difr_min = RHQY(I,J,K)*cv(i,j,k)*XNUE
     &                 / ( DYB(J)*geoxg(i)*geozg(k)*4.d0*ARAD*
     &                 ((DDY1(j)*T(I,J,K)+DDY0(j)*T(i,j-1,K))**3.d0)
     &                 *pr_num)
                  difry(i,j,k) = max(difry(i,j,k),difr_min)
               endif
            ENDDO
         ENDDO
      ENDDO
** Z
** difrz = (1/DZ*r)*lambdaz*c/(rho*(kappa+sig)) - defined on kA
      DO K=1,locNZ+1
         DO J=1,locNY
            DO I=1,locNX
               difrz(I,J,K)= FD*FLIMZ(I,J,K)/ 
     &              ( DZB(K)*geoxh(i)*RHQZ(I,J,K)*
     &              (DDZ1(K)*(xkapR(i,j,k)+sig(I,J,K))
     &              +DDZ0(K)*(xkapR(I,J,K-1)+sig(I,J,K-1))) )
               if(ieta.eq.5) then
                  difr_min = RHQZ(I,J,K)*cv(i,j,k)*XNUE
     &                 / ( DZB(K)*geoxh(i)*4.d0*ARAD*
     &                 ((DDZ1(k)*T(I,J,K)+DDZ0(k)*T(i,j,K-1))**3.d0)
     &                 *pr_num)
                  difrz(i,j,k) = max(difrz(i,j,k),difr_min)
               endif
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DIFF
 
      SUBROUTINE DIFFELBC
**    calculate radiative boundary conditions
      IMPLICIT NONE
ci ibdrad/10=4  -- no calculations done in this routine->open compleatly;
!     used to be:open x; wrap around y and z
      RETURN
      END SUBROUTINE DIFFELBC

      SUBROUTINE NEWTONINAN_HEATING
      USE input_init
      USE fluid_var_init
      USE newton_heat_init
      IMPLICIT NONE
      INTEGER :: I,J,K
      if(modtyp.eq.2.and.modver.eq.6) then 
         CALL CALCULATE_HELD_TEQ
         CALL CALCULATE_HELD_kT
      else
         print *,'Newtonian Teq and kT N/A for this modtyp/ver'
         call clean_stop
      endif
      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX
               T(i,j,k)=T(i,j,k) - 
     %              DELT*(T(i,j,k)-Teq_newton(i,j,k))*kt_newton(i,j,k)
            enddo
         enddo
      enddo
      CALL BOUNDS
      RETURN
      END SUBROUTINE NEWTONINAN_HEATING

      SUBROUTINE CALCSTELLARTAU
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      INTEGER :: I,J,K
      double precision :: dr,totalstellartau
!-- calculate optical depths for both photons
      DR = XXB(2)-XXB(1)
      do k=-1,locNZ+2
         do j=-1,locNY+2
!---  Outside area should be uniformly 0
            if(upperbnd(j,k).lt.locNX) then
               STELLARTAU(upperbnd(j,k)+1:locNX+1,j,k) = 0.d0
            endif
            TOTALstellarTAU = 0.0
            do i=UPPERBND(j,k),1,-1
c--   stellar heating optical depth
               TOTALstellarTAU = TOTALstellarTAU+
     %              (xkapA(i,j,k)*RH(i,j,k)*DR)
               STELLARTAU(i,j,k) = TOTALstellarTAU
            enddo
         enddo
      enddo
!-calculate the tau_star on the A grid
      CALL STELLAR_DTAU
      RETURN
      END SUBROUTINE CALCSTELLARTAU

      SUBROUTINE STELLAR_DTAU
**  Calculates stellar tau at intermediate Gridpoints (A-Grid)
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      IMPLICIT NONE
      INTEGER :: I,J,K
      do k=-1,locNZ+2
         do j=-1,locNY+2
            do i=2,locNX
               dtaus(i,j,k)=(stellartau(i-1,j,k)-stellartau(i+1,j,k))
     $              /2.d0
            enddo
            dtaus(1,j,k) = dtaus(2,j,k) !this is wrong, but tau_star is irrelavent at this depth
            dtaus(locNX+1,j,k) = dtaus(locNX,j,k) !this is wrong, but tau_star=0 here
         enddo
      enddo
      RETURN
      END SUBROUTINE STELLAR_DTAU

**  EXPLICIT HEATING TERMS
      SUBROUTINE STELLARHEATING
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      INTEGER :: I,J,k
      double precision :: Tincident,phi0,W_hubeny
!--- find the opacities
      CALL KAPPA
!--- find the stellar optical depth
      CALL CALCSTELLARTAU
!--- calculate the heating for energy equation
      DO K=1,locNZ
         DO J=1,locNY
            STELLARINPUT(upperbnd(j,k):locNX,j,k) = 0.d0
            DO I=upperbnd(j,k)-1,1,-1
!--   irradiation temperature
!---  standard cos*cos dependance
               if(ncosys.eq.1) then
                  print *,'N/A...'
                  stop
               elseif(ncosys.eq.2) then
                  Tincident = Tirr*((cos(xyb(j))*
     %                 cos(xzb(k)))**0.25)
               endif
               if(ikaptyp.eq.20.or.ikaptyp.eq.21.or.
     $              ikaptyp.eq.22.or.ikaptyp.eq.23.or.
     $              ikaptyp.eq.24) then
                  print *,'the explicit heating in STELLARHEATING might'
                  print *,'not be correct'
                  call clean_stop
               endif
               Tincident=max(Tincident,0.d0)
!--   heating term
               STELLARINPUT(I,J,K) = RH(I,J,K)*
     &              xkapA(I,J,K)*SBCONST*(Tincident**4.d0)*
     &              exp(-STELLARTAU(I,J,K))
!--   update energy equation
               T(I,J,K)=T(I,J,K) + DELT *(
     &              STELLARINPUT(I,J,K)/CV(I,J,K)/RH(I,J,K))
            ENDDO
         ENDDO
      ENDDO
      CALL BOUNDS
      RETURN
      END SUBROUTINE STELLARHEATING

      FUNCTION CALC_TINCIDENT(jindx,kindx)
!-- note: set first_tincident outside do loops when calling this function...
!-- F_INCIDENT=1 : explicit stellar heating
!-- F_INCIDENT=2 : implicit stellar heating w/cos (**)
!-- F_INCIDENT=3 : implicit stellar heating w/cos ramped up over entire run
!-- F_INCIDENT=4 : sph. sym. heating
!-- F_INCIDENT=5 : sph. sym. heating ramped up over entire run
!-- F_INCIDENT=6 : eccentric orbit
!-- F_INCIDENT=7 : non-syncronized spin
!-- F_INCIDENT=8 : implicit stellar heating (=Tstar*sqrt(Rs/a) w/cos and 
!                   min T=Tirr (includes ramps of Tirr OR of Tstar*sqrt(Rs/a)
!----------------------------------------------------------------------
      USE input_init
      USE grid_var_init
      USE fluid_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      INTEGER :: jindx,kindx
      double precision :: calc_tincident,eccentric_anomaly
      double precision :: angle,Pforce
      double precision :: ramp_Tstar_subsol,Tmin_irr,Tmax_irr
      SELECT CASE(F_INCIDENT)
      CASE(0:1) !-- explicit stellar heating
         calc_Tincident= 0.d0
      CASE(2:3) !-- implicit stellar heating w/cos
         angle = (cos(xyb(jindx)))*(cos(xzb(kindx)))
         if(angle.ge.0) then
            calc_Tincident = Tirr*((angle)**0.25)
c            calc_Tincident = Tstar*sqrt(RSTAR*RSOL/semi)*
c     %           ((angle)**0.25)
            calc_Tincident = max(calc_Tincident,0.d0)
         else
            calc_Tincident = 0.d0
         endif
      CASE(4:5) !-- sph. sym. heating
         calc_Tincident = Tirr
      CASE(6:7) !-- eccentric orbit OR non-syncronized spin
         if(first_Tincident) then !--calculate the current radius and sub-stellar point
            Pforce = ((1.d0/Porb)-(1.d0/Prot))**(-1.d0)
            if(f_incident.eq.6.and.hydrotime.gt.ecc_start*DAY) then
               ecc_anomaly = ECCENTRIC_ANOMALY(hydrotime-ecc_start*DAY,
     $              ecc,Porb*DAY)
               radius_orb = semi*(1.d0-ecc*cos(ecc_anomaly))
            else !-non syncronous, but fixed semi-major axis
               radius_orb = semi
            endif
            if(hydrotime.lt.ecc_start*DAY) then
               phi0_ecc = 0.d0
            else
               phi0_ecc = (2.d0*PI/(Pforce*DAY))*
     %              mod((hydrotime-ecc_start*DAY),Pforce*day)
            endif
            first_Tincident = .false.
            if((iprstep.gt.0).and.
     %           ((num_iter/iprstep)*iprstep.eq.num_iter).and.
     %           (myid.eq.0)) then
               write(*,'(A,7(1x,e12.6))')
     &              ' r(AU),phi_0,Pf',radius_orb/AU,phi0_ecc,Pforce,
     $              hydrotime/DAY,ecc_start,
     $              mod((hydrotime-ecc_start*DAY),Pforce*day),
     $              Tstar*sqrt(RSTAR*RSOL/radius_orb)
            endif
         endif
         angle = (cos(xyb(jindx)-phi0_ecc))*(cos(xzb(kindx)))
         if(angle.ge.0) then
            calc_Tincident = Tstar*sqrt(RSTAR*RSOL/radius_orb)*
     %           ((angle)**0.25)
            calc_Tincident = max(calc_Tincident,0.d0)
         else
            calc_Tincident = 0.d0
         endif
      CASE(8) !-- implicit stellar heating w/cos, but minimum given by Tirr
c         Tmax_irr = Tstar*sqrt(RSTAR*RSOL/semi)
         Tmax_irr = Tirr
c         Tmin_irr = Tirr !-for the fixed Tmin_irr case
!-ramp down the minimum irradiation from Tirr to 0 (replaces const. Tmin_irr)
         Tmin_irr = Tirr - Tirr*(1.d0*III/(1.d0*DAMP_STEPS))
         Tmin_irr = max(Tmin_irr,0.d0)
         IF (((III-1)/(IREA/10))*(IREA/10).eq.(III-1)) THEN
            if(myid.eq.0) print *,'Tmin_irr_ramp=',Tmin_irr
         ENDIF
         angle = (cos(xyb(jindx)))*(cos(xzb(kindx)))
         if(angle.ge.0) then
!-ramp up the maximum irradiation from Tirr to Tstar*sqrt(Rs/a)
c            ramp_Tstar_subsol = (Tmax_irr-Tmin_irr)*(1.d0*III/
c     $           (1.d0*DAMP_STEPS)) + Tmin_irr
c            if(ramp_Tstar_subsol.ge.Tmax_irr)
c     $           ramp_Tstar_subsol = Tmax_irr
c            IF (first_Tincident.and.
c     $           ((III-1)/(IREA/10))*(IREA/10).eq.(III-1)) THEN
c               if(myid.eq.0) print *,'Tstar_ramp=',III,ramp_Tstar_subsol
c               first_Tincident = .false.
c            ENDIF
c            calc_Tincident = ramp_Tstar_subsol*((angle)**0.25)
c- use this for minimum ramp down... and comment out ramp_tstar_subsol calc
            calc_Tincident = Tmax_irr*((angle)**0.25)
            calc_Tincident = max(calc_Tincident,Tmin_irr)
         else
            calc_Tincident = Tmin_irr
         endif
      CASE(9:)
         print *,'F_incident not defined:calc_tincident',F_INCIDENT
         stop
      CASE DEFAULT
         print *,'F_incident also not defined:calc_tincident',F_INCIDENT
         stop
      END SELECT
      RETURN
      END FUNCTION CALC_TINCIDENT


      SUBROUTINE STELLARINPUT_HEATING
c!-- Routine to calculate the stellarinput term used in bounds and
c     energy.  This should be used for both the wavelength dependant
c     stellar heating rate and the grey models. This is choosen by
c     selecting ikaptyp
      USE input_init
      USE grid_var_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_var_init
      USE global_constants
      IMPLICIT NONE
      INTEGER :: i,j,k,lb
      double precision :: Tincident,CALC_TINCIDENT
      double precision :: irr_ramp,TIME_START,RAMP_TIME
      double precision :: SPHSYM_INPUT(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      if(f_incident.eq.0) then
         STELLARINPUT = 0.d0
      else
!-recalculate mu_star for eccentric models
         if((f_incident.eq.6.or.f_incident.eq.7).and.
     1        hydrotime.gt.ecc_start*DAY) then
            CALL CALCMUSTAR
         endif
!-calculate wavelength dependent stellar_tau and d(stellar_tau)
         if(ikaptyp.eq.20.or.ikaptyp.eq.21.or.
     $        ikaptyp.eq.22.or.ikaptyp.eq.23.or.ikaptyp.eq.24) then 
            CALL WAVELENGTH_CALCSTELLARTAU
            STELLARINPUT=0.d0
            do k=-1,locNZ+2
               do j=-1,locNY+2
                  do i=1,locNX+2
                     do lb=1,nwave_bins
                        STELLARINPUT(i,j,k)  = STELLARINPUT(i,j,k) +
     %                       BSTAR_BIN(lb)*
     &                       (exp(-stellartau_wave(I,J,K,lb)/
     &                       mu_star(j,k)))*dtaus_wave(i,j,k,lb)
                     enddo
                     STELLARINPUT(i,j,k)=STELLARINPUT(i,j,k)*
     $                    ((RSTAR*RSOL/semi)**(2.0))/dxb(i)
                  enddo
               enddo
            enddo
!-ramp from spherically symmetric to antisymmetric
            TIME_START = 5.d0*DAY
            RAMP_TIME = 10.d0*DAY
            if(f_incident.eq.8.and.
     $           hydrotime.le.(time_start+ramp_time)) then
               SPHSYM_INPUT=0.d0
               do k=-1,locNZ+2
                  do j=-1,locNY+2
                     do i=1,locNX+2
                        do lb=1,nwave_bins
                           SPHSYM_INPUT(i,j,k)  = SPHSYM_INPUT(i,j,k) +
     %                          BSTAR_BIN(lb)*
     &                          (exp(-stellartau_wave(I,J,K,lb)))
     &                          *dtaus_wave(i,j,k,lb)
                        enddo
                        
                        if(hydrotime.lt.time_start) then
                           STELLARINPUT(i,j,k) = SPHSYM_INPUT(i,j,k)*
     $                          ((RSTAR*RSOL/semi)**(2.0))/dxb(i)
                        else
                           SPHSYM_INPUT(i,j,k)=SPHSYM_INPUT(i,j,k)*
     $                          ((RSTAR*RSOL/semi)**(2.0))/dxb(i)
                           STELLARINPUT(i,j,k) = (hydrotime-time_start)*
     $                        ((STELLARINPUT(i,j,k)-SPHSYM_INPUT(i,j,k))
     $                          /RAMP_TIME) + SPHSYM_INPUT(i,j,k)
                        endif
                     enddo
                  enddo            
               enddo
               IF (((III-1)/(IREA/10))*(IREA/10).eq.(III-1)) THEN
                  if(myid.eq.numprocs/2) write(*,'(A,I8,5(1x,e12.6))')
     $                 ' Nightside cooling ramp:',myid,
     $                 TIME_START/DAY,RAMP_TIME/DAY,hydrotime/DAY,
     $                 SPHSYM_INPUT(upperbnd(locNY/2,locNZ/2),
     $                 locNY/2,locNZ/2),
     $                 STELLARINPUT(upperbnd(locNY/2,locNZ/2),
     $                 locNY/2,locNZ/2)
               ENDIF
            endif
         else
            CALL CALCSTELLARTAU
            do k=-1,locNZ+2
               do j=-1,locNY+2
                  Tincident = CALC_TINCIDENT(j,k)
                  do i=1,locNX+2
                     STELLARINPUT(i,j,k)  = SBCONST*
     $                    ((Tincident**4.d0)/mu_star(j,k))*
     &                    (exp(-stellartau(I,J,K)/
     &                    mu_star(j,k)))*dtaus(i,j,k)/dxb(i)
                  enddo
               enddo
            enddo   
         endif
!-perform ramp here
         if(F_INCIDENT.eq.3.or.F_INCIDENT.eq.5) then !--- RAMP-UP
            irr_ramp = (1.d0*III/(1.d0*DAMP_STEPS))
            if(irr_ramp.ge.1.0) irr_ramp = 1.d0
            STELLARINPUT(i,j,k)=STELLARINPUT(i,j,k)*irr_ramp
            IF(((III-1)/(IREA/10))*(IREA/10).eq.(III-1).and.
     $           myid.eq.0) THEN
               print *,'stellarinput %=',irr_ramp
            ENDIF
         endif
      endif
      RETURN
      END SUBROUTINE STELLARINPUT_HEATING


      SUBROUTINE WAVELENGTH_CALCSTELLARTAU
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      IMPLICIT NONE
      include "mpif.h"
      INTEGER :: I,J,K
      double precision :: dr,totalstellartau(nwave_bins)
!-- Zero optical depth everywhere initially
      stellartau_wave = 0.d0
!-- calculate optical depths for both photons
      DR = XXB(2)-XXB(1)
      do k=-1,locNZ+2
         do j=-1,locNY+2
!---  Outside area should be uniformly 0
            TOTALstellarTAU(:) = 0.0
            do i=UPPERBND(j,k),1,-1
c--   stellar heating optical depth
               TOTALstellarTAU(:)=TOTALstellarTAU(:)+
     %              (xkapW(i,j,k,:)*RH(i,j,k)*DR)
               stellartau_wave(I,J,K,:)=TOTALstellarTAU(:)
            enddo
         enddo
      enddo
!-calculate the dtau_star across the A grid
      CALL WAVELENGTH_STELLAR_DTAU
      RETURN
      END SUBROUTINE WAVELENGTH_CALCSTELLARTAU

      SUBROUTINE WAVELENGTH_STELLAR_DTAU
**  Calculates stellar tau at intermediate Gridpoints (A-Grid)
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      IMPLICIT NONE
      INTEGER :: I,J,K
      do k=-1,locNZ+2
         do j=-1,locNY+2
            do i=2,locNX
               dtaus_wave(i,j,k,:)=(stellartau_wave(i-1,j,k,:)-
     $              stellartau_wave(i+1,j,k,:))/2.d0
            enddo
!this is wrong, but tau_star is irrelavent at depth
            dtaus_wave(1,j,k,:) = dtaus_wave(2,j,k,:)
            dtaus_wave(locNX+1,j,k,:) = 0.d0
         enddo
      enddo
      RETURN
      END SUBROUTINE WAVELENGTH_STELLAR_DTAU
