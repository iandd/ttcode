!-
!-- general set of routines to run the Held&Suraez 1994 benchmark calculation
!-

      SUBROUTINE CALCULATE_HELD_TEQ
      USE newton_heat_init
      USE fluid_var_init
      USE grid_var_init
      IMPLICIT NONE
      INTEGER :: I,J,K
      double precision :: deltaT,deltaTheta,P0,kappa,tmp
      deltaT = 60.
      deltaTheta = 10.0
      P0 = 1e6 !=1000mb
      kappa = 2./7.
      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX
               tmp = 315.0 - (deltaT*((sin(xzb(k)))**2)) - 
     #              (deltaTheta*log10(PG(i,j,k)/P0)*((cos(xzb(k)))**2))
               tmp = tmp*((PG(i,j,k)/P0)**kappa)
               Teq_newton(i,j,k) = max(200.0,tmp)
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE CALCULATE_HELD_TEQ


      SUBROUTINE CALCULATE_HELD_KT
      USE fluid_var_init
      USE newton_heat_init
      USE global_constants
      USE grid_var_init
      IMPLICIT NONE
      INTEGER :: I,J,K
      double precision :: kA,kS,sigmaB,P0,sigma,tmp
      kA = (1./40.)/DAY
      kS = 0.25/DAY
      sigmaB = 0.7
      P0 = 1e6 !=1000mb  
      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX
               sigma = PG(i,j,k)/P0
               tmp = (sigma-sigmaB)/(1.0-sigmaB)
               kT_newton(i,j,k) = kA + (kS - kA)*max(0.d0,tmp)*
     $              ((cos(xzb(k)))**4.)
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE CALCULATE_HELD_KT


      FUNCTION CALCULATE_HELD_KV
      USE fluid_var_init
      USE newton_heat_init
      USE global_constants
      USE grid_var_init
      IMPLICIT NONE
      INTEGER :: I,J,K
      double precision :: calculate_held_kv(locNX,locNY,locNZ)
      double precision :: kF,sigmaB,P0,sigma,tmp
      kF = 1.0/DAY
      sigmaB = 0.7
      P0 = 1e6 !=1000mb  
      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX
               sigma = PG(i,j,k)/P0
               tmp = (sigma-sigmaB)/(1.0-sigmaB)
               calculate_held_kv(i,j,k) = kF*max(0.d0,tmp)
            enddo
         enddo
      enddo
      RETURN
      END
