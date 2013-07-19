      SUBROUTINE PRESSG
*   Calculates the Pressure, Cv, mu, etc. of the gas,
*    but not the sound speed
      USE input_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
!-for hs pressure
      USE force_var_init
      USE grid_var_init
      IMPLICIT NONE
      integer :: i,j,k
      SELECT CASE(igastyp)
      CASE(0)
**  ZERO PRESSURE
         PG = 0.d0
      CASE(1)
**  IDEAL GAS 
!-- note for varying mu, i need to change this
         PG=Rgas*T*RH/mu_gas
      CASE(2)
**  ISOTHERMAL (P= CS^2*RH)
         IF(ICSTYP.NE.1.AND.ICSTYP.NE.2) then
            CALL SOUND
         ENDIF
         PG = (CS**2.d0)*RH
      CASE(3)
**  POLYTROPIC 
!-- note that Tdisk0 and sigma0 are only defined for partial disk simulation
         K_poly = RGAS*Tdisk0*(sigma0**(1.d0-gamma))/mu_gas
         PG = K_poly*(RH**gamma)
      CASE(4)
**  NEBULA POLYTROPIC (SET BY OUTER BOUNDARY CONDITION)
         K_poly = RGAS*T_nebula*(RH_nebula**(1.d0-GAMMA))/mu_gas
         PG = K_poly*(RH**(gamma))
      CASE(5:)
         print *,'igastyp undefined',igastyp
         stop
      CASE DEFAULT
         print *,'igastyp undefined',igastyp
         print *,'you may also need to redefine cv'
         stop
      END SELECT
      RETURN
      END SUBROUTINE PRESSG

      SUBROUTINE QDENS
**  Calculates Density at intermediate Gridpoints (A-Grid)
**                and the Boundary Conditions
      USE fluid_var_init
      USE grid_var_init
      IMPLICIT NONE
      integer i,j,k
      DO K=-1,locNZ+2
         DO J=-1,locNY+2
            DO I=0,locNX+2
               RHQX(I,J,K)=DDX1(I)*RH(I,J,K)+DDX0(I)*RH(I-1,J,K)
            ENDDO
         ENDDO
      ENDDO
      DO K=-1,locNZ+2
         DO J=0,locNY+2
            DO I=-1,locNX+2
               RHQY(I,J,K)=DDY1(J)*RH(I,J,K)+DDY0(J)*RH(I,J-1,K)
            ENDDO
         ENDDO
      ENDDO
      DO K=0,locNZ+2
         DO J=-1,locNY+2
            DO I=-1,locNX+2
               RHQZ(I,J,K)=DDZ1(K)*RH(I,J,K)+DDZ0(K)*RH(I,J,K-1)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE QDENS


      SUBROUTINE SOUND
**  Calculates the sound velocity of the gas,
      USE input_init
      USE fluid_var_init
      USE global_constants
      USE grid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: Tgas
      SELECT CASE(icstyp)
      CASE(0)
**  Cs=sqrt(Gam*Pg/Rh)=sqrt(Gam*Rgas*T/XXMY), adiabatic sound speed
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               DO I=-1,locNX+2
                  cs(i,j,k)=sqrt(gamma*pg(i,j,k)/rh(i,j,k))
               ENDDO
            ENDDO
         ENDDO
      CASE(1)
** Constant soundspeed (at Tgas)
         Tgas = 500.d0
         CS = (Rgas*Tgas/mu_gas)**0.5
      CASE(2)
** Constant aspect ratio (H/a)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               DO I=-1,locNX+2
                  CS(I,J,K)=aspect*((GRAV*MSOL/xxb(i))**(0.5))
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
** Locally Isothermal with T-fixed in time
         CS=(gamma*Rgas*T/mu_gas)**0.5
      CASE(4:)
         print *,'icstyp=',icstyp,'not defined in sound'
         stop
      CASE DEFAULT
         print *,'icstyp=',icstyp,'also not defined in sound'
         stop
      END SELECT
      RETURN
      END SUBROUTINE SOUND


      SUBROUTINE CALCINTERNALENERGY(internal)
      USE fluid_var_init
      IMPLICIT NONE
      double precision :: internal(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      integer :: i,j,k
c      internal = RH*T*CV
      DO K=-1,locNZ+2
         DO J=-1,locNY+2
            DO I=0,locNX+2
c               internal(i,j,k) = G(i,j,k)
               internal(i,j,k) = V(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CALCINTERNALENERGY
