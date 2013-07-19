      SUBROUTINE BOUNDS
!---  General MPI SHIFT routines and Boundary Conditions
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE scalar_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: f,n
      integer :: i,j,k
!---  CALL SHIFTING ROUTINE
      if(modtyp.eq.2) then !-- shifting routine for planet mpi grid setup
!--   the standard shifting routine
         do f=1,5
            CALL PLANETSHIFT(f,-1)
         enddo
!--   now correct the polar shifts if necessary
         if(poles) then
            if(proc_coords(3).eq.proc_dims(3)-1) then
               do f=1,5
                  CALL NORTH_POLE_SHIFT(f,-1)
               enddo
            endif
            if(proc_coords(3).eq.0) then
               do f=1,5
                  CALL SOUTH_POLE_SHIFT(f,-1)
               enddo
            endif
         endif
!--  routine for 3D disk w/2D processors
      elseif(modtyp.eq.3.and.modver.eq.3) then
         do f=1,5
            CALL SHIFTVAR2D(f)
         enddo
!--  xz, axi-symmetric shifting routine
      elseif(modtyp.eq.6) then
         do f=1,5
            CALL SHIFTVAR_AXISYMMETRIC(f)
         enddo
!-- shifting routine for everything else
      else
         if(locNZ.eq.1) then
            do f=1,4
               CALL SHIFTVAR2D(f)
            enddo
         else
            do f=1,5
               CALL NEWSHIFTVAR3D(f)
            enddo
         endif
      endif
!---- OVERALL BOUNDARY CONDITIONS      
      SELECT CASE(MODTYP)
      CASE (1)
         CALL BOUNDCARTESIAN
      CASE (2)
         CALL BOUNDPLANET
c         if(iscal.gt.0) CALL SCALAR_BOUNDS
      CASE (3)
         CALL BOUNDDISK
      CASE (4)
         CALL BOUNDDREDGE
      CASE (5)
         CALL BOUNDDISK
      CASE (6)
         CALL BOUNDPACC
      CASE (7:)
         print *,'modtyp=',modtyp,'not yet defined: BOUNDS' 
      CASE DEFAULT
         print *,'modtyp=',modtyp,'also not defined: BOUNDS'
      END SELECT

      if(modtyp.eq.2) then !-- shifting routine for planet mpi grid setup
         do f=1,5
            CALL PLANETSHIFT(f,-1)
         enddo
!--   now correct the polar shifts if necessary
         if(poles) then
            if(proc_coords(3).eq.proc_dims(3)-1) then
               do f=1,5
                  CALL NORTH_POLE_SHIFT(f,-1)
               enddo
            endif
            if(proc_coords(3).eq.0) then
               do f=1,5
                  CALL SOUTH_POLE_SHIFT(f,-1)
               enddo
            endif
         endif
      endif
      RETURN
      END SUBROUTINE BOUNDS

      SUBROUTINE BOUNDPLANET
!-- PLANET BOUNDARY ROUTINE      
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      USE grid_var_init
      USE fluid_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      integer :: shift_coord,shift_dir
      integer :: stat(MPI_STATUS_SIZE)
      integer :: UPPERBND_TYPE
      double precision :: VELRAMP,V_START
      double precision :: V_POLE_AVG(0:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: density_at_mbar
      integer :: local_max_ub_indx,global_max_ub_indx

      if(planet_num.eq.5) then
         UPPERBND_TYPE = 1      !- 1=FIXED VALUE
      else
         UPPERBND_TYPE = 0      !-0=DENSITY BOUNDARY
      endif

      if(modver.ne.1) then
         IF(UPPERBND_TYPE.EQ.0) THEN !--- FIND LOCATION OF UPPER BOUNDARY BASED ON DENSITY, EXCEPT FOR BD MODELS
            if(III.eq.0.and.myid.eq.0) then
               print *,'Utilizing density for upperbnd calculation'
            endif
            DO K=1,locNZ
               DO J=1,locNY
!-----INITIALLY ASSUME THE BOUNDARY IS ABOVE DOMAIN
                  upperbnd(j,k) = locNX+1
                  DO I=1,locNX
                     if(RH(i,j,k).lt.1.0e-9) then
                        upperbnd(j,k) = i
                        GOTO 11
                     endif
                  ENDDO
!-----IF CONDITION IS NOT MET, SET UPPER CELL TO locNX-1
                  upperbnd(j,k) = locNX-1
 11               CONTINUE
               ENDDO
            ENDDO
!--Put an upperlimit to the upperbnd for the Roche overflow simulation
            if(modtyp.eq.2.and.modver.eq.4.and.planet_num.eq.5) then
               call ROCHE_OVERFLOW
            endif
!---  PASS INFO BETWEEN PROCESSORS
            CALL PASSUPPERBND
         ENDIF
         IF (UPPERBND_TYPE.EQ.1) THEN !---  FIX LOCATION OF THE UPPERBND
c            upperbnd(:,:) = 150
            upperbnd(:,:) = 97
c            upperbnd(:,:) = 70
            if(III.eq.0.and.myid.eq.0) then
               print *,'Utilizing fixed upperbnd of',upperbnd(1,1)
            endif
         ENDIF
      endif


!------ NOW DO THE BOUNDARY CONDITIONS
!--   modver = 0,2,3,4,5
      if(modver.eq.0.or.modver.eq.2.or.modver.eq.3.or.
     %     modver.eq.4.or.modver.eq.5) then
         if(BCXMIN.eq.5) then 
            CALL BOUNDXMIN(5)   !-- use to fix D_T and damp V
         elseif(BCXMIN.eq.6) then 
            CALL BOUNDXMIN(6)   !-- use to fix FLUX, slip free
         elseif(BCXMIN.eq.7) then 
            CALL BOUNDXMIN(7)   !-- use to fix D_T
         elseif(BCXMIN.eq.8) then 
            CALL BOUNDXMIN(8)   !-- use to fix D_T, with duphi/dr=0
         elseif(BCXMIN.eq.9) then 
            CALL BOUNDXMIN(9)   !--fixed flux, but free T,dTdr=f(F,T),V=H=0, duphi/dr=0
         elseif(BCXMIN.eq.11) then 
            CALL BOUNDXMIN(11)   !--fixed all values and derivatives at the bottom
         elseif(BCXMIN.eq.12) then 
            CALL BOUNDXMIN(12)
         else
            print *,'bcxmin=',bcxmin,'also not defined: BOUNDPLANET'
            stop
         endif
         CALL BOUNDXMAX(2)
         if(MPIlower.eq.MPI_PROC_NULL) then
            print *,'myid=',myid,'is in ymin'
            stop
         endif
         if(MPIupper.eq.MPI_PROC_NULL) then
            print *,'myid=',myid,'is in ymax'
            stop
         endif
         if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
            call BOUNDZMIN(0)
         endif
         if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL) then
            call BOUNDZMAX(0)
         endif

!- at the pole, the azimuthal velocity -> 0, just like if there really was a pole
         if(poles.and.proc_coords(3).eq.0) then
            DO I=-1,locNX+2
               DO J=0,locNY+2
                  G(I,j,1)=0.d0
                  G(I,j,0)=0.d0
                  G(I,j,-1)=0.d0
               ENDDO
            ENDDO
         endif
         if(poles.and.proc_coords(3).eq.proc_dims(3)-1) then
            DO I=-1,locNX+2
               DO J=0,locNY+2
                  G(I,j,locNZ)=0.d0
                  G(I,j,locNZ+1)=0.d0
                  G(I,j,locNZ+2)=0.d0
               ENDDO
            ENDDO
         endif

!--   modver=1 conditions for browndwarf
      elseif(modver.eq.1) then
!--   set floor for the density to 10^-8
         RH = MAX(RH,10.d0**(-8.d0))
         if(BCXMIN.eq.5) then 
            CALL BOUNDXMIN(5)   !-- use to fix D_T and damp V
         elseif(BCXMIN.eq.6) then 
            CALL BOUNDXMIN(6)   !-- use to fix FLUX and damp V
         elseif(BCXMIN.eq.7) then 
            CALL BOUNDXMIN(7)   !-- use to fix D_T
         elseif(BCXMIN.eq.8) then 
            CALL BOUNDXMIN(8)   !-- use to fix D_T, with duphi/dr=0
         elseif(BCXMIN.eq.9) then 
            CALL BOUNDXMIN(9)   !--fixed flux, but free T, dTdr=f(F,T), V=H=0, duphi/dr=0
         elseif(BCXMIN.eq.11) then 
            CALL BOUNDXMIN(11)   !--fixed all values and derivatives at the bottom
         elseif(BCXMIN.eq.12) then 
            CALL BOUNDXMIN(12)
         else
            print *,'bcxmin=',bcxmin,'not defined: BOUNDPLANET'
         endif
         CALL BOUNDXMAX(1)
         if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
            call BOUNDZMIN(0)
         endif
         if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL) then
            call BOUNDZMAX(0)
         endif
      elseif(modver.ge.6) then
         print *,'modver=',modver,'not yet defined: BOUNDPLANET' 
      else
         print *,'modver=',modver,'also not yet defined: BOUNDPLANET' 
      endif
      RETURN
      END SUBROUTINE BOUNDPLANET

      SUBROUTINE BOUNDCARTESIAN
!-- CONVECTION BOUNDARY ROUTINE      
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      SELECT CASE(MODVER)
      CASE (0)
         if(MPIleft.eq.MPI_PROC_NULL) then
            CALL BOUNDXMIN(0)
         endif
         if(MPIright.eq.MPI_PROC_NULL) then
            CALL BOUNDXMAX(0)
         endif
         if(MPIlower.eq.MPI_PROC_NULL) then
            call BOUNDYMIN(0)
         endif
         if(MPIupper.eq.MPI_PROC_NULL) then
            call BOUNDYMAX(0)
         endif
         if(global_NZ.eq.1) then
            print *,'does this mess up the boundary for the 2d case?'
            call clean_stop
         endif
         if(MPIbelow.eq.MPI_PROC_NULL) then
            call BOUNDZMIN(0)
         endif
         if(MPIabove.eq.MPI_PROC_NULL) then
            call BOUNDZMAX(0)
         endif
      CASE (1) !-- Sod shock wave test
         if(MPIleft.eq.MPI_PROC_NULL) then
            CALL BOUNDXMIN(2) !--open
         endif
         if(MPIright.eq.MPI_PROC_NULL) then
            CALL BOUNDXMAX(7) !--open
         endif
         if(MPIlower.eq.MPI_PROC_NULL.or.MPIupper.eq.MPI_PROC_NULL.or.
     %      MPIbelow.eq.MPI_PROC_NULL.or.MPIabove.eq.MPI_PROC_NULL) then
            print *,'error in the MPI designation:boundcartesian'
         endif
      CASE (2)
         if(MPIleft.eq.MPI_PROC_NULL) then
c            CALL BOUNDXMIN(0)   !-constant dT/dr
            CALL BOUNDXMIN(1)   !-constant T
         endif
         if(MPIright.eq.MPI_PROC_NULL) then
            CALL BOUNDXMAX(0)   !-fixed T
c            CALL BOUNDXMAX(8)   !-zero dT/dx
         endif
         if(MPIlower.eq.MPI_PROC_NULL.or.MPIupper.eq.MPI_PROC_NULL.or.
     $        MPIbelow.eq.MPI_PROC_NULL.or.MPIabove.eq.MPI_PROC_NULL) 
     $        then
            print *,'ERROR in processor boundary: ppbou'
            print *,' should be periodic in the y and z directions'
            call CLEAN_STOP
         endif
      CASE (3:)
         print *,'modver=',modver,'not yet defined: BOUNDCARTESIAN'
         call clean_stop
      CASE DEFAULT
         print *,'modver=',modver,'also not defined: BOUNDCARTESIAN'
         call clean_stop
      END SELECT
      RETURN
      END SUBROUTINE BOUNDCARTESIAN

      SUBROUTINE BOUNDDISK
!-- DISK BOUNDARY ROUTINE      
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      SELECT CASE(MODVER)
      CASE (0)
         if(MPIleft.eq.MPI_PROC_NULL) then
            CALL BOUNDXMIN(3)
         endif
         if(MPIright.eq.MPI_PROC_NULL) then
            CALL BOUNDXMAX(3)
         endif
         if(MPIlower.eq.MPI_PROC_NULL) then
            call BOUNDYMIN(1)
         endif
         if(MPIupper.eq.MPI_PROC_NULL) then
            call BOUNDYMAX(1)
         endif
         H = 0.d0
      CASE (1)
         if(MPIleft.eq.MPI_PROC_NULL) then
            CALL BOUNDXMIN(3)
         endif
         if(MPIright.eq.MPI_PROC_NULL) then
            CALL BOUNDXMAX(3)
         endif
         if(MPIlower.eq.MPI_PROC_NULL) then
            call BOUNDYMIN(1)
            print *,'myid=',myid,'is in ymin:BOUNDDISK'
            stop
         endif
         if(MPIupper.eq.MPI_PROC_NULL) then
            call BOUNDYMAX(1)
            print *,'myid=',myid,'is in ymax:BOUNDDISK'
            stop
         endif
      CASE (2:)
         print *,'modver=',modver,'not yet defined: BOUNDDISK' 
      CASE DEFAULT
         print *,'modver=',modver,'also not defined: BOUNDDISK'
      END SELECT
      RETURN
      END SUBROUTINE BOUNDDISK


      SUBROUTINE BOUNDDREDGE
!-- DREDGE-UP BOUNDARY ROUTINE      
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      if(MPIleft.eq.MPI_PROC_NULL) then
c         SELECT CASE(MODVER)
c         CASE (0:1)
c            CALL BOUNDXMIN(0)   !-- slip lower boundary 
c         CASE (2) 
         CALL BOUNDXMIN(4)      !-- fixed T and RH
c         CASE(3:)
c            print *,'MODVER',MODVER,'not yet defined:BOUNDDREDGE'
c            stop
c         CASE DEFAULT
c            print *,'MODVER',MODVER,'also not defined:BOUNDDREDGE'
c            stop
c         END SELECT
      endif
      if(MPIright.eq.MPI_PROC_NULL) then
         SELECT CASE(MODVER)
         CASE (0) 
            CALL BOUNDXMAX(4) !-- no wind, fixed top temperature and density
         CASE (1:2) 
            CALL BOUNDXMAX(5) !-- disk wind, fixed top temperature and density
         CASE(3:)
            print *,'MODVER',MODVER,'not yet defined:BOUNDDREDGE'
            stop
         CASE DEFAULT
            print *,'MODVER',MODVER,'also not defined:BOUNDDREDGE'
            stop
         END SELECT
      endif
      if(MPIlower.eq.MPI_PROC_NULL) then
         call BOUNDYMIN(1)
         print *,'myid=',myid,'is in ymin:BOUNDDREDGE'
         stop
      endif
      if(MPIupper.eq.MPI_PROC_NULL) then
         call BOUNDYMAX(1)
         print *,'myid=',myid,'is in ymax:BOUNDDREDGE'
         stop
      endif
      RETURN
      END SUBROUTINE BOUNDDREDGE

      SUBROUTINE BOUNDMIGRATION
!-- DISK BOUNDARY ROUTINE      
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      if(MPIleft.eq.MPI_PROC_NULL) then
         CALL BOUNDXMIN(3)
      endif
      if(MPIright.eq.MPI_PROC_NULL) then
         CALL BOUNDXMAX(3)
      endif
      if(MPIlower.eq.MPI_PROC_NULL) then
         print *,'myid=',myid,'is in ymin:BOUNDDISK'
         stop
      endif
      if(MPIupper.eq.MPI_PROC_NULL) then
         print *,'myid=',myid,'is in ymax:BOUNDDISK'
         stop
      endif
      call BOUNDZMIN(0)
      call BOUNDZMAX(0)
      RETURN
      END SUBROUTINE BOUNDMIGRATION

      SUBROUTINE BOUNDPACC
!-- PLANET ACCRETION BOUNDARY ROUTINE
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      USE fluid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      integer :: shift_coord,shift_dir
      integer :: stat(MPI_STATUS_SIZE)
!------ NOW DO THE BOUNDARY CONDITIONS
      if(modver.eq.0) then !--   modver=0
         RH = MAX(RH,RH_nebula) !-- set floor for the density to the nebula density
         if(MPIleft.eq.MPI_PROC_NULL) then
            if(BCXMIN.eq.5) then
               CALL BOUNDXMIN(5) !-- use to fix D_T and damp V
            elseif(BCXMIN.eq.6) then
               CALL BOUNDXMIN(6) !-- use to fix FLUX and damp V
            elseif(BCXMIN.eq.7) then
               CALL BOUNDXMIN(7) !-- use to fix D_T
            elseif(BCXMIN.eq.8) then
               CALL BOUNDXMIN(8) !-- use to fix D_T, with duphi/dr=0
            elseif(BCXMIN.eq.9) then
               CALL BOUNDXMIN(9) !--fixed flux, but free T, dTdr=f(F,T), V=H=0, duphi/dr=0
            elseif(BCXMIN.eq.10) then
               CALL BOUNDXMIN(10) !--RH, T, and H continuous; V=G=0
            elseif(BCXMIN.eq.11) then 
               CALL BOUNDXMIN(11) !--fixed all values and derivatives at the bottom
            else
               print *,'bcxmin=',bcxmin,'not defined: BOUNDPLANET'
            endif
         endif
         if(MPIright.eq.MPI_PROC_NULL) then
            CALL BOUNDXMAX(6) !-- NEBULA CONDITIONS; RH=RH_N; T=T_N; H=0; uphi Cont.; V=outflow 
         endif
         if(MPIbelow.eq.MPI_PROC_NULL) then
            call BOUNDZMIN(2) !-- RH, T, V, and G continous; H=0
         endif
         if(MPIabove.eq.MPI_PROC_NULL) then
            call BOUNDZMAX(2) !-- T=T_NEB; RH=RH_NEB; V=0; G Cont.; H=outflow
         endif
         CALL BOUNDYMIN(2) !- axi-symmetric
         CALL BOUNDYMAX(2) !- axi-symmetric



c         XMAX_B = XMAX-(DXA(1)/2.d0)
c               RH(0,j,K)  = RH_nebula*XMAX_B/ (
c     %              ((xxb(0)**(2.d0))+(xzb(k)**(2.d0)))**(0.5))
c               RH(-1,j,K) = RH_nebula*XMAX_B/ (
c     %              ((xxb(-1)**(2.d0))+(xzb(k)**(2.d0)))**(0.5))
c               T(1,j,K)  = T_nebula*(XMAX_B/ (
c     %              ((xxb(1)**(2.d0))+(xzb(k)**(2.d0)))**(0.5)))**
c     %              (GAMMA-1.d0)
c               T(0,j,K)  = T_nebula*(XMAX_B/ (
c     %              ((xxb(0)**(2.d0))+(xzb(k)**(2.d0)))**(0.5)))**
c     %              (GAMMA-1.d0)
c               T(-1,j,K) = T_nebula*(XMAX_B/ (
c     %              ((xxb(-1)**(2.d0))+(xzb(k)**(2.d0)))**(0.5)))**
c     %              (GAMMA-1.d0)



      elseif(modver.ge.1) then
         print *,'modver=',modver,'not yet defined: BOUNDPACC'
      else
         print *,'modver=',modver,'also not yet defined: BOUNDPACC'
      endif
      RETURN
      END SUBROUTINE BOUNDPACC

      SUBROUTINE PASSUPPERBND
!-- PLANET BOUNDARY ROUTINE      
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      USE fluid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      integer :: shift_coord,shift_dir
      integer :: stat(MPI_STATUS_SIZE)
c     note boundary passing in both y and z -directions chooses y-direction 
c     because y wraps around
!-- info from lowerbelow
      call MPI_SENDRECV(upperbnd(locNY,locNZ),1,MPI_INTEGER,
     %     MPIupperabove,1,upperbnd(0,0),1,MPI_INTEGER,
     %     MPIlowerbelow,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(locNY-1,locNZ-1),1,MPI_INTEGER,
     %     MPIupperabove,1,upperbnd(-1,-1),1,MPI_INTEGER,
     %     MPIlowerbelow,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(locNY,locNZ-1),1,MPI_INTEGER,
     %     MPIupperabove,1,upperbnd(0,-1),1,MPI_INTEGER,
     %     MPIlowerbelow,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(locNY-1,locNZ),1,MPI_INTEGER,
     %     MPIupperabove,1,upperbnd(-1,0),1,MPI_INTEGER,
     %     MPIlowerbelow,1,COMM_CART,stat,ierr)
!-- info from lower
      call MPI_SENDRECV(upperbnd(locNY,1:locNZ),locNZ,MPI_INTEGER,
     %     MPIupper,1,upperbnd(0,1:locNZ),locNZ,MPI_INTEGER,
     %     MPIlower,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(locNY-1,1:locNZ),locNZ,MPI_INTEGER,
     %     MPIupper,1,upperbnd(-1,1:locNZ),locNZ,MPI_INTEGER,
     %     MPIlower,1,COMM_CART,stat,ierr)
!-- info from lowerabove
      call MPI_SENDRECV(upperbnd(locNY,1),1,MPI_INTEGER,
     %     MPIupperbelow,1,upperbnd(0,locNZ+1),1,MPI_INTEGER,
     %     MPIlowerabove,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(locNY-1,2),1,MPI_INTEGER,
     %     MPIupperbelow,1,upperbnd(-1,locNZ+2),1,MPI_INTEGER,
     %     MPIlowerabove,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(locNY-1,1),1,MPI_INTEGER,
     %     MPIupperbelow,1,upperbnd(-1,locNZ+1),1,MPI_INTEGER,
     %     MPIlowerabove,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(locNY,2),1,MPI_INTEGER,
     %     MPIupperbelow,1,upperbnd(0,locNZ+2),1,MPI_INTEGER,
     %     MPIlowerabove,1,COMM_CART,stat,ierr)
!-- info from above
      call MPI_SENDRECV(upperbnd(1:locNY,1),locNY,MPI_INTEGER,
     %     MPIbelow,1,upperbnd(1:locNY,locNZ+1),locNY,MPI_INTEGER,
     %     MPIabove,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(1:locNY,2),locNY,MPI_INTEGER,
     %     MPIbelow,1,upperbnd(1:locNY,locNZ+2),locNY,MPI_INTEGER,
     %     MPIabove,1,COMM_CART,stat,ierr)
!-- info from upperabove
      call MPI_SENDRECV(upperbnd(1,1),1,MPI_INTEGER,
     %     MPIlowerbelow,1,upperbnd(locNY+1,locNZ+1),1,MPI_INTEGER,
     %     MPIupperabove,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(2,2),1,MPI_INTEGER,
     %     MPIlowerbelow,1,upperbnd(locNY+2,locNZ+2),1,MPI_INTEGER,
     %     MPIupperabove,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(1,2),1,MPI_INTEGER,
     %     MPIlowerbelow,1,upperbnd(locNY+1,locNZ+2),1,MPI_INTEGER,
     %     MPIupperabove,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(2,1),1,MPI_INTEGER,
     %     MPIlowerbelow,1,upperbnd(locNY+2,locNZ+1),1,MPI_INTEGER,
     %     MPIupperabove,1,COMM_CART,stat,ierr)
!-- info from upper
      call MPI_SENDRECV(upperbnd(1,1:locNZ),locNZ,MPI_INTEGER,
     %     MPIlower,1,upperbnd(locNY+1,1:locNZ),locNZ,MPI_INTEGER,
     %     MPIupper,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(2,1:locNZ),locNZ,MPI_INTEGER,
     %     MPIlower,1,upperbnd(locNY+2,1:locNZ),locNZ,MPI_INTEGER,
     %     MPIupper,1,COMM_CART,stat,ierr)
!-- info from upperbelow
      call MPI_SENDRECV(upperbnd(1,locNZ),1,MPI_INTEGER,
     %     MPIlowerabove,1,upperbnd(locNY+1,0),1,MPI_INTEGER,
     %     MPIupperbelow,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(2,locNZ-1),1,MPI_INTEGER,
     %     MPIlowerabove,1,upperbnd(locNY+2,-1),1,MPI_INTEGER,
     %     MPIupperbelow,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(1,locNZ-1),1,MPI_INTEGER,
     %     MPIlowerabove,1,upperbnd(locNY+1,-1),1,MPI_INTEGER,
     %     MPIupperbelow,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(2,locNZ),1,MPI_INTEGER,
     %     MPIlowerabove,1,upperbnd(locNY+2,0),1,MPI_INTEGER,
     %     MPIupperbelow,1,COMM_CART,stat,ierr)
!-- info from below
      call MPI_SENDRECV(upperbnd(1:locNY,locNZ),locNY,MPI_INTEGER,
     %     MPIabove,1,upperbnd(1:locNY,0),locNY,MPI_INTEGER,
     %     MPIbelow,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(1:locNY,locNZ-1),locNY,MPI_INTEGER,
     %     MPIabove,1,upperbnd(1:locNY,-1),locNY,MPI_INTEGER,
     %     MPIbelow,1,COMM_CART,stat,ierr)
!-- PROCESSOR BOUNDARIES

!--- set boundary for lowerbelow
c      if(.not.poles.and.MPIlowerbelow.eq.MPI_PROC_NULL) then
      if(MPIlowerbelow.eq.MPI_PROC_NULL) then
         upperbnd(0,0) = upperbnd(0,1)
         upperbnd(0,-1) = upperbnd(0,1)
         upperbnd(-1,0) = upperbnd(-1,1)
         upperbnd(-1,-1) = upperbnd(-1,1)
      endif
!--- set boundary for lower
      if(MPIlower.eq.MPI_PROC_NULL) then
         print *,'MPIlower=',MPIlower,'ERROR in passupperbnd'
         stop
      endif
!--- set boundary for lowerabove
c      if(.not.poles.and.MPIlowerabove.eq.MPI_PROC_NULL) then
      if(MPIlowerabove.eq.MPI_PROC_NULL) then
         upperbnd(0,locNZ+1) = upperbnd(0,locNZ)
         upperbnd(0,locNZ+2) = upperbnd(0,locNZ)
         upperbnd(-1,locNZ+1) = upperbnd(-1,locNZ)
         upperbnd(-1,locNZ+2) = upperbnd(-1,locNZ)
      endif
!--- set boundary for above
c      if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL) then
      if(MPIabove.eq.MPI_PROC_NULL) then
         do j=1,locNY
            upperbnd(j,locNZ+1)=upperbnd(j,locNZ)
            upperbnd(j,locNZ+2)=upperbnd(j,locNZ)
         enddo
      endif
!--- set boundary for upperabove
c      if(.not.poles.and.MPIupperabove.eq.MPI_PROC_NULL) then
      if(MPIupperabove.eq.MPI_PROC_NULL) then
         upperbnd(locNY+1,locNZ+1) = upperbnd(locNY+1,locNZ)
         upperbnd(locNY+1,locNZ+2) = upperbnd(locNY+1,locNZ)
         upperbnd(locNY+2,locNZ+1) = upperbnd(locNY+2,locNZ)
         upperbnd(locNY+2,locNZ+2) = upperbnd(locNY+2,locNZ)
      endif
!--- set boundary for upper
      if(MPIupper.eq.MPI_PROC_NULL) then
         print *,'MPIupper=',MPIupper,'ERROR in passupperbnd'
         stop
      endif
!--- set boundary for upperbelow
c      if(.not.poles.and.MPIupperbelow.eq.MPI_PROC_NULL) then
      if(MPIupperbelow.eq.MPI_PROC_NULL) then
         upperbnd(locNY+1,0) = upperbnd(locNY+1,1)
         upperbnd(locNY+1,-1) = upperbnd(locNY+1,1)
         upperbnd(locNY+2,0) = upperbnd(locNY+2,1)
         upperbnd(locNY+2,-1) = upperbnd(locNY+2,1) 
      endif
!--- set boundary for below
c      if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
      if(MPIbelow.eq.MPI_PROC_NULL) then
         do j=1,locNY
            upperbnd(j,0)=upperbnd(j,1)
            upperbnd(j,-1)=upperbnd(j,1)
         enddo
      endif
c      if(poles.and.proc_coords(3).eq.proc_dims(3)-1) then
c         CALL NORTH_POLE_PASSUPPERBND
c      endif
c      if(poles.and.proc_coords(3).eq.0) then
c         CALL SOUTH_POLE_PASSUPPERBND
c      endif
      RETURN
      END SUBROUTINE PASSUPPERBND

      SUBROUTINE NORTH_POLE_PASSUPPERBND
!-routine to reset the upperbound passing for the north pole
!-- communicate only with MPIabove_pole because above_pole already has
!   ghost cells filled
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      USE fluid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      integer :: shift_coord,shift_dir
      integer :: stat(MPI_STATUS_SIZE)
!-- info from above
      call MPI_SENDRECV(upperbnd(-1:locNY+2,locNZ-2),locNY+4,MPI_INTEGER
     %     ,MPIabove_pole,1,upperbnd(-1:locNY+2,locNZ+1),locNY+4,
     %     MPI_INTEGER,MPIabove_pole,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(-1:locNY+2,locNZ-1),locNY+4,MPI_INTEGER
     %     ,MPIabove_pole,1,upperbnd(-1:locNY+2,locNZ),locNY+4,
     %     MPI_INTEGER,MPIabove_pole,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(-1:locNY+2,locNZ-3),locNY+4,MPI_INTEGER
     %     ,MPIabove_pole,1,upperbnd(-1:locNY+2,locNZ+2),locNY+4,
     %     MPI_INTEGER,MPIabove_pole,1,COMM_CART,stat,ierr)
      RETURN
      END SUBROUTINE NORTH_POLE_PASSUPPERBND

      SUBROUTINE SOUTH_POLE_PASSUPPERBND
!-routine to reset the upperbound passing for the south pole
!-- communicate only with MPIbelow_pole because below_pole already has
!   ghost cells filled
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      USE fluid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      integer :: shift_coord,shift_dir
      integer :: stat(MPI_STATUS_SIZE)
!-- info from MPIbelow_pole
      call MPI_SENDRECV(upperbnd(-1:locNY+2,1),locNY+4,MPI_INTEGER,
     %     MPIbelow_pole,1,upperbnd(-1:locNY+2,0),locNY+4,MPI_INTEGER,
     %     MPIbelow_pole,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(upperbnd(-1:locNY+2,2),locNY+4,MPI_INTEGER,
     %     MPIbelow_pole,1,upperbnd(-1:locNY+2,-1),locNY+4,MPI_INTEGER,
     %     MPIbelow_pole,1,COMM_CART,stat,ierr)
      RETURN
      END SUBROUTINE SOUTH_POLE_PASSUPPERBND


      SUBROUTINE BOUNDXMIN(type)
      USE global_constants
      USE fluid_var_init
      USE input_init
      USE grid_var_init
      USE mpi_var_init
      USE rad_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: type,i,j,k,globalj,globalk
      double precision :: tboundary,bottomDTDR,boundary_ramp
      double precision :: kappaT_bottom
      double precision :: CALC_UNIVERSAL_BOTTOM_DERIV
! type = 0:  closed, fixed bottom T-gradient
! type = 1:  closed, fixed bottom T
! type = 2:  open
! type = 3:  keplarian
! type = 4:  closed, no-slip, fixed bottom T and rho
! type = 5:  fixed T and rh derivatives, slip vel (d(G/r)/dr = d(H/r)/dr = 0)
! type = 6:  Flux set by ipp (yielding fixed T,rh,dT/dr,& dRHdr), slip vel: (d(G)/dr = d(H)/dr = 0)
!               vr=0
! type = 7:  use to fix D_T
! type = 8:  fixed flux, free T, dTdr=f(F,<T>)=UNIVERSAL,Stree-free, inpenetrable:(dG/dr=dH/dr=0,v=0)
! type = 9:  fixed flux, free T, dTdr=f(F,T),Stree-free, inpenetrable:(dG/dr=dH/dr=0,v=0)
! type = 10: fixed bottom RH and T, 1/r in RH and T, zero velocities
! type = 11:  Fixed T,rh,dT/dr,& dRHdr.. derivatives set by ppinit values, V=G=H=0
! type = 12:  Stree free, inpenetrable: (dG/dr=dH/dr=0, v=0), fixed D_T
      SELECT CASE(type)
      CASE (0) !-- closed, fixed bottom dT/dr
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               T(0,j,K)   = T(1,j,k)+D_T*(xxb(0)-xxb(1))
               T(-1,j,K)  = T(1,j,k)+D_T*(xxb(-1)-xxb(1))
               V(1,J,K)   = 0.d0
               V(0,J,K)   = 0.d0
            ENDDO
         ENDDO
! Stree free, dG/dr=dH/dr=0. As in Gilman and Glatz (1981)
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(0,J,K)=G(1,J,K)
               G(-1,J,K)=G(1,J,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(0,J,K) = H(1,j,k)
               H(-1,J,K) = H(1,j,k)
            ENDDO
         ENDDO
      CASE (1) !-- closed, fixed bottom T
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               T(1,j,K)   = Tinitbottom
               T(0,j,K)   = Tinitbottom
               T(-1,j,K)  = Tinitbottom
               V(1,J,K)   = 0.d0
               V(0,J,K)   = 0.d0
            ENDDO
         ENDDO
         DO K=-1,locNZ+2! Stree free, dG/dr=dH/dr=0. As in Gilman and Glatz (1981)
            DO J=0,locNY+2
               G(0,J,K)=G(1,J,K)
               G(-1,J,K)=G(1,J,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(0,J,K) = H(1,j,k)
               H(-1,J,K) = H(1,j,k)
            ENDDO
         ENDDO
      CASE (2) !-- open
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(0,J,K)  = RH(1,J,K)
               RH(-1,J,K) = RH(1,J,K)
               T(0,j,K)   = T(1,J,K)
               T(-1,j,K)  = T(1,J,K)
               V(1,J,K)   = V(2,J,K)
               V(0,J,K)   = V(2,J,K)
               if(j.ne.-1) then
                  G(0,J,K)  = G(1,j,k)
                  G(-1,J,K) = G(1,j,k)
               endif
               if(k.ne.-1) then
                  H(0,J,K)  = H(1,j,k)
                  H(-1,J,K) = H(1,j,k)
               endif
            ENDDO
         ENDDO
      CASE (3) !-- 2DKEPLARIAN
c--- 2D Tanigawa and Watanabe (outflow with zero pressure gradient)
         K=1
         DO J=1,locNY
c-- zero pressure gradient means constant density
cx-------------------
            RH(0,J,k)=RH(1,J,k)
            RH(-1,J,k)=RH(1,J,k)
c-- keep temperature profile constant
            T(0,J,k)=T(1,J,k)*((xxb(0)/xxb(1))**(-0.5))
            T(-1,J,k)=T(1,J,k)*((xxb(-1)/xxb(1))**(-0.5))
c--   outflow, with damping
            V(4,J,K)   = 0.85*V(4,J,K)
            V(3,J,K)   = 0.60*V(3,J,K)
            V(2,J,K)   = 0.35*V(2,J,K)
            V(1,J,K)   = min(V(2,J,K),0.d0)
            V(0,J,K)   = v(1,j,k)
c-- for zero pressure gradient, the flow becomes purly Keplarian
            G(1,j,k) = ((GRAV*MSOL/(xxb(1)**3.d0))**0.5)-OMROT
            G(0,j,k) = ((GRAV*MSOL/(xxb(0)**3.d0))**0.5)-OMROT
            G(-1,j,k) = ((GRAV*MSOL/(xxb(-1)**3.d0))**0.5)-OMROT
         ENDDO
!-- lowerleft processor sets j=-1:0
         if(MPIlower.eq.MPI_PROC_NULL) then 
            DO J=-1,0
               RH(0,J,k)=RH(1,J,k)
               RH(-1,J,k)=RH(1,J,k)
               T(0,J,k)=T(1,J,k)*((xxb(0)/xxb(1))**(-0.5))
               T(-1,J,k)=T(1,J,k)*((xxb(-1)/xxb(1))**(-0.5))
               V(4,J,K)   = 0.85*V(4,J,K)
               V(3,J,K)   = 0.60*V(3,J,K)
               V(2,J,K)   = 0.35*V(2,J,K)
               V(1,J,K)   = min(V(2,J,K),0.d0)
               V(0,J,K)   = v(1,j,k)
               if(j.ne.-1) then
                  G(1,j,k) = ((GRAV*MSOL/(xxb(1)**3.d0))**0.5)-OMROT
                  G(0,j,k) = ((GRAV*MSOL/(xxb(0)**3.d0))**0.5)-OMROT
                  G(-1,j,k) = ((GRAV*MSOL/(xxb(-1)**3.d0))**0.5)-OMROT
               endif
            ENDDO
         endif
!-- upperleft processor sets j=locNY+1,locNY+2
         if(MPIupper.eq.MPI_PROC_NULL) then 
            DO J=locNY+1,locNY+2
               RH(0,J,k)=RH(1,J,k)
               RH(-1,J,k)=RH(1,J,k)
               T(0,J,k)=T(1,J,k)*((xxb(0)/xxb(1))**(-0.5))
               T(-1,J,k)=T(1,J,k)*((xxb(-1)/xxb(1))**(-0.5))
               V(4,J,K)   = 0.85*V(4,J,K)
               V(3,J,K)   = 0.60*V(3,J,K)
               V(2,J,K)   = 0.35*V(2,J,K)
               V(1,J,K)   = min(V(2,J,K),0.d0)
               V(0,J,K)   = v(1,j,k)
               G(1,j,k) = ((GRAV*MSOL/(xxb(1)**3.d0))**0.5)-OMROT
               G(0,j,k) = ((GRAV*MSOL/(xxb(0)**3.d0))**0.5)-OMROT
               G(-1,j,k) = ((GRAV*MSOL/(xxb(-1)**3.d0))**0.5)-OMROT
            ENDDO
         endif
      CASE (4) !-- closed, no-slip, fixed bottom T and rho
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(1,J,K)  = rhinitbottom
               RH(0,J,K)  = rhinitbottom
               RH(-1,J,K) = rhinitbottom
               T(1,j,K)   = TinitBottom 
               T(0,j,K)   = TinitBottom 
               T(-1,j,K)  = TinitBottom
               V(1,J,K)   = 0.d0
               V(0,J,K)   = -V(2,J,K)
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(0,J,K) = 0.d0
               G(-1,J,K) = 0.d0
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(0,J,K) = 0.d0
               H(-1,J,K) = 0.d0
            ENDDO
         ENDDO
!-- fixed T and rh derivatives, damped vr near the bottom
!-- slip vel: (d(G*r*cos(theta))/dr = d(H*r)/dr = 0)
      CASE (5)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(0,J,K)  = RH(1,J,K) - D_RH*(xxb(1)-xxb(0))
               RH(-1,J,K) = RH(1,J,K) - D_RH*(xxb(1)-xxb(-1))
               T(0,j,K)   = T(1,J,K) -  D_T*(xxb(1)-xxb(0))
               T(-1,j,K)  = T(1,J,K) -  D_T*(xxb(1)-xxb(-1))
C-- try slowly (spatially) damping out radial bounce back
               V(4,J,K)   = 0.85*V(4,J,K)
               V(3,J,K)   = 0.60*V(3,J,K)
               V(2,J,K)   = 0.35*V(2,J,K)
               V(1,J,K)   = 0.d0
               V(0,J,K)   = 0.d0
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
c               G(0,J,K)  = G(1,J,K)*xxb(1)/xxb(0)
c               G(-1,J,K) = G(1,J,K)*xxb(1)/xxb(-1)
               G(1,J,K)  = 0.d0
               G(0,J,K)  = 0.d0
               G(-1,J,K) = 0.d0

            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
c               H(0,J,K) = H(1,J,K)*xxb(1)/xxb(0)
c               H(-1,J,K) = H(1,J,K)*xxb(1)/xxb(-1)
               H(1,J,K) = 0.d0
               H(0,J,K) = 0.d0
               H(-1,J,K) = 0.d0
            ENDDO
         ENDDO
!-- fixed T and RH Flux's (set T, RH, dT/dr, and dRHdr)
!-- vr=0, slip vel: (d(G)/dr = d(H)/dr = 0)
      CASE (6)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(1,J,K)  = RHinitbottom
               RH(0,J,K)  = RH(1,J,K) - D_RH*(xxb(1)-xxb(0))
               RH(-1,J,K) = RH(1,J,K) - D_RH*(xxb(1)-xxb(-1))
               T(1,J,K)   = TinitBottom
               bottomDTDR = -3.d0*Fbottom*xkapR(1,j,k)*RHinitbottom/
     $              (4.d0*ARAD*FD*(Tinitbottom**3.d0))
               T(2,j,K)   = T(1,J,K) -  bottomDTDR*(xxb(1)-xxb(2))
               T(0,j,K)   = T(1,J,K) -  bottomDTDR*(xxb(1)-xxb(0))
               T(-1,j,K)  = T(1,J,K) -  bottomDTDR*(xxb(1)-xxb(-1))
c-zero vr
               V(1,J,K)   = 0.d0
               V(0,J,K)   = 0.d0
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(0,J,K)  = G(1,J,K)
               G(-1,J,K) = G(1,J,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(0,J,K) = H(1,J,K)
               H(-1,J,K) = H(1,J,K)
            ENDDO
         ENDDO
!-- fixed T and rh derivatives
!-- G = H = 0
      CASE (7)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(0,J,K)  = RH(1,J,K) - D_RH*(xxb(1)-xxb(0))
               RH(-1,J,K) = RH(1,J,K) - D_RH*(xxb(1)-xxb(-1))
               T(0,j,K)   = T(1,J,K) -  D_T*(xxb(1)-xxb(0))
               T(-1,j,K)  = T(1,J,K) -  D_T*(xxb(1)-xxb(-1))
               V(1,J,K)   = 0.d0
               V(0,J,K)   = -V(2,J,K)
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(1,J,K)  = 0.d0
               G(0,J,K)  = 0.d0
               G(-1,J,K) = 0.d0
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(1,J,K) = 0.d0
               H(0,J,K) = 0.d0
               H(-1,J,K) = 0.d0
            ENDDO
         ENDDO
!--fixed flux, but free T, dTdr=f(F,<T>)=UNIVERSAL
!   Stree-free, inpenetrable:(dG/dr=dH/dr=0,v=0)
      CASE (8)
         bottomDTDR = CALC_UNIVERSAL_BOTTOM_DERIV(.false.)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               T(0,j,K)   = T(1,J,K) -  bottomDTDR*(xxb(1)-xxb(0))
               T(-1,j,K)  = T(1,J,K) -  bottomDTDR*(xxb(1)-xxb(-1))
               V(1,J,K)   = 0.d0
               V(0,J,K)   = 0.d0
            ENDDO
         ENDDO
! Stree free, dG/dr=dH/dr=0. As in Gilman and Glatz (1981)
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(0,J,K)=G(1,J,K)
               G(-1,J,K)=G(1,J,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(0,J,K) = H(1,j,k)
               H(-1,J,K) = H(1,j,k)
            ENDDO
         ENDDO
!--fixed flux, but free T, dTdr=f(F,T),Stree-free, inpenetrable:(dG/dr=dH/dr=0,v=0)
      CASE (9)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               if(ieta.eq.5) then
                  bottomDTDR = -Fbottom/
     $                 (RH(1,j,k)*cv(1,j,k)*(XNUE/PR_NUM))
               else
                  if(III.gt.0) then 
                     bottomDTDR = -3.d0*Fbottom*xkapR(1,j,k)*RH(1,j,k)/
     $                    (4.d0*ARAD*FD*(T(1,j,k)**3.d0))
                  else
                     call BurrowsOPC(T(1,1,1),RH(1,1,1),xkapP(1,1,1),
     $                    xkapR(1,1,1),xkapA(1,1,1),xkapW(1,1,1,:),
     $                    xkapLW(1,1,1,:),local_BB(1,1,1,:),
     $                    dxkapP_dT(1,1,1),dxkapA_dT(1,1,1),1,1,1)                     
                     bottomDTDR = -3.d0*Fbottom*xkapR(1,1,1)*RH(1,1,1)/
     $                    (4.d0*ARAD*FD*(T(1,1,1)**3.d0))
                  endif
               endif
               T(0,j,K)   = T(1,J,K) -  bottomDTDR*(xxb(1)-xxb(0))
               T(-1,j,K)  = T(1,J,K) -  bottomDTDR*(xxb(1)-xxb(-1))
               V(1,J,K)   = 0.d0
               V(0,J,K)   = 0.d0
            ENDDO
         ENDDO
! Stree free, dG/dr=dH/dr=0. As in Gilman and Glatz (1981)
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(0,J,K)=G(1,J,K)
               G(-1,J,K)=G(1,J,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(0,J,K) = H(1,j,k)
               H(-1,J,K) = H(1,j,k)
            ENDDO
         ENDDO
      CASE (10) !--RH, T, and H continuous; V=G=0; 
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(0,j,K)  = RH(1,j,K)
               RH(-1,j,K) = RH(1,j,K)
               T(0,j,K)   = T(1,j,K)
               T(-1,j,K)  = T(1,j,K)
               V(1,J,K)   = 0.d0
               V(0,J,K)   = 0.d0
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(0,J,K)  = 0.d0
               G(-1,J,K) = 0.d0
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(0,J,K)  = H(1,J,K)
               H(-1,J,K) = H(1,J,K)
            ENDDO
         ENDDO
      CASE (11)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(1,J,K)  = RHinitbottom
               RH(0,J,K)  = RH(1,J,K) - D_RH*(xxb(1)-xxb(0))
               RH(-1,J,K) = RH(1,J,K) - D_RH*(xxb(1)-xxb(-1))
               T(1,J,K)  = Tinitbottom
               T(2,j,K)  = T(1,J,K) - D_T*(xxb(1)-xxb(2))
               T(0,J,K)  = T(1,J,K) - D_T*(xxb(1)-xxb(0))
               T(-1,J,K) = T(1,J,K) - D_T*(xxb(1)-xxb(-1))
               V(1,J,K)   = 0.d0
               V(0,J,K)   = 0.d0
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(1,J,K)=0.d0
               G(0,J,K)=0.d0
               G(-1,J,K)=0.d0
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(1,J,K) = 0.d0
               H(0,J,K) = 0.d0
               H(-1,J,K) = 0.d0
            ENDDO
         ENDDO
      CASE (12)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               T(0,j,K)   = T(1,j,k)+D_T*(xxb(0)-xxb(1))
               T(-1,j,K)  = T(1,j,k)+D_T*(xxb(-1)-xxb(1))
               V(1,J,K)   = 0.d0
               V(0,J,K)   = 0.d0
            ENDDO
         ENDDO
! Stree free, dG/dr=dH/dr=0. As in Gilman and Glatz (1981)
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(0,J,K)=G(1,J,K)
               G(-1,J,K)=G(1,J,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(0,J,K) = H(1,j,k)
               H(-1,J,K) = H(1,j,k)
            ENDDO
         ENDDO
      CASE (13:)
         print *,'type',type,'not yet defined: BOUNDXMIN' 
         stop
      CASE DEFAULT
         print *,'type',type,'also not defined: BOUNDXMIN'
         stop
      END SELECT
      RETURN
      END SUBROUTINE BOUNDXMIN

      SUBROUTINE BOUNDXMAX(type)
! type = 0: closed
! type = 1: open outflow
! type = 2: moveable - closed
! type = 3: keplarian
! type = 4: no wind, fixed top temperature and density
! type = 5: disk wind, fixed top temperature and density
! type = 6: nebula conditions, zero velocity
! type = 7: simple open for shock wave test
! type = 8: slip-free, zero T gradient
      USE global_constants
      USE fluid_var_init
      USE force_var_init
      USE input_init
      USE grid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE rad_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: type,i,j,k,min_indx_G,min_indx_H
      INTEGER :: RADBOUNDVER,vdamp_steps
      double precision :: TJupHeat!,phi0
      double precision :: Tincident,CALC_TINCIDENT,irr_ramp,HP
      double precision :: mslope,tau_vdamp,vdamp_ramp,depth_vrdamp
      SELECT CASE(type)
      CASE (0) !-- closed
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               T(locNX,J,K)=Tinittop
               T(locNX+1,J,K)=Tinittop
               T(locNX+2,J,K)=Tinittop
c               V(locNX+1,J,K)=0.d0
c               V(locNX+2,J,K)=0.d0


!-try this for inflow/outflow
               V(locNX+1,J,K)=V(locNX,J,K)
               V(locNX+2,J,K)=V(locNX,J,K)
               RH(locNX+1,J,K) = RH(locNX,J,K)
               RH(locNX+2,J,K) = RH(locNX,J,K)


c               RH(locNX,J,K)   = RHinittop
c               RH(locNX+1,J,K) = RHinittop
c               RH(locNX+2,J,K) = RHinittop
            ENDDO
         ENDDO
! Stree free, dG/dr=dH/dr=0. As in Gilman and Glatz (1981)
         DO K=-1,locNZ+2
            DO J=0,locNY+2
c               G(locNX,J,K)=G(locNX-1,J,K) !added for earth_iso_3D/E
               G(locNX+1,J,K)=G(locNX,J,K)
               G(locNX+2,J,K)=G(locNX,J,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
c               H(locNX,J,K)=H(locNX-1,J,K) !added for earth_iso_3D/E
               H(locNX+1,J,K)=H(locNX,J,K)
               H(locNX+2,J,K)=H(locNX,J,K)
            ENDDO
         ENDDO
      CASE (1) !-- outflow condition
         first_Tincident = .true.
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(locNX,J,K)   = RH(locNX-1,J,K)
               RH(locNX+1,J,K) = RH(locNX-1,J,K)
               RH(locNX+2,J,K) = RH(locNX-1,J,K)
               if(F_INCIDENT.eq.3.or.F_INCIDENT.eq.5) then !--- RAMP-UP
                  irr_ramp = (1.d0*III/(1.d0*DAMP_STEPS))
                  if(irr_ramp.ge.1.0) irr_ramp = 1.d0
               else
                  irr_ramp = 1.0
               endif
               Tincident = CALC_TINCIDENT(j,k)
               TjupHeat = (( (ER(upperbnd(j,k)-1,J,K)/ARAD)+0.25*
     %              ((gam_kap(upperbnd(j,k)-1,j,k)*
     %              Tincident)**4.d0)*irr_ramp*
     %              (exp(-stellartau(upperbnd(j,k)-1,j,k)))/
     $              mu_star(j,k))**(0.25))
               T(locNX,J,K)   = TjupHeat
               T(locNX+1,J,K) = TjupHeat
               T(locNX+2,J,K) = TjupHeat
               V(locNX-1,J,K) = max(V(locNX-1,J,K),0.d0)
               V(locNX,J,K)   = V(locNX-1,J,K)
               V(locNX+1,J,K) = V(locNX-1,J,K)
               V(locNX+2,J,K) = V(locNX-1,J,K)
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
c-- d(uphi)/dr = 0
               if(j.ne.-1) then
                  G(locNX,J,K)=G(locNX-1,J,K)*geoxg(locNX-1)/
     %                 geoxg(locNX)
                  G(locNX+1,J,K)=G(locNX-1,J,K)*geoxg(locNX-1)/
     %                 geoxg(locNX+2)
                  G(locNX+2,J,K)=G(locNX-1,J,K)*geoxg(locNX-1)/
     %                 geoxg(locNX+2)
c                  G(locNX,J,K)   = G(locNX-1,J,K)
c                  G(locNX+1,J,K) = G(locNX-1,J,K)
c                  G(locNX+2,J,K) = G(locNX-1,J,K)
               endif
               if(k.ne.-1) then
c--d(utheta)/dr=0
                  H(locNX,J,K)=H(locNX-1,J,K)*geoxg(locNX-1)/
     %                 geoxg(locNX)
                  H(locNX+1,J,K)=H(locNX-1,J,K)*geoxg(locNX-1)/
     %                 geoxg(locNX+2)
                  H(locNX+2,J,K)=H(locNX-1,J,K)*geoxg(locNX-1)/
     %                 geoxg(locNX+2)
c                  H(locNX,J,K)   = H(locNX-1,J,K)
c                  H(locNX+1,J,K) = H(locNX-1,J,K)
c                  H(locNX+2,J,K) = H(locNX-1,J,K)
               endif
            ENDDO
         ENDDO
      CASE (2) !-- movable (v=0, grad(G,RH,T)=0)
         RADBOUNDVER = IBDRAD
         first_Tincident = .true.
c         DO K=-1,locNZ+2
c            DO J=-1,locNY+2
         DO K=1,locNZ
            DO J=1,locNY
               if(upperbnd(j,k).lt.locNX+1) then
                  IF(RADBOUNDVER.eq.2) then
                     TJupHeat = Tirr*cos(xyb(j)/ymax)
                     TJupHeat = MAX(TJupHeat,100.d0)
                     T(upperbnd(j,k),j,k)= TJupHeat
                  ELSEIF(RADBOUNDVER.eq.3) then
                     if(ncosys.eq.1) then
                        print *,'N/A...'
                        stop
                     elseif(ncosys.eq.2) then
                        TJupHeat = Tirr*((cos(xyb(j))*
     %                       cos(xzb(k)))**0.25)
                     else
                        print *,'bad ncosys:boundxmax',ncosys
                        stop
                     endif
                     TJupHeat = MAX(TJupHeat,100.d0)
c                     TJupHeat = MAX(TJupHeat,500.d0)
                     T(upperbnd(j,k),j,k)= TJupHeat
                  ELSEIF(RADBOUNDVER.eq.4) then ! constant upper temperature
                     TJupHeat = Tirr
                     T(upperbnd(j,k),j,k)= TJupHeat
                  ELSEIF(RADBOUNDVER.eq.5) then ! thermal tide
                     if(ncosys.eq.2) then
c                        phi0 = (2.d0*PI/(thermalP_phi*DAY))*
c     %                       mod((zeit-thermal_t0*DAY),
c     %                       (thermalP_phi*DAY))
c                        TJupHeat = Tirr*((cos(xyb(j)-phi0)*
c     %                       cos(xzb(k)))**0.25)
                        print *,'redo thermal tide with Porb'
                        call clean_stop
                     else
                        print *,'bad ncosys:boundxmax (5)',ncosys
                        stop
                     endif
                     TJupHeat = MAX(TJupHeat,100.d0)
                     T(upperbnd(j,k),j,k)= TJupHeat
                  ELSEIF(RADBOUNDVER.eq.6) then ! temporal pulsations
c                     TJupHeat = Tirr*((cos(xyb(j))* !--standard spatial dependance
c     %                    cos(xzb(k)))**0.25)
c                     TJupHeat = TJupHeat*cos(2.d0*PI*(zeit-thermal_t0)/
c     %                    (0.1*DAY)) !-- temporal pulsations
c                     TJupHeat = MAX(TJupHeat,100.d0)


c                     T(upperbnd(j,k),j,k)= TJupHeat
                     print *,'redo temporal tide ',radboundver
                     call clean_stop
                  ELSEIF(RADBOUNDVER.eq.7) then ! gradT = 0 boundary condition
                     TJupHeat = T(upperbnd(j,k)-1,j,k)
                  ELSEIF(RADBOUNDVER.eq.8) then ! SAME AS 3, BUT DIFFERENT TMIN FOR BD
                     TJupHeat = Tirr*((cos(xyb(j))*
     %                    cos(xzb(k)))**0.25)
                     TJupHeat = MAX(TJupHeat,950.d0)
                  ELSEIF(RADBOUNDVER.eq.9) then
!--T=(ER/a+SI/(4*sigma*rh*kapP))^(0.25)
!--T=(Fx/(2*sigmarad)+SI/(4*sigma*rh*kapP))^(0.25)
                     IF (III.le.1) THEN
!-just the first time, because the opacities aren't defined
                        TJupHeat = T(upperbnd(j,k)-1,j,k)
                     else                        
                        if(irad.eq.5.or.irad.eq.6) then !-try using flux instead of ER (Aug 1, 2012):
                           TJupHeat = ( (flux_x(upperbnd(j,k),J,K)/
     $                          (2.d0*SBCONST))+
     $                          (stellarinput(upperbnd(j,k)-1,j,k)/
     $                          (4.d0*SBCONST*rh(upperbnd(j,k)-1,j,k)*
     $                          xkapP(upperbnd(j,k)-1,j,k))) )**(0.25)
                        else
                           TJupHeat = ( (ER(upperbnd(j,k)-1,J,K)/ARAD) +
     $                          (stellarinput(upperbnd(j,k)-1,j,k)/
     $                          (4.d0*SBCONST*rh(upperbnd(j,k)-1,j,k)*
     $                          xkapP(upperbnd(j,k)-1,j,k))) )**(0.25)
                        endif
c                        if(F_INCIDENT.eq.3.or.F_INCIDENT.eq.5) then !--- RAMP-UP
c                           irr_ramp = (1.d0*III/(1.d0*DAMP_STEPS))
c                           if(irr_ramp.ge.1.0) irr_ramp = 1.d0
c                        else
c                           irr_ramp = 1.0
c                        endif
c                        TJupHeat = Tjupheat*irr_ramp
                     endif
                  ELSEIF(RADBOUNDVER.eq.10) then !--like 9, but with a wavelength dep. source term
                     IF (III.le.1) THEN !-just the first time, because the opacities aren't defined
                        TJupHeat = T(upperbnd(j,k)-1,j,k) 
                     else
                        print *,'I HAVE NOT PUT THE BC IN BOUNDXMAX'
                        call CLEAN_STOP
                     endif
                  ENDIF
                  T(upperbnd(j,k),j,k) = TJupHeat

!-does this get rid of the spikes at upperbnd?
                  V(upperbnd(j,k),j,k)=
     $                 0.99*max(V(upperbnd(j,k),j,k),0.d0)

                  DO I=upperbnd(j,k)+1,locNX+2
                     T(I,J,K)=TJupHeat
!-try putting this back.. is it causing the random spikes
                     RH(i,j,k) = 1.0e-10
                     V(i,j,k)= 0.d0
                  ENDDO
                  IF(upperbnd(j,k).eq.locNX) THEN
                     V(locNX,j,k)= 0.d0
                     RH(locNX,j,k) = 1.0e-10
                     T(locNX,j,k)=TJupHeat
                  ENDIF
               ELSE
                  write(*,'(A,5(I8))') 
     $                 'error in boundxmax',III,upperbnd(j,k),j,k,myid
                  call clean_stop
               endif
            ENDDO
         ENDDO
! Stree free, dG/dr=dH/dr=0
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(upperbnd(j,k)+1:locNX+2,J,K)=G(upperbnd(j,k),J,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(upperbnd(j,k)+1:locNX+2,J,K)=H(upperbnd(j,k),J,K)
            ENDDO
         ENDDO
      CASE (3) !-- KEPLARIAN
         K=1
         DO J=1,locNY
c-- zero pressure gradient means constant density
            RH(locNX+1,J,k)=RH(locNX,J,K)
            RH(locNX+2,J,k)=RH(locNX,J,K)
c--   keep temperature profile constant
            T(locNX+1,J,k)=T(locNX,J,K)*
     %           ((xxb(locNX+1)/xxb(locNX))**(-0.5))
            T(locNX+2,J,k)=T(locNX,J,K)*
     %           ((xxb(locNX+2)/xxb(locNX))**(-0.5))
c--   outflow, with damping
            V(locNX-2,J,K)  = 0.85*V(locNX-2,J,K)
            V(locNX-1,J,K)  = 0.60*V(locNX-1,J,K)
            V(locNX,J,K)    = 0.35*V(locNX,J,K)
            V(locNX+1,J,K)  = max(0.d0,V(locNX,J,K))
            V(locNX+2,J,K)  = V(locNX+1,J,K)
c--   for zero pressure gradient, the flow becomes purly Keplarian
            G(locNX,J,K)=((GRAV*MSOL/(xxb(locNX)**3.d0))**0.5)-OMROT
            G(locNX+1,J,K)=((GRAV*MSOL/(xxb(locNX+1)**3.d0))**0.5)-
     %           OMROT
            G(locNX+2,J,K)=((GRAV*MSOL/(xxb(locNX+2)**3.d0))**0.5)-
     %           OMROT
         ENDDO
!-- lowerright processor sets j=-1:0
         if(MPIlower.eq.MPI_PROC_NULL) then 
            DO J=-1,0
               RH(locNX+1,J,k)=RH(locNX,J,K)
               RH(locNX+2,J,k)=RH(locNX,J,K)
               T(locNX+1,J,k)=T(locNX,J,K)*
     %              ((xxb(locNX+1)/xxb(locNX))**(-0.5))
               T(locNX+2,J,k)=T(locNX,J,K)*
     %              ((xxb(locNX+2)/xxb(locNX))**(-0.5))
               V(locNX-2,J,K)  = 0.85*V(locNX-2,J,K)
               V(locNX-1,J,K)  = 0.60*V(locNX-1,J,K)
               V(locNX,J,K)    = 0.35*V(locNX,J,K)
               V(locNX+1,J,K)  = max(0.d0,V(locNX,J,K))
               V(locNX+2,J,K)  = V(locNX+1,J,K)
               if(j.ne.-1) then
                  G(locNX,J,K)=((GRAV*MSOL/(xxb(locNX)**3.d0))**0.5)
     %                 -OMROT
                  G(locNX+1,J,K)=((GRAV*MSOL/(xxb(locNX+1)**3.d0))**0.5)
     %                 -OMROT
                  G(locNX+2,J,K)=((GRAV*MSOL/(xxb(locNX+2)**3.d0))**0.5)
     %                 -OMROT
               endif
            ENDDO
         endif
!-- upperright processor sets j=locNY+1:locNY+2
         if(MPIupper.eq.MPI_PROC_NULL) then 
            DO J=locNY+1,locNY+2
               RH(locNX+1,J,k)=RH(locNX,J,K)
               RH(locNX+2,J,k)=RH(locNX,J,K)
               T(locNX+1,J,k)=T(locNX,J,K)*
     %              ((xxb(locNX+1)/xxb(locNX))**(-0.5))
               T(locNX+2,J,k)=T(locNX,J,K)*
     %              ((xxb(locNX+2)/xxb(locNX))**(-0.5))
               V(locNX-2,J,K)  = 0.85*V(locNX-2,J,K)
               V(locNX-1,J,K)  = 0.60*V(locNX-1,J,K)
               V(locNX,J,K)    = 0.35*V(locNX,J,K)
               V(locNX+1,J,K)  = max(0.d0,V(locNX,J,K))
               V(locNX+2,J,K)  = V(locNX+1,J,K)
               G(locNX,J,K)=((GRAV*MSOL/(xxb(locNX)**3.d0))**0.5)
     %              -OMROT
               G(locNX+1,J,K)=((GRAV*MSOL/(xxb(locNX+1)**3.d0))**0.5)
     %              -OMROT
               G(locNX+2,J,K)=((GRAV*MSOL/(xxb(locNX+2)**3.d0))**0.5)
     %              -OMROT
            ENDDO
         endif
      CASE (4) !-- NO WIND, FIXED TOP DENSITY AND TEMPERATURE
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(locNX,J,K)=rhinittop
               RH(locNX+1,J,K)=rhinittop
               RH(locNX+2,J,K)=rhinittop
               T(locNX,J,K)=Tinittop
               T(locNX+1,J,K)=Tinittop
               T(locNX+2,J,K)=Tinittop
               V(locNX,J,K)=0.d0
               V(locNX+1,J,K)=0.d0
               V(locNX+2,J,K)=0.d0
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(locNX,J,K)=0.d0
               G(locNX+1,J,K)=0.d0
               G(locNX+2,J,K)=0.d0
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(locNX+1,J,K)=0.d0
               H(locNX+2,J,K)=0.d0
            ENDDO
         ENDDO
      CASE (5) !-- DISK WIND, FIXED TOP DENSITY AND TEMPERATURE
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(locNX,J,K)=rhinittop
               RH(locNX+1,J,K)=rhinittop
               RH(locNX+2,J,K)=rhinittop
               T(locNX,J,K)=Tinittop
               T(locNX+1,J,K)=Tinittop
               T(locNX+2,J,K)=Tinittop
               V(locNX,J,K)=windAMP*COS(XYB(j))
               V(locNX+1,J,K)=windAMP*COS(XYB(j))
               V(locNX+2,J,K)=windAMP*COS(XYB(j))
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(locNX,J,K)=-(windAMP/XXB(locNX))*SIN(XYA(j))
               G(locNX+1,J,K)=-(windAMP/XXB(locNX+1))*SIN(XYA(j))
               G(locNX+2,J,K)=-(windAMP/XXB(locNX+2))*SIN(XYA(j))
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(locNX+1,J,K)=0.d0
               H(locNX+2,J,K)=0.d0
            ENDDO
         ENDDO
      CASE (6) !-- NEBULA CONDITIONS; RH=RH_N; T=T_N; H=0; uphi Continuous; V=outflow 
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(locNX,J,K)   = RH_nebula
               RH(locNX+1,J,K) = RH_nebula
               RH(locNX+2,J,K) = RH_nebula
               T(locNX,J,K)    = T_nebula
               T(locNX+1,J,K)  = T_nebula
               T(locNX+2,J,K)  = T_nebula
               V(locNX,J,K)    = max(0.d0,V(locNX,J,K))
               V(locNX+1,J,K)  = max(0.d0,V(locNX,J,K))
               V(locNX+2,J,K)  = max(0.d0,V(locNX,J,K))
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(locNX+1,J,K) = G(locNX,J,K)*geoxg(locNX+1)/geoxg(locNX)
               G(locNX+2,J,K) = G(locNX,J,K)*geoxg(locNX+1)/geoxg(locNX)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(locNX+1,J,K)=0.d0
               H(locNX+2,J,K)=0.d0
            ENDDO
         ENDDO
      CASE (7) !-- OPEN
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               RH(locNX+1,J,K) = RH(locNX,J,K)
               RH(locNX+2,J,K) = RH(locNX,J,K)
               T(locNX+1,J,K)  = T(locNX,J,K)
               T(locNX+2,J,K)  = T(locNX,J,K)
               V(locNX+1,J,K)  = V(locNX,J,K)
               V(locNX+2,J,K)  = V(locNX,J,K)
               if(j.ne.-1) then
                  G(locNX+1,J,K) = G(locNX,J,K)
                  G(locNX+2,J,K) = G(locNX,J,K)
               endif
               if(k.ne.-1) then
                  H(locNX+1,J,K) = H(locNX,J,K)
                  H(locNX+2,J,K) = H(locNX,J,K)
               endif
            ENDDO
         ENDDO
      CASE (8) !-- slip-free and zero gradient
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               T(locNX,J,K)=T(locNX-1,J,K)
               T(locNX+1,J,K)=T(locNX,J,K)
               T(locNX+2,J,K)=T(locNX,J,K)
               V(locNX+1,J,K)=0.d0
               V(locNX+2,J,K)=0.d0
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO J=0,locNY+2
               G(locNX+1,J,K)=G(locNX,J,K)
               G(locNX+2,J,K)=G(locNX,J,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO J=-1,locNY+2
               H(locNX+1,J,K)=H(locNX,J,K)
               H(locNX+2,J,K)=H(locNX,J,K)
            ENDDO
         ENDDO
      CASE (9:)
         print *,'type',type,'not yet defined: BOUNDXMAX' 
         stop
      CASE DEFAULT
         print *,'type',type,'also not defined: BOUNDXMAX'
         stop
      END SELECT
      RETURN
      END SUBROUTINE BOUNDXMAX

      SUBROUTINE BOUNDYMAX(type)
! type = 0: closed
! type = 1: keplarian
      USE global_constants
      USE fluid_var_init
      USE input_init
      USE grid_var_init
      USE mpi_grid_init
      implicit none
      include "mpif.h"
      integer :: type,i,j,k
      SELECT CASE(type)
      case(0) !-- closed
         DO K=-1,locNZ+2
            DO I=-1,locNX+2
               RH(I,locNY+1,K)=RH(I,locNY,K)
               RH(I,locNY+2,K)=RH(I,locNY,K)
               T(I,locNY+1,K)=T(I,locNY,K)
               T(I,locNY+2,K)=T(I,locNY,K)
               G(I,locNY+1,K)=0.d0
               G(I,locNY+2,K)=-G(I,locNY,K)
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO I=0,locNX+2
               V(I,locNY+1,K)=V(I,locNY,K)
               V(I,locNY+2,K)=V(I,locNY,K)
            ENDDO
         ENDDO 
         DO K=0,locNZ+2
            DO I=-1,locNX+2
               H(I,locNY+1,K)=H(I,locNY,K)
               H(I,locNY+2,K)=H(I,locNY,K)
            ENDDO
         ENDDO
      case(1) !-- Keplarian
         K=1
         DO i=1,locNX
            RH(i,locNY+1,K)=RH(i,locNY,K)
            RH(i,locNY+2,K)=RH(i,locNY,K)
            T(i,locNY+1,K)=T(i,locNY,K)
            T(i,locNY+2,K)=T(i,locNY,K)
            G(I,locNY+1,K) = ((  (1.d0/(RHQY(i,locNY+1,k)*xxb(i)))*
     %           ((PG(i,locNY,K)-PG(i-1,locNY,K))/(XXB(i)-XXB(i-1))) + 
     %           ((OMEGAKEP(I)+OMROT)**2.d0))**0.5d0) - OMROT
            G(i,locNY+2,k) = G(i,locNY+1,k)
            V(i,locNY,k)  =0.d0
            V(i,locNY+1,k)=0.d0
            V(I,locNY+2,K)=0.d0
         ENDDO
!--   lowerleft processor sets i=-1:0
         if(MPIleft.eq.MPI_PROC_NULL) then
            DO i=-1,0
               RH(i,locNY+1,K)=RH(i,locNY,K)
               RH(i,locNY+2,K)=RH(i,locNY,K)
               T(i,locNY+1,K)=T(i,locNY,K)
               T(i,locNY+2,K)=T(i,locNY,K)
               G(i,locNY+1,k) = OMEGAKEP(i)
               G(i,locNY+2,k) = OMEGAKEP(i)
               if(i.ne.-1) then
                  V(i,locNY,k)  =0.d0
                  V(i,locNY+1,k)=0.d0
                  V(I,locNY+2,K)=0.d0
               endif
            ENDDO
         endif
!--   lowerright processor sets i=locNX+1:locNX+2
         if(MPIright.eq.MPI_PROC_NULL) then
            DO i=locNX+1,locNX+2
               RH(i,locNY+1,K)=RH(i,locNY,K)
               RH(i,locNY+2,K)=RH(i,locNY,K)
               T(i,locNY+1,K)=T(i,locNY,K)
               T(i,locNY+2,K)=T(i,locNY,K)
               G(i,locNY+1,k) = OMEGAKEP(i)
               G(i,locNY+2,k) = OMEGAKEP(i)
               V(i,locNY,k)  =0.d0
               V(i,locNY+1,k)=0.d0
               V(I,locNY+2,K)=0.d0
            ENDDO
         endif
      case(2) !-- axi-symmetric
         DO K=-1,locNZ+2
            DO i=-1,locNX+2
               RH(i,locNY+1,K)=RH(i,locNY,K)
               RH(i,locNY+2,K)=RH(i,locNY,K)
               T(i,locNY+1,K)=T(i,locNY,K)
               T(i,locNY+2,K)=T(i,locNY,K)
               G(i,locNY+1,K)=G(i,locNY,K)
               G(i,locNY+2,K)=G(i,locNY,K)
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO i=0,locNX+2
               V(i,locNY+1,K)=V(i,locNY,K)
               V(i,locNY+2,K)=V(i,locNY,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO i=-1,locNX+2
               H(i,locNY+1,K)=H(i,locNY,K)
               H(i,locNY+2,K)=H(i,locNY,K)
            ENDDO
         ENDDO
      case(3:)
         print *,'type=',type,'not yet defined: BOUNDYMAX'
         stop
      CASE DEFAULT
         print *,'type=',type,'also not defined: BOUNDYMAX'
         stop
      END SELECT
      RETURN
      END SUBROUTINE BOUNDYMAX

      SUBROUTINE BOUNDYMIN(type)
! type = 0: closed 
! type = 1: keplarian
      USE global_constants
      USE fluid_var_init
      USE input_init
      USE global_constants
      USE grid_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: type,i,j,k
      SELECT CASE(type)
      CASE(0) !-- closed
         DO K=-1,locNZ+2
            DO I=-1,locNX+2
               RH(I,0,K)=RH(I,1,K) 
               RH(I,-1,K)=RH(I,1,K) 
               T(I,0,K)=T(I,1,K)
               T(I,-1,K)=T(I,1,K)
               G(I,0,K)=-G(I,2,K)
               G(I,1,K)=0.d0
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO I=0,locNX+2
               V(I,0,K)  = V(I,1,K)
               V(I,-1,K) = V(I,1,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO I=-1,locNX+2
               H(I,0,K)  = H(I,1,K)
               H(I,-1,K) = H(I,1,K)
            ENDDO
         ENDDO
      CASE(1) !-- 2D Keplarian
         K=1
         DO I=1,locNX
            RH(I,0,K)=RH(I,1,K) 
            RH(I,-1,K)=RH(I,1,K) 
            T(I,0,K)=T(I,1,K)
            T(I,-1,K)=T(I,1,K)
            G(I,1,K) = ((  (1.d0/(RHQY(i,1,k)*xxb(i)))*
     %           ((PG(i,1,K)-PG(i-1,1,K))/(XXB(i)-XXB(i-1))) + 
     %           ((OMEGAKEP(I)+OMROT)**2.d0))**0.5d0) - OMROT
            G(i,0,k) = G(i,1,k)
            V(i,1,k)=0.0
            V(i,0,k)=0.0
            V(I,-1,K) = 0.d0
         ENDDO         
!-- lowerleft processor sets i=-1:0
         if(MPIleft.eq.MPI_PROC_NULL) then
            DO i=-1,0
               RH(I,0,K)=RH(I,1,K) 
               RH(I,-1,K)=RH(I,1,K) 
               T(I,0,K)=T(I,1,K)
               T(I,-1,K)=T(I,1,K)
               G(i,1,k) = OMEGAKEP(i)
               G(i,0,k) = OMEGAKEP(i)
               if(i.ne.-1) then
                  V(i,1,k)=0.0
                  V(i,0,k)=0.0
                  V(I,-1,K) = 0.d0
               endif
            ENDDO
         endif
!-- lowerright processor sets i=-1:0
         if(MPIright.eq.MPI_PROC_NULL) then
            DO i=locNX+1,locNX+2
               RH(I,0,K)=RH(I,1,K) 
               RH(I,-1,K)=RH(I,1,K) 
               T(I,0,K)=T(I,1,K)
               T(I,-1,K)=T(I,1,K)
               G(i,1,k) = OMEGAKEP(i)
               G(i,0,k) = OMEGAKEP(i)
               V(i,1,k)=0.0
               V(i,0,k)=0.0
               V(I,-1,K) = 0.d0
            ENDDO
         endif
      CASE(2) !-- axi-symmetric
         DO K=-1,locNZ+2
            DO i=-1,locNX+2
               RH(i,0,K) = RH(i,1,K)
               RH(i,-1,K)= RH(i,1,K)
               T(i,0,K)  = T(i,1,K)
               T(i,-1,K) = T(i,1,K)
               G(i,0,K)  = G(i,1,K)
            ENDDO
         ENDDO
         DO K=-1,locNZ+2
            DO i=0,locNX+2
               V(i,0,K)  = V(i,1,K)
               V(i,-1,K) = V(i,1,K)
            ENDDO
         ENDDO
         DO K=0,locNZ+2
            DO i=-1,locNX+2
               H(i,0,K)  = H(i,1,K)
               H(i,-1,K) = H(i,1,K)
            ENDDO
         ENDDO
      CASE(3:)
         print *,'type=',type,'not yet defined: BOUNDYMIN'
         stop
      CASE DEFAULT
         print *,'type=',type,'also not defined: BOUNDYMIN'
         stop
      END SELECT    
      RETURN
      END SUBROUTINE BOUNDYMIN

      SUBROUTINE BOUNDZMIN(type)
! type = 0: closed 
! type = 1: disk
      USE global_constants
      USE fluid_var_init
      USE input_init
      USE grid_var_init
      IMPLICIT NONE
      integer :: type,i,j,k
      SELECT CASE(type)
      CASE(0) !-- closed
         DO I=-1,locNX+2
            DO J=-1,locNY+2
!- ADDED 6/4/08
               RH(I,j,1)=RH(I,j,2)

               RH(I,j,0)=RH(I,j,1)
               RH(I,j,-1)=RH(I,j,1)

!- ADDED 6/3/08
               T(I,j,1)=T(I,j,2)

               T(I,j,0)=T(I,j,1)
               T(I,j,-1)=T(I,j,1)

!- ADDED 6/19/08
               H(I,j,2)=0.d0

               H(I,j,1)=0.d0
               H(I,j,0)=0.d0
            ENDDO
         ENDDO
!--- d(vr)/dz=0
         DO I=0,locNX+2
            DO J=-1,locNY+2
!- ADDED 6/3/08
               V(I,j,1)=0.d0

               V(I,j,0)=V(I,j,1)
               V(I,j,-1)=V(I,j,1)
            ENDDO
         ENDDO
!--- d(vphi)/dz=d(Grcos(theta))/dz=0
         DO I=-1,locNX+2
            DO J=0,locNY+2
c               G(I,j,0)=G(I,j,1)
c               G(I,j,-1)=G(I,j,1)
c-- d(u_phi)/dtheta = d(G*rcos(theta))/dtheta = 0

!- ADDED 6/19/08
               G(I,j,1)=G(I,j,2)*surzb(2)/surzb(1)

               G(I,j,0)=G(I,j,1)*surzb(1)/surzb(0)
               G(I,j,-1)=G(I,j,1)*surzb(1)/cos(xzb(-1))
            ENDDO
         ENDDO
      CASE(1)
         DO J=-1,locNY+2
            DO I=-1,locNX+2
               RH(I,j,0)=RH(I,j,1)
               RH(I,j,-1)=RH(I,j,1)
               T(I,j,0)=T(I,j,1)
               T(I,j,-1)=T(I,j,1)
               H(I,j,1)=0.d0
               H(I,j,0)=0.d0
            ENDDO
         ENDDO
         DO J=-1,locNY+2
            DO I=0,locNX+2
               V(I,j,0)=V(I,j,1)
               V(I,j,-1)=V(I,j,1)
            ENDDO
         ENDDO
         DO J=0,locNY+2
            DO I=-1,locNX+2
*     no-torque, slip condition
               G(I,j,0)=G(I,j,1)
               G(I,j,-1)=G(I,j,1)
            ENDDO
         ENDDO
      CASE(2) !-- RH, T, V, and G continous; H=0
         DO J=-1,locNY+2
            DO I=-1,locNX+2
               RH(I,j,0)  = RH(i,j,1)
               RH(I,j,-1) = RH(i,j,1)
               T(I,j,0)   = T(i,j,1)
               T(I,j,-1)  = T(i,j,1)
               H(I,j,1)   = 0.d0
               H(I,j,0)   = 0.d0
            ENDDO
         ENDDO
         DO J=-1,locNY+2
            DO I=0,locNX+2
               V(I,j,0)  = V(I,j,1)
               V(I,j,-1) = V(I,j,1)
            ENDDO
         ENDDO
         DO J=0,locNY+2
            DO I=-1,locNX+2
               G(I,j,0)  = G(I,j,1)
               G(I,j,-1) = G(I,j,1)
            ENDDO
         ENDDO
      CASE(3:)
         print *,'type=',type,'not yet defined: BOUNDZMIN'
         stop
      CASE DEFAULT
         print *,'type=',type,'also not defined: BOUNDZMIN' 
        stop
      END SELECT    
      RETURN
      END SUBROUTINE BOUNDZMIN

      SUBROUTINE BOUNDZMAX(type)
! type = 0: closed 
! type = 1: disk
      USE global_constants
      USE fluid_var_init
      USE input_init
      USE grid_var_init
      IMPLICIT NONE
      integer :: type,i,j,k
      SELECT CASE(type)
      CASE(0) !-- closed
         DO I=-1,locNX+2
            DO J=-1,locNY+2
!- ADDED 6/4/08
               RH(I,j,locNZ)=RH(I,j,locNZ-1)
               RH(I,j,locNZ+1)=RH(I,j,locNZ)
               RH(I,j,locNZ+2)=RH(I,j,locNZ)

!- ADDED 6/3/08
               T(i,j,locNZ)=T(i,j,locNZ-1)

               T(i,j,locNZ+1)=T(i,j,locNZ)
               T(i,j,locNZ+2)=T(i,j,locNZ)

!- ADDED 6/3/08
               H(I,j,locNZ)=0.d0

               H(I,j,locNZ+1)=0.d0
               H(I,j,locNZ+2)=0.d0
            ENDDO
         ENDDO
         DO I=0,locNX+2
            DO J=-1,locNY+2
!- ADDED 6/3/08
               V(I,j,locNZ)=0.d0

               V(I,j,locNZ+1)=V(I,j,locNZ)
               V(I,j,locNZ+2)=V(I,j,locNZ)
            ENDDO
         ENDDO      
         DO I=-1,locNX+2
            DO J=0,locNY+2
c               G(I,j,locNZ+1)=G(I,j,locNZ)
c               G(I,j,locNZ+2)=G(I,j,locNZ)
c-- d(u_phi)/dtheta = d(G*rcos(theta))/dtheta = 0
!- ADDED 6/19/08
               G(I,j,locNZ)=G(I,j,locNZ-1)*surzb(locNZ-1)/surzb(locNZ)


               G(I,j,locNZ+1)=G(I,j,locNZ)*surzb(locNZ)/surzb(locNZ+1)
               G(I,j,locNZ+2)=G(I,j,locNZ)*surzb(locNZ)/surzb(locNZ+2)
            ENDDO
         ENDDO
      CASE(1)
         DO J=-1,locNY+2
            DO I=-1,locNX+2
               RH(I,j,locNZ+1)=RH(I,j,locNZ)
               RH(I,j,locNZ+2)=RH(I,j,locNZ)
               T(I,j,locNZ+1)=T(I,j,locNZ)
               T(I,j,locNZ+2)=T(I,j,locNZ)
               H(I,j,locNZ+1)=0.d0
               H(I,j,locNZ+2)=0.d0
            ENDDO
         ENDDO
         DO J=-1,locNY+2
            DO I=0,locNX+2
               V(I,j,locNZ+1)=V(I,j,locNZ)
               V(I,j,locNZ+2)=V(I,j,locNZ)
            ENDDO
         ENDDO
         DO J=0,locNY+2
            DO I=-1,locNX+2
*     no-torque, slip condition
               G(I,j,locNZ+1)=G(I,j,locNZ)
               G(I,j,locNZ+2)=G(I,j,locNZ)
            ENDDO
         ENDDO
      CASE (2) !-- T=T_NEB; RH=RH_NEB; V=0; G Cont.; H=outflow
         DO J=-1,locNY+2
            DO I=-1,locNX+2
               RH(I,j,locNZ)  = RH_nebula
               RH(I,j,locNZ+1)= RH_nebula
               RH(I,j,locNZ+2)= RH_nebula
               T(I,j,locNZ)   = T_nebula
               T(I,j,locNZ+1) = T_nebula
               T(I,j,locNZ+2) = T_nebula
               H(I,j,locNZ+1) = max(H(i,k,locNZ),0.d0)
               H(I,j,locNZ+2) = max(H(i,k,locNZ),0.d0)
            ENDDO
         ENDDO
         DO J=-1,locNY+2
            DO I=0,locNX+2
               V(I,j,locNZ+1) = 0.d0
               V(I,j,locNZ+2) = 0.d0
            ENDDO
         ENDDO
         DO J=0,locNY+2
            DO I=-1,locNX+2
*     no-torque, slip condition
               G(I,j,locNZ+1) = G(I,j,locNZ)
               G(I,j,locNZ+2) = G(I,j,locNZ)
            ENDDO
         ENDDO
      CASE(3:)
         print *,'type=',type,'not yet defined: BOUNDZMAX'
         stop
      CASE DEFAULT
         print *,'type=',type,'also not defined: BOUNDZMAX'
         stop
      END SELECT    
      RETURN
      END SUBROUTINE BOUNDZMAX

      SUBROUTINE RADIATIVEBOUND
**    Boundary Conditions the Radiation Energy (THIS ROUTINE IS ONLY USED WITH IRAD=1)
!
!-- NOTE: CHANGES TO RADBOUNDVER=2-4 NEED TO BE MADE TO BOUNDXMAX(2) TOO
!
!--RADBOUNDVER=0 : not defined
!--RADBOUNDVER=1 : fixed bottom dt/dr, periodic sides, fixed top
!--RADBOUNDVER=2-7: fixed bottom T flux (set dTdr and T), open sides
!--RADBOUNDVER=2 : movable/sinsodial (in y) at xmax
!--RADBOUNDVER=3 : movable/sinsodial (in y/z) at xmax, TMIN = 100
!--RADBOUNDVER=4 : movable/constant at xmax
!--RADBOUNDVER=5 : movable/sinsodial (in y(time)/z) at xmax
!--RADBOUNDVER=6 : movable/sinsodial (in t/y/z) at xmax
!--RADBOUNDVER=7 : fixed top dFdr = (stellar heating) (SPHERICALLY SYMETRIC!!) at xmax,open sides
!--RADBOUNDVER=8 : movable/sinsodial (in y/z) at xmax, TMIN for BD
!--RADBOUNDVER=9 : T=(ER/a+1/4*(gam*Tir)^4*exp(-tauA))^(0.25)
      USE input_init
      USE rad_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      INTEGER :: RADBOUNDVER
      double precision :: TJupHeat!,phi0
      double precision :: Tincident,irr_ramp,CALC_TINCIDENT,dERdr
!--- OVERALL BOUNDARY CONDITION.
      RADBOUNDVER = IBDRAD
      if(periods(1)) then
         print *, 'IS periodic x avalible?: radbounds'
         stop
      endif
      SELECT CASE(RADBOUNDVER)
      CASE (0)
         print *,'radboundver',radboundver,'not defined: RADIATIVEBOUND' 
         stop
      CASE (1) !-- fixed bottom, open top and sides
         if(MPIleft.eq.MPI_PROC_NULL) then !** X-BOUNDARY (bottom): fixed dT/dr
            DO K=-1,locNZ+2
               DO J=-1,locNY+2
                  dERdr = 4.d0*((arad*(ER(1,j,k)**3))**(0.25))*D_T
                  ER(-1,J,K)= ER(1,J,K)+dERdr*(xxb(0)-xxb(1))
                  ER(0,J,K) = ER(1,J,K)+dERdr*(xxb(-1)-xxb(1))
               ENDDO
            ENDDO
         endif
         if(MPIright.eq.MPI_PROC_NULL) then !** X-BOUNDARY (top)
            DO K=-1,locNZ+2
               DO J=-1,locNY+2
                  ER(locNX,J,K)   = arad*Tinittop**4
                  ER(locNX+1,J,K) = arad*Tinittop**4
                  ER(locNX+2,J,K) = arad*Tinittop**4
               ENDDO
            ENDDO
         endif
         if(MPIlower.eq.MPI_PROC_NULL.or.MPIupper.eq.MPI_PROC_NULL) then !** Y-BOUNDARY
            print *,'should not be here: radiativebound'
            print *,'MPIlower/upper=',MPIlower,MPIupper,myid
            call clean_stop
         endif
         if(MPIlower.eq.MPI_PROC_NULL.or.MPIupper.eq.MPI_PROC_NULL) then !** Z-BOUNDARY
            print *,'should not be here: radiativebound'
            print *,'MPIabove/below=',MPIabove,MPIbelow,myid
            call clean_stop
         endif
      CASE (2:9) !-- fixed T derivative, movable/sinsodial top and open sides
!----     2 = sinsodial in y
!----     3 = sinsodial in y and z, TMIN = 100
!----     4 = constant at xmax
!----     5 = sinsodial in y(time) and z
!----     6 = temporal pulsing, and sinsodial in y and z
!----     7 = fixed dTdr=0 at top
!----     8 = sinsodial in y and z, TMIN = 950
!-------- xmin condition -- fix the temperature and density gradients
c         if(modver.le.2) then
         if(modver.le.5) then
            DO K=-1,locNZ+2
               DO J=-1,locNY+2
                  if(BCXMIN.eq.6) then
                     T(1,J,K)   = TinitBottom !--fix the entire flux
                  endif
                  T(0,j,K)   = T(1,J,K) -  D_T*(xxb(1)-xxb(0))
                  T(-1,j,K)  = T(1,J,K) -  D_T*(xxb(1)-xxb(-1))
                  ER(1,J,K)  = arad*T(1,j,k)**4
                  ER(0,J,K)  = arad*T(0,j,k)**4
                  ER(-1,J,K) = arad*T(-1,j,k)**4
               ENDDO
            ENDDO
         else
            print *,'modver=',modver,'not avalible in radiativbnd'
            stop
         endif
!-------- xmax condition
         first_Tincident = .true.
         DO K=1,locNZ
            DO J=1,locNY
               IF(upperbnd(j,k).lt.locNX+1) then
                  IF(RADBOUNDVER.eq.2) then
                     TJupHeat = Tirr*cos(xyb(j)/ymax)
                     TJupHeat = MAX(TJupHeat,100.d0)
                  ELSEIF(RADBOUNDVER.eq.3) then
                     if(ncosys.eq.1) then
                        print *,'N/A...'
                        stop
                     elseif(ncosys.eq.2) then
                        TJupHeat = Tirr*((cos(xyb(j))*
     %                       cos(xzb(k)))**0.25)
                     else
                        print *,'bad ncosys:RADIATIVEBND',ncosys
                        stop
                     endif
                     TJupHeat = MAX(TJupHeat,100.d0)
c                     TJupHeat = MAX(TJupHeat,500.d0)
                  ELSEIF(RADBOUNDVER.eq.4) then
                     TJupHeat = Tirr
                  ELSEIF(RADBOUNDVER.eq.5) then ! thermal tide
                     if(ncosys.eq.2) then
c                        phi0 = (2.d0*PI/(thermalP_phi*DAY))*
c     %                       mod((zeit-thermal_t0*DAY),
c     %                       (thermalP_phi*DAY))
c                        TJupHeat = Tirr*((cos(xyb(j)-phi0)*
c     %                       cos(xzb(k)))**0.25)
                        print *,'redo thermal tide with Porb'
                        call clean_stop
                     else
                        print *,'bad ncosys:boundxmax (5)',ncosys
                        stop
                     endif
                     TJupHeat = MAX(TJupHeat,100.d0)
                  ELSEIF(RADBOUNDVER.eq.6) then ! temporal pulsations
c                     TJupHeat = Tirr*((cos(xyb(j))* !--standard spatial dependance
c     %                    cos(xzb(k)))**0.25)
c                     TJupHeat = TJupHeat*cos(2.d0*PI*(zeit-thermal_t0)/
c     %                    (0.1*DAY))  !-- temporal pulsations
c                     TJupHeat = MAX(TJupHeat,100.d0)
                     print *,'redo temporal tide '
                     call clean_stop
                  ELSEIF(RADBOUNDVER.eq.7) then ! gradT = 0 boundary condition
                     TJupHeat = T(upperbnd(j,k)-1,j,k)
                  ELSEIF(RADBOUNDVER.eq.8) then ! SAME AS 3, BUT DIFFERENCT TMIN FOR BD
                     TJupHeat = Tirr*((cos(xyb(j))*
     %                    cos(xzb(k)))**0.25)
                     TJupHeat = MAX(TJupHeat,950.d0)
                  ELSEIF(RADBOUNDVER.eq.9) then !--T=(ER/a+1/4*(gam*Tir)^4*exp(-tauA))^(0.25)
                     IF (III.le.1) THEN
                        TJupHeat = T(upperbnd(j,k)-1,j,k) !-just the first time, because the opacities aren't defined
                     else
c                        TJupHeat = ( (ER(upperbnd(j,k)-1,J,K)/ARAD) +
                        TJupHeat = ( (ER(upperbnd(j,k),J,K)/ARAD) +
     $                       (stellarinput(upperbnd(j,k)-1,j,k)/
     $                       (4.d0*SBCONST*rh(upperbnd(j,k)-1,j,k)*
     $                       xkapP(upperbnd(j,k)-1,j,k))) )**(0.25)
c-old grey method
c                        Tincident = CALC_TINCIDENT(j,k)
c                        if(F_INCIDENT.eq.3.or.F_INCIDENT.eq.5) then !--- RAMP-UP
c                           irr_ramp = (1.d0*III/(1.d0*DAMP_STEPS))
c                           if(irr_ramp.ge.1.0) irr_ramp = 1.d0
c                        else
c                           irr_ramp = 1.0
c                        endif
c                        TJupHeat = ( (ER(upperbnd(j,k)-1,J,K)/ARAD) +
c     %                       irr_ramp*
c     $                       0.25*(dtaus(upperbnd(j,k)-1,j,k)/
c     $                       (rh(upperbnd(j,k)-1,j,k)*
c     $                       xkapP(upperbnd(j,k)-1,j,k)*
c     $                       dxb(upperbnd(j,k)-1)))*
c     $                       ((Tincident**4.d0)/mu_star(j,k))*
c     $                       exp(-stellartau(upperbnd(j,k),j,k)/
c     $                       mu_star(j,k)))**(0.25)
                     endif
                  ENDIF
                  DO i=upperbnd(j,k),locNX+2
                     ER(I,J,K)=arad*TJupHeat**4.d0
                  ENDDO
               ELSE
                  print *,'error in radiativebnd',upperbnd(j,k)
                  stop
               ENDIF
            ENDDO
         ENDDO
!----------ymin condition
         if(MPIlower.eq.MPI_PROC_NULL) then
            print *,'didnt'
            stop
            DO K=-1,locNZ+2
               DO I=-1,locNX+2
                  ER(I,-1,K)=ER(I,1,K)
                  ER(I,0,K)=ER(I,1,K)
               ENDDO
            ENDDO
         endif
!---------- ymax condition
         if(MPIupper.eq.MPI_PROC_NULL) then
            DO K=-1,locNZ+2
               DO I=-1,locNX+2
                  ER(I,locNY+1,K)=ER(I,locNY,K)
                  ER(I,locNY+2,K)=ER(I,locNY,K)
               ENDDO
            ENDDO
         endif
!--------- zmin condition
         if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
            DO J=-1,locNY+2
               DO I=-1,locNX+2
                  ER(I,J,-1)=ER(I,J,1)
                  ER(I,J,0)=ER(I,J,1)
               ENDDO
            ENDDO
         endif
!-------- zmax condition
         if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL) then
            DO J=-1,locNY+2
               DO I=-1,locNX+2
                  ER(I,J,locNZ+1)=ER(I,J,locNZ)
                  ER(I,J,locNZ+2)=ER(I,J,locNZ)
               ENDDO
            ENDDO
         endif
      CASE (10:)
         print *,'radboundver',radboundver,'not defined: RADIATIVEBOUND' 
         stop
      CASE DEFAULT
         print *,'radboundver',radboundver,'not defined:RADIATIVEBOUND'
         stop
      END SELECT
!---PASS ER DATA BETWEEN PROCESSORS
      if(modtyp.eq.2) then !-- shifting routine for planet mpi grid setup
         CALL PLANETSHIFT(6,-1)
!--   now correct the polar shifts if necessary
         if(poles) then
            if(proc_coords(3).eq.proc_dims(3)-1) then
               CALL NORTH_POLE_SHIFT(6,-1)
            endif
            if(proc_coords(3).eq.0) then
               CALL SOUTH_POLE_SHIFT(6,-1)
            endif
         endif
      else  !-- shifting routine for everything else
         IF(locNZ.eq.1) THEN
            CALL SHIFTVAR2D(7)
         ENDIF 
         IF(locNZ.ne.1) THEN
            CALL NEWSHIFTVAR3D(6)
c            CALL SHIFTVAR3D(6)
         ENDIF
      endif
      RETURN
      END SUBROUTINE RADIATIVEBOUND

     
      SUBROUTINE PLANETMASSSINK     
**    Planet Boundary Mass Sink
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: cosb,Rtwo,rplanet,phiplanet,om
      double precision :: rsink,tsink,rhill,locMacc,global_Macc,area
      rplanet = rp*AU
      phiplanet = phip*PI
      rhill = rplanet*((mstar2/3)**0.33)
      rsink = 0.5*rhill
      tsink = 1.d0
      locMacc = 0.d0
      k=1
      DO J=1,locNY
         cosb=rplanet*cos(XYB(J)-phiplanet)
         DO I=1,locNX
            Rtwo=dsqrt(rplanet*rplanet + XXA(I)*XXA(I)
     &           -2.d0*XXA(I)*cosb)
            if(Rtwo.le.rsink) then
               om = (GRAV*MSOL/(xxb(i)**3.d0))**0.5
               area = (XXA(I+1)-XXA(I))*XXB(I)*(XYA(J+1)-XYA(J))
               locMacc = locMacc + (RH(i,j,k)*delt*om/tsink)*area
               RH(i,j,k) = RH(i,j,k)*(1.d0-delt*om/tsink)
            endif
         ENDDO
      ENDDO
      call MPI_ALLREDUCE(locMacc,global_Macc,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      Maccreated = Maccreated + global_Macc
      RETURN
      END SUBROUTINE PLANETMASSSINK


      SUBROUTINE CALCGAPMASS
**    TOTAL GAS IN GAP REGION
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: rplanet,seperatrix
      double precision :: rhill,locMgap,area,global_Mgap
      rplanet = rp*AU
      rhill = rplanet*((mstar2/3)**0.33)
      seperatrix = dsqrt(12.d0)*rhill
      k=1
      locMgap = 0.d0
      DO J=1,locNY
         DO I=1,locNX
            if(abs(XXA(I)-rplanet).le.seperatrix) then
               area = (XXA(I+1)-XXA(I))*XXB(I)*(XYA(J+1)-XYA(J))
               locMgap = locMgap + RH(i,j,k)*area
            endif
         ENDDO
      ENDDO
      call MPI_ALLREDUCE(locMgap,global_Mgap,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      Mgap = global_Mgap
      RETURN
      END SUBROUTINE CALCGAPMASS


!     Calculate a universal derivative (dT/dr or dQdr) at the bottom by
!     using average rh,T, and cv. Call with: 
!     CALC_UNIVERSAL_BOTTOM_DERIV(.false.) for dT/dr
!     CALC_UNIVERSAL_BOTTOM_DERIV(.true) for dQ/dr
      FUNCTION CALC_UNIVERSAL_BOTTOM_DERIV(calc_DQDR)
      USE global_constants
      USE fluid_var_init
      USE input_init
      USE grid_var_init
      USE mpi_var_init
      USE rad_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      double precision :: calc_universal_bottom_deriv,bottomDTDR
      logical :: calc_DQDR
      double precision :: local_avg_RH,avg_RH
      double precision :: local_avg_cv,avg_cv
      double precision :: local_avg_T,avg_T
      double precision :: local_avg_difrx,avg_difrx
      double precision :: local_avg_kapR,avg_kapR
!-calculate local averages
      local_avg_RH = sum(RH(1,1:locNY,1:locNZ))/(1.d0*locNY*locNZ)
      local_avg_cv = sum(cv(1,1:locNY,1:locNZ))/(1.d0*locNY*locNZ)
      local_avg_T = sum(T(1,1:locNY,1:locNZ))/(1.d0*locNY*locNZ)
      if(III.eq.0) then !-for first calc of kappa isn't defined.. do it manually
         call BurrowsOPC(T(1,1,1),RH(1,1,1),xkapP(1,1,1),
     $        local_avg_kapR,xkapA(1,1,1),xkapW(1,1,1,:),
     $        xkapLW(1,1,1,:),local_BB(1,1,1,:),
     $        dxkapP_dT(1,1,1),dxkapA_dT(1,1,1),1,1,1)
      else
         local_avg_kapR = sum(xkapR(1,1:locNY,1:locNZ))/
     $        (1.d0*locNY*locNZ)
      endif
!-sum all local averages
      call MPI_ALLREDUCE(local_avg_RH,avg_RH,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(local_avg_cv,avg_cv,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(local_avg_T,avg_T,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(local_avg_kapR,avg_kapR,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!-find overall averages
      avg_RH=avg_RH/(1.d0*numprocs)
      avg_cv=avg_cv/(1.d0*numprocs)
      avg_T=avg_T/(1.d0*numprocs)
      avg_kapR=avg_kapR/(1.d0*numprocs)
!-put it all together for dT/dr via radiative gradient
      bottomDTDR =-3.d0*Fbottom*avg_kapR*avg_RH/
     $     (4.d0*ARAD*FD*(avg_T**3.d0))
      if(ieta.eq.5) then        !-fixed PR and viscosity at the bottom
         bottomDTDR = -Fbottom/(avg_RH*avg_cv*(XNUE/PR_NUM))
      endif
!- augment by rh*cv for dQ/dr if necessary
      if(calc_dQdr) then
         calc_universal_bottom_deriv = avg_RH*avg_CV*bottomDTDR
      else 
         calc_universal_bottom_deriv = bottomDTDR
      endif
      RETURN
      END FUNCTION CALC_UNIVERSAL_BOTTOM_DERIV
