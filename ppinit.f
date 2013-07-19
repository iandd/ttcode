      SUBROUTINE localinit
      USE fluid_var_init
      USE force_var_init
      USE global_var_init
      USE grid_var_init
      USE input_init
      USE rad_var_init
      USE global_constants
      USE mpi_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      INTEGER :: i,j,k,f,loci,locj,lock
      double precision :: period_start,npoly
      double precision :: xkapW_bottom(nwave_bins)
      double precision :: xkapLW_bottom(nwave_bins)
      double precision :: local_BB_bottom(nwave_bins)
      double precision :: xkapP_bottom,xkapA_bottom
      double precision :: dxkapP_dT_bottom,dxkapA_dT_bottom
      double precision :: SETPLANETXGRAV,tmp_gravx
      double precision :: Fbot_test,kappaT_bottom
      double precision :: HP
!---  Initilize local copies of the variables
      do k=1,locNZ
         lock = k+locNZ*proc_coords(3)
         if(lock.gt.global_NZ.or.lock.lt.1) then
            print *,'ERROR IN INIT-k'
            call clean_stop
         endif
         do j=1,locNY
            locj = j+locNY*proc_coords(2)
            if(locj.gt.global_NY.or.locj.lt.1) then
               print *,'ERROR IN INIT-j'
               call clean_stop
            endif
            do i=1,locNX
               loci = i+locNX*proc_coords(1)
               if(loci.gt.global_NX.or.loci.lt.1) then
                  print *,'ERROR IN INIT-i'
                  call clean_stop
               endif
               T(i,j,k)  = global_T(loci,locj,lock)
               ER(i,j,k) = global_ER(loci,locj,lock)
               RH(i,j,k) = global_RH(loci,locj,lock)
               V(i,j,k)  = global_V(loci,locj,lock)
               G(i,j,k)  = global_G(loci,locj,lock)
               H(i,j,k)  = global_H(loci,locj,lock)
            enddo
         enddo
      enddo
!---set Flux through the bottom (via ipp value)
      FBottom = Fbot
!---set all of H to zero for two-dimensional runs
      if(locNZ.eq.1) then
         H = 0.d0
      endif
!-- un-movable upper boundary for the overflowing BD model
      if(MODTYP.EQ.2.and.MODVER.EQ.1) THEN
         upperbnd = locNX
         if(myid.eq.0) print *,'initupperbnd',upperbnd(2,2)
      endif
!-- un-movable upper boundary for non irradiation models
      if(MODTYP.NE.2) THEN
         upperbnd = locNX
      ENDIF
!-- call boundary routine
      CALL BOUNDS
      if(modtyp.eq.2) then
         CALL PLANETSHIFT(6,-1) !---Shift ER for the initialization
!--   now correct the polar shifts if necessary
         if(poles) then
            if(proc_coords(3).eq.proc_dims(3)-1) 
     $           CALL NORTH_POLE_SHIFT(6,-1)
            if(proc_coords(3).eq.0) CALL SOUTH_POLE_SHIFT(6,-1)
         endif
      else
         CALL NEWSHIFTVAR3D(6)
      endif
!--FOR THESE BCXMIN, RH IS NOT SET IN BOUNDXMIN. SET HERE FOR Q-CALCULATION
c      if(bcxmin.eq.0.or.bcxmin.eq.8.or.bcxmin.eq.9.or.bcxmin.eq.12) then
c         DO K=-1,locNZ+2
c            DO J=-1,locNY+2
c               RH(0,j,K)  = RH(1,j,K)
c               RH(-1,j,K) = RH(1,j,K)
c            ENDDO
c         ENDDO
c      endif
      CALL BOUNDS
!-- calculate the diffusion constants for boundary conditions (especially bcxmin=8/9)
      if(irad.eq.1.or.irad.eq.2.or.irad.eq.3.or.irad.eq.4.or.irad.eq.6) 
     $     CALL DIFF
!--- calculate the stellar optical depth for initial boundary calls
!    the heating functions call the appropriate stellar_tau function
      if(irad.gt.0.and.ikaptyp.ne.0) then
         CALL KAPPA
         if(F_INCIDENT.eq.1) then
            if(FDELT.gt.0.0) then
               print *,'not sure you can have'
               print *,'fdelt=0 and explicit heating'
               print *,'stopping in ppinit'
               call clean_stop
            endif
            CALL STELLARHEATING
         else
            CALL STELLARINPUT_HEATING 
         endif
!--- recalculate bounds with correct opacities
         CALL BOUNDS
      endif

!---  set the initial boundary temperatures and densities
!---  only do this for initial runs. Otherwise it is read in from AKTTD
      IF (MOD(IREA,10).EQ.0) THEN
         Tinittop     = global_T(global_nx,1,1)
         rhinittop    = global_RH(global_NX,1,1)
         TinitBottom  = global_T(1,1,1)
         rhinitbottom = global_RH(1,1,1)
      ENDIF
      if(myid.eq.0) then
         print *,'LOCALINIT SET:'
         print *,' Bottom D_RH, D_T =',D_RH,D_T
         print *,' Tinitbottom/top',TinitBottom,Tinittop
         print *,' rhobottom/top',rhinitbottom,rhinittop
         write(*,'(1x,A,1x,1(1x,e12.6))') 'Fbottom',
     $        Fbottom
         print *,'---------------------------'
      endif

!-- set the bottom density and temperature derivatives
      if(modtyp.eq.1) then !cartesian box model (changes need to be made to initcartesian too)
         SELECT CASE(MODVER)
         CASE (0) ! simple convecting box
            call fgrav
            npoly = 1.d0
            D_T = -ABS(gravx(1,1,1))*mu_gas/(RGAS*(npoly+1.))
            D_RH = 0.d0
         CASE (1)               !-- SOD Shockwave test
            D_T = 0.d0
            D_RH = 0.d0
         CASE (2) !--- conv. zone / isothermal earth model
            D_T = (global_T(2,1,1)-global_T(1,1,1))/(xxb(2)-xxb(1))
            D_RH = 0.d0
         CASE (3:)
            print *,'MODVER not defined in ppinit, cartesian'
         CASE DEFAULT
            print *,'MODVER also not defined in ppinit, cartesian'
         END SELECT

      elseif(modtyp.eq.2) then
!-- values are from peters/adams input models
         SELECT CASE(MODVER)
         CASE (0)
!---  these are for 60_cells_peter.dat
            D_RH = -2.1870361*(10.**(-12.0))
            D_T = -1.3611882*(10.**(-6.0))
         CASE (1) !--- these are for 100_1808_poly.dat
            D_RH  =  -2.7934540*(10.**(-9.0))
            D_T   = -0.00097568917
         CASE (2) !--- these are for 60_cells_deep.dat
            D_RH  =  -2.3748328*(10.**(-11.0))
            D_T = -4.0182995*(10.**(-6.0))
            CALL PRINTPLANETNAME(PLANET_NUM)
         CASE (3)
c            if(global_NX.eq.60) then !---  these are for 60_cells_burrows.dat
c               D_RH =  -1.9417530*(10.**(-10.0))
c               D_T  =  -6.2500764*(10.**(-6.0))
c            endif
c            if(global_NX.eq.100) then  
c               if(planet_num.eq.1) then !---  these are for 100_HD209_Hansen.dat
c                  D_RH =   -9.6439106*(10.**(-11.))
c                  D_T =   -3.4225022*(10.**(-06.0))
c               endif
c               if(planet_num.eq.2) then !---  these are for 100_HD189_Hansen.dat
c                  call BurrowsOPC(Tinitbottom,RHinitbottom,xkapP_bottom,
c     $                 xkapRinitbottom,xkapA_bottom,xkapW_bottom(:),
c     $                 xkapLW_bottom(:),local_BB_bottom(:),
c     $                 dxkapP_dT_bottom,dxkapA_dT_bottom,1,1,1)
c                  D_RH = 0.d0
c                  D_T =  -3.d0*Fbottom*xkapRinitbottom*RHinitbottom/
c     $                 (4.d0*ARAD*FD*(Tinitbottom**3.d0))
c               endif
c            endif
            D_RH = 0.d0
            call BurrowsOPC(Tinitbottom,RHinitbottom,xkapP_bottom,
     $           xkapRinitbottom,xkapA_bottom,xkapW_bottom(:),
     $           xkapLW_bottom(:),local_BB_bottom(:),
     $           dxkapP_dT_bottom,dxkapA_dT_bottom,1,1,1)
            D_T =  -3.d0*Fbottom*xkapRinitbottom*RHinitbottom/
     $           (4.d0*ARAD*FD*(Tinitbottom**3.d0))
            if(myid.eq.0) then
               print *,'Setting D_T with radiative flux and'
               print *,' Fbottom=',Fbottom
               print *,' D_T',D_T
               write(*,'(A,1x,3(1x,e12.6))') 
     $              '  factors: kapR,rh,T=',
     $              xkapRinitbottom,RHinitbottom,Tinitbottom
            endif
         CASE (4)
            if(planet_num.eq.1) then !      planet_num=1 : HD209458b
c               if(global_NX.eq.100) then !---  these are for 100_cells_209458_poly.dat
c                  D_RH =  -1.1134852*(10.**(-10.0))
c                  D_T  =  -2.4122013*(10.**(-5.0))
c               endif
               if(global_NX.eq.100) then !---  these are for 100_cells_209458_poly_MAY.dat and test_subsolar209.dat
                  D_RH =  -2.1976300*(10.d0**(-10.d0))
                  D_T  =  -5.0305445*(10.d0**(-6.d0))
               endif
c              if(global_NX.eq.100) then !- for 100_160_64_cells_209458b_poly.dat
c                  D_RH =  -2.1976300*(10.d0**(-10.d0))
c                  D_T  =  -5.0305445*(10.d0**(-6.d0))
c               endif
            elseif(planet_num.eq.2) then !      planet_num=2 : HD189733b
               if(global_NX.eq.100) then !---  these are for 100_cells_189733b_poly.dat
                  D_RH =  -4.1081340*(10.**(-10.0))
                  D_T  =  -4.6779189*(10.**(-5.0))
               endif
            elseif(planet_num.eq.5) then !      planet_num=5 : WASP-12b
               if(global_NX.eq.100) then !---  these are for 100log_wasp12_RLfill_radzone.dat
                  D_RH =   -5.5412161e-12
                  D_T =   -4.1841993e-06
               endif
            else
               print *,'No initiliation file defined for planet_num',
     $              planet_num
               call clean_stop
            endif
         CASE(5)
!-- set the initial guess to burrows values
            if(global_NX.eq.60) then !---  these are for 60_cells_burrows.dat
               D_RH =  -1.9417530*(10.**(-10.0))
               D_T  =  -6.2500764*(10.**(-6.0))
            endif
            if(global_NX.eq.100) then !---  these are for 100_cells_burrows.dat
               D_RH =  -1.9756314*(10.**(-10.0))
               D_T  =  -5.9942251*(10.**(-6.0))
            endif
         CASE (6:)
            print *,'MODVER is not defined:init',modver
         CASE DEFAULT
            print *,'MODVER is also not defined:init',modver
         END SELECT
      else
         print *,'PROBABLY NEED TO SET D_T for modtyp=',modtyp
         call clean_stop
      endif
      if(myid.eq.0) then
         print *,' Bottom D_RH, D_T =',D_RH,D_T
      endif

!-- calculate the rotation rate
      CALL ROTATION_INITILIZATION
      RETURN
      END SUBROUTINE localinit

      SUBROUTINE GLOBALINIT
!--- Initilize variables
      USE input_init
      IMPLICIT NONE
      SELECT CASE(MODTYP)
      CASE (1)
         CALL INITCARTESIAN
      CASE (2)
         CALL INITPLANETMODEL
      CASE (3)
         CALL INITDISK
      CASE (4)
         CALL INITDREDGEUP
      CASE (5)
         CALL INITDISK
      CASE (6)
         CALL INITPACC
      CASE (7:)
         print *,'MODTYP not defined'
      CASE DEFAULT
         print *,'MODTYP also not defined'
      END SELECT
      RETURN
      END SUBROUTINE GLOBALINIT

      SUBROUTINE INITPLANETMODEL
      USE global_var_init
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE global_constants
      USE grid_var_init
      IMPLICIT NONE
      integer :: i,j,k,ipeter
      double precision :: r_test(global_NX)
      double precision :: Tincident4,v1,v2,v3,v4,kapP,kapA
      double precision :: tau_star,dr,testr,testnr,testnphi,gam_init
!----- initial velocities
      global_V = 0.d0
      global_G = 0.d0
      global_H = 0.d0

!----- initial temperature and density
      SELECT CASE(MODVER)
      CASE (0) !--planet
         if(global_NX.eq.60) then
            open(40,file="60_cells_peter.dat",status="old") !-- low resolution
            if(myid.eq.0) print *,'reading:60_cells_peter.dat'
         endif
         if(global_NX.eq.160) then
            open(40,file="160_cells_peter.dat",status="old") !-- high resolution
            if(myid.eq.0) print *,'reading:160_cells_peter.dat'
         endif
      CASE (1) !--brown dwarf
         open(40,file="100_1808_poly.dat",status="old")
         if(myid.eq.0) print *,'reading:100_1808_poly.dat'
      CASE (2) !--deep planet model
         if(global_NX.eq.60) then
            open(40,file="60_cells_deep.dat",status="old")      
            if(myid.eq.0) print *,'reading:60_cells_deep.dat'
         endif
         if(global_NX.eq.160) then
            open(40,file="160_cells_deep.dat",status="old")      
            if(myid.eq.0) print *,'reading:160_cells_deep.dat'
         endif
      CASE (3) !--burrows planet model
         if(global_NX.eq.60) then
            open(40,file="60_cells_burrows.dat",status="old")      
            if(myid.eq.0) print *,'reading:60_cells_burrows.dat'
         endif
         if(global_NX.eq.100) then
            if(planet_num.eq.1) then
               open(40,file="100_HD209_Hansen.dat",status="old")      
               if(myid.eq.0) print *,'reading:100_HD209_Hansen.dat'
            endif
            if(planet_num.eq.2) then
               open(40,file="100_HD189_Hansen.dat",status="old")      
               if(myid.eq.0) print *,'reading:100_HD189_Hansen.dat'
            endif
            if(planet_num.eq.5) then
c               open(40,file="100_wasp12_radzone.dat",status="old")      
c               if(myid.eq.0) print *,'100_wasp12_radzone.dat'
               open(40,file="100_wasp12_PROPOSAL.dat",status="old")
               if(myid.eq.0) print *,'100_wasp12_PROPOSAL.dat'
            endif
         endif
         if(global_NX.eq.200) then
c            open(40,file="200_HD209_Hansen_deep.dat",status="old")      
c            if(myid.eq.0) print *,'reading:200_HD209_Hansen_deep.dat'
            open(40,file="200_HD209_Hansen.dat",status="old")      
            if(myid.eq.0) print *,'reading:200_HD209_Hansen.dat'
         endif
      CASE (4)
         if(planet_num.eq.1) then !-HD209458
            if(global_NX.eq.100) then
c               open(40,file="100_cells_209458b_poly.dat",status="old")
c               if(myid.eq.0) 
c     $              print *,'reading:100_cells_209458b_poly.dat'
c               open(40,file="100_cells_209458b_poly_MAY.dat",
c     $              status="old")
c               if(myid.eq.0) print *,
c     $              'reading:100_cells_209458b_poly_MAY.dat'
               open(40,file="test_subsolar209.dat",status="old")
               if(myid.eq.0) print *,
     $              'reading:test_subsolar209.dat'
c     open(40,file="100_160_64_cells_209458b_poly.dat"
c     $              ,status="old")
c     if(myid.eq.0) print *,
c     $              'reading:100_160_64_cells_209458b_poly.dat'
            endif  
         elseif(planet_num.eq.2) then !-HD189733b
            if(global_NX.eq.100) then
               open(40,file="100_cells_189733b_poly.dat",status="old")
               if(myid.eq.0) 
     $              print *,'reading:100_cells_189733b_poly.dat'
            endif
         elseif(planet_num.eq.5) then !-wasp12b
            if(global_NX.eq.100) then
c-this is the RL filling model of only the radiative zone
c               open(40,file="wasp12_RLfill_radzone_100.dat",
c     $              status="old")
c               if(myid.eq.0) 
c     $              print *,'reading:wasp12_RLfill_radzone_100.dat'               
c               open(40,file="100_wasp12_RLfill_radzone.dat",
c     $              status="old")
               if(myid.eq.0) 
     $              print *,'100_wasp12_RLfill_radzone.dat'
               open(40,file="100log_wasp12_RLfill_radzone.dat",
     $              status="old")
               if(myid.eq.0) 
     $              print *,'100log_wasp12_RLfill_radzone.dat'
            endif
         endif
         if(global_NX.eq.200) then
            print *,'NO INITILIZATION FILE'
            call clean_stop
         endif
      CASE (5)
!--use burrows planet model as initial guess
         if(global_NX.eq.60) then
            open(40,file="60_cells_burrows.dat",status="old")      
            if(myid.eq.0) print *,'reading:60_cells_burrows.dat'
         endif
         if(global_NX.eq.100) then
            open(40,file="100_cells_burrows.dat",status="old")      
            if(myid.eq.0) print *,'reading:100_cells_burrows.dat'
         endif
      CASE (6)
         if(global_NX.ne.20) then
            print *,'ERROR: Held/Suarez model should have nx=20'
            call clean_stop
         endif
         open(40,file="Held_Suarez.dat.dat",status="old")      
         if(myid.eq.0) print *,'reading:Held_Suarez.dat'
!--use the Held Suarez initial model
      CASE (7:)
         print *,'MODVER is not defined:initplanet',modver
      CASE DEFAULT
         print *,'MODVER is also not defined:initplanet',modver
      END SELECT

!-READ IN THE INFORMATION FROM THE FILE

!--extra read in info for brown-dwarf model
      if(MODVER.eq.1) then
         read(40,*) Tinitbottom
         read(40,*) RHinitbottom
         read(40,*) D_T
         read(40,*) D_RH
         if(myid.eq.0) then
            print*,'test',Tinitbottom,RHinitbottom,D_T,D_RH
         endif
         global_upperbnd = global_NX
      endif

!--extra read in info for Held_Suarez.dat
      if(MODVER.eq.6) then
         read(40,*) Tinitbottom
         read(40,*) RHinitbottom
         global_upperbnd = global_NX
      endif

      read(40,*) ipeter
      if(ipeter.gt.global_NX) THEN
         print *,'problem with model read in',ipeter,global_NX
         call clean_stop
      ENDIF
!------the initial upper bound=upperbnd (ipeter here)
      global_upperbnd = ipeter
      DO I=1,ipeter
         read(40,*) global_T(i,1,1),global_RH(i,1,1),r_test(i)
      ENDDO
      DO K=1,global_NZ
         DO J=1,global_NY
            DO I=1,ipeter
               global_T(i,j,k) = global_T(i,1,1)
               global_RH(i,j,k) = global_RH(i,1,1)
            ENDDO
            if(ipeter+1.lt.global_NX) then 
               DO I=ipeter+1,global_NX
                  global_T(i,j,k) = global_T(ipeter,1,1)
                  global_RH(i,j,k) = 10.d0**(-10.d0)
c                  global_RH(i,j,k) = global_RH(ipeter,1,1)
               ENDDO
            endif
         ENDDO
      ENDDO

!-CLOSE THE INITIALIZATION FILE      
      close(40)

!--- Initial ER distribution. Assume optically thick everywhere
c      global_ER = ARAD*((global_T)**(4.0))
c      global_ER = max(global_ER,0.3)
c-use radiation steadystate as the initial condition(paper 2, eq18)
c Here gamma=0.5
      gam_init=0.5
      global_ER = ( (4.d0*SBCONST*(global_T)**(4.0)) - 
     $     (gam_init**4)*(SBCONST*(Tirr**(4.0))) )/FD


!-PRINTOUT THE INTIAL STATE AT SUBSTELLAR POINT TO THE SCREEN
      if(myid.eq.0) then
         print *, 'ipeter',ipeter
         DO I=1,global_NX
            write(*,'(I8,1x,3(1x,e12.6))') 
     &           i,global_T(i,1,1),global_RH(i,1,1),xxb(i)
         ENDDO
      endif

      RETURN
      END SUBROUTINE INITPLANETMODEL


      SUBROUTINE ROTATION_INITILIZATION
      USE fluid_var_init
      USE input_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      double precision :: period_start
!--- Start rotation
      SELECT CASE(IROTATE)
      CASE (0)
         OMROT = 0.d0
      CASE (1) !-- ignore secondary mass in omega calculation
         OMROT = 2.d0*PI/(PROT*DAY)
      CASE (2) !-- include both masses
         if(mstar2.eq.0.d0) then
            OMROT = (GRAV*MSOL*mstar1/((rp*AU)**3.d0))**0.5
            if(myid.eq.0) then
               write(*,'(A,e12.3)') 'OMROT=',OMROT
            endif
         else
            if(myid.eq.0) then
               print *,'Including Secondary in OMEGA calc'
               write(*,'(A,e12.3)') 'M_1 (MSOL)=',mstar1
               write(*,'(A,e12.3)') 'M_2 (MSOL)=',mstar2
               write(*,'(A,e12.3)') 'a (AU)=',rp
            endif
            OMROT = (GRAV*MSOL*(mstar1+mstar2)/((rp*AU)**3.d0))**0.5
         endif
      CASE (3)
         rot_start = ramp_start*DAY
         rot_full  = ramp_full*DAY
         if(hydrotime.lt.rot_start) then
            OMROT = 0.d0
         elseif (hydrotime.ge.rot_start.and.
     %           hydrotime.lt.rot_full) then
            OMROT = (2.d0*PI/(PROT*DAY))*
     %           ((hydrotime-rot_start)/(rot_full-rot_start))
         elseif(hydrotime.ge.rot_full) then
            OMROT = 2.d0*PI/(PROT*DAY)
         endif
         if(myid.eq.0) then
            print *,'irotate=3 information:'
            print *,'rot_start  =',rot_start
            print *,'rot_full   =',rot_full
            print *,'init_time   =',zeit
            print *,'init_hydrotime   =',hydrotime
            print *,'omrot_init =',omrot
         endif
      CASE (4) !(MODIFY period_start here AND in pp.f)
         period_START = 3.52
         OMROT = 2.d0*PI/(period_START*DAY)
         if(myid.eq.0) then
            print *,'irotate=4 information:'
            print *,'period_start  =',period_start
         endif
      CASE (5:)
         print *,'IROTATE not defined',irotate
      CASE DEFAULT
         print *,'IROTATE also not defined'
      END SELECT
      RETURN
      END SUBROUTINE ROTATION_INITILIZATION


      SUBROUTINE INITCARTESIAN
      USE global_var_init
      USE grid_var_init
      USE input_init
      USE fluid_var_init
      USE force_var_init
      USE mpi_var_init
      USE global_constants
      IMPLICIT NONE
      integer :: i,j,k,type,rand_size
      double precision :: P_val,rho_val,npoly
      double precision :: Ttop,rhtop,rrr,vpert_mag
      integer,allocatable,dimension(:) :: k_rand
      double precision :: r_test(global_NX)
      if(ncosys.ne.0) then
         print *,'This mode is only catesian',ncosys
         call clean_stop
      endif 
      SELECT CASE(MODVER)
      CASE (0) !--convecting box
         type = 0
!---  Initilization for a convective region
!-----initial velocities
         global_V = 0.d0
         global_G = 0.d0
         global_H = 0.d0
!----- x-bottom
         if(type.eq.0) then !-any changes to D_T should be made in localinit too
            call fgrav
            npoly = 1.d0
            D_T = -ABS(gravx(1,1,1))*mu_gas/(RGAS*(npoly+1.))
            Ttop = 200.0
            rhtop = 1.0
c            Ttop = 2000.0
c            rhtop = 2e-8
            if(myid.eq.0) print *,'D_T=',D_T
            if(myid.eq.0) print *,'npoly=',npoly
            if(myid.eq.0) print *,'Ttop=',Ttop
            if(myid.eq.0) print *,'rhtop=',rhtop
            DO i=global_NX,1,-1
!-from brummel 2002 polytropic EOS
               global_T(i,1,1)=D_T*(xxb(i)-xxb(global_nx)) + Ttop
               global_RH(i,1,1)=rhtop*((global_T(i,1,1)/Ttop)**(npoly))
            ENDDO
!-set random seed to myid value
            call random_seed(size=rand_size)                                                         
            Allocate(k_rand(rand_size))
            k_rand = myid
            call random_seed(put=k_rand(1:rand_size) )
!-random vertical velocity distrubuted as +/- vpert_mag
            vpert_mag = 0.1
            DO k=1,global_NZ
               DO J=1,global_NY
                  DO I=1,global_NX
                     global_T(i,j,k) = global_T(i,1,1)
                     global_RH(i,j,k) = global_RH(i,1,1)
                     call random_number(rrr)
                     global_V(i,j,k) = 2.0*vpert_mag*rrr-vpert_mag
                     if(myid.eq.0.and.k.eq.1.and.j.eq.1) then
                        write(*,'(A,I8,3(1x,e12.3))') 
     $                       'global init(i,T,rh,Vx)',
     $                       i,global_T(i,1,1),global_rh(i,1,1),
     $                       global_V(i,1,1)
                     endif
                  ENDDO
               ENDDO
            ENDDO
!------y-bottom
         elseif(type.eq.1) then
            DO j=1,global_NY
               global_T(1,j,1)=j*(500.d0-1000.d0)/
     %              (1.d0*global_NY)+1000.d0
               global_RH(1,j,1)=j*(0.1-0.7)/(1.d0*global_NY)+0.7
            ENDDO
            DO k=1,global_NZ
               DO J=1,global_NY
                  DO I=1,global_NX
                     global_T(i,j,k) = global_T(1,j,1)
                     global_RH(i,j,k) = global_RH(1,j,1)
                  ENDDO
               ENDDO
            ENDDO
         endif
         TinitBottom = global_T(1,1,1)
         Tinittop    = global_T(global_nx,1,1)
         rhinitbottom = global_RH(1,1,1)
         rhinittop   = global_RH(global_NX,1,1)
c         DO k=1,global_NZ
c            DO J=1,global_NY
c               global_upperbnd(j,k) = locNX
c            ENDDO
c         ENDDO
         global_upperbnd = locNX
      CASE (1) !-- SOD Shockwave test
         global_V = 0.d0
         global_G = 0.d0
         global_H = 0.d0
         DO I=1,global_NX
            if(i.lt.global_NX/2) then
               P_val = 1.0
               rho_val = 1.0
            else
               P_val = 0.1
               rho_val = 0.125
            endif
            global_T(i,1,1) = mu_gas*P_val/(RGAS*rho_val)
            global_RH(i,1,1) = rho_val
            if(myid.eq.0) write(*,'(A,I8,2(1x,e12.3))') 'global init',i,
     %           global_T(i,1,1),global_rh(i,1,1)
         enddo
         DO k=1,global_NZ
            DO J=1,global_NY
               DO I=1,global_NX
                  global_T(i,j,k) = global_T(i,1,1)
                  global_RH(i,j,k) = global_RH(i,1,1)
               ENDDO
               global_upperbnd(j,k) = locNX
            ENDDO
         ENDDO
      CASE (2)


!---------------------------------------------------------------------------
!---------SIMPLE BOX CONVECTION (SAME AS MODVER=0) -------------------------
!---------------------------------------------------------------------------
!---  Initilization for a convective region
!-----initial velocities
c         global_V = 0.d0
c         global_G = 0.d0
c         global_H = 0.d0
c         call fgrav
c         npoly = 1.d0
c         D_T = -ABS(gravx(1,1,1))*mu_gas/(RGAS*(npoly+1.))
c         Ttop = 200.0
c         rhtop = 1.0
c         if(myid.eq.0) print *,'D_T=',D_T
c         if(myid.eq.0) print *,'npoly=',npoly
c         if(myid.eq.0) print *,'Ttop=',Ttop
c         if(myid.eq.0) print *,'rhtop=',rhtop
c         DO i=global_NX,1,-1
c!-from brummel 2002 polytropic EOS
c            global_T(i,1,1)=D_T*(xxb(i)-xxb(global_nx)) + Ttop
c            global_RH(i,1,1)=rhtop*((global_T(i,1,1)/Ttop)**(npoly))
c         ENDDO
c!-set random seed to myid value
c         call random_seed(size=rand_size)                                                         
c         Allocate(k_rand(rand_size))
c         k_rand = myid
c         call random_seed(put=k_rand(1:rand_size) )
c!-random vertical velocity distrubuted as +/- vpert_mag
c         vpert_mag = 0.1
c         DO k=1,global_NZ
c            DO J=1,global_NY
c               DO I=1,global_NX
c                  global_T(i,j,k) = global_T(i,1,1)
c                  global_RH(i,j,k) = global_RH(i,1,1)
c                  call random_number(rrr)
c                  global_V(i,j,k) = 2.0*vpert_mag*rrr-vpert_mag
c                  if(myid.eq.0.and.k.eq.1.and.j.eq.1) then
c                     write(*,'(A,I8,3(1x,e12.3))') 
c     $                    'global init(i,T,rh,Vx)',
c     $                    i,global_T(i,1,1),global_rh(i,1,1),
c     $                    global_V(i,1,1)
c                  endif
c               ENDDO
c            ENDDO
c         ENDDO

!---------------------------------------------------------------------------
!---------EARTH CONVECTION -------------------------------------------------
!---------------------------------------------------------------------------

!---  Initilization for Earth atmosphere
         if(global_NX.eq.300) then
            open(40,file="300_cells_earth_dry_isothermal.dat",
     $           status="old")
            if(myid.eq.0) 
     $           print *,'300_cells_earth_dry_isothermal.dat'
         endif

         if(global_NX.eq.100) then
            open(40,file="100_cells_earth_dry_isothermal.dat",
     $           status="old")
            if(myid.eq.0) 
     $           print *,'100_cells_earth_dry_isothermal.dat'
         endif

         if(global_NX.eq.50) then
            open(40,file="50_cells_earth_dry_isothermal.dat",
     $           status="old")
            if(myid.eq.0) 
     $           print *,'50_cells_earth_dry_isothermal.dat'
         endif

         DO I=1,global_nx
            read(40,*) global_T(i,1,1),global_RH(i,1,1),r_test(i)
         ENDDO
         D_T = (global_T(2,1,1)-global_T(1,1,1))/(xxb(2)-xxb(1))
         if(myid.eq.0) print *,'D_T=',D_T
!-----initial velocities
         global_V = 0.d0
         global_G = 0.d0
         global_H = 0.d0

!-set random seed to myid value
         call random_seed(size=rand_size)                                                         
         Allocate(k_rand(rand_size))
         k_rand = myid
         call random_seed(put=k_rand(1:rand_size) )
!-random vertical velocity distrubuted as +/- vpert_mag
         vpert_mag = 0.0001
         DO k=1,global_NZ
            DO J=1,global_NY
               DO I=1,global_NX
                  global_T(i,j,k) = global_T(i,1,1)
                  global_RH(i,j,k) = global_RH(i,1,1)
                  call random_number(rrr)
                  global_V(i,j,k) = 2.0*vpert_mag*rrr-vpert_mag
                  if(myid.eq.0.and.k.eq.1.and.j.eq.1) then
                     write(*,'(A,I8,3(1x,e12.3))') 
     $                    'global init(i,T,rh,Vx)',
     $                    i,global_T(i,1,1),global_rh(i,1,1),
     $                    global_V(i,1,1)
                  endif
               ENDDO
            ENDDO
         ENDDO
         
!---------------------------------------------------------------------------
!-------------STUFF FOR BOTH------------------------------------------------
!---------------------------------------------------------------------------
         TinitBottom = global_T(1,1,1)
         Tinittop    = global_T(global_nx,1,1)
         rhinitbottom = global_RH(1,1,1)
         rhinittop   = global_RH(global_NX,1,1)
         global_upperbnd = locNX
      CASE (3:)
         print *,'MODVER is not defined:initcartesian',modver
      CASE DEFAULT
         print *,'MODVER is also not defined:initcartesian',modver
      END SELECT      

      RETURN
      END SUBROUTINE INITCARTESIAN


      SUBROUTINE INITDISK
      USE global_var_init
      USE input_init
      USE fluid_var_init
      USE global_constants
      USE mpi_var_init
      USE grid_var_init
      IMPLICIT NONE
      integer :: i,j,k,type
      double precision :: r
!--- Initilization for a convective region
!----- initial velocities (G IS SET TO KEPLARIAN IN KEPLARINIT)
      global_V = 0.d0
      global_G = 0.d0
      global_H = 0.d0
      !--- DISK PARAMETERS
      sigma0 = 950.d0
      ndisk = -1.5d0
      Tdisk0 = 280.d0
      if(myid.eq.0) write(*,'(A,e12.3)') 'Sigma_o =',sigma0
      if(myid.eq.0) write(*,'(A,e12.3)')
     %     'n for surface density =',ndisk
      if(myid.eq.0) write(*,'(A,e12.3)') 'Tdisk0 =',Tdisk0
      if(myid.eq.0)  print *,'initial model: r(AU), rh(r), T(r):'
      DO i=1,global_NX
         r = xmin+(DXB(1)/2.d0) + (xmax-xmin-DXB(1))*(i-1)/(global_NX-1)
         global_RH(i,1,1) = sigma0*((r/AU)**(ndisk))
         global_T(i,1,1)  = Tdisk0*((r/AU)**(-0.5))
         if(myid.eq.0)  write(*,'(I8,e12.3,e12.3,e12.3)') 
     $        i,r/AU,global_RH(i,1,1),global_T(i,1,1)
      ENDDO
      DO k=1,global_NZ
         DO J=1,global_NY
            DO I=1,global_NX
               global_T(i,j,k) = global_T(i,1,1)
               global_RH(i,j,k) = global_RH(i,1,1)
            ENDDO
         ENDDO
      ENDDO
      TinitBottom = global_T(1,1,1)
      Tinittop    = global_T(global_nx,1,1)
      rhinitbottom = global_RH(1,1,1)
      rhinittop   = global_RH(global_NX,1,1)
      if(myid.eq.0) then
         print *,'init rh',rhinitbottom,rhinittop
      endif
      DO k=1,global_NZ
         DO J=1,global_NY
            global_upperbnd(j,k) = locNX
         ENDDO
      ENDDO
!-- not used, but set ER
      global_ER = ARAD*((global_T)**(4.0))
      RETURN
      END SUBROUTINE INITDISK


      SUBROUTINE INITDREDGEUP
      USE global_var_init
      USE input_init
      USE fluid_var_init
      USE global_constants
      USE mpi_var_init
      USE grid_var_init
      IMPLICIT NONE
      integer :: i,j,k,type
      double precision :: r,npoly,rinf,rhoinf,Rbondi,Kpoly
      double precision :: t0,t1,t2,t3
      logical :: shulin,polytrope,analytic
      
      shulin=.false.
      polytrope=.false.
      analytic=.true.

      if(shulin.and.polytrope) then
         print *,'choose initial model for dredgeup'
         call clean_stop
      endif

!--- Initial velocities
      global_V = 0.d0
      global_G = 0.d0
      global_H = 0.d0


!-- POLYTROPIC INITIAL CONDITIONS
      if(polytrope) then
         if(myid.eq.0) then
            print *,'reading: 120_cells_polytrope.dat'
         endif
         open(40,file="120_cells_polytrope.dat",status="old")
         DO I=1,global_NX
            read(40,*) global_T(i,1,1),global_RH(i,1,1)
         ENDDO
         close(40)
      endif

!-- POLYTROPIC INITIAL CONDITIONS
      if(analytic) then
         if(myid.eq.0) then
            print *,'reading: 120_cells_analytic_deep.dat'
         endif
         open(40,file="120_cells_analytic_deep.dat",status="old")
         read(40,*) K_poly
         read(40,*) n_poly
         DO I=1,global_NX
            read(40,*) global_T(i,1,1),global_RH(i,1,1)
         ENDDO
         close(40)
      endif

!--- READ SHULIN'S FILE
      if(shulin) then
         open(40,file="120_cells_Shulin.dat",status="old")
         DO I=1,global_NX
            read(40,*) global_T(i,1,1),global_RH(i,1,1)
         ENDDO
         close(40)
      endif

!---  PRINT OUT THE INTIAL STATE
      if(myid.eq.0) then
         do i=1,global_NX
            print *,'initialt/rh',i,global_T(i,1,1),global_RH(i,1,1)
         enddo
      endif

!---  SAVE INITIAL BOUNDARY INFOMATION
      TinitBottom = global_T(1,1,1)
      Tinittop    = global_T(global_nx,1,1)
      rhinitbottom = global_RH(1,1,1)
      rhinittop   = global_RH(global_NX,1,1)
      
!---  PRINT INITIAL BOUNDARY INFOMATION
      if(myid.eq.0) print *,'Bottom Temperature =',TinitBottom
      if(myid.eq.0) print *,'Top Temperature    =',Tinittop
      if(myid.eq.0) print *,'Wind Density (top) =',rhinittop
      if(myid.eq.0) print *,'Wind Amplitude =',windAMP

!---  FILL IN THE REST OF THE ARRAY ASSUMING SPHERICAL SYMMETRY
      DO K=1,global_NZ
         DO J=1,global_NY
            DO I=1,global_NX
               global_T(i,j,k) = global_T(i,1,1)
               global_RH(i,j,k) = global_RH(i,1,1)
            ENDDO
         ENDDO
      ENDDO

!-------set upper boundaries to the top
      DO k=1,global_NZ
         DO J=1,global_NY
            global_upperbnd(j,k) = locNX
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE INITDREDGEUP
      
      SUBROUTINE INITPACC
      USE global_var_init
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE global_constants
      USE grid_var_init
      IMPLICIT NONE
      integer :: i,j,k
      double precision :: x_dist,z_dist,dist,max_dist
!----- initial velocities
      global_V = 0.d0
      global_G = 0.d0
      global_H = 0.d0

!---- set exterior nebula conditions here
      T_nebula = 100.0
      RH_nebula = (10.d0)**(-8.d0)
      
      SELECT CASE(MODVER)
      CASE (0) !--1/r profiles for density, polytrope EOS determined T
         max_dist = XMAX - (dxb(1)/2.d0)
         DO K=1,global_NZ
            DO J=1,global_NY
               DO I=1,global_NX
                  x_dist = (XMIN + (dxb(1)/2.d0)) + 
     %                 dxb(1)*(i-1)
                  z_dist = (ZMIN + (dzb(1)/2.d0)) + 
     %                 dzb(1)*(k-1)
                  dist=((x_dist**(2.d0))+(z_dist**(2.d0)))**(0.5)
                  global_T(i,j,k)  = T_nebula*
     %                 (max_dist/dist)**(GAMMA-1.d0)
                  global_RH(i,j,k) = RH_nebula*(max_dist/dist)
                  if(myid.eq.0.and.k.eq.1.and.j.eq.1) then
                     write(*,'(I8,1x,4(1x,e12.6))') 
     1                    i,dist,x_dist,global_T(i,1,1),global_RH(i,1,1)
                  endif
               ENDDO
            ENDDO
         ENDDO
!-- set floor for T and RH
         global_T  = MAX(global_T,T_nebula)
         global_RH = MAX(global_RH,RH_nebula)
      CASE (1:)
         print *,'MODVER is not defined:initpacc',modver
      CASE DEFAULT
         print *,'MODVER is also not defined:initpacc',modver
      END SELECT      
!---Initial ER distribution
      global_ER = ARAD*((global_T)**(4.0))

!---Printout initial model
      if(myid.eq.0) then
         DO I=1,global_NX
            write(*,'(I8,1x,3(1x,e12.6))') 
     1           i,global_T(i,1,1),global_RH(i,1,1)
         ENDDO
         print *,'now at k=10....'
         DO I=1,global_NX
            write(*,'(I8,1x,3(1x,e12.6))') 
     1           i,global_T(i,1,10),global_RH(i,1,10)
         ENDDO
      endif

!---fix upperbnd=locNX+1
      global_upperbnd = global_NX+1

      RETURN
      END SUBROUTINE INITPACC

      SUBROUTINE GITTER
**  Generates the Gridsystem
      USE input_init
      USE grid_var_init
      USE fluid_var_init
      USE mpi_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: YMAX1,YMIN1
      double precision :: ZMAX1,ZMIN1
      double precision :: DELX,DELY,DELZ
      double precision :: DRAN3,XXAN3,DYAN3,XYAN3
      double precision :: DZAN3,XZAN3
      double precision :: prop,xxa1
      INTEGER :: IERROR
      logical :: full2pi
      integer,dimension(locNX+1) :: int_array
      double precision ,dimension(locNX+1) :: rescaled
      double precision ,dimension(2) :: scaleFactor
      double precision :: maxRange,minRange,vectorMax,vectorMin
      if (ncosys.eq.0) then
**  If Cartesian, Input length units!!
         YMAX1=YMAX
         YMIN1=YMIN
         ZMAX1=ZMAX
         ZMIN1=ZMIN
      else if (ncosys.eq.1) then
**  If Cylindrical: Y-Coordinates in Degrees, Z-Coordinate in length units
        YMAX1=YMAX*2.d0*pi
        YMIN1=YMIN*2.d0*pi
        ZMAX1=ZMAX
        ZMIN1=ZMIN
      else if (ncosys.eq.2) then
**  If Spherical Polars: Y and Z-Coordinates in Degrees!!
         YMAX1=YMAX*2.d0*PI
         YMIN1=YMIN*2.d0*PI
         ZMAX1=(ZMAX/2.d0)*pi
         ZMIN1=(ZMIN/2.d0)*pi
      else
         write(*,*) ' Coordinate System NOT possible: Check NCOSYS'
         print *,'NCOSYS=',NCOSYS
         call clean_STOP
      end if

**  Initialization of X-Grid
      IF(int(DXMIN).eq.0) THEN !-linear grid
         DELX=(XMAX-XMIN)/global_NX
         DXA=DELX
         DXA2=DXA*DXA
         prop=1.d0
         XXA(1)=XMIN+(proc_coords(1)*(XMAX-XMIN)/proc_dims(1))
         XXA(0)=XXA(1)-DXA(1)
         do i=2,locNX+2
            XXA(I)=XXA(I-1)+DXA(I)
         enddo
         do i=0,locNX+1
            XXB(I)=(XXA(I)+XXA(I+1))/2.d0
         enddo
         DRAN3=DXA(locNX+2)*PROP
         XXAN3=XXA(locNX+2)+DRAN3
         XXB(locNX+2)=(XXA(locNX+2)+XXAN3)/2.d0
         DRAN3=DXA(1)/PROP
         XXAN3=XXA(0)-DRAN3
         XXB(-1)=(XXA(0)+XXAN3)/2.d0
      ELSEIF(int(DXMIN).eq.1) then 
         if(proc_dims(1).ne.1) then
            print *,'The log-grid setup only works with proc_dim(1)=1'
            call clean_stop
         endif
         vectorMin= 1.d0
         vectorMax= locNX+1
         int_array=(/(i,i=int(vectorMin),int(vectorMax),1)/)         
         maxRange= DLOG10(XMAX)
         minRange= DLOG10(XMIN)
         scaleFactor = [((minRange*vectorMax)-(maxRange*vectorMin))/
     $        (vectorMax-vectorMin), (maxRange-minRange)/
     $        (vectorMax-vectorMin)]
         rescaled= int_array * scaleFactor(2) + scaleFactor(1)         
         XXA(1:locNX+1) = 10.d0**(rescaled)
         XXA(0) = XXA(1)-(XXA(2)-XXA(1))
         XXA(LOCNX+2) = XXA(LOCNX+1)+(XXA(LOCNX+1)-XXA(locNX))
         do i=0,locNX+1
            XXB(I)=(XXA(I)+XXA(I+1))/2.d0
         enddo
         XXB(locNX+2)=XXB(locNX+1)+(XXB(locNX+1)-XXB(locNX))
         XXB(-1)=XXB(0)-(XXB(1)-XXB(0))
!-directly calculate dxa values (dxb is done below)
         do i=1,locNX+2
            dxa(i) = xxa(i)-xxa(i-1)
         enddo
         DXA2=DXA*DXA
!-write out grid to screen
         if(myid.eq.0) then
            write(*,'(A)') 'LOG X-Grid: i,xxa(i),xxb(i)'
            do i=0,locNX+2
               write(*,'(I8,2(1x,e12.3))') i,xxa(i),xxb(i)
            enddo
         endif
      else
         print *,'DXMIN=',DXMIN,'is not available'
         call clean_stop
      endif

**  Initialization of Y-Grid
      DELY=(YMAX1-YMIN1)/global_NY
      DYA=DELY
      DYA2=DYA*DYA
      prop=1.d0

      XYA(1)=YMIN1+(proc_coords(2)*(YMAX1-YMIN1)/proc_dims(2))
      XYA(0)=XYA(1)-DYA(1)
      DO J=2,locNY+2
         XYA(J)=XYA(J-1)+DYA(J)
      ENDDO
      DO J=0,locNY+1
         XYB(J)=(XYA(J)+XYA(J+1))/2.d0
      ENDDO
      if(PROP.ne.1.d0) then
         print *,'problem in gitter'
         print *,'also note hard coded dr in photosphere'
         call clean_stop
      endif
      DYAN3=DYA(locNY+2)*PROP
      XYAN3=XYA(locNY+2)+DYAN3
      XYB(locNY+2)=(XYA(locNY+2)+XYAN3)/2.d0
      DYAN3=DYA(1)/PROP
      XYAN3=XYA(0)-DYAN3
      XYB(-1)=(XYA(0)+XYAN3)/2.d0

**  Initialization of Z-Grid 
      DELZ=(ZMAX1-ZMIN1)/global_NZ
      DZA=DELZ
      DZA2=DZA*DZA
      prop=1.d0

      XZA(1)=ZMIN1+(proc_coords(3)*(ZMAX1-ZMIN1)/proc_dims(3))
      XZA(0)=XZA(1)-DZA(1)
      DO K=2,locNZ+2
         XZA(K)=XZA(K-1)+DZA(K)
      ENDDO
      DO K=0,locNZ+1
         XZB(K)=(XZA(K)+XZA(K+1))/2.d0
      ENDDO
      DZAN3=DZA(locNZ+2)*PROP
      XZAN3=XZA(locNZ+2)+DZAN3
      XZB(locNZ+2)=(XZA(locNZ+2)+XZAN3)/2.d0
      DZAN3=DZA(1)/PROP
      XZAN3=XZA(0)-DZAN3
      XZB(-1)=(XZA(0)+XZAN3)/2.d0

!--   Calculate the delta's for the B-grids
      DO I=0,locNX+2
         DXB(I)=XXB(I)-XXB(I-1)
      ENDDO
      DO J=0,locNY+2
         DYB(J)=XYB(J)-XYB(J-1)
      ENDDO
      DO K=0,locNZ+2
         DZB(K)=XZB(K)-XZB(K-1)
      ENDDO
!-- the squared delta_B's
      DXB2 = DXB*DXB
      DYB2 = DYB*DYB
      DZB2 = DZB*DZB
!
!   Calculation of the Additional Grids
!

      IF(DXMIN.eq.0) THEN !-hard code 1/2 for standard grid
         DO I=0,locNX+2
            DDX0(I)=0.5d0
            DDX1(I)=0.5d0
         ENDDO
      ELSE
         DO I=0,locNX+2
            DDX0(I)=(XXB(I)-XXA(I))/DXB(I)
            DDX1(I)=(XXA(I)-XXB(I-1))/DXB(I)
         ENDDO
      ENDIF
      IF(DYMIN.eq.0) THEN !-hard code 1/2 for standard grid
         DO J=0,locNY+2
            DDY0(J)=0.5d0
            DDY1(J)=0.5d0
         ENDDO
      ELSE
         DO J=0,locNY+2
            DDY0(J)=(XYB(J)-XYA(J))/DYB(J)
            DDY1(J)=(XYA(J)-XYB(J-1))/DYB(J)
         ENDDO
      ENDIF
      IF(DYMIN.eq.0) THEN !-hard code 1/2 for standard grid
         DO J=1,locNY+2
            DCY0(J)=0.5d0
            DCY1(J)=0.5d0
         ENDDO
      ELSE
         DO J=1,locNY+2
            DCY0(J)=(XYA(J)-XYB(J-1))/DYA(J)
            DCY1(J)=(XYB(J-1)-XYA(J-1))/DYA(J)
         ENDDO
      ENDIF
      IF(DZMIN.eq.0) THEN       !-hard code 1/2 for standard grid
         DO K=0,locNZ+2
            DDZ0(k)=0.5d0
            DDZ1(k)=0.5d0
         ENDDO
      ELSE
         DO K=0,locNZ+2
            DDZ0(k)=(XZB(k)-XZA(k))/DZB(k)
            DDZ1(k)=(XZA(k)-XZB(k-1))/DZB(k)
         ENDDO
      ENDIF
** add Keplarian ones for disk problem
      DO I=1,locNX+2
         DDX0kep(I)=DDX0(I)*XXB(I-1)**(1.5)
         DDX1kep(I)=DDX1(I)*XXB(I)**(1.5)
      ENDDO

**    in order to avoid possible divide-by-zero errors
**    replace zero with a small positive number
      xxa1=xxa(1)
      if (xxa1.eq.0.) xxa1=1.e-20
**  Geometry dependent Surface and Volume Elements
      IF (NCOSYS.EQ.0) THEN
**  Cartesian:  x, y and z
**     VISC requires:
**      SUR must be defined from 0 to N+2
**      VOLB must be defined from 0 to N+1
         DO I=0,locNX+1
            VOLXB(I)=XXA(I+1)-XXA(I)
         ENDDO
         DO I=0,locNX+2
            SURXA(I)=1.d0
            CCXA(I)=0.d0
            VOLXA(I)=XXB(I)-XXB(I-1)
            SURXB(I)=1.d0
            CCXB(I)=0.d0
            geoxga(I)=1.d0
            geoxha(I)=1.d0
         ENDDO
         DO I=-1,locNX+2
            geoxg(I)=1.d0
            geoxh(I)=1.d0
         ENDDO
         DO J=0,locNY+1
            VOLYB(J)=XYA(J+1)-XYA(J)
         ENDDO
         DO J=0,locNY+2
            SURYA(J)=1.d0
            CCYA(J)=0.d0
            VOLYA(J)=XYB(J)-XYB(J-1)
            SURYB(J)=1.d0
            CCYB(J)=0.d0
         ENDDO
         DO K=0,locNZ+1
            VOLZB(K)=XZA(K+1)-XZA(K)
         ENDDO
         DO K=0,locNZ+2
            SURZA(K)=1.d0
            CCZA(K)=0.d0
            VOLZA(K)=XZB(K)-XZB(K-1)
            SURZB(K)=1.d0
            CCZB(K)=0.d0
            geozga(k)=1.d0
         ENDDO
         DO K=-1,locNZ+2
            geozg(k)=1.d0
         ENDDO

      ELSE IF (NCOSYS.EQ.1) THEN
**  Cylindrical:  x<>r, y<>phi,  and  z<>z 
        DO I=0,locNX+1
          VOLXB(I)=(XXA(I+1)**2-XXA(I)**2)/2.d0
        ENDDO
        DO I=0,locNX+2
           SURXA(I)=XXA(I)
           CCXA(I)=1.
           VOLXA(I)=(XXB(I)**2-XXB(I-1)**2)/2.d0
           SURXB(I)=XXB(I)
           CCXB(I)=1.
           geoxga(I)=XXA(I)
           geoxha(I)=1.d0
        ENDDO
        geoxga(1)=xxa1
        DO I=-1,locNX+2
           geoxg(I)=XXB(I)
           geoxh(I)=1.d0
        ENDDO
        DO J=0,locNY+1
          VOLYB(J)=XYA(J+1)-XYA(J)
        ENDDO
        DO J=0,locNY+2
           SURYA(J)=1.d0
           CCYA(J)=0.d0
           VOLYA(J)=XYB(J)-XYB(J-1)
           SURYB(J)=1.d0
           CCYB(J)=0.d0
        ENDDO
        DO K=0,locNZ+1
           VOLZB(K)=XZA(K+1)-XZA(K)
        ENDDO
        DO K=0,locNZ+2
           SURZA(K)=1.d0
           CCZA(K)=0.d0
           VOLZA(K)=XZB(K)-XZB(K-1)
           SURZB(K)=1.d0
           CCZB(K)=0.d0
           geozga(k)=1.d0
        ENDDO
        DO K=-1,locNZ+2
           geozg(k)=1.d0
        ENDDO
      ELSEIF (NCOSYS.EQ.2) THEN
**  Spherical:  x<>r, y<>Phi, and  z<>Zeta
         DO I=0,locNX+1
! d(r^3)/3 = r^2*dr
            VOLXB(I)=(XXA(I+1)**3-XXA(I)**3)/3.
         ENDDO
         DO I=0,locNX+2
            SURXA(I)=XXA(I)*XXA(I)
            CCXA(I)=XXA(I)
            VOLXA(I)=(XXB(I)**3-XXB(I-1)**3)/3.
            SURXB(I)=XXB(I)*XXB(I)
            CCXB(I)=XXB(I)
            geoxga(I)=xxa(i)
            geoxha(I)=xxa(i)
         ENDDO         
         geoxga(1)=xxa1
         geoxha(1)=xxa1
         DO I=-1,locNX+2
            geoxg(I)=xxb(i)
            geoxh(I)=xxb(i)
         ENDDO
         DO J=0,locNY+1
            VOLYB(J)=XYA(J+1)-XYA(J)
         ENDDO
         DO J=0,locNY+2
            SURYA(J)=1.d0
            CCYA(J)=0.d0
            VOLYA(J)=XYB(J)-XYB(J-1)
            SURYB(J)=1.d0
            CCYB(J)=0.d0
         ENDDO
         DO K=0,locNZ+1
!-d(sin) = cos(theta)*dtheta
            VOLZB(K)=sin(XZA(K+1))-sin(XZA(K))
         ENDDO
         DO K=0,locNZ+2
            SURZA(K)=cos(xza(K))
            CCZA(K)=sin(xza(K))/2.d0
            VOLZA(K)=abs(sin(XZB(K))-sin(XZB(K-1)))
            if(VOLZA(K).eq.0.) then
               print *,'VOLUME TERM=0: gitter'
               call clean_stop
            endif
ci this is to make the volume term at the poles ne 0
c            if(k.eq.1.or.k.eq.locNZ+1) then
c               VOLZA(K)=abs(sin(XZB((locNZ/2)+1))-sin(XZB(locNZ/2)))
cci changed this 7/19/05    VOLZA(K)=VOLZA((locNZ/2)+1)
c            endif
            SURZB(K)=cos(xzb(K))
            CCZB(K)=sin(xzb(K))/2.d0
            geozga(k)=cos(xza(k))
         ENDDO
c         print *,'VOLUMES=',volza(locNZ+1),volza(locNZ),volza(locNZ-1)
         DO K=-1,locNZ+2
            geozg(k)=cos(xzb(k))
         ENDDO
*     special conditions so viscosities don't blow up at poles
         if (geozga(locNZ+1).eq.0.) geozga(locNZ+1)=1.e-10
         if (geozga(1).eq.0.) geozga(1)=1.e-10
      ELSE
         print *,'this coordinate system is not defined:gitter'
         call clean_stop
      ENDIF
      RETURN
      END SUBROUTINE GITTER

      SUBROUTINE GASINIT
!--- sets gas parameters
      USE input_init
      USE fluid_var_init
      USE global_constants
      IMPLICIT NONE
      xmue = mu_gas
      cv = Rgas/(mu_gas*(gamma-1.d0))
      RETURN
      END SUBROUTINE GASINIT

      SUBROUTINE RADINIT
      USE fluid_var_init
      IMPLICIT NONE
      STELLARINPUT = 0.d0
      RETURN
      END SUBROUTINE RADINIT

      SUBROUTINE ROTATIONINIT
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      IMPLICIT NONE
      if(MODTYP.NE.2) then !-- initialize the planet with zero velocity in the rotating frame
         if(myid.eq.0) then
            Print *,'Figure out about G for MODTYP=',MODTYP,
     %           ':ROTATIONINIT'
            call clean_stop
         endif
         print *, 'Subtracting OMROT from G'
         G=G-OMROT
      endif
      RETURN
      END SUBROUTINE ROTATIONINIT

      SUBROUTINE initEdep
      USE fluid_var_init
      IMPLICIT NONE
      EdepFor = 0.d0
      EdepAdv = 0.d0
      EdepRad = 0.d0
      INTERNAL1 = 0.d0
      INTERNAL1 = 0.d0
      RETURN
      END SUBROUTINE initEdep

      SUBROUTINE VISCINIT
      USE input_init
      USE force_var_init
      USE fluid_var_init
      IMPLICIT NONE
      if(IVIS.eq.0) then
         DIVTX = 0.d0
         DIVTY = 0.d0
         DIVTZ = 0.d0
      else
         print *,'IVIS=',ivis,' is not Avalible'
         call clean_stop
      endif
      RETURN
      END SUBROUTINE VISCINIT 

      SUBROUTINE KEPLARINIT
!--- ROUTINE TO CALCULATE THE KEPLARIAN VELOCITY AND ORBITAL FREQUENCY
      USE input_init
      USE fluid_var_init
      USE force_var_init
      USE global_constants
      USE grid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      integer :: stat(MPI_STATUS_SIZE)
      double precision :: P,Q,oldG
      DO I=-1,locNX+2
        OMEGAKEP(i) = ((GRAV*mstar1*MSOL/(xxb(i)**3.d0))**0.5)-OMROT
      ENDDO
!-- new run ( old run , just set omegakep)
      IF (MOD(IREA,10).EQ.0) THEN
         Maccreated = 0.d0
         IF(IGASTYP.EQ.0) THEN
            if(myid.eq.0) then
               print *,'Velocity initilized WITHOUT pressure support'
            endif
            DO K=-1,locNZ+2
               DO J=0,locNY+2
                  DO I=0,locNX+1
                     G(I,J,K) = OMEGAKEP(I)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            CALL BOUNDS
            CALL SOUND
            CALL PRESSG
            if(myid.eq.0) then
               do i=1,locNX
                  write(*,'(A,4(1x,e12.6))') 'test press:',
     %                 PG(i,1,1),RH(i,1,1),T(i,1,1),CS(i,1,1)
               enddo
            endif
            CALL QDENS
            CALL FGRAV
            if(myid.eq.0) then
               print *,'Velocity initilized WITH pressure support'
            endif
            K=1
            J=1
            DO I=1,locNX
               G(I,J,K) = ((  (1.d0/(RHQY(i,j,k)*xxb(i)))*
     %              ((PG(i,J,K)-PG(i-1,J,K))/(XXB(i)-XXB(i-1))) + 
     %              ((OMEGAKEP(I)+OMROT)**2.d0))**0.5d0) - OMROT
               if(myid.eq.0) then
                  write(*,'(I8,6(1x,e12.6))')
     %                 i,G(i,j,k),xxb(i),omegakep(i),
     %                 (PG(I,J,K)-PG(I-1,J,K))/
     %                 (xxb(i)*rhqy(i,j,k)*(XXB(I)-XXB(I-1))),pg(i,j,k),
     %                 (OMEGAKEP(I)+OMROT)**2.d0
               endif   
            ENDDO
c--   fill in the rest of the G array
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     G(I,J,K) = G(I,1,1)
                  ENDDO
               ENDDO
            ENDDO
            CALL BOUNDS
         ENDIF 
      ENDIF
      RETURN
      END SUBROUTINE KEPLARINIT

      SUBROUTINE READRESTART
      USE input_init
      USE global_var_init
      USE fluid_var_init
      USE mpi_var_init
      USE global_constants
      IMPLICIT NONE
      logical :: use_binary
      integer :: i,j,k,dummy_global_NX,dummy_global_NY,dummy_global_NZ
      integer :: dummy_locNX,dummy_locNY,dummy_locNZ,dummy_numprocs
      double precision :: dummy_XXB(global_NX)
      double precision :: dummy_XYB(global_NY)
      double precision :: dummy_XZB(global_NZ)
      double precision :: dummy_PG(global_NX,global_NY,global_NZ)
      double precision :: dummy_CS(global_NX,global_NY,global_NZ)
!-internal variable for reading in an ascii file (AKTTD_ASCII) rather then
!  the standard binary file (AKTTD)
      use_binary = .true.

      if(use_binary) then
         if(myid.eq.0) print *,'DOING A RESTART: READING AKTTD'
         OPEN (12,FILE='AKTTD',STATUS='OLD',FORM='UNFORMATTED')
         rewind(12)
         READ (12) dummy_global_NX
         READ (12) dummy_global_NY
         READ (12) dummy_global_NZ         
         READ (12) dummy_locNX,dummy_locNY,dummy_locNZ,dummy_numprocs
         READ (12) ZEIT,DELT,NUM_ITER
         READ (12) (dummy_XXB(I),I=1,global_NX)
         READ (12) (dummy_XYB(J),J=1,global_NY)
         READ (12) (dummy_XZB(K),K=1,global_NZ)
         READ (12) (((global_V(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12) (((global_G(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12) (((global_H(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12) (((global_RH(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12) (((global_T(I,J,K) ,I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12) (((global_ER(I,J,K) ,I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
!---  use this this read in old files
c     if(myid.eq.0) then
c     print *,'------------------------------------------------'
c     print *,'-----  SETTING ER=AT^4 on readRESTART-----------'
c     print *,'------------------------------------------------'
c     endif
c     global_ER = ARAD*((global_T)**(4.0))
         READ (12) (((dummy_PG(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12) (((dummy_CS(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12) ((global_upperbnd(j,k),J=1,global_NY),K=1,global_NZ)
         READ (12) TinitBottom,Tinittop,rhinittop,rhinitbottom,
     %        hydrotime,Maccreated
      ELSE
         if(myid.eq.0) print *,'DOING A RESTART: READING AKTTD_ASCII'
         OPEN (12,FILE='AKTTD_ASCII',STATUS='OLD')
         rewind(12)
         READ (12,*) dummy_global_NX
         READ (12,*) dummy_global_NY
         READ (12,*) dummy_global_NZ
         READ (12,*) dummy_locNX,dummy_locNY,dummy_locNZ,dummy_numprocs
         READ (12,*) ZEIT,DELT,NUM_ITER
         READ (12,*) (dummy_XXB(I),I=1,global_NX)
         READ (12,*) (dummy_XYB(J),J=1,global_NY)
         READ (12,*) (dummy_XZB(K),K=1,global_NZ)
         READ (12,*) (((global_V(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12,*) (((global_G(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12,*) (((global_H(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12,*) (((global_RH(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12,*) (((global_T(I,J,K) ,I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12,*) (((global_ER(I,J,K) ,I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12,*) (((dummy_PG(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12,*) (((dummy_CS(I,J,K),I=1,global_NX),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12,*) ((global_upperbnd(j,k),J=1,global_NY),
     %        K=1,global_NZ)
         READ (12,*) TinitBottom,Tinittop,rhinittop,rhinitbottom,
     %        hydrotime,Maccreated
      ENDIF
      CLOSE(12)
      if(myid.eq.0) then
         print *,'Restart Info:'
         print *,' global NX,NY,NZ=',dummy_global_NX,dummy_global_NY,
     %        dummy_global_NZ
         print *,' Starting ZEIT/NUM_ITER=',ZEIT,NUM_ITER
         print *,' Starting hydrotime=',hydrotime
         print *,' Tinitbottom/top',TinitBottom,Tinittop
         print *,' rhobottom/top',rhinitbottom,rhinittop
         print *,'---------------------------'
      endif
      RETURN
      END SUBROUTINE READRESTART
