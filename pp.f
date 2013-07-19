      PROGRAM convectingbox
!
!    PARALLEL COMPUTATION OF A CONVECTING BOX, ETC.
!
!    DATE     :  8.30.05
!     AUTHOR   :  IAN DOBBS-DIXON
!
!   MODTYP =  1  :  CARTESIAN BOX
!     MODVER = 0 :    CONVECTING BOX
!     MODVER = 1 :    SOD SHOCKWAVE TEST
!     MODVER = 2 :    EARTH CONVECTION
!   MODTYP =  2  :  IRRADIATED PLANET/BD
!     MODVER = 0 :    PETERS PLANET MODEL
!     MODVER = 1 :    PETERS 0.05M_SUN BROWNDWARF MODEL
!     MODVER = 2 :    PETERS PLANET MODEL WITH DEEPER INNER BOUNDARY
!     MODVER = 3 :    ADAM BURROWS MODEL/Hansen-neil
!     MODVER = 4 :    MY POLYTROPE MODELS (OR NEILS) MATCHED TO HANSEN SOLUTION
!     MODVER = 5 :    INITIAL MATCHED SEMI STRUC CALCULATION (WITH PLANET_NUM)
!     MODVER = 6 :    HELD & SUAREZ 1994 BENCHMARK CALCULATION
!      IBDRAD = 1 : fixed bottom dt/dr,periodic sides, fixed top
!      IBDRAD = 2-7: fixed bottom T flux (set dTdr and T), open sides
!      IBDRAD = 2: TJupHeat varies in phi-direction
!      IBDRAD = 3: TJupHeat varies in phi & theta directions, TMIN =100
!      IBDRAD = 4: TJupHeat is constant
!      IBDRAD = 5: TJupHeat varies in phi & theta directions and moves
!      IBDRAD = 6: TJupHeat varies in phi & theta directions and pulses in time
!      IBDRAD = 7: gradT = 0 boundary condition
!      IBDRAD = 8: TJupHeat varies in phi & theta directions, TMIN =950
!      IBDRAD = 9: TJupHeat=(ER/a+1/4*(gam*Tir)^4*exp(-tauA))^(0.25)
!       F_INCIDENT = 1: anti-symetric explicit stellar heating
!       F_INCIDENT = 2: anti-symetric implicit stellar heating
!       F_INCIDENT = 3: anti-symetric implicit stellar heating ramped up over damp_steps
!       F_INCIDENT = 4: sph. symetric implicit stellar heating
!       F_INCIDENT = 5: sph. symetric implicit stellar heating ramped up over damp_steps
!       F_INCIDENT = 6: eccentric implicit stellar heating
!       F_INCIDENT = 7: non-syncronized spin, zero ecc
!       F_INCIDENT = 8: ramp from symmetric to anti-symetric implicit stellar heating
!      planet_num=1 : HD209458b
!      planet_num=2 : HD189733b
!      planet_num=3 : Tres-1
!      planet_num=4 : HD149026
!      planet_num=5 : WASP-12
!      planet_num=6 : HD17156b
!      planet_num=7 : HAT-P-2
!      planet_num=8 : HAT-P-7
!      ISCAL=1: density testcase. SCAL=RH at init, read-restart not implemented
!      ISCAL=2: diff testcase. read-restart not implemented
!      ISCAL=3: XCO
!   MODTYP =  3  :  DISK WITH MASS ACCRETION
!     MODVER = 0 :    PARTIAL DISK
!     MODVER = 1 :    FULL DISK,2D
!     MODVER = 2 :    FULL DISK,3D
!   MODTYP =  4  :  INERTIAL DREDGE UP
!     MODVER = 0 :    NO WIND
!     MODVER = 1 :    WIND
!     MODVER = 2 :    WIND + NO-SLIP
!   MODTYP =  5  :  MIGRATION
!     MODVER = 0 :    STATIONARY PLANET
!     MODVER = 1 :    MOVING PLANET
!   MODTYP = 6   :   PLANETARY ACCRETION MODEL
!
!--- SET MODULES
      USE input_init
      USE input_alloc
      USE mpi_var_init
      USE mpi_var_alloc
      USE global_var_init
      USE global_var_alloc
      USE fluid_var_init
      USE fluid_var_alloc
      USE force_var_init
      USE force_var_alloc
      USE rad_var_init
      USE rad_var_alloc
      USE grid_var_init
      USE grid_var_alloc
c      USE particle_var_init
c      USE particle_var_alloc
      USE deltaER_var_init
      USE deltaER_var_alloc
      USE global_constants
      USE mpi_grid_init
      USE hypre_var_init
      USE hypre_var_alloc
      USE sor_var_init
      USE sor_var_alloc
      USE relax_var_init
      USE relax_var_alloc
      USE scalar_var_init
      USE scalar_var_alloc
      USE dif_var_init
      USE dif_var_alloc
      USE saumon_var_init
      USE saumon_var_alloc
      USE artvis_var_init
      USE artvis_var_alloc
      USE newton_heat_init
      USE newton_heat_alloc
      IMPLICIT NONE

!--- Include files
      include "mpif.h"
      include 'bHYPRE_SStructVariable.inc'

!--- VARIABLE DECLERATIONS
      integer :: i,j,k,FILE_RSRT,ioerror,s
      double precision :: debug1,debug2,debug3,debug4
      double precision :: debug5,debug6,debug7,debug8
      double precision :: a1,Pforce,eccentric_anomaly
      double precision :: period_start,period_ramp,rot_ramp,flux_ramp
      double precision :: vel_tmp1,vel_tmp2
      double precision :: G_tmp1,G_tmp2,H_tmp1,H_tmp2
      double precision :: rh_tmp1,rh_tmp2
      double precision :: T_tmp1,T_tmp2
      integer :: rot_STEPS      
      character proc_num*3
      character print_num*3
      character filename*12
      character stopfile*8
      character endname*9
      character tmpfname*11

!--- MPI INITILIZATION
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

!--- READ IN THE IPP CONTROL FILE
      call readinitfile
!--- ALLOCATE THE VARIOUS COMPONENTS
      call allocatempivar(MPI_STATUS_SIZE)
      call allocatefluidvar !THIS HAS TO BE HERE: locNX/NY DECLERATIONS
      call allocateforcevar
      call allocateradvar      
      call allocategridvar
      call allocateglobalvar
c      call allocateparticlevar
      call allocatedeltaERvar
      call allocatehyprevar
      call allocatesorvar
      call allocaterelaxvar
      if(iscal.gt.0.or.idif.gt.0) call allocatescalarvar
      if(idif.gt.0) call allocatedifvar
      if(modtyp.eq.2.and.modver.eq.5) call allocate_saumon_eos
      if(iartvis.gt.0) call allocateartvisvar
      if(irad.eq.4) call allocatenewtonheatvar

!--- HYPRE COMMUNICATOR
      if(irad.eq.2) then
         mpi_comm = MPI_COMM_WORLD
         call bHYPRE_MPICommunicator_CreateF_f( mpi_comm,bHYPRE_mpicomm,
     1        except)
      endif

!--- SET UP MPI GRID AND CARTESIAN_COMMUNICATOR
      CALL MPI_CART_CREATE(MPI_COMM_WORLD,proc_ndims,proc_dims,periods,
     %     reorder,COMM_CART,ierr)

!--- FUNCTION TO GET EVERYONE TOGETHER - USE ANYWHERE, BUT include "mpif.h"
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!--- FIND THE COODINATES OF PROCESSOR
      call MPI_CART_COORDS(COMM_CART,myid,proc_ndims,proc_coords,ierr)
c      PRINT *,MYID,proc_coords(1),proc_coords(2),proc_coords(3)
c      stop

!--- CREATE CHARACTER STRING OF ID NUMBER
      write(proc_num,'(i3.3)') myid

!--- FIND THE NEIGHBORS FOR EACH PROCESSOR
      IF(MODTYP.eq.1) then
         CALL SETUPNEW3DMPIGRID
      ELSEIF(MODTYP.eq.2) then
         call SETUP3DMPIGRID
      ELSEIF(MODTYP.eq.3.or.MODTYP.eq.4.or.MODTYP.eq.5) then
         call SETUP2DMPIGRID
      ELSEIF(MODTYP.eq.6) then
         call SETUP_AXISYMMETRIC
      ELSE
         print *,'incorrect MODTYP for setup of proc:pp.f'
         stop
      ENDIF


!--- SETUP THE HYPRE GRID AND STENCIL
      if(irad.eq.2) then
         call SETUPHYPREGRID
         call SETUPHYPRESTENCIL
         if(myid.eq.0) print *,'Done seting up hypre'
      endif

!-- set debug flag to 1
      DEBUGFLAG = 1

!-- set the iteration count for this run to 0 (for initilization)
      III=0

!--- GENERATE THE PHYSICAL GRID SYSTEM
      call GITTER

      if(myid.eq.0) write(*,'(A,3(1x,e12.3))') 'DXA,DYA,DZA=',
     %     DXA(1),DYA(1),DZA(1)

!--- CALCULATE SEMI-MAJORAXIS
      semi = (GRAV*MSOL*(mstar1+mstar2)*((Porb*DAY)**2.d0)/
     $     (4.0*PI*PI))**(1.d0/3.d0)
      IF(MODTYP.EQ.2.AND.MYID.EQ.0) THEN
         PRINT *,'semi-major axis (AU)',semi/AU
      ENDIF

!--- CALCULATE CURRENT ORBITAL RADIUS AND PHI0_ECC FOR FIRST FGRAV CALL
      IF(IGRAV.EQ.5.OR.IGRAV.EQ.6) THEN ! FOR ROCHE POTENTIAL
         if(f_incident.eq.6.and.hydrotime.gt.ecc_start*DAY) then
            ecc_anomaly = ECCENTRIC_ANOMALY(hydrotime-ecc_start*DAY,
     $           ecc,Porb*DAY)
            radius_orb = semi*(1.d0-ecc*cos(ecc_anomaly))
         else
            radius_orb = semi
         endif
         if((f_incident.eq.6.or.f_incident.eq.7).and.
     1        hydrotime.gt.ecc_start*DAY) then
            Pforce = ((1.d0/Porb)-(1.d0/Prot))**(-1.d0)
            phi0_ecc = (2.d0*PI/(Pforce*DAY))*
     %           mod((hydrotime-ecc_start*DAY),Pforce*day)
         else
            phi0_ecc = 0.d0
         endif
         if(myid.eq.0) then
            print *,'initial radius (AU),phi0_ecc=',radius_orb/AU,
     %           phi0_ecc
         endif
      ENDIF
      
!-printout maximum irradiation temperature
      IF(MODTYP.EQ.2.AND.MYID.EQ.0) THEN
         print *,'Maximum irradiation temperature=',
     $        Tstar*sqrt(RSTAR*RSOL/semi)
      ENDIF

!--- INITILIZE THE CV, AND MU OF GAS
      CALL GASINIT

!-- CALCULATE THE ANGLE FROM THE NORMAL TO THE STAR (mu_gas)
      CALL CALCMUSTAR

!------ READ IN AN OPACITY TABLE IF NEEDED
      if(myid.eq.0) print *,'ikapyp =',ikaptyp
      if(myid.eq.0) print *,'zkapmult =',zkapmult
      if(irad.gt.0.and.irad.ne.4.and.ikaptyp.gt.0) then
         SELECT CASE(ikaptyp)
         CASE (1)
            xkapR=ZKAPMULT
         CASE(2)
            CALL READALEXANDEROPCTABLE
         CASE(4:5)
            CALL GASINIT !-- I need mu
            CALL READFREEDMANOPCTABLE
         CASE(6:7)
            CALL READFREEDMANOPCTABLE
            CALL READCOLUMNDEPTHOPCTABLE
         CASE(9:12)
            CALL READBURROWSOPCTABLE
         CASE(14)
            CALL READBURROWSOPCTABLE
         CASE(15:19)
            CALL READELIZAOPCTABLE
         CASE(20:25)
            CALL READBURROWSOPCTABLE
         CASE(26:)
            print *,'ikaptyp=',ikaptyp,'is not avalible:pp.f'
            call clean_stop
         CASE DEFAULT
            print *,'also not defined: KAPPA (pp.f)',ikaptyp
            call clean_stop
         END SELECT
      endif


      if(myid.eq.0) then
         print *,'---------------------------'
         print *,'Simulation Info:'
         write(*,'(A,3(I7))') '  global NX,NY,NZ=',
     &        global_NX,global_NY,global_NZ
         write(*,'(A,3(I7))') '  locNX/NY/NZ=',locNX,locNY,locNZ
         write(*,'(A,I8)') '  numproc=',numprocs
         print *,'---------------------------'
      endif
      


!--- INFORMATION FOR A NEW RUN
      IF (MOD(IREA,10).EQ.0) THEN
         num_iter = 0
         DELT=DTMAX
         zeit=0.
!------ INITILIZE GLOBAL_ ON ALL PROCESSORS
         CALL GLOBALINIT
c         print *,'second2 ',proc_num,' ',myid
!------ INITILIZE VARIABLES LOCALLY
         CALL LOCALINIT
c         print *,'second2B ',proc_num,' ',myid
         write(proc_num,'(i3.3)') myid
c         print *,'second2B2 ',proc_num,' ',myid
!------ FOR MODTYP=2 AND MODVER=5 ITERATE TO INITIAL STRUCTURE
!-----    USING LOADED MODEL AS A GUESS
         if(modtyp.eq.2.and.modver.eq.5) then
            CALL CALC_INITIAL_MODEL
         endif
!------ ZERO ENERGY DEPOSITION ARRAYS
         CALL initEdep
!------ INITILIZE THE PASSIVE SCALAR(S)
         if(iscal.gt.0) then 
            CALL SCALARINIT
         endif
!------ ZERO HYDRODYNAMICAL TIME
         hydrotime = 0.d0
!------ INITILIZE THE TRACER POSITIONS
c         if(nparticles.gt.0) then
c            if(myid.eq.0) print *,'numparticle=',nparticles
c            CALL INITPARTICLES
c            if(myid.eq.0) print *,'PARTICLES INIT DONE'
c         endif
         
!--- READ IN INFORMATION FOR A RESTART (AKTTD)
      ELSEIF (mod(IREA,10).EQ.1) THEN
         CALL READRESTART
!------ INITILIZE VARIABLES LOCALLY
c         print *,'third ',proc_num,' ',myid
         CALL LOCALINIT
c         print *,'third2 ',proc_num,' ',myid
         write(proc_num,'(i3.3)') myid
c         print *,'third3 ',proc_num,' ',myid
!------ FOR MODTYP=2 AND MODVER=5 ITERATE TO INITIAL STRUCTURE
!-----    USING LOADED MODEL AS A GUESS
         if(modtyp.eq.2.and.modver.eq.5) then
            CALL CALC_INITIAL_MODEL
         endif

!------ INITILIZE THE PASSIVE SCALAR(S)
         if(iscal.gt.0) CALL SCALARINIT

!-Replace above line with the following to start tracers for the
!- first time from a previously run purely hydro model
c         if(iscal.eq.3) then
c            if(icstyp.eq.1.or.icstyp.eq.2) CALL SOUND
c            CALL PRESSG
c            CALL CO_SCALAR_init
c            print *,'calling calc_carbon_mass'
c            call CALC_CARBON_MASS
c         endif

!------ READ IN LAST TRACER FILE (for now INITILIZE THE TRACER POSITIONS)
c         if(nparticles.gt.0) then
c            if(myid.eq.0) print *,'numparticle=',nparticles
c            CALL INITPARTICLES
c            if(myid.eq.0) print *,'PARTICLES RESTART DONE'
c         endif

      ELSE
         print *,'irea',irea,' is undefined'
         STOP
      END IF


!--- OPEN M-ACCRETION FILE
      IF(MODTYP.EQ.3) THEN
         IF (MOD(IREA,10).EQ.0) THEN
            OPEN(20,FILE='MACC.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         ELSE
            OPEN(20,FILE='MACC.dat',STATUS='OLD',FORM='FORMATTED',
     %           POSITION='APPEND')
         ENDIF
      ENDIF
      
!--- OPEN PLANT LOCATION FILE
      IF(MODTYP.EQ.5) THEN
         if(myid.eq.0) print *,'FILE1'
         IF (MOD(IREA,10).EQ.0) THEN
            if(myid.eq.0) print *,'FILE2'
            OPEN(30,FILE='PLANET.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         ELSE
            OPEN(30,FILE='PLANET.dat',STATUS='OLD',FORM='FORMATTED',
     %           POSITION='APPEND')
         END IF
      ENDIF

!--- OPEN MASS TRANSFER FILE
      IF(MODTYP.EQ.2.and.((MODVER.eq.1).or.
     %     (modver.eq.5.and.planet_num.eq.5).or.
     %     (modver.eq.4.and.planet_num.eq.5))) THEN
         OPEN(40,FILE='MASSTRANFER.dat',STATUS='UNKNOWN',
     %        FORM='FORMATTED')
!-- zero the mass transfer
         Mtransfer = 0.d0
      ENDIF

!--- TEST THE PARTICLE FINDING ROUTINE
c      if(myid.eq.0) print *,'STARTING PARTICLE INTERPOLATION'
c      CALL INTERPOLATEVELOCITIES

!------FILL OPACITY ARRAYS FOR FIRST PRINTOUT
      if(irad.gt.0.and.irad.ne.4.and.ikaptyp.gt.0) then
         CALL KAPPA
      endif

!--- SET THE PRINT COUNTER. ZERO OR FILESTART
!------- FILESTART should be the last GLOB/SNAP file in previous run
      IF (MOD(IREA,10).EQ.0) THEN !- new run
         num_print = 0
         FILE_RSRT = 0
      ELSE                      !- restart
         num_print = FILESTART
         FILE_RSRT = 1
      END IF

!--- CALCULATE KEPLARIAN VELOCITY
      IF(MODTYP.EQ.3.OR.MODTYP.EQ.5) CALL KEPLARINIT
!--- FIND BOUNDARY CONDITIONS
      CALL BOUNDS
!--- INITILIZE THE ANGULAR VELOCITY FOR NEW RUNS(BUT NOT FOR MODTYP=3 & 5)
      IF (MOD(IREA,10).EQ.0) THEN
         IF (IROTATE.NE.0.and.MODTYP.NE.3.and.MODTYP.NE.5) then
            CALL ROTATIONINIT
         ENDIF
      endif
!--- PRINT OUT THE ROTATION FREQUENCY
      if(myid.eq.0) write(*,'(A,e12.3)') 'OMROT=',OMROT
!--- PRINT OUT HEATING VERSION
      if(myid.eq.0) then
         write(*,'(A,I5)') 'F_INCIDENT=',F_INCIDENT
         if(F_INCIDENT.eq.6.or.F_INCIDENT.eq.7) 
     $        write(*,'(1x,A,4(1x,e12.3))')
     $        'Po(d),Pr(d),ecc,ecc_start=',Porb,Prot,ecc,ecc_start
      endif
!--- PRINT OUT THE ARTIFICIAL VISCOSITY COEFFICIENT AND INITILIZE DT_ARTVIS
      if(iartvis.gt.0) then 
         if(myid.eq.0) write(*,'(A,e12.3)') 'Artvisc_coeff= ',
     %        artvisc_coeff
         dt_artvis = DTMAX
      endif
!---  INITILIZE CS (BEFORE 1st PRESSG/SCHRITT CALL FOR (LOCALLY) ISOTHERMAL EOS)
      if(icstyp.eq.1.or.icstyp.eq.2) CALL SOUND
!---  CALCULATE PRESSURE
      CALL PRESSG
!---  INITILIZE CS (GENERAL CALL)
      CALL SOUND
!---  ZERO VISCOUS FORCE AND DISSIPATION IF NOT NEEDED
      IF(IVIS.EQ.0) THEN
         DIVTX = 0.d0
         DIVTY = 0.d0
         DIVTZ = 0.d0
      ENDIF
      IF(IDISS.EQ.0) THEN
         DISFN = 0.d0
      ENDIF
!---  CALCULATE DENSITY AT MIDPOINTS
      CALL QDENS
!---  CALCULATE THE TIMESTEP
      CALL SCHRITT
!--- CALCULATE GRAVITY FORCE
      CALL FGRAV
!--- FILL THE RADIATION COEFFICIENTS FOR THE FIRST PRINTOUT
      if(irad.gt.0.and.irad.ne.7) CALL DIFF

!--- READ IN RELAXATION FILE      
c      CALL READRELAX

!--- CALL THE MISC_CHECK FILES TO CHECK CONSITANCY BETWEEN VARIABLES
      CALL MISC_CHECKS

!------------------------------------------------------------
!---  START OF ITERATIONS
!------------------------------------------------------------
      startingtime = ZEIT
      if (myid.eq.0.and.iprstep.gt.0) then
         CALL SCREEN_TIME_PRINTOUT
      endif

!--- FIND MAX VELOCITIES INITIALLY
      CALL MAXVELOCITIES


      DO III=1,ISTEPMAX
c         if (myid.eq.0) print *,'step number=',III,delt,'------------'
         
!------ DAMP VELOCITIES AT THE POLES ONLY (USE INC_POLES=2)
         if(inc_poles.eq.2) 
     $        CALL POLEVELOCITYDAMPING(III,locNX,locNY,locNZ,V,G,H,
     $        hydrotime)

!---  ARTIFICAL GAS COOLING ROUTINE
c         CALL GASCOOL

!--- SLOW RAMP UP OF ROTATION
         if(irotate.eq.3) then
            if(hydrotime.lt.rot_start) then
               OMROT = 0.d0
            elseif (hydrotime.ge.rot_start.and.
     %              hydrotime.lt.rot_full) then
               OMROT = (2.d0*PI/(PROT*DAY))*
     %              ((hydrotime-rot_start)/(rot_full-rot_start))
            elseif(hydrotime.ge.rot_full) then
               OMROT = 2.d0*PI/(PROT*DAY)
            endif
         endif

!--- GRADUAL CHANGE IN ROTATION (MODIFY period_start here AND in localinit)
         if(irotate.eq.4) then
            rot_STEPS = 400000
            rot_ramp = (1.d0*III/(1.d0*rot_STEPS))
            if(rot_ramp.ge.1.0) rot_ramp = 1.d0
            period_START = 3.52
            period_ramp = period_start - (period_start-PROT)*rot_ramp
            OMROT = 2.d0*PI/(period_ramp*DAY)
            IF(((III-1)/(IREA/10))*(IREA/10).eq.(III-1).and.
     $           myid.eq.0) THEN
               write(*,'(A,2(1x,e12.6),1x,I8)') 
     $              'rotation ramp %=',period_ramp,OMROT,III
            ENDIF
         endif

!--- GRADUAL CHANGE IN FLUX FROM INTERIOR
         if(iFbotramp.eq.1) then
            flux_ramp = (1.d0*III/(1.d0*Framp_steps))
            if(flux_ramp.ge.1.0) flux_ramp = 1.d0
c            Fbottom = Fbot - (Fbot-Fbot_final)*flux_ramp
            Fbottom = 10.d0**(LOG10(Fbot) - 
     %           (LOG10(Fbot)-LOG10(Fbot_final))*flux_ramp)

            IF(((III-1)/(IREA/10))*(IREA/10).eq.(III-1).and.
     $        myid.eq.0) THEN
               write(*,'(A,2(1x,e12.6),1x,I8)') 
     $              'interior flux ramp %=',flux_ramp,Fbottom,III
            ENDIF
         endif

!--- PERURBATIONS AZIMUTHAL VELOCITY AT PHOTOSPHERE
         if(iperturb.eq.1.and.III.eq.perturb_iter) then
            if(myid.eq.0) print *,'Perturbing velocity. iii=',iii
            CALL PERTURBATIONS
         endif         

!---- LIMIT ON VELOCITIES
         if (abs(IVEL).gt.0) then
            CALL VARLIMIT (IVEL,V,-VELMAX,VELMAX)
            if(ivel.lt.3) then
               if(ncosys.eq.0) then
                  CALL VARLIMIT (IVEL,G,-VELMAX,VELMAX)
                  CALL VARLIMIT (IVEL,H,-VELMAX,VELMAX)
               elseif(ncosys.eq.2) then
                  CALL SPHERICALVELLIMIT
               else
                  print *,'ncosys not defined:pp.f'
                  stop
               endif 
            endif
         endif

!---- LIMIT ON TEMPERATURE
         IF (ABS(ITEMP).GT.0) THEN
            CALL VARLIMIT(itemp,T,tempmin,tempmax)
         ENDIF

!---- PRINTOUT SNAP FILE
         IF (MOD(IREA,10).LE.1) THEN
            IF ((IREA/10).ne.0 .and. 
     &           ((III-1)/(IREA/10))*(IREA/10).eq.(III-1)) THEN
               write(print_num,'(i3.3)') num_print
               filename = 'SNAP_'//TRIM(ADJUSTL(proc_num))//'.'//
     $              print_num
               IF (FILE_RSRT.EQ.1) THEN !-designate first restart file
                  filename = 'RSRT_'//TRIM(ADJUSTL(proc_num))//'.'//
     $                 print_num
                  FILE_RSRT = 0
               ENDIF
               OPEN (12,FILE=filename,STATUS='UNKNOWN',
     $              FORM='UNFORMATTED')
c               OPEN (12,FILE=filename,STATUS='UNKNOWN',
c     $              FORM='UNFORMATTED',CONVERT='LITTLE_ENDIAN')
               CALL BOUNDS
               REWIND(12)
               CALL PRINTOUT(1)
               CLOSE(12)
!---   scalar printout routine
               if(iscal.gt.0) then
                  filename = 'SCAL_'//TRIM(ADJUSTL(proc_num))//'.'//
     $                 print_num
                  OPEN (12,FILE=filename,STATUS='UNKNOWN',
     $                 FORM='UNFORMATTED')
                  CALL PRINTSCALAR
                  CLOSE(12)
               endif
               num_print=num_print+1
            ENDIF
!---- END OF SNAP PRINTOUT PROCEDURES
!---- PRINTOUT INFO ON SETING TIMESTEP
c            CALL CHECK_MAX_DT
         ENDIF
         
!---- WRITE OUT MASS ACCRETION FILE
         IF(MODTYP.EQ.3) THEN
            IF(MOD(num_iter,1000).EQ.0) THEN
               call CALCGAPMASS
               if(myid.eq.0) write(20,'(I8,1x,3(1x,e12.6))') 
     %              num_iter,zeit,Maccreated,Mgap
            ENDIF
         ENDIF

!---- WRITE OUT MASS ACCRETION FILE
c         IF(MODTYP.EQ.5) THEN
c            IF(MOD(num_iter,100).EQ.0) THEN
c               if(myid.eq.0) write(30,'(I8,1x,3(1x,e12.6))') 
c     %              num_iter,zeit,rp,phip
c            ENDIF
c         ENDIF

!---- WRITE OUT MASS TRANSFER FILE
c         IF(MODTYP.EQ.2.and.((MODVER.eq.1).or.
c     %        (modver.eq.5.and.planet_num.eq.5).or.
c     %        (modver.eq.4.and.planet_num.eq.5)).and.
c     %        MOD(num_iter,100).EQ.0) THEN            
         IF(MODTYP.EQ.2.and.modver.eq.4.and.planet_num.eq.5.and.
     $        MOD(num_iter,100).EQ.0) THEN
            if(myid.eq.0) then
               write(40,'(I8,1x,2(1x,e12.6))') 
     %              num_iter,zeit,Mtransfer
               write(*,'(A,I8,1x,2(1x,e12.6))') 'III,time,Mtransfer',
     %              num_iter,zeit,Mtransfer
!--   rezero mtransfer
            endif
            Mtransfer = 0.d0
         ENDIF

!--- ADVANCE TIME AND ITERATION COUNTS
         num_iter=num_iter+1
         ZEIT=ZEIT+DELT
         if(IADV.gt.0.0) then
            hydrotime=hydrotime+DELT
         endif

c         if(myid.eq.0) print *,'bf adv'
c         print *,'bf adv',myid

c         if (myid.eq.19) 
c     %        write(*,'(A,I8,4(1x,e12.6))') ,'BF ADV',III,
c     %        rh(97,6,6),rh(98,6,6),rh(99,6,6),rh(100,6,6),
c     %        V(99,0,0)

**  Advection
         CALL CALCINTERNALENERGY(INTERNAL1)
         IF (IADV.EQ.1) THEN
            if(iscal.gt.0) then !-- do scalar first so it sees same velocities as rh
               CALL SCALAR_BOUNDS !- Boundary Conditions
               do s=1,nscalar
                  pass_sc(:,:,:) = scalar(s,:,:,:)
                  CALL ADVECTSCALAR3D(s)
                  scalar(s,:,:,:) = pass_sc(:,:,:)
               enddo
               CALL SCALAR_BOUNDS !- Boundary Conditions
               IF(iscal.eq.3) then
                  RH_PRE_ADV = RH
               endif
            endif
            if(locNZ.eq.1) then
               CALL TESTADVECT
               !CALL ADVECT2D
            elseif(locNZ.ne.1) then
               CALL ADVECT3D
            else
               print *,'something is incorrect:pp.f - advect'
               call clean_stop
            endif
            IF (ABS(IRHO).GT.0) THEN
               CALL VARLIMIT(irho,RH,rhmin,1.d10)
            ENDIF
            IF (ABS(ITEMP).GT.0) THEN
               CALL VARLIMIT(itemp,T,tempmin,tempmax)
            ENDIF
!- Make adjustment for scalar being a mass fraction rather then a density
            IF(iscal.eq.3) THEN
               do s=1,nscalar
                  pass_sc(:,:,:) = scalar(s,:,:,:)
                  CALL SCALAR_ADJUST_XCO(s)
                  scalar(s,:,:,:) = pass_sc(:,:,:)
               enddo
               CALL SCALAR_BOUNDS !- Boundary Conditions
            ENDIF
         ENDIF
         CALL CALCINTERNALENERGY(INTERNAL2)
         EdepAdv = INTERNAL2-INTERNAL1
            
c         if(myid.eq.0) print *,'bf for'
c         print *,'bf for',myid
c         if (myid.eq.19) 
c     %        write(*,'(A,I8,4(1x,e12.6))') ,'BF FOR',III,
c     %        rh(97,6,6),rh(98,6,6),rh(99,6,6),rh(100,6,6)

**  Force
         IF (IFOR.EQ.1) THEN
            IF(locNZ.EQ.1) THEN
               CALL FORGAS2D
            ELSEIF(locNZ.NE.1) THEN
               CALL FORGAS3D
            ELSE
               print *,'something is incorrect:pp.f - force'
               stop
            ENDIF
            IF (ABS(ITEMP).GT.0) THEN
               CALL VARLIMIT(itemp,T,tempmin,tempmax)
            ENDIF 
         ENDIF

**  Explicit Stellar Heating
         if(F_INCIDENT.eq.1.and.FDELT.gt.0.0) then
            CALL STELLARHEATING
         endif
         CALL CALCINTERNALENERGY(INTERNAL1)
         EdepFor = INTERNAL1 - INTERNAL2

c         if(myid.eq.0) print *,'bf rad'
c         print *,'bf rad',myid

c         if (myid.eq.19) 
c     %        write(*,'(A,I8,4(1x,e12.6))') ,'BF RAD',III,
c     %        rh(97,6,6),rh(98,6,6),rh(99,6,6),rh(100,6,6)

**  Radiation
         IF (IRAD.EQ.1) THEN
            CALL RADDIF
         ELSEIF(IRAD.eq.2.or.IRAD.eq.3) THEN
            CALL COUPLED_ENERGY_EQS
            CALL BOUNDS
c            CALL LAPLACETEST !-simple test routine
         ELSEIF(IRAD.eq.4) THEN
            CALL NEWTONINAN_HEATING
         ELSEIF(IRAD.eq.5.or.IRAD.eq.6) THEN
            CALL WAVELENGTH_RADIAL_HEATING
         ELSEIF(IRAD.EQ.7) THEN
            CALL GREY_TWOSTREAM
         ELSEIF(IRAD.GE.8) THEN
            print *,'IRAD=',IRAD,'not avalible'
            stop
         END IF
         CALL CALCINTERNALENERGY(INTERNAL2)
         EdepRad = INTERNAL2-INTERNAL1
c         if(myid.eq.0) print *,'af rad'


**  CO relaxation scheme
         if(iscal.eq.3) then
            if(icstyp.eq.1.or.icstyp.eq.2) CALL SOUND
            CALL PRESSG
            call CO_relaxation
            CALL SCALAR_BOUNDS
         endif

**  Diffusion
         if(idif.eq.1) then
            CALL SCALAR_BOUNDS  !- Boundary Conditions
            do s=1,nscalar
               pass_sc(:,:,:) = scalar(s,:,:,:)
               CALL DIFFUSION(s)
               scalar(s,:,:,:) = pass_sc(:,:,:)
            enddo
            CALL SCALAR_BOUNDS  !- Boundary Conditions
         endif
         if(idif.eq.2) then
            CALL BOUNDS
            pass_sc(:,:,:)= T(:,:,:)
            CALL DIFFUSION(s)
            T(:,:,:)= pass_sc(:,:,:)
            CALL BOUNDS
         endif

**  Mass Sink
         IF(MODTYP.eq.3.and.mstar2.gt.0.0) then
            CALL PLANETMASSSINK
         endif

**  set the boundary values, calc PG & maybe T (for printout only) 
         CALL BOUNDS
         CALL PRESSG
         CALL SCHRITT

**    SCREEN PRINTOUT
         if((iprstep.gt.0).and.
     %        ((num_iter/iprstep)*iprstep.eq.num_iter)) then
            if(myid.eq.0) then
               CALL SCREEN_TIME_PRINTOUT
               if(irad.eq.3) CALL MATRIX_ITERATIONS_PRINTOUT                  
               if(irotate.eq.3) then
                  print *,'omrot=',omrot
               endif

            endif
            call MAXVELOCITIES
         endif

**     ARTIFICIAL STOPPING ROUTINE (DUMP INFO AND MPI FINALIZE)
c         open(666,file="STOP",status="old",iostat=ioerror)
c         if(ioerror.eq.0) then
c            close(666)
c            print *,'ARTIFICALLY STOPPING PROGRAM. MYID=',myid
c            stopfile = 'STOP_'//TRIM(ADJUSTL(proc_num))
c            OPEN (12,FILE=stopfile,STATUS='UNKNOWN',
c     $                 FORM='UNFORMATTED')
c            REWIND(12)
c            CALL PRINTOUT(1)
c            CLOSE(12)
c            if(irad.eq.2) then
c               call bHYPRE_SStructGrid_deleteRef_f( grid, except )
c               call bHYPRE_SStructGraph_deleteRef_f( graph, except )
c               call bHYPRE_SStructStencil_deleteRef_f( stencil, except )
c               call bHYPRE_MPICommunicator_deleteRef_f( bHYPRE_mpicomm, 
c     $              except )
c            endif
c            call MPI_FINALIZE(ierr)
c            stop
c         endif

      ENDDO
**  END OF ITERATIONS

c      CLOSE(321)

**  FINAL BINARY DATA --> FILE 'AKTTD'
      endname = 'AKTTD_'//TRIM(ADJUSTL(proc_num))
      OPEN(12,FILE=endname,STATUS='UNKNOWN',FORM='UNFORMATTED')
      rewind(12)
      CALL PRINTOUT(1)
      CLOSE (12)

      if(myid.eq.0) then
         print *,'Total Simulated Time=     ',zeit
         print *,'Total Simulated Hydrotime=',hydrotime
      endif

      IF(MODTYP.EQ.3) CLOSE(20)
      IF(MODTYP.EQ.5) CLOSE(30)
c      IF(MODTYP.EQ.2.and.MODVER.eq.1) CLOSE(40)
      IF(MODTYP.EQ.2.and.((MODVER.eq.1).or.
     %     (modver.eq.5.and.planet_num.eq.5).or.
     %     (modver.eq.4.and.planet_num.eq.5))) CLOSE(40)
!-----WRAP UP HYPRE
      if(irad.eq.2) then
         call bHYPRE_SStructGrid_deleteRef_f( grid, except )
         call bHYPRE_SStructGraph_deleteRef_f( graph, except )
         call bHYPRE_SStructStencil_deleteRef_f( stencil, except )
         call bHYPRE_MPICommunicator_deleteRef_f(bHYPRE_mpicomm,except)
      endif
!----- WRAP UP MPI
      call MPI_FINALIZE(ierr)
      STOP
      END PROGRAM convectingbox
