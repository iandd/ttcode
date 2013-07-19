!---- INPUT VARIABLES
      MODULE input_init
      IMPLICIT NONE
      SAVE
      integer :: ISTEPMAX,IREA,iprstep,FILESTART,ietot,isoth,isoden
      integer :: IADV,IFOR,IGRAV,IVIS,IRAD,INUC,ISCAL,IDIF
      integer :: IRHO,ITEMP,IVEL,IVISTYP,IDISS,IROTATE,IETA,IARTVIS
      integer :: IGASTYP,ICSTYP,IKAPTYP,ISIGTYP,DAMP_STEPS
      integer :: IBDRAD,IFLIMTYP,F_INCIDENT,BCXMIN
      integer :: NCOSYS,MODTYP,MODVER,PLANET_NUM
      integer :: ITMAX,iperturb,perturb_typ,perturb_iter
      integer :: NUMDIM,global_NX,global_NY,global_NZ
      integer :: nparticles=3,iFbotramp,Framp_steps
      integer :: inc_poles
      double precision :: FDELT,DTMAX
      double precision :: RHMIN
      double precision :: TEMPMAX,TEMPMIN
      double precision :: VELMAX,VELDAMP
      double precision :: XNUE,ALPHA,GAMMA,mu_gas,aspect
      DOUBLE PRECISION :: Pr_num
      double precision :: ZKAPMULT,gamma_OPC
      double precision :: PROT,ramp_start,ramp_full,windAMP
      double precision :: SORPARAM,EPSMAT
      double precision :: Porb,ecc,ecc_start
      double precision :: STARTYP,RSTAR,TSTAR,MSTAR1,TIRR,RH0
      double precision :: rstar2,mstar2,rp,phip,psoft,Fbot,Fbot_final
      double precision :: XMIN,XMAX,DXMIN
      double precision :: YMIN,YMAX,DYMIN
      double precision :: ZMIN,ZMAX,DZMIN
      double precision :: dbl_tmp(6)
      integer :: int_tmp(2)
      logical :: poles
      END MODULE input_init

      MODULE input_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE readinitfile
          USE input_init
          IMPLICIT NONE
          integer :: i
          OPEN(15,FILE='ipp',FORM='FORMATTED',STATUS='OLD')
          rewind(15)
          READ(15,*) ISTEPMAX,IREA,FDELT,DTMAX
          READ(15,*) iprstep,FILESTART
          READ(15,*) ietot,isoth,isoden
          READ(15,*) IADV,IFOR,IGRAV,IVIS,IRAD,INUC,ISCAL,IDIF
          READ(15,*) IRHO,RHMIN
          READ(15,*) ITEMP,TEMPMAX,TEMPMIN
          READ(15,*) IVEL,VELMAX,VELDAMP,DAMP_STEPS
          READ(15,*) IVISTYP,XNUE,IDISS,ALPHA,IETA,IARTVIS,PR_NUM
          READ(15,*) IGASTYP,GAMMA,mu_gas,ICSTYP,aspect
          READ(15,*) IKAPTYP,ISIGTYP,ZKAPMULT,gamma_OPC
          READ(15,*) IROTATE,PROT,ramp_start,ramp_full,windAMP
          READ(15,*) IBDRAD,IFLIMTYP,F_INCIDENT,BCXMIN,Fbot
          READ(15,*) SORPARAM,ITMAX,EPSMAT
          READ(15,*) !- More Switches (some open)
          READ(15,*) iperturb,perturb_typ,perturb_iter
          READ(15,*) Porb,ecc,ecc_start
          READ(15,*) iFbotramp,Framp_steps,Fbot_final
          READ(15,*) int_tmp(1),dbl_tmp(3),dbl_tmp(2),dbl_tmp(1)
          READ(15,*) !- Fixed Modelparameters
          READ(15,*) NCOSYS,inc_poles
          READ(15,*) MODTYP,MODVER,PLANET_NUM
          READ(15,*) STARTYP,RSTAR,TSTAR,MSTAR1,TIRR,RH0
          READ(15,*) rstar2,mstar2,rp,phip,psoft
          READ(15,*) NUMDIM
          READ(15,*) global_NX,XMIN,XMAX,DXMIN
          READ(15,*) global_NY,YMIN,YMAX,DYMIN
          READ(15,*) global_NZ,ZMIN,ZMAX,DZMIN
          CLOSE(15)
!--set the logic variable for including poles
          if(inc_poles.ge.1) then 
             POLES = .TRUE.
          else 
             POLES = .FALSE.
          endif
      END SUBROUTINE readinitfile
      END MODULE input_alloc
!----------------------
!
!----------------------
!
!---- MPI GRID VARIABLES
      MODULE mpi_grid_init
      IMPLICIT NONE
      SAVE
      integer :: MPIupper,MPIlower
      integer :: MPIright,MPIleft
      integer :: MPIabove,MPIbelow
!-- y-x plane
      integer :: MPIupperright,MPIupperleft
      integer :: MPIlowerright,MPIlowerleft
!-- y-z plane
      integer :: MPIupperabove,MPIupperbelow
      integer :: MPIlowerabove,MPIlowerbelow
!-- z-x plane
      integer :: MPIaboveright,MPIbelowright
      integer :: MPIaboveleft,MPIbelowleft 
!-- right corners
      integer :: MPIupperaboveright,MPIupperbelowright
      integer :: MPIloweraboveright,MPIlowerbelowright
!-- left corners
      integer :: MPIupperaboveleft,MPIupperbelowleft
      integer :: MPIloweraboveleft,MPIlowerbelowleft
!-- polar variables
      integer :: MPIabove_pole,MPIbelow_pole
      END MODULE mpi_grid_init
!----------------------
!
!---- MPI VARIABLES
      MODULE mpi_var_init
      IMPLICIT NONE
      SAVE
      integer :: proc_ndims,root=0
      integer,ALLOCATABLE,DIMENSION(:) :: proc_dims,proc_coords
      integer,ALLOCATABLE,DIMENSION(:) :: neig_coord
      logical,ALLOCATABLE,DIMENSION(:) :: periods
      integer,ALLOCATABLE,DIMENSION(:) :: proc_stat
      integer :: ierr,myid,numprocs,irc
      integer :: comm_cart
      logical :: reorder=.false.
      integer :: tag,yourrank
      integer :: dest,sendtag,source,recvtag
      integer :: shift_xpos_source,shift_xpos_dest
      integer :: shift_ypos_source,shift_ypos_dest
      integer :: shift_zpos_source,shift_zpos_dest
      END MODULE mpi_var_init

      MODULE mpi_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatempivar(MPI_STATUS_SIZE)
         USE mpi_var_init
         USE input_init
         IMPLICIT NONE
         integer i,MPI_STATUS_SIZE
         proc_ndims = numdim
         ALLOCATE(proc_dims(proc_ndims),proc_coords(proc_ndims))
         ALLOCATE(periods(proc_ndims))
         ALLOCATE(proc_stat(MPI_STATUS_SIZE))         
         if(mod(numprocs,4).ne.0) then
            if(myid.eq.0) then
               print *,'There may be a diff. between the number of'
               print *,'  proc. you think your using and the actual #'
               print *,'  This has to do with columbia assigning'
               print *,'  processors in blocks of 4'
               print *,'  Manually set $NCPUS in submission script'
               print *,'  numproc=',numprocs,myid
            endif
         endif
         SELECT CASE(MODTYP)
         CASE (1)
            proc_dims(1) = 1
            proc_dims(2) = INT((numprocs)**(1.d0/2.d0))
            proc_dims(3) = INT((numprocs)**(1.d0/2.d0))
            periods(1) = .false.
            periods(2) = .true.
            periods(3) = .true.
            if(myid.eq.0) then
               print *, 'info:',proc_dims(1),proc_dims(2),proc_dims(3)
               print *, 'info:',periods(1),periods(2),periods(3)
            endif
         CASE (2) !- planet model is a 2D decomposition
            proc_dims(1) = 1
            proc_dims(2) = INT((numprocs)**(1.d0/2.d0))
            proc_dims(3) = INT((numprocs)**(1.d0/2.d0))
            periods(1) = .false.
            periods(2) = .true.
            periods(3) = .false.
            if(myid.eq.0) then
               print *, 'info:',proc_dims(1),proc_dims(2),proc_dims(3)
               print *, 'info:',periods(1),periods(2),periods(3)
            endif
         CASE (3)
            proc_dims(1) = INT((numprocs)**(1.d0/2.d0))
            proc_dims(2) = INT((numprocs)**(1.d0/2.d0))
            proc_dims(3) = 1
            periods(1) = .false.
            periods(3) = .false.
            SELECT CASE(MODVER)
            CASE (0)
               periods(2) = .false.
            CASE (1:2)
               periods(2) = .true.
            CASE DEFAULT
               print *,'MODTYP not defined:mpi_alloc=',MODTYP
               stop
            END SELECT
         CASE (4)
            proc_dims(1) = INT((numprocs)**(1.d0/2.d0))
            proc_dims(2) = INT((numprocs)**(1.d0/2.d0))
            proc_dims(3) = 1
            periods(1) = .false.
            periods(2) = .true.
            periods(3) = .false.
         CASE (5)
            proc_dims(1) = INT((numprocs)**(1.d0/2.d0))
            proc_dims(2) = INT((numprocs)**(1.d0/2.d0))
            proc_dims(3) = 1
            periods(1) = .false.
            periods(2) = .true.
            periods(3) = .false.
         CASE (6) !- planet accretion model: 2D cylindrical in r-z
            SELECT CASE(MODVER)
            CASE (0)
               proc_dims(1) = INT((numprocs)**(1.d0/2.d0))
               proc_dims(2) = 1
               proc_dims(3) = INT((numprocs)**(1.d0/2.d0))
               periods(1) = .false.
               periods(2) = .false.
               periods(3) = .false.
            CASE (1:)
               print *,'MODTYP not defined:mpi_alloc=',MODTYP
            CASE DEFAULT
               print *,'MODTYP also not defined:mpi_alloc=',MODTYP
               stop
            END SELECT
            if(myid.eq.0) then
               print *, 'info:',proc_dims(1),proc_dims(2),proc_dims(3)
               print *, 'info:',periods(1),periods(2),periods(3)
            endif
         CASE (7:)
            print *,'unknown NCOSYS=',NCOSYS
            stop
         CASE DEFAULT
            print *,'also not defined NCOSYS=',NCOSYS
            stop
         END SELECT

      END SUBROUTINE allocatempivar
      END MODULE mpi_var_alloc
!---------------------------
!
!----- GLOBAL VARIABLES
      MODULE global_var_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_T
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_ER
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_RH
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_V,global_G
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: global_H
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: global_upperbnd
      END MODULE global_var_init

      MODULE global_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocateglobalvar
         USE global_var_init
         USE input_init
         IMPLICIT NONE
         ALLOCATE(global_T(global_NX,global_NY,global_NZ))
         ALLOCATE(global_ER(global_NX,global_NY,global_NZ))
         ALLOCATE(global_RH(global_NX,global_NY,global_NZ))
         ALLOCATE(global_V(global_NX,global_NY,global_NZ))
         ALLOCATE(global_G(global_NX,global_NY,global_NZ))
         ALLOCATE(global_H(global_NX,global_NY,global_NZ))
         ALLOCATE(global_upperbnd(global_NY,global_NZ))
      END SUBROUTINE allocateglobalvar
      END MODULE global_var_alloc
!---------------------
!
!---- FLUID VARIABLES
      MODULE fluid_var_init
      IMPLICIT NONE
      SAVE
      INTEGER :: locNX,locNY,locNZ
      INTEGER :: NTOP,num_iter,III
      INTEGER :: DEBUGFLAG
      integer :: num_print
      INTEGER,ALLOCATABLE,DIMENSION(:,:)  :: upperbnd
      DOUBLE PRECISION :: ZEIT,DELT,hydrotime,startingtime
      DOUBLE PRECISION :: TinitBottom,Tinittop,rhinittop,rhinitbottom
      DOUBLE PRECISION :: FBottom,xkapRinitBottom,T_nebula,RH_nebula
      DOUBLE PRECISION :: OMROT,rot_start,rot_full
      DOUBLE PRECISION :: k_poly,n_poly
      double precision :: ecc_anomaly,semi,radius_orb,phi0_ecc
      logical :: first_Tincident
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:),TARGET  :: T,RH
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: ER_OLD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:),TARGET  :: ER
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:),TARGET  :: V,G,H
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:),TARGET  :: GD
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: PG,CS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: CV,XMUE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: RHQX,RHQY,RHQZ
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: EdepFor,EdepAdv
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: EdepRad
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: INTERNAL1
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: INTERNAL2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: TAU,STELLARINPUT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: COOLING
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: PHOTOSPHERE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: OMEGAKEP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: TEST_VADV
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: CONV_FLUX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: TEST_ROCHE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: ARTVIS_T,ARTVIS_V
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: ARTVIS_G,ARTVIS_H
      DOUBLE PRECISION :: oldetot,oldpe,oldke,oldie
      DOUBLE PRECISION :: sigma0,ndisk,Tdisk0,Maccreated,Mgap,Mtransfer
      DOUBLE PRECISION :: D_RH,D_T
      END MODULE fluid_var_init

      MODULE fluid_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatefluidvar
          USE fluid_var_init
          USE input_init
          USE mpi_var_init
          IMPLICIT NONE
          locNX=global_NX/proc_dims(1)
          locNY=global_NY/proc_dims(2)
          locNZ=global_NZ/proc_dims(3)
          NTOP = MAX(locNX,locNY,locNZ) + 2
          ALLOCATE(upperbnd(-1:locNY+2,-1:locNZ+2))
          ALLOCATE(T(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(ER(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(ER_old(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(RH(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(V(0:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(G(-1:locNX+2,0:locNY+2,-1:locNZ+2))
          ALLOCATE(H(-1:locNX+2,-1:locNY+2,0:locNZ+2))
          ALLOCATE(GD(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(PG(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(CS(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(CV(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(XMUE(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(RHQX(0:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(RHQY(-1:locNX+2,0:locNY+2,-1:locNZ+2))
          ALLOCATE(RHQZ(-1:locNX+2,-1:locNY+2,0:locNZ+2))
          ALLOCATE(EdepFor(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(EdepAdv(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(EdepRad(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(INTERNAL1(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(INTERNAL2(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(TAU(locNX,locNY,locNZ))
          ALLOCATE(STELLARINPUT(locNX,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(COOLING(locNX,locNY,locNZ))
          ALLOCATE(PHOTOSPHERE(locNY,locNZ))
          ALLOCATE(OMEGAKEP(-1:locNX+2))
          ALLOCATE(test_vadv(3,0:locNX+1,0:locNY+1,0:locNZ+1))
          ALLOCATE(conv_flux(3,0:locNX+1,0:locNY+1,0:locNZ+1))
          ALLOCATE(test_roche(3,-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(ARTVIS_T(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(ARTVIS_V(0:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(ARTVIS_G(-1:locNX+2,0:locNY+2,-1:locNZ+2))
          ALLOCATE(ARTVIS_H(-1:locNX+2,-1:locNY+2,0:locNZ+2))
      END SUBROUTINE allocatefluidvar
      END MODULE fluid_var_alloc
!---------------------
!
!---- SCALAR VARIABLES
      MODULE scalar_var_init
      IMPLICIT NONE
      SAVE
      INTEGER :: nscalar
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:),TARGET  :: SCALAR
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:),  TARGET  :: PASS_SC
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: RH_PRE_ADV
      END MODULE scalar_var_init

      MODULE scalar_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatescalarvar
          USE scalar_var_init
          USE fluid_var_init
          USE input_init
          IMPLICIT NONE
!--density test case
          if(iscal.eq.1) nscalar=1
!--diffusion test case
          if(iscal.eq.2) nscalar=1
!--diffusion test case
          if(iscal.eq.3) nscalar=1
          ALLOCATE(SCALAR(nscalar,-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(PASS_SC(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(RH_PRE_ADV(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
      END SUBROUTINE allocatescalarvar
      END MODULE scalar_var_alloc
!---------------------------
!
!----- ARTIFICAL VISCOSITY VARIABLES
      MODULE artvis_var_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION :: dt_artvis
      DOUBLE PRECISION, PARAMETER :: artvisc_coeff = 3.0d0
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: QviscX,QviscY
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: QviscZ
      END MODULE artvis_var_init

      MODULE artvis_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocateartvisvar
         USE artvis_var_init
         USE fluid_var_init
         IMPLICIT NONE
         ALLOCATE(QviscX(0:locNX,locNY,locNZ))
         ALLOCATE(QviscY(locNX,0:locNY,locNZ))
         ALLOCATE(QviscZ(locNX,locNY,0:locNZ))
      END SUBROUTINE allocateartvisvar
      END MODULE artvis_var_alloc
!---------------------
!
!---- FORCE VARIABLES
      MODULE force_var_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: centcorx,centcory
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: centcorz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: gravx,gravy,gravz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: divtx,divty,divtz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: DISFN,eta
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: visc_flux
      DOUBLE PRECISION :: VISC_VAL,MSTAR1_RAMP
      END MODULE force_var_init

      MODULE force_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocateforcevar
          USE force_var_init
          USE fluid_var_init
          USE input_init
          USE mpi_var_init
          IMPLICIT NONE
          ALLOCATE(centcorx(locNX+1,locNY+1,locNZ+1))
          ALLOCATE(centcory(locNX+1,locNY+1,locNZ+1))
          ALLOCATE(centcorz(locNX+1,locNY+1,locNZ+1))
          ALLOCATE(gravx(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(gravy(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(gravz(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(divtx(locNX+1,0:locNY+1,0:locNZ+1))
          ALLOCATE(divty(0:locNX+1,locNY+1,0:locNZ+1))
          ALLOCATE(divtz(0:locNX+1,0:locNY+1,locNZ+1))
          ALLOCATE(DISFN(0:locNX+1,0:locNY+1,0:locNZ+1))
          ALLOCATE(eta(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(visc_flux(3,locNX,locNY,locNZ))
      END SUBROUTINE allocateforcevar
      END MODULE force_var_alloc
!---------------------
!
!---- RADIATION VARIABLES
      MODULE rad_var_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: B,RHS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: AX,CX
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: AY,CY
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: AZ,CZ
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: FLIMX,FLIMY,FLIMZ
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: difrx,difry,difrz
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: flux_x,flux_y
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: flux_z
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: flux_r_up
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: flux_r_down
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: flux_r_up_grey
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: flux_r_down_grey
      double precision,ALLOCATABLE,DIMENSION(:,:,:) :: tau_grey
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: xkapR,xkapP,xkapA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:) :: xkapW,xkapLW
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: gam_kap
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: dxkapP_dT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: dxkapA_dT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: sig
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: C1,C2,C3
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: D1,D2,D3
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: E1,E2,E3
      DOUBLE PRECISION :: global_SORPARAM,global_RHSNORM,global_RNORM
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:),TARGET  :: delta_ER
      double precision,ALLOCATABLE,DIMENSION(:,:,:) :: STELLARTAU
      double precision,ALLOCATABLE,DIMENSION(:,:,:) :: DTAUS
      double precision,ALLOCATABLE,DIMENSION(:,:,:,:) :: STELLARTAU_WAVE
      double precision,ALLOCATABLE,DIMENSION(:,:,:,:) :: DTAUS_WAVE
      double precision,ALLOCATABLE,DIMENSION(:,:,:,:) :: TAU_LW
      double precision,ALLOCATABLE,DIMENSION(:,:,:,:) :: DTAU_LW
      double precision,ALLOCATABLE,DIMENSION(:) :: BSTAR_BIN 
      double precision,ALLOCATABLE,DIMENSION(:,:,:,:) :: local_BB
      double precision,ALLOCATABLE,DIMENSION(:,:) :: mu_star
      double precision,ALLOCATABLE,DIMENSION(:,:,:) :: saveRRX
      double precision :: conv_parm
      INTEGER :: energy_iter
      INTEGER :: nwave_bins = 30 !-change this in opac_module too
      END MODULE rad_var_init

      MODULE rad_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocateradvar
          USE rad_var_init
          USE fluid_var_init
          IMPLICIT NONE
          ALLOCATE(B(locNX,locNY,locNZ))
          ALLOCATE(RHS(locNX,locNY,locNZ))
          ALLOCATE(AX(locNX,locNY,locNZ))
          ALLOCATE(CX(locNX,locNY,locNZ))
          ALLOCATE(AY(locNX,locNY,locNZ))
          ALLOCATE(CY(locNX,locNY,locNZ))
          ALLOCATE(AZ(locNX,locNY,locNZ))
          ALLOCATE(CZ(locNX,locNY,locNZ))
          ALLOCATE(FLIMX(locNX+1,locNY,locNZ))
          ALLOCATE(FLIMY(locNX,locNY+1,locNZ))
          ALLOCATE(FLIMZ(locNX,locNY,locNZ+1))
          ALLOCATE(difrx(locNX+1,locNY,locNZ))
          ALLOCATE(difry(locNX,locNY+1,locNZ))
          ALLOCATE(difrz(locNX,locNY,locNZ+1))
          ALLOCATE(flux_X(locNX+1,locNY,locNZ))
          ALLOCATE(flux_Y(locNX,locNY+1,locNZ))
          ALLOCATE(flux_Z(locNX,locNY,locNZ+1))
          ALLOCATE(flux_r_up(locNX+1,locNY,locNZ,nwave_bins))
          ALLOCATE(flux_r_down(locNX+1,locNY,locNZ,nwave_bins))
          ALLOCATE(flux_r_up_grey(locNX+1,locNY,locNZ))
          ALLOCATE(flux_r_down_grey(locNX+1,locNY,locNZ))
          ALLOCATE(tau_grey(locNX+1,locNY,locNZ))
          ALLOCATE(xkapR(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(xkapP(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(xkapA(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(xkapW(-1:locNX+2,-1:locNY+2,-1:locNZ+2,nwave_bins))
          ALLOCATE(xkapLW(-1:locNX+2,-1:locNY+2,-1:locNZ+2,nwave_bins))
          ALLOCATE(gam_kap(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(dxkapP_dT(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(dxkapA_dT(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(sig(0:locNX+1,0:locNY+1,0:locNZ+1))
          ALLOCATE(C1(locNX,locNY,locNZ))
          ALLOCATE(C2(locNX,locNY,locNZ))
          ALLOCATE(C3(locNX,locNY,locNZ))
          ALLOCATE(D1(locNX,locNY,locNZ))
          ALLOCATE(D2(locNX,locNY,locNZ))
          ALLOCATE(D3(locNX,locNY,locNZ))
          ALLOCATE(E1(locNX,locNY,locNZ))
          ALLOCATE(E2(locNX,locNY,locNZ))
          ALLOCATE(E3(locNX,locNY,locNZ))
          ALLOCATE(delta_ER(-1:locNX+2,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(STELLARTAU(locNX+1,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(DTAUS(locNX+1,-1:locNY+2,-1:locNZ+2))
          ALLOCATE(STELLARTAU_WAVE(locNX+1,-1:locNY+2,-1:locNZ+2,
     %         nwave_bins))
          ALLOCATE(TAU_LW(0:locNX+1,1:locNY,1:locNZ,nwave_bins))
          ALLOCATE(DTAU_LW(0:locNX+1,1:locNY,1:locNZ,nwave_bins))
          ALLOCATE(DTAUS_WAVE(locNX+1,-1:locNY+2,-1:locNZ+2,nwave_bins))
          ALLOCATE(BSTAR_BIN(nwave_bins))
          ALLOCATE(local_BB(-1:locNX+2,-1:locNY+2,-1:locNZ+2,
     %         nwave_bins))
          ALLOCATE(mu_star(-1:locNY+2,-1:locNZ+2))
          ALLOCATE(saveRRX(locNX,locNY,locNZ))
      END SUBROUTINE allocateradvar
      END MODULE rad_var_alloc
!---------------------------
!
!----- NEWTONIAN HEATING VARIABLES
      MODULE newton_heat_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: Teq_newton
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: kT_newton
      END MODULE newton_heat_init

      MODULE newton_heat_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatenewtonheatvar
         USE newton_heat_init
         USE fluid_var_init
         IMPLICIT NONE
         ALLOCATE(Teq_newton(locNX,locNY,locNZ))
         ALLOCATE(kT_newton(locNX,locNY,locNZ))
      END SUBROUTINE allocatenewtonheatvar
      END MODULE newton_heat_alloc
!--------------------
!
!---- GRID VARIABLES
      MODULE grid_var_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XXA,XXB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DXA,DXB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XYA,XYB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DYA,DYB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XZA,XZB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DZA,DZB
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DXA2,DXB2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DYA2,DYB2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DZA2,DZB2
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DDX0,DDX1
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DDY0,DDY1
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DDZ0,DDZ1
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DCY0,DCY1
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: DDX0kep,DDX1kep

      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: surxa,surxb,volxa
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: volxb,ccxa,ccxb
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: surya,suryb,volya
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: volyb,ccya,ccyb
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: surza,surzb,volza
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: volzb,ccza,cczb
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: geoxga,geoxha
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: geoxg,geoxh
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: geozg,geozga
      END MODULE grid_var_init

      MODULE grid_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocategridvar
          USE grid_var_init
          USE fluid_var_init
          IMPLICIT NONE
! grids
          ALLOCATE(XXA(0:locNX+2),XXB(-1:locNX+2))
          ALLOCATE(DXA(locNX+2),DXB(0:locNX+2))
          ALLOCATE(DXA2(locNX+2),DXB2(0:locNX+2))
          ALLOCATE(XYA(0:locNY+2),XYB(-1:locNY+2))
          ALLOCATE(DYA(locNY+2),DYB(0:locNY+2))
          ALLOCATE(DYA2(locNY+2),DYB2(0:locNY+2))
          ALLOCATE(XZA(0:locNZ+2),XZB(-1:locNZ+2))
          ALLOCATE(DZA(locNZ+2),DZB(0:locNZ+2))
          ALLOCATE(DZA2(locNZ+2),DZB2(0:locNZ+2))
          ALLOCATE(DDX0(0:locNX+2),DDX1(0:locNX+2))
          ALLOCATE(DDY0(0:locNY+2),DDY1(0:locNY+2))
          ALLOCATE(DDZ0(0:locNZ+2),DDZ1(0:locNZ+2))
          ALLOCATE(DCY0(1:locNY+2),DCY1(1:locNY+2))
          ALLOCATE(DDX0kep(locNX+2))
          ALLOCATE(DDX1kep(locNX+2))

! volumes
          ALLOCATE(surxa(0:locNX+2),surxb(0:locNX+2),volxa(0:locNX+2))
          ALLOCATE(volxb(0:locNX+2),ccxa(0:locNX+2),ccxb(0:locNX+2))
          ALLOCATE(surya(0:locNY+2),suryb(0:locNY+2),volya(0:locNY+2))
          ALLOCATE(volyb(0:locNY+2),ccya(0:locNY+2),ccyb(0:locNY+2))
          ALLOCATE(surza(0:locNZ+2),surzb(0:locNZ+2),volza(0:locNZ+2))
          ALLOCATE(volzb(0:locNZ+2),ccza(0:locNZ+2),cczb(0:locNZ+2))
! angular velocities correction 
          ALLOCATE(geoxga(0:locNX+2),geoxha(0:locNX+2))
          ALLOCATE(geoxg(-1:locNX+2),geoxh(-1:locNX+2))
          ALLOCATE(geozg(-1:locNZ+2),geozga(0:locNZ+2))
      END SUBROUTINE allocategridvar
      END MODULE grid_var_alloc
!---------------------
!
!---- OPACITY VARIABLES
      MODULE opc_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION :: ZKAP(14,60)
      DOUBLE PRECISION :: TVAL(60),RVAL(14)
!-- new files from Richard (3/19/08)
      double precision :: FreedmanXKAP_plnk(42,18)
      double precision :: FreedmanXKAP_ross(42,18)
      double precision :: FreedmanXKAP_abs(42,18)
      double precision :: FreedmanTVAL(42),FreedmanPVAL(18)
      double precision :: FreedmanCVAL(8)
      double precision :: ColumnXKAP_abs(42,18,8)
!-- Burrows values (5/8/08)
      double precision :: BurrowsT(50),BurrowsRH(50)
      double precision :: BurrowsXKAP_plnk(50,50)
      double precision :: BurrowsXKAP_ross(50,50)
      double precision :: BurrowsXKAP_abs(50,50)
      double precision :: BurrowsXKAP_waveopac(50,50,30)
      double precision :: BurrowsXKAP_localwave(50,50,30)
      double precision :: localBB_wavelgth(50,50,30)
!-- Eliza Kemptons values (4/2011)
      double precision :: ElizaT(39),ElizaPG(57)
      double precision :: ElizaCH4(15)
      double precision :: ElizaXKAP_plnk_woC(39,57)
      double precision :: ElizaXKAP_plnk_CO(39,57)
      double precision :: ElizaXKAP_plnk_CH4(39,57)
      double precision :: ElizaXKAP_abs_woC(39,57)
      double precision :: ElizaXKAP_abs_CO(39,57)
      double precision :: ElizaXKAP_abs_CH4(39,57)
      double precision :: ElizaXKAP_ross(39,57,15)
!-- values for Alexander and Polack opacities
      DATA TVAL/ 
     +     2.000,2.097,2.176,2.243,2.301,2.352,2.398,2.439,
     +     2.477,2.512,2.544,2.574,2.602,2.628,2.653,2.677,
     +     2.699,2.740,2.778,2.813,2.845,2.875,2.903,2.929,
     +     2.954,2.978,3.000,3.021,3.041,3.061,3.079,3.097,
     +     3.114,3.130,3.146,3.161,3.176,3.204,3.230,3.255,
     +     3.279,3.301,3.350,3.400,3.450,3.500,3.550,3.600,
     +     3.650,3.700,3.800,3.900,4.000,4.079,4.176,4.301,
     +     4.477,4.699,4.845,5.000/
      DATA RVAL/-12.d0,-11.d0,-10.d0,-9.d0,-8.d0,-7.d0,-6.d0,-5.d0,
     +     -4.d0,-3.d0,-2.d0,-1.d0,0.d0,1.d0/
!-- values for Freedman et all opacities
      DATA FreedmanTVAL/
     %     75.d0,100.d0,110.d0,120.d0,130.d0,140.d0,150.d0,160.d0,
     %     170.d0,180.d0,190.d0,200.d0,210.d0,220.d0,230.d0,240.d0,
     %     250.d0,275.d0,300.d0,400.d0,500.d0,575.d0,650.d0,725.d0,
     %     800.d0,900.d0,1000.d0,1100.d0,1200.d0,1300.d0,1400.d0,
     %     1500.d0,1600.d0,1700.d0,1800.d0,1900.d0,2000.d0,2300.d0,
     %     2600.d0,3000.d0,3500.d0,4000.d0/
      DATA FreedmanPVAL/        !--   values in mbars  
     %     0.001d0,0.003d0,0.01d0,0.03d0,0.1d0,0.3d0,1.d0,3.d0,10.d0,
     %     30.d0,100.d0,300.d0,1000.d0,3000.d0,10000.d0,30000.d0,
     %     100000.d0,300000.d0/
      DATA FreedmanCVAL/        !--   values in cm^2/mol
     %     -26.d0,-25.d0,-24.d0,-23.d0,-22.d0,-21.d0,-20.d0,-19.d0/
!-- T/RH values for Burrows opacities
c      DATA BurrowsRH/(non-log table)
c     $     1.00102e-12,1.60163e-12,2.56260e-12,4.10014e-12,6.56020e-12,
c     $     1.04963e-11,1.67940e-11,2.68703e-11,4.29923e-11,6.87875e-11,
c     $     1.10059e-10,1.76095e-10,2.81750e-10,4.50799e-10,7.21276e-10,
c     $     1.15404e-09,1.84645e-09,2.95431e-09,4.72688e-09,7.56298e-09,
c     $     1.21007e-08,1.93611e-08,3.09776e-08,4.95641e-08,7.93022e-08,
c     $     1.26883e-07,2.03012e-07,3.24818e-07,5.19707e-07,8.31529e-07,
c     $     1.33044e-06,2.12870e-06,3.40590e-06,5.44943e-06,8.71906e-06,
c     $     1.39504e-05,2.23206e-05,3.57129e-05,5.69122e-05,9.10593e-05,
c     $     0.000145694,0.000233110,0.000372602,0.000596161,0.000953854,
c     $     0.00152616,0.00244185,0.00390694,0.00625109,0.0100017/
c      DATA BurrowsT/ (non-log table)
c     $     49.9988,54.9267,60.3403,66.2874,72.8207,79.9979,87.8824,
c     $     96.5441,106.059,116.513,127.996,140.611,154.470,169.695,
c     $     186.420,204.793,224.977,247.151,271.510,298.270,327.668,
c     $     359.963,395.440,434.415,477.231,524.266,575.938,632.702,
c     $     695.061,763.566,838.823,920.576,1011.31,1110.98,1220.48,
c     $     1340.77,1472.92,1618.09,1777.57,1952.76,2145.23,2356.66,
c     $     2588.93,2844.09,3124.41,3432.35,3770.64,4142.28,4550.54,
c     $     4999.04/
      DATA BurrowsRH/
     $     -11.9996, -11.7954, -11.5913, -11.3872, -11.1831,
     $     -10.9790, -10.7748, -10.5707, -10.3666, -10.1625,
     $     -9.95837, -9.75425, -9.55014, -9.34602, -9.14190,
     $     -8.93778, -8.73366, -8.52954, -8.32543, -8.12131,
     $     -7.91719, -7.71307, -7.50895, -7.30483, -7.10071,
     $     -6.89660, -6.69248, -6.48836, -6.28424, -6.08012,
     $     -5.87600, -5.67189, -5.46777, -5.26365, -5.05953,
     $     -4.85541, -4.65129, -4.44718, -4.24479, -4.04068,
     $     -3.83656, -3.63244, -3.42875, -3.22464, -3.02052,
     $     -2.81640, -2.61228, -2.40816, -2.20404, -1.99993/
      DATA BurrowsT/
     $     1.69896,1.73978,1.78061,1.82143,1.86225,
     $     1.90308,1.94390,1.98473,2.02555,2.06637,
     $     2.10720,2.14802,2.18884,2.22967,2.27049,
     $     2.31132,2.35214,2.39296,2.43379,2.47461,
     $     2.51543,2.55626,2.59708,2.63790,2.67873,
     $     2.71955,2.76038,2.80120,2.84202,2.88285,
     $     2.92367,2.96406,3.00488,3.04571,3.08653,
     $     3.12735,3.16818,3.20900,3.24983,3.29065,
     $     3.33147,3.37230,3.41312,3.45394,3.49477,
     $     3.53559,3.57642,3.61724,3.65806,3.69889/

      DATA ElizaT/
     $     100.d0,150.d0,200.d0,250.d0,300.d0,350.d0,400.d0,450.d0,
     $     500.d0,550.d0,600.d0,650.d0,700.d0,750.d0,800.d0,850.d0,
     $     900.d0,950.d0,1000.d0,1100.d0,1200.d0,1300.d0,1400.d0,
     $     1500.d0,1600.d0,1700.d0,1800.d0,1900.d0,2000.d0,2100.0,
     $     2200.0,2300.0,2400.0,2500.0,2600.0,2700.0,2800.0,2900.0,
     $     3000.0/
      DATA ElizaPG/
     $     1.0000000E-04, 1.7782790E-04, 3.1622782E-04, 5.6234130E-04,
     $     1.0000000E-03, 1.7782790E-03, 3.1622779E-03, 5.6234132E-03,
     $     1.0000000E-02, 1.7782789E-02, 3.1622782E-02, 5.6234129E-02,
     $     1.0000000E-02, 1.7782791E-01, 3.1622779E-01, 5.6234127E-01,
     $     1.0000000E+00, 1.7782789E+00, 3.1622779E+00, 5.6234126E+00,
     $     1.0000000E+01, 1.7782789E+01, 3.1622780E+01, 5.6234131E+01,
     $     1.0000000E+02, 1.7782790E+02, 3.1622778E+02, 5.6234131E+02,
     $     1.0000000E+03, 1.7782789E+03, 3.1622781E+03, 5.6234131E+03,
     $     1.0000000E+04, 1.7782791E+04, 3.1622781E+04, 5.6234133E+04,
     $     1.0000000E+05, 1.7782789E+05, 3.1622778E+05, 5.6234131E+05,
     $     1.0000000E+06, 1.7782790E+06, 3.1622780E+06, 5.6234130E+06,
     $     1.0000000E+07, 1.7782790E+07, 3.1622780E+07, 5.6234128E+07,
     $     1.0000000E+08, 1.7782790E+08, 3.1622781E+08, 5.6234125E+08,
     $     1.0000000E+09, 1.7782790E+09, 3.1622781E+09, 5.6234132E+09,
     $     1.0000000E+10/
      DATA ElizaCH4/
     $     0.d0,0.01d0,0.05d0,0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,
     $     0.7d0,0.8d0,0.9d0,0.95d0,0.99d0,1.d0/
      END MODULE opc_init
!---------------------
!
!---- OPLIN PARAMETERS
      MODULE oplin_constants
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION, PARAMETER :: power1=4.44444444 E - 2
      DOUBLE PRECISION, PARAMETER :: power2=2.381 E - 2
      DOUBLE PRECISION, PARAMETER :: power3=2.267 E - 1
      DOUBLE PRECISION, PARAMETER :: t234=1.6 E 3
      DOUBLE PRECISION, PARAMETER :: t456=5.7 E 3
      DOUBLE PRECISION, PARAMETER :: t678=2.28 E 6
      DOUBLE PRECISION, PARAMETER :: ak1=2.0 E - 4
      DOUBLE PRECISION, PARAMETER :: ak2=2.0 E 16
      DOUBLE PRECISION, PARAMETER :: ak3=5.0 E - 3
      DOUBLE PRECISION, PARAMETER :: bk3=50.D0
      DOUBLE PRECISION, PARAMETER :: bk4=2.0 E - 2
      DOUBLE PRECISION, PARAMETER :: bk5=2.0 E 4
      DOUBLE PRECISION, PARAMETER :: bk6=1.D4
      DOUBLE PRECISION, PARAMETER :: bk7=1.5D10
      DOUBLE PRECISION, PARAMETER :: bk8=0.348
      END MODULE oplin_constants
!---------------------
!
!---- TRACER PARTICLE VARIABLES
      MODULE particle_var_init
      IMPLICIT NONE
      SAVE
      INTEGER,ALLOCATABLE,DIMENSION(:) :: part_proc_id
      INTEGER,ALLOCATABLE,DIMENSION(:) :: send_to_upper,send_to_lower
      INTEGER,ALLOCATABLE,DIMENSION(:) :: send_to_above,send_to_below
      INTEGER,ALLOCATABLE,DIMENSION(:) :: recieve_from_upper
      INTEGER,ALLOCATABLE,DIMENSION(:) :: recieve_from_lower
      INTEGER,ALLOCATABLE,DIMENSION(:) :: recieve_from_above
      INTEGER,ALLOCATABLE,DIMENSION(:) :: recieve_from_below
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: part_posx,part_posy
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: part_posz
c      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: part_vx
c      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: part_vy
c      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: part_vz
      END MODULE particle_var_init

      MODULE particle_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocateparticlevar
          USE particle_var_init
          USE input_init
          IMPLICIT NONE
          ALLOCATE(part_proc_id(nparticles))
          ALLOCATE(send_to_upper(nparticles))
          ALLOCATE(send_to_lower(nparticles))
          ALLOCATE(send_to_above(nparticles))
          ALLOCATE(send_to_below(nparticles))
          ALLOCATE(recieve_from_upper(nparticles))
          ALLOCATE(recieve_from_lower(nparticles))
          ALLOCATE(recieve_from_above(nparticles))
          ALLOCATE(recieve_from_below(nparticles))
          ALLOCATE(part_posx(nparticles))
          ALLOCATE(part_posy(nparticles))
          ALLOCATE(part_posz(nparticles))
c          ALLOCATE(part_vx(nparticles))
c          ALLOCATE(part_vy(nparticles))
c          ALLOCATE(part_vz(nparticles))
      END SUBROUTINE allocateparticlevar
      END MODULE particle_var_alloc
!---------------------
!
!---- ENERGY VARIABLES
      MODULE deltaER_var_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: RHS_ER
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: EVAL1,EVAL2,EVAL3
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: EVAL4,EVAL5,EVAL6
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: EVAL7,EVAL8
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: EVAL_1_6
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: Jacob_eta
!---
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: FQ
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: FR
      END MODULE deltaER_var_init

      MODULE deltaER_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatedeltaERvar
          USE deltaER_var_init
          USE fluid_var_init
          IMPLICIT NONE
          ALLOCATE(RHS_ER(locNX,locNY,locNZ))
          ALLOCATE(EVAL1(locNX,locNY,locNZ))
          ALLOCATE(EVAL2(locNX,locNY,locNZ))
          ALLOCATE(EVAL3(locNX,locNY,locNZ))
          ALLOCATE(EVAL4(locNX,locNY,locNZ))
          ALLOCATE(EVAL5(locNX,locNY,locNZ))
          ALLOCATE(EVAL6(locNX,locNY,locNZ))
          ALLOCATE(EVAL7(locNX,locNY,locNZ))
          ALLOCATE(EVAL8(locNX,locNY,locNZ))
          ALLOCATE(EVAL_1_6(locNX,locNY,locNZ))
          ALLOCATE(Jacob_eta(locNX,locNY,locNZ))
!-----
          ALLOCATE(FQ(locNX,locNY,locNZ))
          ALLOCATE(FR(locNX,locNY,locNZ))
!---
      END SUBROUTINE allocatedeltaERvar
      END MODULE deltaER_var_alloc
!---------------------
!
!---- SOR VARIABLES
      MODULE sor_var_init
      IMPLICIT NONE
      SAVE
      INTEGER :: SOR_ITER
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: SOR_A
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: SOR_B
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: SOR_C
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: SOR_D
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: SOR_E
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: SOR_F
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: SOR_G
      END MODULE sor_var_init

      MODULE sor_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatesorvar
          USE sor_var_init
          USE fluid_var_init
          IMPLICIT NONE
          ALLOCATE(SOR_A(locNX,locNY,locNZ))
          ALLOCATE(SOR_B(locNX,locNY,locNZ))
          ALLOCATE(SOR_C(locNX,locNY,locNZ))
          ALLOCATE(SOR_D(locNX,locNY,locNZ))
          ALLOCATE(SOR_E(locNX,locNY,locNZ))
          ALLOCATE(SOR_F(locNX,locNY,locNZ))
          ALLOCATE(SOR_G(locNX,locNY,locNZ))
      END SUBROUTINE allocatesorvar
      END MODULE sor_var_alloc
!---------------------
!
!---- DIFFUSION VARIABLES
      MODULE dif_var_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: RHS_DIF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: DVAL1,DVAL2,DVAL3
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: DVAL4,DVAL5,DVAL6
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: DVAL_1_6
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: difrx_DIF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: difry_DIF
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: difrz_DIF
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: DIFF_PTR
      END MODULE dif_var_init

      MODULE dif_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatedifvar
          USE dif_var_init
          USE fluid_var_init
          IMPLICIT NONE
          ALLOCATE(RHS_DIF(locNX,locNY,locNZ))
          ALLOCATE(difrx_DIF(locNX+1,locNY,locNZ))
          ALLOCATE(difry_DIF(locNX,locNY+1,locNZ))
          ALLOCATE(difrz_DIF(locNX,locNY,locNZ+1))
          ALLOCATE(DVAL1(locNX,locNY,locNZ))
          ALLOCATE(DVAL2(locNX,locNY,locNZ))
          ALLOCATE(DVAL3(locNX,locNY,locNZ))
          ALLOCATE(DVAL4(locNX,locNY,locNZ))
          ALLOCATE(DVAL5(locNX,locNY,locNZ))
          ALLOCATE(DVAL6(locNX,locNY,locNZ))
          ALLOCATE(DVAL_1_6(locNX,locNY,locNZ))
      END SUBROUTINE allocatedifvar
      END MODULE dif_var_alloc
!---------------------
!
!---- HYPRE VARIABLES
      MODULE hypre_var_init
      IMPLICIT NONE
      SAVE
!-- Main Hypre pointers
      integer*8 bHYPRE_mpicomm
      integer*8 mpi_comm,grid,stencil,graph
      integer*8 H_matrix,x_vec,b_vec
      integer hypre_ndim,hypre_nparts,hypre_part,hypre_var
!-- Cast Pointers
      integer*8 dummy,sA,vb,vx,opA
!-- Solver and preconditioner Pointers
      integer *8 PCGsolver,SMGprecond,precond
      integer*8 except
      integer ilower(3),iupper(3)
      integer offsets(3,7),sentry
      integer HYPRE_STRUCT,object_type
      parameter( HYPRE_STRUCT =  1111 )
      integer stencil_indices(7),nentries, nvalues
      integer neighbor_stencil_indices(1)
      double precision tol
      integer vartypes(1)
      integer MAX_BOUNDARY_VALUES,MAX_BOUNDARY_STENCILS
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: H_values
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: b_values
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: x_values
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: boundary_vals
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: boundary_matrix_vals
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: neighbor_matrix_vals
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: neighbor_b_vals
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: answer
      END MODULE hypre_var_init

      MODULE hypre_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatehyprevar
         USE hypre_var_init
         USE fluid_var_init
         IMPLICIT NONE
         allocate(H_values(locNX*locNY*locNZ*7))
         ALLOCATE(x_values(locNX*locNY*locNZ))
         ALLOCATE(b_values(locNX*locNY*locNZ))
         MAX_BOUNDARY_VALUES=MAX(locNX*locNY,locNX*locNZ,locNY*locNZ)
         MAX_BOUNDARY_STENCILS = MAX_BOUNDARY_VALUES*7
         ALLOCATE(boundary_matrix_vals(MAX_BOUNDARY_STENCILS))
         ALLOCATE(boundary_vals(MAX_BOUNDARY_VALUES))
         ALLOCATE(neighbor_matrix_vals(MAX_BOUNDARY_VALUES))
         ALLOCATE(neighbor_b_vals(MAX_BOUNDARY_VALUES))
         ALLOCATE(answer(locNX*locNY*locNZ))
      END SUBROUTINE allocatehyprevar
      END MODULE hypre_var_alloc
!---------------------
!
!---- ARTIFICAL COOLING VARIABLES
      MODULE relax_var_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: Trelax
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: ERrelax
      END MODULE relax_var_init

      MODULE relax_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocaterelaxvar
          USE relax_var_init
          USE input_init
          USE fluid_var_init
          IMPLICIT NONE
          ALLOCATE(Trelax(locNX,locNY,locNZ))
          ALLOCATE(ERrelax(locNX,locNY,locNZ))
      END SUBROUTINE allocaterelaxvar
      END MODULE relax_var_alloc
!---------------------
!
!---- SAUMON EOS VARIABLES
      MODULE saumon_var_init
      IMPLICIT NONE
      SAVE
      integer :: NTMP,NPRS
      integer, parameter :: ngmax=10000
      DOUBLE PRECISION :: XX,YY
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: HTABLE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: HETABLE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: HTLOG
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: HETLOG
      END MODULE saumon_var_init

      MODULE saumon_var_alloc
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocate_saumon_eos
         USE saumon_var_init
         XX = 0.7
         YY = 0.3
         NTMP=63
         NPRS=76
         ALLOCATE(HTABLE(NTMP,NPRS,11))
         ALLOCATE(HETABLE(NTMP,NPRS,11))
         ALLOCATE(HTLOG(NTMP))
         ALLOCATE(HETLOG(NTMP))
      END SUBROUTINE allocate_saumon_eos
      END MODULE saumon_var_alloc
!---------------------
!
!---- PLANETARY VARIABLES
!      MODULE planet_var
!      IMPLICIT NONE
!      SAVE
!      DOUBLE PRECISION :: rplanet,phiplanet
!      END MODULE planet_var
!----
!      MODULE rkdumb_path
!      USE nrtype REAL(SP), DIMENSION(:), ALLOCATABLE:: xx
!      REAL(SP), DIMENSION(:,:), ALLOCATABLE :: y
!      END MODULE rkdumb_path
!---------------------
!
!---- GLOBAL CONSTANTS
      MODULE global_constants
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793D0
      DOUBLE PRECISION, PARAMETER :: MSOL   = 1.989  E  33
      DOUBLE PRECISION, PARAMETER :: RSOL   = 6.9598  E  10
      DOUBLE PRECISION, PARAMETER :: MJUP   = 1.8987  E  30
      DOUBLE PRECISION, PARAMETER :: RJUP   = 7.1492 E 9
      DOUBLE PRECISION, PARAMETER :: GRAV   = 6.672  E - 8
      DOUBLE PRECISION, PARAMETER :: FD     = 2.99792458 E 10
      DOUBLE PRECISION, PARAMETER :: ARAD   = 7.56596 E -15
      DOUBLE PRECISION, PARAMETER :: RGAS   = 8.315  E 7
      DOUBLE PRECISION, PARAMETER :: DAY    = 8.64  E 4
      DOUBLE PRECISION, PARAMETER :: AU     = 1.495979 E 13
      DOUBLE PRECISION, PARAMETER :: MEARTH = 5.976  E  27
      DOUBLE PRECISION, PARAMETER :: SBCONST = 5.6705 E -5
      DOUBLE PRECISION, PARAMETER :: MPROTON = 1.672621 E -24
      DOUBLE PRECISION, PARAMETER :: MHYRODGEN = 1.67357 E -24
      DOUBLE PRECISION, PARAMETER :: MHELIUM   = 6.646442 E -24
      DOUBLE PRECISION, PARAMETER :: KBOLTZ = 1.3806504 E -16
      DOUBLE PRECISION, PARAMETER :: PLANCK = 6.62606885 E -27
      DOUBLE PRECISION, PARAMETER :: TINY = 1.0 E -30
      END MODULE global_constants

