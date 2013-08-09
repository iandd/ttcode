      SUBROUTINE SCHRITT
**  Calculates the Magnitude of the Time Step
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE force_var_init
      USE grid_var_init
      USE artvis_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: V1,V2,V3,ZC1,ZC2,ZC3
      double precision :: localDT
      double precision :: XNU,ZCF
**  Check for negative values of rho, T, Er or Pg
      IF(IGASTYP.eq.0) then
         if(MINVAL(RH).le.0.d0.or.MINVAL(T).le.0.d0) then
            print *,'non-positive value of rh/t'
            print *,'III,num_iter',III,num_iter
            print *,'min rh/T',MINVAL(RH),MINVAL(T)
            print *,'(i,j,k) of RH_MIN',minloc(RH)
            print *,'(i,j,k) of T_MIN',minloc(T)
            stop
         endif
      else
         if(MINVAL(RH).le.0.d0.or.MINVAL(T).le.0.d0) then
            do k=1,locNZ
               do j=1,locNY
                  do i=1,locNX
                     if(RH(i,j,k).le.0.d0.or.T(i,j,k).le.0.d0.or.
     $                    PG(i,j,k).le.0.d0) then
                        print *,'non-positive value of rh/pg/t'
                        print *,'III,num_iter,myid',III,num_iter,myid
                        print *,'min rh/T/P',RH(i,j,k),T(i,j,k),
     $                       PG(i,j,k),myid
                        print *,'V,G,H,myid',V(i,j,k),G(i,j,k),H(i,j,k),
     $                       myid
                        print *,'(i,j,k,myid) of error',i,j,k,myid
                        print *,'(upperbnd,myid) of error',upperbnd(j,k)
     $                       ,myid
                        call clean_stop
                     endif
                  enddo
               enddo
            enddo
         endif
      endif      
**  if fdelt=0. >> delt=dtmax
      IF (FDELT.EQ.0.0) THEN
         DELT=DTMAX
      ELSE
         localDT=DELT*100.d0/FDELT
**  SOUND VELOCITY
         IF(ICSTYP.NE.1.AND.ICSTYP.NE.2) then
            CALL SOUND
         ENDIF
         localDT=DTMAX
!-- COURANT CONDITION IN ALL THREE DIMENTIONS
         if(locNZ.ne.1) then
            DO K=1,locNZ
               DO J=1,locNY
                  DO I=1,locNX
                     V1=ABS(V(I+1,J,K))
                     V2=ABS(G(I,J+1,K))
                     V3=ABS(H(I,J,K+1))
                     ZC1=DXA(I+1)/(V1+CS(i,j,k))
                     ZC2=DYA(J+1)/(V2+CS(i,j,k)/geoxg(i)/geozg(k))
                     ZC3=DZA(K+1)/(V3+CS(i,j,k)/geoxh(i))
                     localDT=MIN(ZC1,ZC2,ZC3,localDT)
                     IF (localDT.LE.0.d0) THEN
                        write (*,*) 'local dt negative',localDT,zc1,zc2,
     %                       zc3,V1,v2,v3,cs(i,j,k),T(i,j,k),
     %                       dxa(i+1),dya(j+1),dza(i+1),i,j,k,myid
                        STOP
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
!-- COURANT CONDITION IN TWO DIMENTIONS
         elseif(locNZ.eq.1) then
            K=1
            DO J=1,locNY
               DO I=1,locNX
                  V1=ABS(V(I+1,J,K))
                  V2=ABS(G(I,J+1,K))
                  ZC1=DXA(I+1)/(V1+CS(i,j,k))
                  ZC2=DYA(J+1)/(V2+CS(i,j,k)/geoxg(i)/geozg(k))
                  localDT=MIN(ZC1,ZC2,localDT)
                  IF (localDT.LE.0.d0) THEN
                     write (*,*) 'local dt negative',localDT,zc1,zc2,
     %                    V1,v2,cs(i,j,k),T(i,j,k),i,j,k,myid
                     STOP
                  ENDIF
               ENDDO
            ENDDO
         endif
!-- VISCOUS TIMESTEP RESTRICTION
!----- 3D VISCOUS TIMESTEP
         IF(IVIS.GT.0) THEN
            if(locNZ.ne.1) then
               DO K=1,locNZ
                  DO J=1,locNY
                     DO I=1,locNX
                        XNU= ETA(I,J,K) / RH(I,J,K)
                        ZCF= (1.d0/8.d0)/ XNU /
     &                       (DXA(I+1)**(-2.d0) +
     &                       (DYA(J+1)*geoxg(i)*geozg(k))**(-2.d0) +
     &                       (DZA(K+1)*geoxh(i))**(-2.d0))
                        localDT=MIN(ZCF,localDT)
                     ENDDO
                  ENDDO
               ENDDO
!----- 2D VISCOUS TIMESTEP
            elseif(locNZ.eq.1) then
               K=1
               print *,'check 2d visc'
               call clean_stop
            endif
         ENDIF

!----- ARTIFICIAL VISCOSITY TIMESTEP
c         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c         if(localDT.gt.dt_artvis) then
c         write(*,'(A,2(I8,1x),2(1x,e12.3))') 'dt_ARTVISC:',iii,myid,
c     %        localDT,dt_artvis
c         endif
         if(iartvis.gt.0) 
     %        localDT=MIN(dt_artvis,localDT)
         
!--- CHECK MAX TIMESTEP
         localDT=  MIN (localDT,DTMAX) * FDELT
!--- find the globally lowest timestep
         call MPI_ALLREDUCE(localDT,DELT,1,MPI_DOUBLE_PRECISION,MPI_MIN,
     %        MPI_COMM_WORLD,ierr)
      ENDIF
      RETURN
      END SUBROUTINE SCHRITT


      SUBROUTINE VARLIMIT(ilimit,feld,varmin,varmax)
**  Imposes limits on the physical variables
      USE fluid_var_init
      IMPLICIT NONE
      integer :: i,j,k,ilimit
      double precision :: feld(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: varmin,varmax
      if (locNY.eq.1) then
         print *,'not sure what to do here...:VARLIMIT'
         stop
      endif
      feld = min(feld,varmax)
      feld = max(feld,varmin)
      RETURN
      END SUBROUTINE VARLIMIT

      SUBROUTINE CALCMUSTAR
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      IMPLICIT NONE
      INTEGER :: I,J,K
      if(F_INCIDENT.eq.2.or.F_INCIDENT.eq.3.or.F_INCIDENT.eq.8) then
         do k=-1,locNZ+2
            do j=-1,locNY+2
               mu_star(j,k) = max((cos(xyb(j)))*(cos(xzb(k))),
     $              10.**(-7.))
            enddo
         enddo
      elseif(F_INCIDENT.eq.6.or.F_INCIDENT.eq.7) then
         mu_star(j,k) = max((cos(xyb(j)-phi0_ecc))*(cos(xzb(k))),
     $        10.**(-7.))
         print *,'fix/check mu_star for F_inc=',F_incident
         call clean_stop
      else  !spherically symmetric
         mu_star = 1.d0 
      endif
      RETURN
      END SUBROUTINE CALCMUSTAR

      SUBROUTINE VELOCITYDAMPING(VINT,GINT,HINT)
**  Uniformily damps velocity
      USE input_init
      USE mpi_var_init
      USE global_constants
      USE fluid_var_init
      USE grid_var_init
      IMPLICIT NONE
      integer :: i,j,k
      DOUBLE PRECISION :: VINT(locNX,locNY,locNZ)
      DOUBLE PRECISION :: GINT(locNX,locNY,locNZ)
      DOUBLE PRECISION :: HINT(locNX,locNY,locNZ)
      DOUBLE PRECISION :: velramp,veldamp_start,DAMP_TIME
      double precision :: eta_sponge_start,sponge_coeff,eta_sponge,Rw
      IF(VELDAMP.LT.1.d0) THEN
!----- ala Burrket
         VELRAMP = VELDAMP
         VINT=VELRAMP*VINT
         GINT=VELRAMP*GINT
         HINT=VELRAMP*HINT
      ELSEIF(VELDAMP.GE.1.d0.and.VELDAMP.LT.2.d0) THEN
!-----RAMP-DOWN DAMPING (From VELDAMP to 1.0)
         VELRAMP = (1.d0-(VELDAMP-1.d0))*(1.d0*III/
     $        (1.d0*DAMP_STEPS))+(VELDAMP-1.d0)
c         velramp  = velramp + 0.307
         if(velramp.ge.1.0) velramp = 1.d0
         VINT=VELRAMP*VINT
         GINT=VELRAMP*GINT
         HINT=VELRAMP*HINT
      ELSEIF(VELDAMP.GE.2.d0.and.VELDAMP.LT.3.d0) THEN
!-----RAMP-DOWN DAMPING (From VELDAMP to 1.0) USING HYDROTIME
         veldamp_start = 0.d0*DAY
         DAMP_TIME = 20.d0*DAY
         if(hydrotime.ge.veldamp_start) then
            VELRAMP = (1.d0-(VELDAMP-2.d0))*((hydrotime-veldamp_start)/
     $           (DAMP_TIME))+(VELDAMP-2.d0)
            if(velramp.ge.1.0) velramp = 1.d0
         else
            velramp=0.d0
         ENDIF
         VINT=VELRAMP*VINT
         GINT=VELRAMP*GINT
         HINT=VELRAMP*HINT
      ELSEIF(VELDAMP.EQ.4.d0) then
!-----SPONGE LAYER AT THE TOP ON RADIAL VELOCITY
         eta_sponge_start = 0.80
         sponge_coeff = 0.8
         do k=1,locNZ
            do j=1,locNY
               do i=1,upperbnd(j,k)
                  eta_sponge = (xxa(i)-xxa(1))/
     $                 (xxa(upperbnd(j,k))-xxa(1))
                  if(eta_sponge.gt.eta_sponge_start) then
                     Rw=sponge_coeff*(sin((eta_sponge-eta_sponge_start)/
     $                    (1.d0-eta_sponge_start)))**2.d0
c                     if(myid.eq.0.and.III.eq.3.and.j.eq.locNY/2.and.
c     $                    k.eq.locNZ/2) then
c                        write(*,'(A,2(1x,I8),4(1x,e13.6))') 'SPONGE',
c     $                       i,upperbnd(j,k),eta_sponge,Rw,
c     $                       VINT(i,j,k),VINT(i,j,k)/(1.d0+Rw*DELT)
c                     endif
                     VINT(i,j,k)=VINT(i,j,k)/(1.d0+Rw*DELT)
                  endif
               enddo
            enddo
         enddo

      ENDIF
      if(veldamp.ne.4.and.((III-1)/(IREA/10))*(IREA/10).eq.(III-1)) THEN
         if(myid.eq.0) print *,'velramp=',velramp
      ENDIF
      RETURN
      END SUBROUTINE VELOCITYDAMPING

      SUBROUTINE POLEVELOCITYDAMPING(III,locNX,locNY,locNZ,VINT,
     $     GINT,HINT,hydrotime)
**  limit velocity at the polar grid cells
      USE input_init
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,j,k,III,locNX,locNY,locNZ
      DOUBLE PRECISION :: VINT(0:locNX+2,-1:locNY+2,-1:locNZ+2)
      DOUBLE PRECISION :: GINT(-1:locNX+2,0:locNY+2,-1:locNZ+2)
      DOUBLE PRECISION :: HINT(-1:locNX+2,-1:locNY+2,0:locNZ+2)
      DOUBLE PRECISION :: VELRAMP,POLEDAMP,hydrotime
!----- RAMP-DOWN DAMPING (From POLEDAMP to 1.0)                                                     
      POLEDAMP = 0.01
      VELRAMP = (1.d0-(POLEDAMP))*(1.d0*III/
     $     (1.d0*DAMP_STEPS))+(POLEDAMP)
      if(VELRAMP.GE.1.0) VELRAMP = 1.d0
      IF (((III-1)/(IREA/10))*(IREA/10).eq.(III-1)) THEN
         if(myid.eq.0) print *,'polar velramp=',velramp,III,hydrotime
      ENDIF
      if(proc_coords(3).eq.proc_dims(3)-1) then !-north-pole
         do j=-1,locNY+2
            do i=-1,locNX+2
               if(i.ne.-1) then
                  VINT(i,j,locNZ)   = VELRAMP*VINT(i,j,locNZ)
                  VINT(i,j,locNZ+1) = VELRAMP*VINT(i,j,locNZ+1)
                  VINT(i,j,locNZ+2) = VELRAMP*VINT(i,j,locNZ+2)
               endif
               if(j.ne.-1) then
                  GINT(i,j,locNZ)   = VELRAMP*GINT(i,j,locNZ)
                  GINT(i,j,locNZ+1) = VELRAMP*GINT(i,j,locNZ+1)
                  GINT(i,j,locNZ+2) = VELRAMP*GINT(i,j,locNZ+2)
               endif
               HINT(i,j,locNZ)   = VELRAMP*HINT(i,j,locNZ)
               HINT(i,j,locNZ+1) = VELRAMP*HINT(i,j,locNZ+1)
               HINT(i,j,locNZ+2) = VELRAMP*HINT(i,j,locNZ+2)
            enddo
         enddo
      endif
      if(proc_coords(3).eq.0) then !-south-pole
         do j=-1,locNY+2
            do i=-1,locNX+2
               if(i.ne.-1) then
                  VINT(i,j,-1) = VELRAMP*VINT(i,j,-1)
                  VINT(i,j,0)  = VELRAMP*VINT(i,j,0)
c                  VINT(i,j,1)  = VELRAMP*VINT(i,j,1)
               endif
               if(j.ne.-1) then
                  GINT(i,j,-1) = VELRAMP*GINT(i,j,-1)
                  GINT(i,j,0)  = VELRAMP*GINT(i,j,0)
c                  GINT(i,j,1)  = VELRAMP*GINT(i,j,1)
               endif
               HINT(i,j,0) = VELRAMP*HINT(i,j,0)
               HINT(i,j,1) = VELRAMP*HINT(i,j,1)
c               HINT(i,j,2) = VELRAMP*HINT(i,j,2)
            enddo
         enddo
      endif      
      RETURN
      END SUBROUTINE POLEVELOCITYDAMPING


      SUBROUTINE SPHERICALVELLIMIT
!---------------------------------------------------
**  Imposes limits on G and H for NCOSYS=2
!
!-- ivel=1 : set the ceiling to be velmax (from ipp) for ALL velocities
!-- ivel=2 : ramp ceiling from VELDAMP_START to VELMAX over DAMP_TIME for ALL velocities
!-- ivel=3 : set the celing on V_R ONLY to be velmax (from ipp)
!-- ivel=4 : ramp celing on V_R ONLY from VELDAMP_START to VELMAX in DAMP_TIME
!
!---------------------------------------------------
      USE fluid_var_init
      USE input_init
      USE grid_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,j,k
      DOUBLE PRECISION :: VARMIN,VARMAX
      DOUBLE PRECISION :: TIME_START,DAMP_TIME,VELDAMP_START
      if(ivel.eq.1.or.ivel.eq.3) then 
         VARMIN = -VELMAX
         VARMAX = VELMAX
      elseif(ivel.eq.2.or.ivel.eq.4) then
!-----RAMP-UP CELING FOR VELOCITY FROM VELDAMP_START TO VELMAX
         TIME_START = 0.d0*DAY
c         TIME_START = 380.d0*DAY
c         DAMP_TIME = 25.d0*DAY
         DAMP_TIME = 25.d0*DAY
         VELDAMP_START = 0.d0
c         VELDAMP_START =5e4
         if(hydrotime.ge.time_start.and.hydrotime.le.
     $        (time_start+damp_time)) then
            VARMAX = (VELMAX-VELDAMP_START)*(sin((PI/2.d0)*
     $           (hydrotime-TIME_START)/DAMP_TIME)**(4.0)) +
     $           VELDAMP_START
c            VARMAX = ((hydrotime-TIME_START)*VELMAX/DAMP_TIME)
            VARMIN = -VARMAX
         elseif(hydrotime.gt.(time_start+damp_time)) then
            VARMAX=VELMAX
            VARMIN=-VELMAX
         else
            VARMAX=VELDAMP_START
            VARMIN=-VELDAMP_START
         ENDIF
         IF (((III-1)/(IREA/10))*(IREA/10).eq.(III-1)) THEN
            if(myid.eq.0) print *,'velocity ceiling=',
     $           VARMAX,hydrotime/DAY
         ENDIF
      else
         print *,'ivel=',ivel,'is not available'
      endif
      if(ivel.eq.1.or.ivel.eq.2) then
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               DO I=-1,locNX+2
                  if(j.ne.-1) then
                     G(I,J,K)=min(G(I,J,K),VARMAX/(XXB(I)*COS(XZB(K))))
                     G(I,J,K)=max(G(I,J,K),VARMIN/(XXB(I)*COS(XZB(K))))
!     -make sure the polar G=0 conditions hold
                     if(poles.and.proc_coords(3).eq.0.and.k.le.1) 
     $                    G(I,j,k)=0.d0
                     if(poles.and.proc_coords(3).eq.proc_dims(3)-1.and.
     $                    k.ge.locNZ) G(I,j,k)=0.d0
                  endif
                  if(k.ne.-1) then
                     H(I,J,K)=min(H(I,J,K),VARMAX/XXB(I))
                     H(I,J,K)=max(H(I,J,K),VARMIN/XXB(I))
                  endif
               ENDDO
            ENDDO
         ENDDO
      endif
      V = min(V,VARMAX)
      V = max(V,VARMIN)
      RETURN
      END SUBROUTINE SPHERICALVELLIMIT
      
      SUBROUTINE MAXVELOCITIES
**  Finds the maximum velocitys
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE grid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: MAXX,MAXY,MAXZ,MAXCS
      double precision :: global_MAXX,global_MAXY,global_MAXz
      double precision :: global_MAXCS,global_MINX
      integer :: maxi,maxj,maxk
      IF(ICSTYP.NE.1.AND.ICSTYP.NE.2) then
         CALL SOUND
      ENDIF
      maxX  = MAXVAL(ABS(V))          !- these ones I can do as a matrix always
      maxCS = MAXVAL(ABS(CS))
      if(ncosys.ne.0) then
         MAXY = 0.d0
         MAXZ = 0.d0
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  if(ncosys.eq.1) then
                     maxY = max(ABS(G(i,j,k)*xxb(i)),maxY)
                  elseif(ncosys.eq.2) then
                     if(ABS(H(i,j,k)*xxb(i)).gt.maxZ) maxk = k
                     maxY = max(ABS(G(i,j,k)*xxb(i)*cos(XZB(k))),maxY)
                     maxZ = max(ABS(H(i,j,k)*xxb(i)),maxZ)
                  else
                     print *,'Unknown coordinate sys:MAXVELOCITIES',
     $                    ncosys
                     stop
                  endif
               ENDDO
            ENDDO
         ENDDO
      else
         maxY = MAXVAL(ABS(G))
         maxZ = MAXVAL(ABS(H))
      endif
      call MPI_ALLREDUCE(maxX,global_maxX,1,
     %     MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(maxY,global_maxY,1,
     %     MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(maxZ,global_maxZ,1,
     %     MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)         
      call MPI_ALLREDUCE(maxCS,global_maxCS,1,
     %     MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(maxX,global_MINX,1,
     %     MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
c         write(*,'(A,4(1x,e12.3))')'MAX V,G,H,CS',global_MAXX,
c     %        global_MAXY,global_MAXZ,global_MAXCS
         write(*,'(A,4(1x,e10.3))') ' Vmax/min,|Gmax|,|CS|',global_MAXX,
     %        global_MINX,global_MAXY,global_MAXCS
c         PRINT *,'MAX V,G,H,CS',global_MAXX,global_MAXY,global_MAXZ,
c     %        global_MAXCS
c         PRINT *,'MAX V,CS',global_MAXX,global_MAXCS
      endif
      RETURN
      END SUBROUTINE MAXVELOCITIES

      SUBROUTINE MAXVAR
**  Finds the maximum of variable defined below
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      IMPLICIT NONE
      include "mpif.h"
      double precision :: MAXVAR1,MAXVAR2,MINVAR1,MINVAR2
      double precision :: global_MAXVAR1,global_MAXVAR2
      double precision :: global_MINVAR1,global_MINVAR2
      double precision :: VAR1(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: VAR2(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
!- Define variable here, make sure deceration of VAR1/2 agrees dimensionally
      VAR1 = RH
      VAR2 = T
      MAXVAR1 = MAXVAL(VAR1)
      MAXVAR2 = MAXVAL(VAR2)
      MINVAR1 = MINVAL(VAR1)
      MINVAR2 = MINVAL(VAR2)
      call MPI_ALLREDUCE(maxVAR1,global_maxVAR1,1,
     %     MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(maxVAR2,global_maxVAR2,1,
     %     MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(minVAR1,global_minVAR1,1,
     %     MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(minVAR2,global_minVAR2,1,
     %     MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
         PRINT *,'MIN/MAX VAR1',global_MINVAR1,global_MAXVAR1
         PRINT *,'MIN/MAX VAR2',global_MINVAR2,global_MAXVAR2
      endif
      RETURN
      END SUBROUTINE MAXVAR

      SUBROUTINE CALCE
      USE global_constants
      USE grid_var_init
      USE fluid_var_init
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: pe,ke,ie
      double precision :: global_pe,global_ke,global_ie,global_etot
      double precision :: delpe,delke,delie,deletot,totcomp
      double precision :: vphi,vr
      pe = 0.d0
      ke = 0.d0
      ie = 0.d0
      global_pe = 0.d0
      global_ke = 0.d0
      global_ie = 0.d0
      global_etot = 0.d0
      delpe  = 0.d0
      delke  = 0.d0
      delie = 0.d0
      deletot = 0.d0
      totcomp = 0.d0

      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX
               pe   = pe - GRAV*MSOL*RH(I,J,K)/xxb(i)
               vphi = 0.5*(G(i,j,k)+G(i,j+1,k))*xxb(i)
               vr   = 0.5*(V(i,j,k)+V(i+1,j,k))
               ke   = ke + 0.5*RH(I,J,K)*((vphi**2.d0)+(vr**2.d0))
               ie   = ie+ RH(I,J,K)*T(I,J,K)*CV(I,J,K)
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(pe,global_pe,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ke,global_ke,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ie,global_ie,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      global_etot = global_pe+global_ke+global_ie
      delpe   = global_pe-oldpe
      delke   = global_ke-oldke
      delie   = global_ie-oldie
      deletot = global_etot-oldetot
      totcomp = delpe+delke+delie
      if(myid.eq.0) then
c         print *,'DE',delpe,delke,delie,deletot,totcomp
         write(*,'(A,2(e35.24))') 'DE',oldetot,global_etot,deletot
c         print *,'DE',oldetot,global_etot,deletot
      endif
      oldpe = global_pe
      oldke = global_ke
      oldie = global_ie
      oldetot = global_etot           
      RETURN
      END SUBROUTINE CALCE

      SUBROUTINE FINDPHOTOSPHERE
      USE global_constants
      USE grid_var_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: dr,totaltau
      DR = XXB(2)-XXB(1)
      do k=1,locNZ
        do j=1,locNY
           TOTALTAU = 0.0
           do i=UPPERBND(j,k),1,-1
              TOTALTAU = TOTALTAU+(xkapP(i,j,k)+SIG(i,j,k))*RH(i,j,k)*DR
              TAU(i,j,k) = TOTALTAU
           enddo         
           if(upperbnd(j,k).lt.locNX) then
              do i=UPPERBND(j,k),locNX
                 TAU(i,j,k) = 0.0
              enddo
           endif
        enddo
      enddo
      do k=1,locNZ
        do j=1,locNY
           photosphere(j,k) = upperbnd(j,k)
           do i=0,upperbnd(j,k)
              IF(tau(i,j,k).LE.(2./3.)) THEN
                 photosphere(j,k) = i
                 GOTO 100
              ENDIF
           enddo
 100       continue
        enddo
      enddo
      RETURN
      END SUBROUTINE FINDPHOTOSPHERE


      SUBROUTINE PERTURBATIONS
      USE global_constants
      USE grid_var_init
      USE fluid_var_init
      USE mpi_var_init
      USE input_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: AMP,phi_0,delta_phi,theta_0
      double precision :: delta_theta,gauss_fac1,gauss_fac2
!---  PERURBATIONS AZIMUTHAL VELOCITY AT PHOTOSPHERE
      SELECT CASE(perturb_typ)
      CASE (0) !---  LOCATE THE PHOTOSPHERE
         CALL FINDPHOTOSPHERE
!--- PERTURB AZIMUTHAL VELOCITY
         do k=1,locNZ
            do j=1,locNY 
               phi_0 = PI/4.0
               delta_phi = 1.d0/20.d0
               theta_0 = 0.0
               delta_theta = 1.d0/20.d0
               gauss_fac1 = ((xya(i)-phi_0)/delta_phi)**2.0
               gauss_fac2 = ((xzb(i)-theta_0)/delta_theta)**2.0
               AMP = G(i,j,k)/2.0
               G(photosphere(j,k),j,k)=G(i,j,k)+AMP*exp(gauss_fac1)*
     %              exp(gauss_fac2)
            enddo
         enddo
!---  FIND BOUNDARY CONDITIONS
         CALL BOUNDS
      CASE (1:)
         print *,'perturb_typ not defined',perturb_typ
      CASE DEFAULT
         print *,'perturb_typ also not defined'
      END SELECT
      RETURN
      END SUBROUTINE PERTURBATIONS



      SUBROUTINE CHECK_MAX_DT
**  Routine to check(and print) what is setting the maximum timestep
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE force_var_init
      USE grid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: V1,V2,V3,ZC1,ZC2,ZC3
      double precision :: localDT,checkDT
      double precision :: XNU,ZCF
      double precision :: min_Vdt,min_Gdt,min_Hdt,min_viscdt
**  if fdelt=0. >> delt=dtmax
      IF (FDELT.EQ.0.0) THEN
         write(*,'(A)') 'Timestep set by fdelt',fdelt
      ELSE
         min_Vdt    = 100000.0
         min_Gdt    = 100000.0
         min_Hdt    = 100000.0
         min_viscdt = 100000.0
         localDT=DELT*100.d0/FDELT
**  SOUND VELOCITY
         IF(ICSTYP.NE.1.AND.ICSTYP.NE.2) then
            CALL SOUND
         ENDIF
         localDT=DTMAX
!-- COURANT CONDITION IN ALL THREE DIMENTIONS
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,locNX
                  V1=ABS(V(I+1,J,K))
                  V2=ABS(G(I,J+1,K))
                  V3=ABS(H(I,J,K+1))
!--x-velocity
                  ZC1=DXA(I+1)/(V1+CS(i,j,k))
                  min_Vdt = MIN(min_Vdt,zc1)
!--y-velocity
                  ZC2=DYA(J+1)/(V2+CS(i,j,k)/geoxg(i)/geozg(k))
                  min_Gdt = MIN(min_Gdt,zc2)
!--z-velocity
                  ZC3=DZA(K+1)/(V3+CS(i,j,k)/geoxh(i))
                  min_Hdt = MIN(min_Hdt,zc3)
!--calculate min
                  localDT=MIN(ZC1,ZC2,ZC3,localDT)
                  IF(IVIS.GT.0) THEN
                     XNU= ETA(I,J,K) / RH(I,J,K)
                     ZCF= (1.d0/8.d0)/ XNU /
     &                    (DXA(I+1)**(-2.d0) +
     &                    (DYA(J+1)*geoxg(i)*geozg(k))**(-2.d0) +
     &                    (DZA(K+1)*geoxh(i))**(-2.d0))
                     min_viscdt = MIN(min_viscdt,zcf)
                     localDT=MIN(ZCF,localDT)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
!--- CHECK MAX TIMESTEP
         localDT=  MIN (localDT,DTMAX) * FDELT
!--- find the globally lowest timestep
         call MPI_ALLREDUCE(localDT,checkDT,1,MPI_DOUBLE_PRECISION,
     %        MPI_MIN,MPI_COMM_WORLD,ierr)
         
         write(*,'(A,I8,5(1x,e13.6))') 
     &        'timestep check',myid,min_Vdt*FDELT,min_Gdt*FDELT,
     &        min_Hdt*FDELT,min_viscdt*FDELT,checkDT

      ENDIF
      RETURN
      END SUBROUTINE CHECK_MAX_DT


      SUBROUTINE CLEAN_STOP
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      USE hypre_var_init
      USE hypre_var_alloc
      IMPLICIT NONE
      include "mpif.h"
      include 'bHYPRE_SStructVariable.inc'
      character stopfile*8
      character proc_num*3
      print *,'EXCECUTING CLEAN_STOP MYID=',myid
      write(proc_num,'(i3.3)') myid
      stopfile = 'CRSH_'//TRIM(ADJUSTL(proc_num))
      OPEN (12,FILE=stopfile,STATUS='UNKNOWN',
     $     FORM='UNFORMATTED')
      REWIND(12)
      CALL PRINTOUT(1)
      CLOSE(12)
      if(irad.eq.2) then
         call bHYPRE_SStructGrid_deleteRef_f( grid, except )
         call bHYPRE_SStructGraph_deleteRef_f( graph, except )
         call bHYPRE_SStructStencil_deleteRef_f( stencil, except )
         call bHYPRE_MPICommunicator_deleteRef_f(bHYPRE_mpicomm,except)
      endif
      call MPI_FINALIZE(ierr)
      STOP
      RETURN
      END SUBROUTINE CLEAN_STOP


      FUNCTION ECCENTRIC_ANOMALY(time,ecc,period)
!      returns E=eccentric anomaly as a function of
!       time and eccentricity e. The period is used to set the accuracy.
      USE global_constants
      USE mpi_var_init
      implicit none
      integer :: i
      double precision :: eccentric_anomaly,time,ecc,period
      double precision :: M,E,dE
!-mean anomaly
      M = (2.0*PI/period)*(time)
!initial guess for eccentric anomaly
      E=(ecc*sin(M)/(1.0-sin(M+ecc)+sin(M)))+M
      if (ecc.gt.0.95) then
         print *,'something should change for this high ecc'
         print *,'see andrews initial kepler.pro IDL routine'
         stop
      endif
! using Newtons method
      dE=E
      do i=1,100
         dE=-(E-M-ecc*sin(E))/(1-ecc*cos(E))
         E=E+dE
c         write(*,'(A,i8,3(1x,e12.3))') 'E=',i,E,dE,(0.001/period)
         if(dE.le.(0.001/period)) then
            GOTO 100
         endif
      enddo
      if(myid.eq.0) then
         print *,'Eccentric anomoly overflow:',dE,(0.001/period),
     $        time,E,M
      endif
 100  CONTINUE
!-- return E
      eccentric_anomaly = E
      RETURN
      END

      SUBROUTINE ROCHE_OVERFLOW
      USE input_init
      USE mpi_var_init
      USE mpi_grid_init
      USE fluid_var_init
      USE force_var_init
      USE grid_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k,overflow_index
      double precision :: loc_Mloss,global_Mloss,dV,HP,RH_reset
      loc_Mloss = 0.d0 !-set all to zero initially so I can still do a allreduce
      overflow_index = 92
      DO K=1,locNZ
         DO J=1,locNY
            if(upperbnd(j,k).gt.overflow_index) then
               DO I=overflow_index+1,locNX
                  dV = (xxb(i)**(2.d0))*cos(xzb(k))*dxa(i)*dya(j)*dza(k)
!calculate the density ppbou WILL set 
                  HP = -RGAS*T(overflow_index,j,k)/
     $                 (xmue(overflow_index,j,k)*gravx(i,j,k))
                  RH_reset = RH(overflow_index,j,k)*
     $                 exp(-(xxb(i)-xxb(overflow_index))/HP)
                  if(RH(i,j,k).gt.RH_reset) then
                     loc_Mloss = loc_Mloss + (RH(i,j,k)-RH_reset)*dV
                  endif
               enddo
!-set upperbnd to the max index
               upperbnd(j,k) = overflow_index
            endif
         enddo
      enddo
      call MPI_ALLREDUCE(loc_Mloss,global_Mloss,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      Mtransfer = Mtransfer + global_Mloss
      RETURN
      END SUBROUTINE ROCHE_OVERFLOW


      SUBROUTINE MISC_CHECKS
!- general routine to perform all the various misc.
! checks on consistancies that would otherwise have
! to be performed in the middle of loops, etc.
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
      IMPLICIT NONE
      include "mpif.h"
      include 'bHYPRE_SStructVariable.inc'

!--ppenergy doesn't define boundaries in the y-direction
      IF(MPIupper.eq.MPI_PROC_NULL.or.MPIlower.eq.MPI_PROC_NULL) then
         print *,'update y-proc boundaries for ppenergy'
         call clean_stop
      ENDIF

!-Geoff's ietot variable removed
      if (ietot.ge.1) then
         print *,'removed from ttadv'
         call clean_stop
      endif
!-Eliza's opacities only available for modtyp=2
      if((IKAPTYP.ge.15.and.ikaptyp.le.19).and.modtyp.ne.2) then
         print *,'Elizas opacities only available for modtyp=2'
         call clean_stop
      endif
!-wavelength dependent source term only available for ikaptyp=20
      if(ibdrad.eq.10.and.(IKAPTYP.ne.20.or.ikaptyp.ne.21)) then
         print *,'wavelength dependent source term only available'
         print *,'for ikaptyp=20'
         call clean_stop
      endif
      RETURN
      END SUBROUTINE MISC_CHECKS
