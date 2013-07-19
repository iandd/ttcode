      SUBROUTINE KAPPA  !---  Opacity
!--- ikaptyp = 1 : Constant Opacity = ZKAPMULT (from ipp)
!--- ikaptyp = 2 : Opacity tables from Pollack (1985) + AF
!--- ikaptyp = 3 : Dougs opacity functions
!--- ikaptyp = 4 : Freedman(2007) opacities (updated 4/2008) normal table
!--- ikaptyp = 5 : Freedman(2007) opacities (updated 4/2008) without TIVO
!--- ikaptyp = 6 : Freedman(2007) opacities (updated 4/2008) normal table, column opacities
!--- ikaptyp = 7 : Freedman(2007) opacities (updated 4/2008) without TIVO, column opacities
!--- ikaptyp = 8 : H- bound-free absorption opacity
!--- ikaptyp = 9 : Burrows opacities (calc. 5/2008)  with TIVO
!--- ikaptyp = 10 : Burrows opacities (calc. 5/2008) without TIVO
!--- ikaptyp = 11 : Burrows opacities (calc. 5/2008) without TIVO + extra absorber
!--- ikaptyp = 12 : Burrows opacities (calc. 5/2008) without TIVO for xkapR,xkapP, constant gamma
!--- ikaptyp = 13 : constant kapP=zkapmult, constant gamma
!--- ikaptyp = 14 : Burrows opacities (calc. 5/2008) without TIVO for xkapP, xkapA=xkapR=gamma^4*xkapP, constant gamma
!--- ikaptyp = 15 : Eliza Kempton's opacities with CO/CH4 dependence (calc. 4/2011)
!--- ikaptyp = 16 : Eliza Kempton's opacities with CO/CH4 dependence (calc. 4/2011) + extra absorb
!--- ikaptyp = 17 : Eliza Kempton's solar + purly CO opacities + extra absorb
!--- ikaptyp = 18 : Eliza Kempton's solar + purly CH4 opacities + extra absorb
!--- ikaptyp = 19 : Eliza Kempton's equilibrium opacities + extra absorb
!--- ikaptyp = 20 : Burrows opacities with wavelength dependance (updated 6/2011) with TIVO
!--- ikaptyp = 21 : Burrows opacities with wavelength dependance (updated 6/2011) without TIVO
!--- ikaptyp = 22 : Burrows opacities full rad. transfer sol. (updated 6/2011) with TIVO
!--- ikaptyp = 23 : Burrows opacities, full rad. transfer sol. (updated 6/2011) without TIVO
!--- ikaptyp = 24 : Burrows opacities, full rad. transfer sol. without TIVO + extra abs in stellarinput_heating opacity
!--- ikaptyp = 25 : Burrows opacities, full rad. transfer sol. without TIVO + Rayleigh SCat + extra abs for all opac
      USE rad_var_init
      USE input_init
      USE fluid_var_init
      USE scalar_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE grid_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: deltaR,testA
      double precision :: equil_XCO,XCOeq
      SELECT CASE(ikaptyp)
      CASE(1)
!--- ikaptyp = 1 : Constant Opacity = ZKAPMULT (from ipp)
         xkapR = ZKAPMULT
         xkapP = ZKAPMULT
         xkapA = ZKAPMULT
         dxkapP_dT = 0.d0
         dxkapA_dT = 0.d0         
      CASE(2)
!--- ikaptyp = 2 : Opacity tables from Pollack et al (1985) and
!---     Alexander & Ferguson (1994). The Alexander tables don't go
!---     low enough.. coupled to a stellar interior table on the
!---     high end (from Peter)
         if(MODTYP.eq.2) then
            DO k=0,locNZ+1
               DO j=0,locNY+1
                  DO i=0,upperbnd(j,k)
                     call opcPAF(T(i,j,k),RH(i,j,k),xkapR(i,j,k))
                  ENDDO
!---   Outside area should be uniform
                  DO i=upperbnd(j,k)+1,locNX+1
                     xkapR(i,j,k) = 1.d0
                  ENDDO
               ENDDO
            ENDDO
         else
            DO k=0,locNZ+1
               DO j=0,locNY+1
                  DO i=0,locNX+1
                     call opcPAF(T(i,j,k),RH(i,j,k),xkapR(i,j,k))
                  ENDDO
               ENDDO
            ENDDO
         endif
      CASE(3)
!--- ikatpyp = 3 : Dougs opacity functions
!---     smoothed version of the tables
         DO k=0,locNZ+1
            DO j=0,locNY+1
               DO i=0,locNX+1
                  call OPLIN0(T(i,j,k),RH(i,j,k),xkapR(i,j,k))
               ENDDO
            ENDDO
         ENDDO
         xkapP = xkapR
         xkapA = xkapR
         dxkapP_dT = 0.d0
         dxkapA_dT = 0.d0
      CASE(4:5)
!---  Freedman, Marley, & Lodders (2007) opacities (updated March 2008)
!--- ikaptyp = 4 : Freedman(2007) opacities (updated 4/2008) normal table
!--- ikaptyp = 5 : Freedman(2007) opacities (updated 4/2008) without TiVO
!--     these opacities are in PG vs T space, so call pressure function
         CALL PRESSG
         if(MODTYP.eq.2) then
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,upperbnd(j,k)
                     call FreedmanOPC(T(i,j,k),PG(i,j,k),xkapP(i,j,k),
     $                    xkapR(i,j,k),xkapA(i,j,k),
     $                    dxkapP_dT(i,j,k),dxkapA_dT(i,j,k))
                  ENDDO
!---   Outside area should be uniform
                  DO i=upperbnd(j,k)+1,locNX+2
                     xkapR(i,j,k) = xkapR(upperbnd(j,k),j,k)
                     xkapP(i,j,k) = xkapP(upperbnd(j,k),j,k)
                     xkapA(i,j,k) = xkapA(upperbnd(j,k),j,k)
                     dxkapP_dT(i,j,k) = dxkapP_dT(upperbnd(j,k),j,k)
                     dxkapA_dT(i,j,k) = dxkapA_dT(upperbnd(j,k),j,k)
                  ENDDO
               ENDDO
            ENDDO
         else
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,locNX+2
                     call FreedmanOPC(T(i,j,k),PG(i,j,k),xkapP(i,j,k),
     $                    xkapR(i,j,k),xkapA(i,j,k),
     $                    dxkapP_dT(i,j,k),dxkapA_dT(i,j,k))
                  ENDDO
               ENDDO
            ENDDO
         endif
      CASE(6:7)
!--- ikaptyp = 6 : Freedman(2007) opacities (updated 4/2008) normal table, column opacities
!--- ikaptyp = 7 : Freedman(2007) opacities (updated 4/2008) without TIVO, column opacities
         CALL PRESSG
         if(MODTYP.eq.2) then
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,upperbnd(j,k)

c--- testing
                     call FreedmanOPC(T(i,j,k),PG(i,j,k),xkapP(i,j,k),
     $                    xkapR(i,j,k),testA,
     $                    dxkapP_dT(i,j,k),dxkapA_dT(i,j,k))
c--- testing
                     deltaR = xxb(upperbnd(j,k))-xxb(i)
                     call columnOPC(T(i,j,k),PG(i,j,k),RH(i,j,k),
     $                    XMUE(i,j,k),deltaR,xkapP(i,j,k),xkapR(i,j,k),
     $                    xkapA(i,j,k),dxkapP_dT(i,j,k),
     $                    dxkapA_dT(i,j,k))                     
c                     if(myid.eq.3.and.k.eq.5.and.j.eq.5.and.
c     $                    num_iter.lt.10) then
c                        write(*,'(A,I8,1x,3(1x,e12.6))') 
c     $                       'test',i,xkapA(i,j,k),testA
c     $                       -23.9526937 + LOG10(xmue(i,j,k)) -
c     $                       LOG10(RH(i,j,k)) - 
c     $                       LOG10(deltaR+1.0)
c                     endif



                  ENDDO
!---   Outside area should be uniform
                  DO i=upperbnd(j,k)+1,locNX+2
                     xkapR(i,j,k) = xkapR(upperbnd(j,k),j,k)
                     xkapP(i,j,k) = xkapP(upperbnd(j,k),j,k)
                     xkapA(i,j,k) = xkapA(upperbnd(j,k),j,k)
                     dxkapP_dT(i,j,k) = dxkapP_dT(upperbnd(j,k),j,k)
                     dxkapA_dT(i,j,k) = dxkapA_dT(upperbnd(j,k),j,k)
                  ENDDO
               ENDDO
            ENDDO
         else
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,locNX+2
                     deltaR = xxb(upperbnd(j,k))-xxb(i)
                     call columnOPC(T(i,j,k),PG(i,j,k),RH(i,j,k),deltaR,
     $                    xkapP(i,j,k),xkapR(i,j,k),xkapA(i,j,k),
     $                    dxkapP_dT(i,j,k),dxkapA_dT(i,j,k))
                  ENDDO
               ENDDO
            ENDDO
         endif
      CASE(8) 
!--- ikaptyp=6: H- bound-free absorption opacity
!---   function taken from Andrews 'Stellar Interiors' book
         DO k=-1,locNZ+2
            DO j=1,locNY+2
               DO i=-1,upperbnd(j,k)
                  xkapR(I,J,K)=0.95d0*(RH(i,j,k)**(0.5))*
     %                 ((T(i,j,k)/2500.d0)**(9.0))
                  dxkapP_DT(I,J,K)=0.00342*(RH(i,j,k)**(0.5))*
     %                 ((T(i,j,k)/2500.d0)**(8.0))
                  dxkapA_DT(I,J,K)=0.00342*(RH(i,j,k)**(0.5))*
     %                 ((T(i,j,k)/2500.d0)**(8.0))
               ENDDO
!---  Outside area should be uniform
               DO i=upperbnd(j,k)+1,locNX+2
                  xkapR(i,j,k) = xkapR(upperbnd(j,k),j,k)
                  dxkapP_DT(I,J,K)=dxkapP_DT(upperbnd(j,k),j,k)
                  dxkapA_DT(I,J,K)=dxkapA_DT(upperbnd(j,k),j,k)
               ENDDO
            ENDDO
         ENDDO
         xkapP = xkapR
         xkapA = xkapR
      CASE(9:12)
!--- ikaptyp = 9 : Burrows opacities (calc. 5/2008) normal table
!--- ikaptyp = 10 : Burrows opacities (calc. 5/2008) without TIVO
!--- ikaptyp = 11 : Burrows opacities (calc. 5/2008) without TIVO + extra absorber
!--- ikaptyp = 12 : Burrows opacities (calc. 5/2008) without TIVO for xkapR, constant gamma
         if(MODTYP.eq.2) then
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,upperbnd(j,k)
                     call BurrowsOPC(T(i,j,k),RH(i,j,k),xkapP(i,j,k),
     $                    xkapR(i,j,k),xkapA(i,j,k),xkapW(i,j,k,:),
     $                    xkapLW(i,j,k,:),local_BB(i,j,k,:),
     $                    dxkapP_dT(i,j,k),dxkapA_dT(i,j,k),i,j,k)
                  ENDDO
!---   Outside area should be uniform
                  DO i=upperbnd(j,k)+1,locNX+2
                     xkapR(i,j,k) = xkapR(upperbnd(j,k),j,k)
                     xkapP(i,j,k) = xkapP(upperbnd(j,k),j,k)
                     xkapA(i,j,k) = xkapA(upperbnd(j,k),j,k)
                     dxkapP_dT(i,j,k) = dxkapP_dT(upperbnd(j,k),j,k)
                     dxkapA_dT(i,j,k) = dxkapA_dT(upperbnd(j,k),j,k)
                  ENDDO
               ENDDO
            ENDDO
         else
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,locNX+2
                     call BurrowsOPC(T(i,j,k),RH(i,j,k),xkapP(i,j,k),
     $                    xkapR(i,j,k),xkapA(i,j,k),xkapW(i,j,k,:),
     $                    xkapLW(i,j,k,:),local_BB(i,j,k,:),
     $                    dxkapP_dT(i,j,k),dxkapA_dT(i,j,k),i,j,k)
                  ENDDO
               ENDDO
            ENDDO
         endif
      CASE(13)
!--- ikaptyp = 13 : constant kapP=zkapmult, constant gamma
         xkapR = ZKAPMULT
         xkapP = ZKAPMULT
         xkapA = (gamma_OPC**(4.d0))*ZKAPMULT
         dxkapP_dT = 0.d0
         dxkapA_dT = 0.d0         
      CASE(14)
!--- ikaptyp = 14 : Burrows opacities (calc. 5/2008) without TIVO for xkapP, xkapA=xkapR=gamma^4*xkapP, constant gamma
         if(MODTYP.eq.2) then
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,upperbnd(j,k)
                     call BurrowsOPC(T(i,j,k),RH(i,j,k),xkapP(i,j,k),
     $                    xkapR(i,j,k),xkapA(i,j,k),xkapW(i,j,k,:),
     $                    xkapLW(i,j,k,:),local_BB(i,j,k,:),
     $                    dxkapP_dT(i,j,k),dxkapA_dT(i,j,k),i,j,k)
                  ENDDO
!---   Outside area should be uniform
                  DO i=upperbnd(j,k)+1,locNX+2
                     xkapR(i,j,k) = xkapR(upperbnd(j,k),j,k)
                     xkapP(i,j,k) = xkapP(upperbnd(j,k),j,k)
                     xkapA(i,j,k) = xkapA(upperbnd(j,k),j,k)
                     dxkapP_dT(i,j,k) = dxkapP_dT(upperbnd(j,k),j,k)
                     dxkapA_dT(i,j,k) = dxkapA_dT(upperbnd(j,k),j,k)
                  ENDDO
               ENDDO
            ENDDO
         else
            DO k=-1,locNZ+2
               DO j=-1,locNY+2
                  DO i=-1,locNX+2
                     call BurrowsOPC(T(i,j,k),RH(i,j,k),xkapP(i,j,k),
     $                    xkapR(i,j,k),xkapA(i,j,k),xkapW(i,j,k,:),
     $                    xkapLW(i,j,k,:),local_BB(i,j,k,:),
     $                    dxkapP_dT(i,j,k),dxkapA_dT(i,j,k),i,j,k)
                  ENDDO
               ENDDO
            ENDDO
         endif
      CASE(15:18)
!---  ikaptyp = 15 : Eliza Kempton's opacities with CO/CH4 dependence (calc. 4/2011)
!---  ikaptyp = 16 : Eliza Kempton's opacities with CO/CH4 dependence (calc. 4/2011) + extra absorb
!---  ikaptyp = 17 : Eliza Kempton's purly CO opacities + extra absorb
!---  ikaptyp = 18 : Eliza Kempton's purly CH4 opacities + extra absorb
!--     these opacities are in PG vs T space, so call pressure function
         CALL PRESSG
         DO k=-1,locNZ+2
            DO j=-1,locNY+2
               DO i=-1,upperbnd(j,k)
                  call ElizasOPC(T(i,j,k),PG(i,j,k),scalar(1,i,j,k),
     $                 xkapP(i,j,k),xkapR(i,j,k),xkapA(i,j,k),
     $                 dxkapP_dT(i,j,k),dxkapA_dT(i,j,k),i,j,k)
               ENDDO
!---  Outside area should be uniform
               DO i=upperbnd(j,k)+1,locNX+2
                  xkapR(i,j,k) = xkapR(upperbnd(j,k),j,k)
                  xkapP(i,j,k) = xkapP(upperbnd(j,k),j,k)
                  xkapA(i,j,k) = xkapA(upperbnd(j,k),j,k)
                  dxkapP_dT(i,j,k) = dxkapP_dT(upperbnd(j,k),j,k)
                  dxkapA_dT(i,j,k) = dxkapA_dT(upperbnd(j,k),j,k)
               ENDDO
            ENDDO
         ENDDO
      CASE(19)
!---  ikaptyp = 19 : Eliza Kempton's equilibrium chemistry opacities
         CALL PRESSG
         DO k=-1,locNZ+2
            DO j=-1,locNY+2
               DO i=-1,upperbnd(j,k)
!- calculate the equil value
                  equil_XCO = XCOeq(i,j,k)
                  call ElizasOPC(T(i,j,k),PG(i,j,k),equil_XCO,
     $                 xkapP(i,j,k),xkapR(i,j,k),xkapA(i,j,k),
     $                 dxkapP_dT(i,j,k),dxkapA_dT(i,j,k),i,j,k)
               ENDDO
!---  Outside area should be uniform
               DO i=upperbnd(j,k)+1,locNX+2
                  xkapR(i,j,k) = xkapR(upperbnd(j,k),j,k)
                  xkapP(i,j,k) = xkapP(upperbnd(j,k),j,k)
                  xkapA(i,j,k) = xkapA(upperbnd(j,k),j,k)
                  dxkapP_dT(i,j,k) = dxkapP_dT(upperbnd(j,k),j,k)
                  dxkapA_dT(i,j,k) = dxkapA_dT(upperbnd(j,k),j,k)
               ENDDO
            ENDDO
         ENDDO
      CASE(20:25)
!--- ikaptyp = 21 : Burrows opacities with wavelength dependance (updated 6/2011) without TIVO
!--- ikaptyp = 22 : Burrows opacities full rad. transfer sol. (updated 6/2011) with TIVO
!--- ikaptyp = 23 : Burrows opacities, full rad. transfer sol. (updated 6/2011) without TIVO
!--- ikaptyp = 24 : Burrows opacities, full rad. transfer sol. (updated 6/2011) without TIVO + extra abs in stellarinput_heating opacity
!--- ikaptyp = 25 : Burrows opacities, full rad. transfer sol. without TIVO + Rayleigh SCat + extra abs for all opac
         DO k=-1,locNZ+2
            DO j=-1,locNY+2
               DO i=-1,upperbnd(j,k)
                  call BurrowsOPC(T(i,j,k),RH(i,j,k),xkapP(i,j,k),
     $                 xkapR(i,j,k),xkapA(i,j,k),xkapW(i,j,k,:),
     $                 xkapLW(i,j,k,:),local_BB(i,j,k,:),
     $                 dxkapP_dT(i,j,k),dxkapA_dT(i,j,k),i,j,k)
               ENDDO
!---  Outside area should be uniform
c               DO i=upperbnd(j,k)+1,locNX+2
               DO i=upperbnd(j,k),locNX+2
                  xkapR(i,j,k) = xkapR(upperbnd(j,k)-1,j,k)
                  xkapP(i,j,k) = xkapP(upperbnd(j,k)-1,j,k)
                  xkapA(i,j,k) = xkapA(upperbnd(j,k)-1,j,k)
                  xkapW(i,j,k,:) = xkapW(upperbnd(j,k)-1,j,k,:)
                  if(ikaptyp.eq.22.or.ikaptyp.eq.23.or.
     $                 ikaptyp.eq.24.or.ikaptyp.eq.25)
     $                 xkapLW(i,j,k,:) = xkapLW(upperbnd(j,k)-1,j,k,:)
                  dxkapP_dT(i,j,k) = dxkapP_dT(upperbnd(j,k)-1,j,k)
                  dxkapA_dT(i,j,k) = dxkapA_dT(upperbnd(j,k)-1,j,k)
               ENDDO
            ENDDO
         ENDDO
      CASE(26:)
         print *,'ikaptyp=',ikaptyp,'is not avalible:KAPPA'
         stop
      CASE DEFAULT
         print *,'also not defined: KAPPA (pprphy.f)',ikaptyp
         stop
      END SELECT
!-- calculate the spatially dependent ratio of opacities      
      gam_kap = (xkapA/xkapP)**(0.25)
      RETURN
      END SUBROUTINE KAPPA

      SUBROUTINE FLIM
**    FLUX-LIMITER
      USE input_init
      USE rad_var_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      USE mpi_var_init
      USE mpi_grid_init
      USE force_var_init
!      USE hypre_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: stat(MPI_STATUS_SIZE)
      integer :: i,j,k,shift_coord,shift_dir,xmax_index
      double precision,parameter :: xdmin=1.0d-50
      double precision :: GRADX,GRADY,GRADZ,GRAD,RRX,RRY,RRZ
      double precision :: u(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: HP(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: tau_HP(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
**  B.C. for Variables (RH,T), excluding Rad. Energy 
**  all these routines must be called before SUBKAP and DIFF
** this is the only place where BOUNDS is called during the rad substep **
      if(irad.eq.1.or.irad.eq.4) then !---ORIGINAL SOR SOLVER
         CALL BOUNDS
**  B.C. for Rad. Energy 
         CALL RADIATIVEBOUND
**  Opacities and scattering
         CALL KAPPA
         CALL SUBSIG
**  RHQX,Y,Z are needed for DIFF and for FLIM
         CALL QDENS
      endif
!-- only calculate Opacities and scattering on first iteration
      if((irad.eq.2.or.irad.eq.3).and.energy_iter.le.1) then !---COUPLED ENERGY SOLVERS
c      if(irad.eq.2.or.irad.eq.3) then !---COUPLED ENERGY SOLVERS (DOES THIS HELP INIT. CONVERGNCE)
         CALL KAPPA
         CALL SUBSIG
**  RHQX,Y,Z are needed for DIFF and for FLIM
         CALL QDENS
      endif
!--- CALCULATE THE FLUX LIMITER 
      SELECT CASE(IFLIMTYP)
      CASE(1)
!-- CONSTANT FLUX-LIMITER (=1./3., OPTICALLY THICK)
         FLIMX=1.d0/3.d0
         FLIMY=1.d0/3.d0
         FLIMZ=1.d0/3.d0
      CASE(2)
!-- LEVERMORE & POMRANING (1981)
**    Flux-Limiter from Pomraning
!-- x-direction
         DO K=1,locNZ
            DO J=1,locNY
               if(MPIright.eq.MPI_PROC_NULL) then
                  xmax_index = upperbnd(J,K)-1 
               else
                  xmax_index = locNX+1
               endif
               DO I=1,xmax_index
                  GRADX=(ER(I,J,K)-ER(I-1,J,K))/DXB(I)
                  GRADY=(DDX1(I)*(ER(I,J+1,K)-ER(I,J-1,K))
     &                 +DDX0(I)*(ER(I-1,J+1,K)-ER(I-1,J-1,K)))
     &                 /(geoxga(i)*geozg(k)*(xyb(j+1)-xyb(j-1)))
                  GRADZ=(DDX1(I)*(ER(I,J,K+1)-ER(I,J,K-1))
     &                 +DDX0(I)*(ER(I-1,J,K+1)-ER(I-1,J,K-1)))
     &                 /(geoxha(i)*(xzb(k+1)-xzb(k-1)))
                  GRAD= SQRT(GRADX*GRADX+GRADY*GRADY+GRADZ*GRADZ)
     &                 +xdmin
                  RRX=2.d0*GRAD/( RHQX(I,J,K)*
     &                 (DDX1(I)*(xkapR(I,J,K)+sig(i,j,k))
     &                 +DDX0(I)*(xkapR(i-1,j,k)+sig(i,j,k)))
     &                 *(ER(I,J,K)+ER(I-1,J,K)) )

c impose a boundary condition on RXX
c                  if(i.ge.upperbnd(j,k)-10) then
c                     RRX=(((10.**7.)-saveRRX(upperbnd(j,k)-10,j,k))/
c     $                    10.d0)*(1.d0*i)+saveRRX(upperbnd(j,k)-10,j,k)
c                  endif

c-- took the limits out Nov 8,2007 they seemed to not give correct optically thin limit
c                 if(RRX.gt.10.**7.) then
c                     RRX=10.**7.
c                  endif
                  FLIMX(I,J,K)=(2.d0+RRX)/(6.d0+3.d0*RRX+RRX*RRX)
                  saveRRX(i,j,k) = RRX
               ENDDO
               if(MPIright.eq.MPI_PROC_NULL) then
                  DO i=upperbnd(j,k),locNX+1
                     FLIMX(I,J,K)=FLIMX(upperbnd(j,k)-1,J,K)
                     saveRRX(i,j,k) = saveRRX(upperbnd(j,k)-1,j,k)
                  ENDDO
               endif
            ENDDO
         ENDDO

         debugflag=0


!-- y-direction
         DO K=1,locNZ
            DO J=1,locNY+1
               DO I=1,xmax_index
                  GRADX=(DDY1(J)*(ER(I+1,J,K)-ER(I-1,J,K))
     &                 +DDY0(J)*(ER(I+1,J-1,K)-ER(I-1,J-1,K)))
     &                 /(xxb(i+1)-xxb(i-1))
                  GRADY=(ER(I,J,K)-ER(I,J-1,K))
     &                 /(geoxg(i)*geozg(k)*DYB(J))
                  GRADZ=(DDY1(J)*(ER(I,J,K+1)-ER(I,J,K-1))
     &                 +DDY0(J)*(ER(I,J-1,K+1)-ER(I,J-1,K-1)))
     &                 /((xzb(k+1)-xzb(k-1))*geoxh(i))
                  GRAD= SQRT(GRADX*GRADX+GRADY*GRADY+GRADZ*GRADZ)
     &                 +xdmin
                  RRY = 2.d0 * GRAD
     &                 /( RHQY(I,J,K)*
     &                 (DDY1(J)*(xkapR(i,j,k)+sig(i,j,k))
     &                 +DDY0(J)*(xkapR(I,J-1,K)+sig(i,j,k)))
     &                 *(ER(I,J,K)+ER(I,J-1,K)) )
c                  if(RRY.gt.10.**7.) then
c                     RRY=10.**7.
c                  endif
                  FLIMY(I,J,K)=(2.d0+RRY)/(6.d0+3.d0*RRY+RRY*RRY)
               ENDDO
               if(MPIright.eq.MPI_PROC_NULL) then
                  DO i=upperbnd(j,k),locNX
                     FLIMY(I,J,K)=FLIMY(upperbnd(j,k)-1,J,K)
                  ENDDO
               endif
            ENDDO
         ENDDO
!--   z-direction
         DO K=1,locNZ+1
            DO J=1,locNY
               DO I=1,xmax_index
                  GRADX=(DDZ1(K)*(ER(I+1,J,K)-ER(I-1,J,K))
     &                 +DDZ0(K)*(ER(I+1,J,K-1)-ER(I-1,J,K-1)))
     &                 /(xxb(i+1)-xxb(i-1))
                  GRADY=(DDZ1(K)*(ER(I,J+1,K)-ER(I,J-1,K))
     &                 +DDZ0(K)*(ER(I,J+1,K-1)-ER(I,J-1,K-1)))
     &                 /((xyb(j+1)-xyb(j-1))*geoxg(i)*geozga(k))
                  GRADZ=(ER(I,J,K)-ER(I,J,K-1))
     &                 /(geoxh(i)*DZB(K))
                  GRAD= SQRT(GRADX*GRADX+GRADY*GRADY+GRADZ*GRADZ)
     &                 +xdmin
                  RRZ= 2.d0 * GRAD
     &                 /( RHQZ(I,J,K)*
     &                 (DDZ1(K)*(xkapR(i,j,k)+sig(i,j,k))
     &                 +DDZ0(K)*(xkapR(I,J,K-1)+sig(i,j,k)))
     &                 *(ER(I,J,K)+ER(I,J,K-1)) )
c                  if(RRZ.gt.10.**7.) then
c                     RRZ=10.**7.
c                  endif
                  FLIMZ(I,J,K)=(2.d0+RRZ)/(6.d0+3.d0*RRZ+RRZ*RRZ)
               ENDDO
               if(MPIright.eq.MPI_PROC_NULL) then
                  DO i=upperbnd(j,k),locNX
                     FLIMZ(I,J,K)=FLIMZ(upperbnd(j,k)-1,J,K)
                  ENDDO
               endif
            ENDDO
         ENDDO
      CASE(3)
         HP = PG/(RH*abs(gravx))
         tau_HP = xkapR*RH*HP
         DO K=1,locNZ
            DO J=1,locNY
               if(MPIright.eq.MPI_PROC_NULL) then
                  xmax_index = upperbnd(J,K)-1 
               else
                  xmax_index = locNX+1
               endif
               DO I=1,xmax_index
                  GRADX=(DDZ1(K)*(ER(I+1,J,K)-ER(I-1,J,K))
     &                 +DDZ0(K)*(ER(I+1,J,K-1)-ER(I-1,J,K-1)))
     &                 /(xxb(i+1)-xxb(i-1))
                  flimx(i,j,k) = ((1.d0-exp(-tau_HP(i,j,k)))*1.d0/3.d0)+
     $                 (exp(-tau_HP(i,j,k))*xkapR(i,j,k)*
     $                 rh(i,j,k)*ER(i,j,k)/abs(gradx))
                  if(num_iter.le.1.and.myid.eq.0.and.j.eq.1.and.k.eq.1)
     $                 then
                     write(*,'(I8,10(1x,e12.4))') i,flimx(i,j,k),
     $                    tau_hp(i,j,k),hp(i,j,k),PG(i,j,k),
     $                    gravx(i,j,k),abs(gravx(i,j,k)),
     $                    gradx,xkapR(i,j,k),rh(i,j,k),er(i,j,k)
                  endif
               enddo
               if(MPIright.eq.MPI_PROC_NULL) then
                  DO i=upperbnd(j,k),locNX+1
                     FLIMX(I,J,K)=FLIMX(upperbnd(j,k)-1,J,K)
                  ENDDO
               endif
            enddo
         enddo
         flimy = 1.d0/3.d0
         flimz = 1.d0/3.d0
      CASE(4:)
         print *,'iflimtyp=',iflimtyp,' is not avalible'
         stop
      CASE DEFAULT
         print *,'also not defined'
         stop
      END SELECT 
      RETURN
      END SUBROUTINE FLIM

      SUBROUTINE SUBSIG
**  Scattering Coefficient
      USE rad_var_init
      USE input_init
      IMPLICIT NONE
      if (isigtyp.eq.0) then
!----- NO Scattering
         sig=0.d0
      else if (isigtyp.eq.1) then
!----- Constant Thomson Scattering (=0.4)
         sig=0.4d0
      else
         print *,'isigtyp=', isigtyp,'is not defined'
         stop
      end if
      RETURN
      END SUBROUTINE SUBSIG

      SUBROUTINE OPCPAF(TM1,RHO1,BKAP)
!---   Routine to interpolate in Pollack et al (1985) and
!---     Alexander & Ferguson (1994) Tables
!---  NOTE:  Minimum Density     = 1.e-12
!---         Maximum Density     = 1.e+ 1
!---         Minimum Temperature = 1.e+ 2
!---         Maximum Temperature = 1.e+ 5
!-----------------
      USE opc_init
      USE fluid_var_init
      IMPLICIT NONE
      integer :: K,IT,IT1,IW,IW1
      double precision :: TM1,RHO1,BKAP,TLOG,RHOL
      double precision :: WS1,WS,WS2,WS3,Z00,Z10,Z11,Z01,WS4,WS5
      IF(TM1.LE. 0.) then
         print *,' NEG. T. IN AF OPC=',TM1
         STOP
      endif
      TLOG = LOG10(TM1)
      RHOL = LOG10(RHO1)
      IF(TLOG .LT. 2.d0)   TLOG = 2.0000001
      IF(TLOG .GT. 5.d0)   TLOG = 4.9999999
      IF(RHOL .LT. -12.d0) RHOL = -11.999999
      IF(RHOL .GT. 1.d0)   RHOL = 0.9999999
      DO K = 1,60
         IF(TVAL(K) .GT. TLOG) GO TO 15
      ENDDO
      IF(K.EQ.61) THEN
         PRINT *,'Overflow in opacity temperature array:OPCPAF'
         print *,'iteration',num_iter,III
         PRINT *,'log(T)= ',TLOG
         CALL CLEAN_STOP
      ENDIF
15    IT = K
      IT1 = K-1
      IF(IT1.eq.0) THEN
         PRINT *,'underflow in opacity temperature array:OPCPAF'
         PRINT *,'log(T)= ',TLOG
         CALL CLEAN_STOP
      ENDIF
      WS1 = TVAL(IT1) - TLOG
      WS = TVAL(IT1) - TVAL(IT)
      WS1 = WS1 / WS
      WS = 1.d0 - WS1
      DO K = 1,14
         IF(RVAL(K) .GT. RHOL) GO TO 25
      ENDDO
      IF(K.EQ.15) THEN
         PRINT *,'Overflow in opacity density array:OPCPAF'
         print *,'iteration',num_iter,III
         PRINT *,'log(rho)= ',rhol
         CALL CLEAN_STOP
      ENDIF
25    IW = K
      IW1 = K-1
      IF(IW1.EQ.0) THEN
         PRINT *,'underflow in opacity density array:OPCPAF'
         PRINT *,'log(rho)= ',rhol
         STOP
      ENDIF
      WS2 = RHOL - RVAL(IW1)
      WS3 = RVAL(IW) - RVAL(IW1)
      WS2 = WS2 / WS3
      WS3 = 1. - WS2
      Z00 = ZKAP(IW,IT)
      Z10 = ZKAP(IW1,IT)
      Z11 = ZKAP(IW1,IT1)
      Z01 = ZKAP(IW,IT1)
c      IF(Z00.EQ.0..OR.Z01.EQ.0..OR.Z11.EQ.0..OR.Z10.EQ.0.) then
c         print *,'negative ZKAP',Z00,Z01,Z10,Z11
c         stop
c      ENDIF
      WS4 = WS2*(WS1*Z00 + WS*Z01)
      WS5 = WS3*(WS1*Z10 + WS*Z11)
      BKAP = 10.d0**(WS4+WS5)
      IF(BKAP.LE.0.d0) THEN
         print *,'negative/zero BKAP',BKAP,TM1,RHO1
         stop
      ENDIF
      RETURN
      END SUBROUTINE OPCPAF

      SUBROUTINE READALEXANDEROPCTABLE
      USE opc_init
      USE input_init
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      double precision :: min,max,out
      integer :: i,j
      integer :: iT,iRH
      if(F_INCIDENT.gt.0) then
         print *,'YOU CANNOT USE ALEXANDER OPC AND STELLAR HEATING'
      ENDIF
      OPEN(UNIT=9,FILE='opctbl',STATUS='OLD')
      if(myid.eq.0) print *, 'Reading Alexander Opacity Table'
 46   FORMAT(7F7.3)
      DO iT=1,60
         READ(9,46) (ZKAP(iRH,iT),iRH=1,7)
         READ(9,46) (ZKAP(iRH,iT),iRH=8,14)
      ENDDO
      CLOSE(9)
!--- MULTIPLY THE OPACITY TABLE BY ZKAPMULT - BUT TABLE IS LOG..
      ZKAP = ZKAP+LOG10(ZKAPMULT)
      min = 100.d0
      max = -100.d0
      do i=1,14
         do j=1,60
            if(zkap(i,j).lt.min) min = zkap(i,j)
            if(zkap(i,j).gt.max) max = zkap(i,j)
         enddo
      enddo
c      if(myid.eq.0) print *,'Minimum Opacity = ',10.d0**min
c      if(myid.eq.0) print *,'Maximum Opacity = ',10.d0**max
c      call OPCPAF(Tstar,10.0**(-11.0),out)
c      if(myid.eq.0) print *,'Outside Opacity = ',10.d0**out
      RETURN
      END SUBROUTINE READALEXANDEROPCTABLE

      SUBROUTINE OPLIN0(temp,rho,opacty)
**  Doug's Opacity Function for Ice, grains, .... (orig)
      USE input_init
      USE oplin_constants
      IMPLICIT NONE
      double precision :: temp,rho,opacty
      double precision :: ts4,ts42,ts44,ts48
      double precision :: rho13,rho23
      double precision :: o1,o2,o3,o4,o5,o6,o7,o8
      double precision :: o1an,o2an,o3an,o4an,o6an,o7an,o8an
      double precision :: t2,t4,t8,t10
!--------------
!       OPACITIES UP TO 1-1,000,000 K
!	The opacity table is divided into eight regions
!	region 1 is dominated by ice grains 
!	region 2 is where ice grains melts
!	region 3 is dominated by metal grains
!	region 4 is where metal grains melts
!	region 5 is dominated by molecules
!	region 6 is dominated by H-
!	region 7 is dominated by Kramer's law
!	region 8 is dominated by electron scattering  (taken out)
!
!	the opacity for these eight regions are obtained from fitting
!       Bodemheimer's table to power laws of T and rho
!	A smooth function is introduced to bring about the transition 
!	between each region
!--------------
      if(F_INCIDENT.gt.0) then
         print *,'YOU CANNOT USE OPLIN AND STELLAR HEATING'
      ENDIF
c   test T against (T_23 * T_34 * T_34)**0.333333333
      if(temp.gt.t234*rho**power1)then
c   to avoid overflow
         ts4=1.d-4*temp
         rho13=rho**0.333333333
         rho23=rho13*rho13
         ts42=ts4*ts4
         ts44=ts42*ts42
         ts48=ts44*ts44
c   test T against (T_45 * T_56)**0.5
         if (temp.gt.t456*rho**power2) goto 10
c   disjoint opacity laws for 3, 4, and 5.
         o3=bk3*ts4
         o4=bk4*rho23/(ts48*ts4)
         o5=bk5*rho23*ts42*ts4
c   parameters used for smoothing
         o4an=o4**4
         o3an=o3**4
c   smoothed and continuous opacity law for regions 3, 4, and 5.
         opacty=((o4an*o3an/(o4an+o3an))+(o5/
     &        (1+6.561e-5/ts48))**4)**0.25
      else
c   different powers of temperature
         t2=temp*temp
         t4=t2*t2
         t8=t4*t4
         t10=t8*t2
c   disjoint opacity laws
         o1=ak1*t2
         o2=ak2*temp/t8
         o3=ak3*temp
c   parameters used for smoothing
         o1an=o1*o1
         o2an=o2*o2
c   smoothed and continuous opacity law for regions 1, 2, and 3.
         opacty=((o1an*o2an/(o1an+o2an))**2+(o3/
     &        (1+1.e22/t10))**4)**0.25
      endif
      RETURN
 10   if(temp.lt.t678*rho**power3)then
c   disjoint opacity laws for 5, 6, and 7.
         o5=bk5*rho23*ts42*ts4
         o6=bk6*rho13*ts44*ts44*ts42
         o7=bk7*rho/(ts42*sqrt(ts4))
c   parameters used for smoothing
         o6an=o6*o6
         o7an=o7*o7
c   smoothed and continuous opacity law for regions 5, 6, and 7.
         opacty=((o6an*o7an/(o6an+o7an))**2+(o5/
     &        (1+(ts4/(1.1*rho**0.04762))**10))**4)**0.25
      else
c   disjoint opacity laws for 7 and 8.
         o7=bk7*rho/(ts42*sqrt(ts4))
         o8=bk8	
**  No Scattering here (in Sigma)
         o8=0.0	
c   parameters used for smoothing
         o7an=o7*o7
         o8an=o8*o8
c   smoothed and continuous opacity law for regions 7 and 8.
         opacty=(o7an*o7an+o8an*o8an)**0.25
      endif
      RETURN
      END SUBROUTINE OPLIN0

      SUBROUTINE FreedmanOPC(T_IN,PG_IN,BKAP_P,BKAP_R,BKAP_A,
     $     BDXKAPP_DT,BDXKAPA_DT)
!---   Routine to interpolate in Freedman March 19,2008 Opacites
!---     Similar to Freedman, Marley, & Lodders (2007) opacities Tables
!---      but with more P/T points
!---  NOTE:  Minimum Pressure     = 0.001 mbars
!---         Maximum Pressure     = 300000 mbars
!---         Minimum Temperature  = 75.0 K
!---         Maximum Temperature  = 4000.0 K
!-----------------
      USE opc_init
      USE mpi_var_init
      USE fluid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: K,IT,IT1,IW,IW1
      double precision :: T_IN,PG_IN
      double precision :: TM1,PG1,BKAP_P,BKAP_R,BKAP_A
      double precision :: BDXKAPP_DT,BDXKAPA_DT
      double precision :: WS1,WS,WS2,WS3
      double precision :: R00,R10,R11,R01,R_WS4,R_WS5
      double precision :: P00,P10,P11,P01,P_WS4,P_WS5
      double precision :: A00,A10,A11,A01,A_WS4,A_WS5
      character proc_num*3
      character filename*12
      IF(T_IN.LE. 0.) then
         print *,' NEG. T. IN FREEDMAN OPC=',T_IN,myid
         STOP
      endif
!-- PUT PRESSURE IN MBARS FOR TABLE COMPARISON (PG1 ORIGINALLY IN DYN/CM^2)
      PG1 = PG_IN*(0.001d0)
      TM1 = T_IN
      IF(TM1 .LT. 75.d0)   TM1 = 75.001d0
      IF(TM1 .GT. 4000.d0) TM1 = 3999.9d0
      IF(PG1 .LT. 0.001d0) PG1 = 0.001001d0
      IF(PG1 .GT. 300000.d0)   PG1 = 299999.9d0
      DO K = 1,42
         IF(FreedmanTVAL(K) .GT. TM1) GO TO 15
      ENDDO
      IF(K.EQ.43) THEN
         PRINT *,'Overflow in opacity temperature array'
         print *,'iteration',num_iter,III
         PRINT *,'Freedman opacities',TM1,myid
         CALL CLEAN_STOP
      ENDIF
15    IT = K
      IT1 = K-1
      IF(IT1.eq.0) THEN
         PRINT *,'underflow in opacity temperature array'
         PRINT *,'Freedman opacities',TM1
         CALL CLEAN_STOP
      ENDIF
      WS1 = FreedmanTVAL(IT1) - TM1
      WS =  FreedmanTVAL(IT1) - FreedmanTVAL(IT)
      WS1 = WS1 / WS
      WS = 1.d0 - WS1
      DO K = 1,18
         IF(FreedmanPVAL(K) .GT. PG1) GO TO 25
      ENDDO
      IF(K.EQ.19) THEN
         PRINT *,'Overflow in opacity pressure array',PG1
         print *,'iteration',num_iter,III
         CALL CLEAN_STOP
      ENDIF
25    IW = K
      IW1 = K-1
      IF(IW1.EQ.0) THEN
         PRINT *,'underflow in opacity pressure array'
         STOP
      ENDIF
      WS2 = PG1 - FreedmanPVAL(IW1)
      WS3 = FreedmanPVAL(IW) - FreedmanPVAL(IW1)
      WS2 = WS2 / WS3
      WS3 = 1. - WS2

!-- ROSSLAND MEAN OPACITIES
      R00 = FreedmanXKAP_ross(IT,IW)
      R10 = FreedmanXKAP_ross(IT,IW1)
      R11 = FreedmanXKAP_ross(IT1,IW1)
      R01 = FreedmanXKAP_ross(IT1,IW)      
      R_WS4 = WS2*(WS1*R00 + WS*R01)
      R_WS5 = WS3*(WS1*R10 + WS*R11)
      BKAP_R = 10.d0**(R_WS4+R_WS5)

!-- PLANCK MEAN OPACITIES
      P00 = FreedmanXKAP_plnk(IT,IW)
      P10 = FreedmanXKAP_plnk(IT,IW1)
      P11 = FreedmanXKAP_plnk(IT1,IW1)
      P01 = FreedmanXKAP_plnk(IT1,IW)
      P_WS4 = WS2*(WS1*P00 + WS*P01)
      P_WS5 = WS3*(WS1*P10 + WS*P11)
      BKAP_P = 10.d0**(P_WS4+P_WS5)

!-- ABSORPTION OPACITIES
      A00 = FreedmanXKAP_abs(IT,IW)
      A10 = FreedmanXKAP_abs(IT,IW1)
      A11 = FreedmanXKAP_abs(IT1,IW1)
      A01 = FreedmanXKAP_abs(IT1,IW)
      A_WS4 = WS2*(WS1*A00 + WS*A01)
      A_WS5 = WS3*(WS1*A10 + WS*A11)
      BKAP_A = 10.d0**(A_WS4+A_WS5)


!-- CLUDGE FOR NOW
      BDXKAPP_DT = ((10.d0**(P00))-(10.d0**(P01)))/
     $     (FreedmanTVAL(IT)-FreedmanTVAL(IT1))

      BDXKAPA_DT = ((10.d0**(A00))-(10.d0**(A01)))/
     $     (FreedmanTVAL(IT)-FreedmanTVAL(IT1))
      


      if(R00.eq.-7.d0) print *,'R00=-7.d0',TM1,PG1
      if(R10.eq.-7.d0) print *,'R10=-7.d0',TM1,PG1
      if(R11.eq.-7.d0) print *,'R11=-7.d0',TM1,PG1
      if(R01.eq.-7.d0) print *,'R01=-7.d0',TM1,PG1

      IF(BKAP_R.LE.0.d0.OR.BKAP_P.LE.0.d0.OR.BKAP_A.LE.0.d0) THEN
         print *,'neg/zero BKAP_R,BKAP_P,BKAP_A',BKAP_R,BKAP_P,BKAP_A
         PRINT *,'TM1,PG1=',TM1,PG1
         print *,'tval=',FreedmanTVAL(IT),FreedmanTVAL(IT1)
         print *,'pval=',FreedmanPVAL(IW),FreedmanPVAL(IW1)
         print *,'init Ws',WS1,WS,WS2,WS3
         print *,'ZXX1=',R00,R10,R11,R01
         print *,'ZXX2=',P00,P10,P11,P01
         print *,'ZXX3=',A00,A10,A11,A01
         print *,'final Ws',R_WS4,R_WS5
         print *,'final Ws',P_WS4,P_WS5
         print *,'final Ws',A_WS4,A_WS5
         print *, IT,IT1,IW,IW1

         write(proc_num,'(i3.3)') myid
         filename = 'SNAP_'//TRIM(ADJUSTL(proc_num))//'.'//
     $        '100'
         OPEN (12,FILE=filename,STATUS='UNKNOWN',
     $        FORM='UNFORMATTED')
         REWIND(12)
         CALL PRINTOUT(1)
         CLOSE(12)
         stop
      ENDIF
      RETURN
      END SUBROUTINE FreedmanOPC

      SUBROUTINE READFREEDMANOPCTABLE
      USE opc_init
      USE mpi_var_init
      USE input_init
      IMPLICIT NONE
      include "mpif.h"
      character text*20
      integer :: i,j
      integer :: Tindex,Pindex,N_entries
      double precision :: readT,readP,readkapR,readkapP
      double precision :: readkapABS,Tabsorption_read
!--- ikaptyp = 4 : Freedman(2007) opacities (updated 4/2008) normal table
!--- ikaptyp = 5 : Freedman(2007) opacities (updated 4/2008) without TIVO
!--- ikaptyp = 6 : Freedman(2007) opacities (updated 4/2008) normal table, column opacities
!--- ikaptyp = 7 : Freedman(2007) opacities (updated 4/2008) without TIVO, column opacities
      if(ikaptyp.eq.4.or.ikaptyp.eq.6) then
         open(12,file="Freedman6117.dat",status='old')
         if(myid.eq.0) print *,'reading Freedman6117.dat'
      elseif(ikaptyp.eq.5.or.ikaptyp.eq.7) then
         open(12,file="Freedman6117_noTIVO.dat",status='old')
         if(myid.eq.0) print *,'reading Freedman6117_noTiVO.dat'
      else
         print *,'ikaptyp=',ikaptyp,'not avail in readfreedmanopctable' 
         stop
      endif
      N_entries=736
!-- Initially set everything to -7 to account for missing values
      FreedmanXKAP_plnk = -7.d0
      FreedmanXKAP_ross = -7.d0
      FreedmanXKAP_abs = -7.d0
!-- Skip header
      DO i=1,10
         READ(12,*) text
      ENDDO
      READ(12,*) Tabsorption_read
      if(myid.eq.0) then
         print *,'Absorption opacities calculated with T_star=',
     %        Tabsorption_read
      endif
!-- read in actual values
      DO i=1,N_entries
         READ(12,*) readT,readP,readkapP,readkapR,readkapABS
         do j=1,42
            if(FreedmanTVAL(j).eq.readT) goto 100
         enddo
 100     CONTINUE
         Tindex = j
         do j=1,18
            if(FreedmanPVAL(j).eq.readP) goto 101
         enddo
 101     CONTINUE
         Pindex = j
         FreedmanXKAP_plnk(Tindex,Pindex) = readkapP
         FreedmanXKAP_ross(Tindex,Pindex) = readkapR
         FreedmanXKAP_abs(Tindex,Pindex)  = readkapABS
      enddo
      CLOSE(12)
      RETURN
      END SUBROUTINE READFREEDMANOPCTABLE


      SUBROUTINE BurrowsOPC(T_IN,RH_IN,BKAP_P,BKAP_R,BKAP_A,
     $     BKAP_W,BKAP_LW,BKAP_LBB,BDXKAPP_DT,BDXKAPA_DT,
     $     i_indx,j_indx,k_indx)
!---   Routine to interpolate in Burrows opacities tables
!---  NOTE:  Minimum density     = 1.00102e-12 g/cc
!---         Maximum density     = 0.0100017 g/cc
!---         Minimum Temperature  = 49.9988
!---         Maximum Temperature  = 4999.04 K
!-----------------
      USE input_init
      USE opc_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: K,IT,IT1,IW,IW1
      double precision :: T_IN,RH_IN
      double precision :: TM1,RH1,BKAP_P,BKAP_R,BKAP_A
      double precision :: BDXKAPP_DT,BDXKAPA_DT
      double precision :: WS1,WS,WS2,WS3
      double precision :: R00,R10,R11,R01,R_WS4,R_WS5
      double precision :: P00,P10,P11,P01,P_WS4,P_WS5
      double precision :: A00,A10,A11,A01,A_WS4,A_WS5
      double precision :: W00,W10,W11,W01,W_WS4,W_WS5
      double precision :: LW00,LW10,LW11,LW01,LW_WS4,LW_WS5
      double precision :: LBB00,LBB10,LBB11,LBB01,LBB_WS4,LBB_WS5 
      double precision :: BKAP_W(nwave_bins),BKAP_LW(nwave_bins)
      double precision :: BKAP_LBB(nwave_bins)

      integer :: i_indx,j_indx,k_indx,l
      character proc_num*3
      character filename*12
      IF(T_IN.LE. 0.) then
         print *,' NEG. T. IN BURROWS OPC=',T_IN,myid
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         print *,'iteration',III
         call clean_STOP
      endif
      RH1 = LOG10(RH_IN)
      TM1 = LOG10(T_IN)
!-If values fall outside table, reset to inside table
      IF(TM1 .LT. 1.69896)   TM1 = 1.69897
      IF(TM1 .GT. 3.69889) TM1 = 3.69888
      IF(RH1 .LT. -11.9996) RH1 = -11.9995
      IF(RH1 .GT. -1.99993) RH1 = -1.99994
      DO K = 1,50
         IF(BurrowsT(K) .GT. TM1) GO TO 15
      ENDDO
      IF(K.EQ.51) THEN
         PRINT *,'Overflow in opacity temperature array'
         print *,'iteration',num_iter,III
         PRINT *,'Burrows opacities: T/myid=',TM1,myid
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd(j,k)',upperbnd(j_indx,k_indx)
         print *,'upperbnd(j-1,k)',upperbnd(j_indx-1,k_indx)
         print *,'upperbnd(j+1,k)',upperbnd(j_indx+1,k_indx)
         print *,'upperbnd(j,k-1)',upperbnd(j_indx,k_indx-1)
         print *,'upperbnd(j,k+1)',upperbnd(j_indx,k_indx+1)
         print *,'v',v(i_indx,j_indx,k_indx)
         print *,'g',g(i_indx,j_indx,k_indx)
         print *,'h',h(i_indx,j_indx,k_indx)
         print *,'T',T(i_indx,j_indx,k_indx)
         print *,'rh',rh(i_indx,j_indx,k_indx)
         print *,'energy_iter',energy_iter
         CALL CLEAN_STOP
      ENDIF
15    IT = K
      IT1 = K-1
      IF(IT1.eq.0) THEN
         PRINT *,'underflow in opacity temperature array'
         PRINT *,'Burrows opacities',TM1
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'iteration',num_iter,III
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         CALL CLEAN_STOP
      ENDIF
      WS1 = BurrowsT(IT1) - TM1
      WS =  BurrowsT(IT1) - BurrowsT(IT)
      WS1 = WS1 / WS
      WS = 1.d0 - WS1
      DO K = 1,50
         IF(BurrowsRH(K) .GT. RH1) GO TO 25
      ENDDO
      IF(K.EQ.51) THEN
         PRINT *,'Overflow in burrows opacity pressure array',RH1,myid
         print *,'iteration',num_iter,III
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         CALL CLEAN_STOP
      ENDIF
25    IW = K
      IW1 = K-1
      IF(IW1.EQ.0) THEN
         PRINT *,'underflow in opacity pressure array'
         print *,'iteration',num_iter,III
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         CALL CLEAN_STOP
      ENDIF
      WS2 = RH1 - BurrowsRH(IW1)
      WS3 = BurrowsRH(IW) - BurrowsRH(IW1)
      WS2 = WS2 / WS3
      WS3 = 1. - WS2

!-- ROSSLAND MEAN OPACITIES
      R00 = BurrowsXKAP_ross(IT,IW)
      R10 = BurrowsXKAP_ross(IT,IW1)
      R11 = BurrowsXKAP_ross(IT1,IW1)
      R01 = BurrowsXKAP_ross(IT1,IW)      
      R_WS4 = WS2*(WS1*R00 + WS*R01)
      R_WS5 = WS3*(WS1*R10 + WS*R11)
      BKAP_R = 10.d0**(R_WS4+R_WS5)

!-- PLANCK MEAN OPACITIES
      P00 = BurrowsXKAP_plnk(IT,IW)
      P10 = BurrowsXKAP_plnk(IT,IW1)
      P11 = BurrowsXKAP_plnk(IT1,IW1)
      P01 = BurrowsXKAP_plnk(IT1,IW)
      P_WS4 = WS2*(WS1*P00 + WS*P01)
      P_WS5 = WS3*(WS1*P10 + WS*P11)
      BKAP_P = 10.d0**(P_WS4+P_WS5)

!-- ABSORPTION OPACITIES
      A00 = BurrowsXKAP_abs(IT,IW)
      A10 = BurrowsXKAP_abs(IT,IW1)
      A11 = BurrowsXKAP_abs(IT1,IW1)
      A01 = BurrowsXKAP_abs(IT1,IW)
      A_WS4 = WS2*(WS1*A00 + WS*A01)
      A_WS5 = WS3*(WS1*A10 + WS*A11)
      BKAP_A = 10.d0**(A_WS4+A_WS5)

      if(ikaptyp.eq.20.or.ikaptyp.eq.21.or.
     $     ikaptyp.eq.22.or.ikaptyp.eq.23.or.
     $     ikaptyp.eq.24.or.ikaptyp.eq.25) then
         do l=1,nwave_bins
            W00 = BurrowsXKAP_waveopac(IT,IW,l)
            W10 = BurrowsXKAP_waveopac(IT,IW1,l)
            W11 = BurrowsXKAP_waveopac(IT1,IW1,l)
            W01 = BurrowsXKAP_waveopac(IT1,IW,l)
            W_WS4 = WS2*(WS1*W00 + WS*W01)
            W_WS5 = WS3*(WS1*W10 + WS*W11)
            BKAP_W(l) = 10.d0**(W_WS4+W_WS5)
         enddo
      else
         BKAP_W = -7
      endif

      if(ikaptyp.eq.22.or.ikaptyp.eq.23.or.
     $     ikaptyp.eq.24.or.ikaptyp.eq.25) then
         do l=1,nwave_bins
            LW00 = BurrowsXKAP_localwave(IT,IW,l)
            LW10 = BurrowsXKAP_localwave(IT,IW1,l)
            LW11 = BurrowsXKAP_localwave(IT1,IW1,l)
            LW01 = BurrowsXKAP_localwave(IT1,IW,l)
            LW_WS4 = WS2*(WS1*LW00 + WS*LW01)
            LW_WS5 = WS3*(WS1*LW10 + WS*LW11)
            BKAP_LW(l) = 10.d0**(LW_WS4+LW_WS5)
         enddo
      else
         BKAP_LW = -7
      endif


!-note I should fix this because the source is independent of 
!-  density... ie lbb00=lbb10 and lbb11=lb01
      if(ikaptyp.eq.22.or.ikaptyp.eq.23.or.
     $     ikaptyp.eq.24.or.ikaptyp.eq.25) then
         do l=1,nwave_bins
            LBB00 = localBB_wavelgth(IT,IW,l)
            LBB10 = localBB_wavelgth(IT,IW1,l)
            LBB11 = localBB_wavelgth(IT1,IW1,l)
            LBB01 = localBB_wavelgth(IT1,IW,l)
            LBB_WS4 = WS2*(WS1*LBB00 + WS*LBB01)
            LBB_WS5 = WS3*(WS1*LBB10 + WS*LBB11)
            BKAP_LBB(l) = (LBB_WS4+LBB_WS5)
         enddo
      else
         BKAP_LBB = -7
      endif



!-- CLUDGE FOR NOW
      BDXKAPP_DT = ((10.d0**(P00))-(10.d0**(P01)))/
     $     ((10.d0**BurrowsT(IT))-(10.d0**BurrowsT(IT1)))
      BDXKAPA_DT = ((10.d0**(A00))-(10.d0**(A01)))/
     $     ((10.d0**BurrowsT(IT))-(10.d0**BurrowsT(IT1)))



      if(R00.eq.-7.d0) print *,'R00=-7.d0',TM1,RH1
      if(R10.eq.-7.d0) print *,'R10=-7.d0',TM1,RH1
      if(R11.eq.-7.d0) print *,'R11=-7.d0',TM1,RH1
      if(R01.eq.-7.d0) print *,'R01=-7.d0',TM1,RH1

      IF(BKAP_R.LE.0.d0.OR.BKAP_P.LE.0.d0.OR.BKAP_A.LE.0.d0) THEN
         print *,'in burrows opc,,,'
         print *,'neg/zero BKAP_R,BKAP_P,BKAP_A',BKAP_R,BKAP_P,BKAP_A
         PRINT *,'TM1,RH1=',TM1,RH1
         print *,'tval=',BurrowsT(IT),BurrowsT(IT1)
         print *,'rhval=',BurrowsRH(IW),BurrowsRH(IW1)
         print *,'init Ws',WS1,WS,WS2,WS3
         print *,'ZXX=',R00,R10,R11,R01
         print *,'ZXX=',P00,P10,P11,P01
         print *,'ZXX=',A00,A10,A11,A01
         print *,'final Ws',R_WS4,R_WS5
         print *,'final Ws',P_WS4,P_WS5
         print *,'final Ws',A_WS4,A_WS5
         print *, IT,IT1,IW,IW1
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         write(proc_num,'(i3.3)') myid
         filename = 'SNAP_'//TRIM(ADJUSTL(proc_num))//'.'//
     $        '100'
         OPEN (12,FILE=filename,STATUS='UNKNOWN',
     $        FORM='UNFORMATTED')
         REWIND(12)
         CALL PRINTOUT(1)
         CLOSE(12)
         stop
      ENDIF
      RETURN
      END SUBROUTINE BurrowsOPC

      SUBROUTINE READBURROWSOPCTABLE
      USE opc_init
      USE mpi_var_init
      USE rad_var_init
      USE input_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      character text*20
      integer :: i,j,l
      integer :: Tindex,RHindex
      double precision :: readT,readRH,readkapR,readkapP
      double precision :: readkapABS,Tabsorption_read
      double precision :: PG_test,kap_ext,test_val,SETPLANETKAPEXTRA
      double precision :: read1,read2,read3,read4,read5,read6,read7
      double precision :: read8,read9,read10,read11,read12,read13,read14
      double precision :: read15,read16,read17,read18,read19,read20
      double precision :: read21,read22,read23,read24,read25,read26
      double precision :: read27,read28,read29,read30
      double precision :: read31,read32,read33,read34,read35,read36
      double precision :: read37,read38,read39,read40,read41,read42
      double precision :: read43,read44,read45,read46,read47,read48
      double precision :: read49,read50,read51,read52,read53,read54
      double precision :: read55,read56,read57,read58,read59,read60
      double precision :: read61,read62,read63,read64,read65,read66
      double precision :: read67,read68,read69,read70,read71,read72
      double precision :: read73,read74,read75,read76,read77,read78
      double precision :: read79,read80,read81,read82,read83,read84
      double precision :: read85,read86,read87,read88,read89,read90
      double precision :: bin_tstar
      integer :: nbin_test
!--- ikaptyp = 9 : Burrows opacities (updated 5/2008) with TIVO
!--- ikaptyp = 10 : Burrows opacities (updated 5/2008) without TIVO
!--- ikaptyp = 11 : Burrows opacities (updated 5/2008) without TIVO, with artificial opc
!--- ikaptyp = 12 : Burrows opacities (calc. 5/2008) without TIVO for xkapR, constant gamma
!--- ikaptyp = 14 : Burrows opacities (calc. 5/2008) without TIVO for xkapP, xkapA=xkapR=gamma^4*xkapP, constant gamma
!--- ikaptyp = 20 : Burrows opacities (updated 5/2008) with TIVO
!--- ikaptyp = 21 : Burrows opacities (updated 5/2008) without TIVO
!--- ikaptyp = 22 : Burrows opacities full radial radiative transfer solution(updated 6/2011) with TIVO
!--- ikaptyp = 23 : Burrows opacities, full radial radiative transfer solution(updated 6/2011) without TIVO
!--- ikaptyp = 24 : Burrows opacities, full rad. transfer sol. (updated 6/2011) without TIVO + extra abs in stellarinput_heating opacity
!--- ikaptyp = 25 : Burrows opacities, full rad. transfer sol. without TIVO + Rayleigh SCat + extra abs for all opac
      if(ikaptyp.eq.9.or.ikaptyp.eq.20) then
         call OPENOPACITYFILE_WITH_TIOVO
      elseif(ikaptyp.eq.10.or.ikaptyp.eq.11.or.ikaptyp.eq.12.or.
     %        ikaptyp.eq.14.or.ikaptyp.eq.21.or.ikaptyp.eq.23.or.
     $        ikaptyp.eq.24.or.ikaptyp.eq.25) then
         CALL OPENOPACITYFILE_WOTIOVO
      else
         print *,'ikaptyp=',ikaptyp,'not avail in readburrowsopctable' 
         call clean_stop
      endif
!-- Initially set everything to -7 to account for missing values
      BurrowsXKAP_plnk = -7.d0
      BurrowsXKAP_ross = -7.d0
      BurrowsXKAP_abs = -7.d0
      BurrowsXKAP_waveopac = -7.d0
!-- Skip header
      DO i=1,10
         READ(12,*) text
      ENDDO
      READ(12,*) Tabsorption_read
      if(myid.eq.0) then
         print *,'Absorption opacities calculated with T_star=',
     %        Tabsorption_read
      endif

!-- read in actual values
      DO Tindex=1,50
         DO RHindex=1,50
            if(ikaptyp.eq.20.or.ikaptyp.eq.21.or.
     $           ikaptyp.eq.22.or.ikaptyp.eq.23.or.ikaptyp.eq.24.or.
     $           ikaptyp.eq.25) then
               if(ikaptyp.eq.20.or.ikaptyp.eq.21) then
                  READ(12,*) readT,readRH,readkapP,readkapR,readkapABS,
     $                 read1,read2,read3,read4,read5,read6,read7,
     $                 read8,read9,read10,read11,read12,read13,read14,
     $                 read15,read16,read17,read18,read19,read20,
     $                 read21,read22,read23,read24,read25,read26,
     $                 read27,read28,read29,read30
               else
                  READ(12,*) readT,readRH,readkapP,readkapR,readkapABS,
     $                 read1,read2,read3,read4,read5,read6,read7,
     $                 read8,read9,read10,read11,read12,read13,read14,
     $                 read15,read16,read17,read18,read19,read20,
     $                 read21,read22,read23,read24,read25,read26,
     $                 read27,read28,read29,read30,read31,read32,read33,
     $                 read34,read35,read36,read37,read38,read39,read40,
     $                 read41,read42,read43,read44,read45,read46,read47,
     $                 read48,read49,read50,read51,read52,read53,read54,
     $                 read55,read56,read57,read58,read59,read60,read61,
     $                 read62,read63,read64,read65,read66,read67,read68,
     $                 read69,read70,read71,read72,read73,read74,read75,
     $                 read76,read77,read78,read79,read80,read81,read82,
     $                 read83,read84,read85,read86,read87,read88,read89,
     $                 read90
               endif             
               BurrowsXKAP_plnk(Tindex,RHindex) = readkapP
               BurrowsXKAP_ross(Tindex,RHindex) = readkapR
               BurrowsXKAP_abs(Tindex,RHindex)  = readkapABS
               BurrowsXKAP_waveopac(Tindex,RHindex,1)  = read1
               BurrowsXKAP_waveopac(Tindex,RHindex,2)  = read2
               BurrowsXKAP_waveopac(Tindex,RHindex,3)  = read3
               BurrowsXKAP_waveopac(Tindex,RHindex,4)  = read4
               BurrowsXKAP_waveopac(Tindex,RHindex,5)  = read5
               BurrowsXKAP_waveopac(Tindex,RHindex,6)  = read6
               BurrowsXKAP_waveopac(Tindex,RHindex,7)  = read7
               BurrowsXKAP_waveopac(Tindex,RHindex,8)  = read8
               BurrowsXKAP_waveopac(Tindex,RHindex,9)  = read9
               BurrowsXKAP_waveopac(Tindex,RHindex,10)  = read10
               BurrowsXKAP_waveopac(Tindex,RHindex,11)  = read11
               BurrowsXKAP_waveopac(Tindex,RHindex,12)  = read12
               BurrowsXKAP_waveopac(Tindex,RHindex,13)  = read13
               BurrowsXKAP_waveopac(Tindex,RHindex,14)  = read14
               BurrowsXKAP_waveopac(Tindex,RHindex,15)  = read15
               BurrowsXKAP_waveopac(Tindex,RHindex,16)  = read16
               BurrowsXKAP_waveopac(Tindex,RHindex,17)  = read17
               BurrowsXKAP_waveopac(Tindex,RHindex,18)  = read18
               BurrowsXKAP_waveopac(Tindex,RHindex,19)  = read19
               BurrowsXKAP_waveopac(Tindex,RHindex,20)  = read20
               BurrowsXKAP_waveopac(Tindex,RHindex,21)  = read21
               BurrowsXKAP_waveopac(Tindex,RHindex,22)  = read22
               BurrowsXKAP_waveopac(Tindex,RHindex,23)  = read23
               BurrowsXKAP_waveopac(Tindex,RHindex,24)  = read24
               BurrowsXKAP_waveopac(Tindex,RHindex,25)  = read25
               BurrowsXKAP_waveopac(Tindex,RHindex,26)  = read26
               BurrowsXKAP_waveopac(Tindex,RHindex,27)  = read27
               BurrowsXKAP_waveopac(Tindex,RHindex,28)  = read28
               BurrowsXKAP_waveopac(Tindex,RHindex,29)  = read29
               BurrowsXKAP_waveopac(Tindex,RHindex,30)  = read30
!- NOTE: The bins are in order of decreasing wavelength (ie the first
!- wavelength bin is the longest wavelength)
               if(ikaptyp.eq.22.or.ikaptyp.eq.23.or.
     $              ikaptyp.eq.24.or.ikaptyp.eq.25) then
                  BurrowsXKAP_localwave(Tindex,RHindex,1)  = read31
                  BurrowsXKAP_localwave(Tindex,RHindex,2)  = read32
                  BurrowsXKAP_localwave(Tindex,RHindex,3)  = read33
                  BurrowsXKAP_localwave(Tindex,RHindex,4)  = read34
                  BurrowsXKAP_localwave(Tindex,RHindex,5)  = read35
                  BurrowsXKAP_localwave(Tindex,RHindex,6)  = read36
                  BurrowsXKAP_localwave(Tindex,RHindex,7)  = read37
                  BurrowsXKAP_localwave(Tindex,RHindex,8)  = read38
                  BurrowsXKAP_localwave(Tindex,RHindex,9)  = read39
                  BurrowsXKAP_localwave(Tindex,RHindex,10)  = read40
                  BurrowsXKAP_localwave(Tindex,RHindex,11)  = read41
                  BurrowsXKAP_localwave(Tindex,RHindex,12)  = read42
                  BurrowsXKAP_localwave(Tindex,RHindex,13)  = read43
                  BurrowsXKAP_localwave(Tindex,RHindex,14)  = read44
                  BurrowsXKAP_localwave(Tindex,RHindex,15)  = read45
                  BurrowsXKAP_localwave(Tindex,RHindex,16)  = read46
                  BurrowsXKAP_localwave(Tindex,RHindex,17)  = read47
                  BurrowsXKAP_localwave(Tindex,RHindex,18)  = read48
                  BurrowsXKAP_localwave(Tindex,RHindex,19)  = read49
                  BurrowsXKAP_localwave(Tindex,RHindex,20)  = read50
                  BurrowsXKAP_localwave(Tindex,RHindex,21)  = read51
                  BurrowsXKAP_localwave(Tindex,RHindex,22)  = read52
                  BurrowsXKAP_localwave(Tindex,RHindex,23)  = read53
                  BurrowsXKAP_localwave(Tindex,RHindex,24)  = read54
                  BurrowsXKAP_localwave(Tindex,RHindex,25)  = read55
                  BurrowsXKAP_localwave(Tindex,RHindex,26)  = read56
                  BurrowsXKAP_localwave(Tindex,RHindex,27)  = read57
                  BurrowsXKAP_localwave(Tindex,RHindex,28)  = read58
                  BurrowsXKAP_localwave(Tindex,RHindex,29)  = read59
                  BurrowsXKAP_localwave(Tindex,RHindex,30)  = read60
                  localBB_wavelgth(Tindex,RHindex,1)  = read61
                  localBB_wavelgth(Tindex,RHindex,2)  = read62
                  localBB_wavelgth(Tindex,RHindex,3)  = read63
                  localBB_wavelgth(Tindex,RHindex,4)  = read64
                  localBB_wavelgth(Tindex,RHindex,5)  = read65
                  localBB_wavelgth(Tindex,RHindex,6)  = read66
                  localBB_wavelgth(Tindex,RHindex,7)  = read67
                  localBB_wavelgth(Tindex,RHindex,8)  = read68
                  localBB_wavelgth(Tindex,RHindex,9)  = read69
                  localBB_wavelgth(Tindex,RHindex,10)  = read70
                  localBB_wavelgth(Tindex,RHindex,11)  = read71
                  localBB_wavelgth(Tindex,RHindex,12)  = read72
                  localBB_wavelgth(Tindex,RHindex,13)  = read73
                  localBB_wavelgth(Tindex,RHindex,14)  = read74
                  localBB_wavelgth(Tindex,RHindex,15)  = read75
                  localBB_wavelgth(Tindex,RHindex,16)  = read76
                  localBB_wavelgth(Tindex,RHindex,17)  = read77
                  localBB_wavelgth(Tindex,RHindex,18)  = read78
                  localBB_wavelgth(Tindex,RHindex,19)  = read79
                  localBB_wavelgth(Tindex,RHindex,20)  = read80
                  localBB_wavelgth(Tindex,RHindex,21)  = read81
                  localBB_wavelgth(Tindex,RHindex,22)  = read82
                  localBB_wavelgth(Tindex,RHindex,23)  = read83
                  localBB_wavelgth(Tindex,RHindex,24)  = read84
                  localBB_wavelgth(Tindex,RHindex,25)  = read85
                  localBB_wavelgth(Tindex,RHindex,26)  = read86
                  localBB_wavelgth(Tindex,RHindex,27)  = read87
                  localBB_wavelgth(Tindex,RHindex,28)  = read88
                  localBB_wavelgth(Tindex,RHindex,29)  = read89
                  localBB_wavelgth(Tindex,RHindex,30)  = read90
               endif
            else
               READ(12,*) readT,readRH,readkapP,readkapR,readkapABS
               BurrowsXKAP_plnk(Tindex,RHindex) = readkapP
               BurrowsXKAP_ross(Tindex,RHindex) = readkapR
               BurrowsXKAP_abs(Tindex,RHindex)  = readkapABS
            endif
            
            if(ikaptyp.eq.12.or.ikaptyp.eq.14) then !-- reset absorption opacity with constant gamma
               BurrowsXKAP_abs(Tindex,RHindex) = readkapP + 
     %              4.d0*log10(gamma_OPC)
            endif
            if(ikaptyp.eq.14) then  !-- reset rosseland opacity to be same as planck
               BurrowsXKAP_ross(Tindex,RHindex) = readkapP
            endif
            
!-- add in an extra absorber for ikaptyp=11, table is in log..
            if(ikaptyp.eq.11) then
               kap_ext = SETPLANETKAPEXTRA(planet_num)
               if(Tindex.eq.1.and.RHindex.eq.1.and.myid.eq.0) then
                  print *,'Adding extra opacity to xkapA'
                  print *,' kappa_extra =',kap_ext,'cm^2/g'
                  print *,' For P <0.03 bars'
               endif
               PG_test = (Rgas*readT*readRH/mu_gas)*(10.d0**(-6.d0))
               if(PG_test.le.0.03) then !-this was used in paper2
                  BurrowsXKAP_abs(Tindex,RHindex) = 
     %                 (10.d0**BurrowsXKAP_abs(Tindex,RHindex))+kap_ext
                  BurrowsXKAP_abs(Tindex,RHindex) = 
     %                 LOG10(BurrowsXKAP_abs(Tindex,RHindex))
               endif
            endif

!-- add in an extra absorber for ikaptyp=20 and 24 (table is in log..)
            if(ikaptyp.eq.20.or.ikaptyp.eq.24) then
               kap_ext = SETPLANETKAPEXTRA(planet_num)
               if(Tindex.eq.1.and.RHindex.eq.1.and.myid.eq.0) then
                  print *,'Adding extra opacity to xkapA'
                  print *,' kappa_extra =',kap_ext,'cm^2/g'
                  print *,' For All pressures'
               endif
               BurrowsXKAP_waveopac(Tindex,RHindex,:) = 
     $              (10.d0**BurrowsXKAP_waveopac(Tindex,RHindex,:))
     $              +kap_ext
               BurrowsXKAP_waveopac(Tindex,RHindex,:) = 
     $              LOG10(BurrowsXKAP_waveopac(Tindex,RHindex,:))
            endif

c-Extra opacity on top of the Rayleigh Scattering
            if(ikaptyp.eq.25) then
               kap_ext = SETPLANETKAPEXTRA(planet_num)
               if(Tindex.eq.1.and.RHindex.eq.1.and.myid.eq.0) then
                  print *,'Adding extra opacity to all opacities'
                  print *,' kappa_extra =',kap_ext,'cm^2/g'
                  print *,' For All pressures'
               endif
c-opacity for incident stellar irradiation
               BurrowsXKAP_waveopac(Tindex,RHindex,:) = 
     $              (10.d0**BurrowsXKAP_waveopac(Tindex,RHindex,:))
     $              +kap_ext
               BurrowsXKAP_waveopac(Tindex,RHindex,:) = 
     $              LOG10(BurrowsXKAP_waveopac(Tindex,RHindex,:))
c-opacity for the local, re-radiated two-stream radiative transfer
               BurrowsXKAP_localwave(Tindex,RHindex,:) = 
     $              (10.d0**BurrowsXKAP_localwave(Tindex,RHindex,:))
     $              +kap_ext
               BurrowsXKAP_localwave(Tindex,RHindex,:) = 
     $              LOG10(BurrowsXKAP_localwave(Tindex,RHindex,:))
            endif

         enddo
      enddo
      CLOSE(12)

      if(ikaptyp.eq.12.and.myid.eq.0) then
         print *,'setting kapA as constant ratio of kapP:'
         print *,' gamma_OPC = (kappaA/kappaP)^(0.25)=',gamma_OPC
      endif
      if(ikaptyp.eq.14.and.myid.eq.0) then
         print *,'setting kapA -AND- kapR as constant ratio of kapP:'
         print *,' gamma_OPC = (kappaA/kappaP)^(0.25)=',gamma_OPC
      endif

!-read in the integrated pi*bstar*dn bin values      
      if(ikaptyp.eq.20.or.ikaptyp.eq.21.or.
     $     ikaptyp.eq.22.or.ikaptyp.eq.23.or.ikaptyp.eq.24.or.
     $     ikaptyp.eq.25) then
         if(planet_num.eq.1) then !HD209458
            open(12,file="BSTAR_WAVEBINS_TS6070.dat",status='old')
            if(myid.eq.0) print *,'reading BSTAR_WAVEBINS_TS6070.dat'
         elseif(planet_num.eq.2) then !--HD189733b
            open(12,file="BSTAR_WAVEBINS_TS5040.dat",status='old')
            if(myid.eq.0) print *,'reading BSTAR_WAVEBINS_TS5040.dat'
         elseif(planet_num.eq.5) then !--WASP-12 (approx and large R)
            open(12,file="BSTAR_WAVEBINS_TS6300.dat",status='old')
            if(myid.eq.0) print *,'reading BSTAR_WAVEBINS_TS6300.dat'
         else
            print *,'Need to create a BSTAR_WAVEBIN file for this Tstar'
            call clean_stop
         endif
         DO i=1,3
            READ(12,*) text
         ENDDO
         READ(12,*) bin_tstar
         READ(12,*) nbin_test
         if(nbin_test.ne.nwave_bins) then
            print *,'Some problem with the number of bins in file'
            call clean_stop
         endif
         do i=1,nwave_bins
            READ(12,*) read1,read2
            BSTAR_BIN(i) = read2
         enddo
         close(12)
      endif
      RETURN
      END SUBROUTINE READBURROWSOPCTABLE

      SUBROUTINE columnOPC(T_IN,PG_IN,RH_IN,MU_IN,deltaR_IN,
     $     BKAP_P,BKAP_R,BKAP_A,BDXKAPP_DT,BDXKAPA_DT)
!---   Routine to interpolate in Freedman April 2008 Opacites WITH 
!---   the COLUMN MODIFICATION interpolation. Local Rosseland and Planck 
!---   opacites are calculated using the standard table
!---     Similar to Freedman, Marley, & Lodders (2007) opacities Tables
!---      but with more P/T points
!---  NOTE:  Minimum Pressure     = 0.001 mbars
!---         Maximum Pressure     = 300000 mbars
!---         Minimum Temperature  = 75.0 K
!---         Maximum Temperature  = 4000.0 K
!---         Minimum kappa cutoff = -26.0
!---         Maximum kappa cutoff = -19.0
!-----------------
      USE opc_init
      USE mpi_var_init
      USE global_constants
      USE fluid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: K,IT,IT1,IW,IW1,IC,IC1
      double precision :: T_IN,PG_IN,RH_IN,MU_IN,deltaR_IN
      double precision :: TM1,PG1,KAPCUTOFF,BKAP_P,BKAP_R,BKAP_A
      double precision :: BDXKAPP_DT,BDXKAPA_DT
      double precision :: WS1,WS,WS2,WS3
      double precision :: CS1,CS0
      double precision :: R00,R10,R11,R01,R_WS4,R_WS5
      double precision :: P00,P10,P11,P01,P_WS4,P_WS5
      double precision :: A000,A010,A011,A001,A_0
      double precision :: A100,A110,A111,A101,A_1
      double precision :: A_WS4,A_WS5
      character proc_num*3
      character filename*12
      IF(T_IN.LE. 0.) then
         print *,' NEG. T. IN FREEDMAN OPC=',T_IN,myid
         STOP
      endif
!-- PUT PRESSURE IN MBARS FOR TABLE COMPARISON (PG1 ORIGINALLY IN DYN/CM^2)
      PG1 = PG_IN*(0.001d0)
      TM1 = T_IN
!-- CALCULATE THE OPACITY CUTOFF (+1.0 is for x=upperbnd-> no column depth)
c      kapcutoff = 2.d0*MU_IN*MPROTON/(3.d0*RH_IN*(deltaR_IN+1.0))
c      kapcutoff = LOG10(kapcutoff)
      kapcutoff = -23.9526937 + LOG10(MU_IN) - LOG10(RH_IN) - 
     $     LOG10(deltaR_IN+1.0)
!-- COLUMN OPACITY LIMITS
      IF(kapcutoff .LT. -26.d0) kapcutoff=-25.9999d0
      IF(kapcutoff .GT. -19.d0) kapcutoff=-19.0001d0
!-- P/T LIMITS
      IF(TM1 .LT. 75.d0)   TM1 = 75.001d0
      IF(TM1 .GT. 4000.d0) TM1 = 3999.9d0
      IF(PG1 .LT. 0.001d0) PG1 = 0.001001d0
      IF(PG1 .GT. 300000.d0)   PG1 = 299999.9d0
!---- temperature interpolation
      DO K = 1,42
         IF(FreedmanTVAL(K) .GT. TM1) GO TO 15
      ENDDO
      IF(K.EQ.42) THEN
         PRINT *,'Overflow in opacity temperature array'
         print *,'iteration',num_iter,III
         PRINT *,'Freedman opacities',TM1,myid
         STOP
      ENDIF
15    IT = K
      IT1 = K-1
      IF(IT1.eq.0) THEN
         PRINT *,'underflow in opacity temperature array'
         PRINT *,'Freedman opacities',TM1
         print *,'iteration',num_iter,III
         STOP
      ENDIF
      WS1 = FreedmanTVAL(IT1) - TM1
      WS =  FreedmanTVAL(IT1) - FreedmanTVAL(IT)
      WS1 = WS1 / WS
      WS = 1.d0 - WS1
!---- pressure interpolation
      DO K = 1,18
         IF(FreedmanPVAL(K) .GT. PG1) GO TO 25
      ENDDO
      IF(K.EQ.19) THEN
         PRINT *,'Overflow in opacity pressure array',PG1
         print *,'iteration',num_iter,III
         STOP
      ENDIF
25    IW = K
      IW1 = K-1
      IF(IW1.EQ.0) THEN
         PRINT *,'underflow in opacity pressure array'
         STOP
      ENDIF
      WS2 = PG1 - FreedmanPVAL(IW1)
      WS3 = FreedmanPVAL(IW) - FreedmanPVAL(IW1)
      WS2 = WS2 / WS3
      WS3 = 1. - WS2

!---- column depth interpolation
      DO K = 1,8
         IF(FreedmanCVAL(K) .GT. kapcutoff) GO TO 35
      ENDDO
      IF(K.EQ.9) THEN
         PRINT *,'Overflow in opacity kappa array',kapcutoff
         print *,'iteration',num_iter,III
         STOP
      ENDIF
35    IC = K
      IC1 = K-1
      IF(IC1.EQ.0) THEN
         PRINT *,'underflow in opacity kappa array',kapcutoff
         print *,'iteration',num_iter,III
         STOP
      ENDIF
      CS1 = kapcutoff - FreedmanCVAL(IC1)
      CS0 = FreedmanCVAL(IC) - FreedmanPVAL(IC1)
      CS1 = CS1 / CS0
      CS0 = 1. - CS1


!-- ROSSLAND MEAN OPACITIES
      R00 = FreedmanXKAP_ross(IT,IW)
      R10 = FreedmanXKAP_ross(IT,IW1)
      R11 = FreedmanXKAP_ross(IT1,IW1)
      R01 = FreedmanXKAP_ross(IT1,IW)      
      R_WS4 = WS2*(WS1*R00 + WS*R01)
      R_WS5 = WS3*(WS1*R10 + WS*R11)
      BKAP_R = 10.d0**(R_WS4+R_WS5)

!-- PLANCK MEAN OPACITIES
      P00 = FreedmanXKAP_plnk(IT,IW)
      P10 = FreedmanXKAP_plnk(IT,IW1)
      P11 = FreedmanXKAP_plnk(IT1,IW1)
      P01 = FreedmanXKAP_plnk(IT1,IW)
      P_WS4 = WS2*(WS1*P00 + WS*P01)
      P_WS5 = WS3*(WS1*P10 + WS*P11)
      BKAP_P = 10.d0**(P_WS4+P_WS5)

!-- ABSORPTION OPACITIES
      A000 = ColumnXKAP_abs(IT,IW,IC)
      A010 = ColumnXKAP_abs(IT,IW1,IC)
      A011 = ColumnXKAP_abs(IT1,IW1,IC)
      A001 = ColumnXKAP_abs(IT1,IW,IC)
      A_WS4 = WS2*(WS1*A000 + WS*A001)
      A_WS5 = WS3*(WS1*A010 + WS*A011)
      A_0   = CS1*(A_WS4+A_WS5)
      A100 = ColumnXKAP_abs(IT,IW,IC1)
      A110 = ColumnXKAP_abs(IT,IW1,IC1)
      A111 = ColumnXKAP_abs(IT1,IW1,IC1)
      A101 = ColumnXKAP_abs(IT1,IW,IC1)
      A_WS4 = WS2*(WS1*A100 + WS*A101)
      A_WS5 = WS3*(WS1*A110 + WS*A111)
      A_1   = CS0*(A_WS4+A_WS5)
      BKAP_A = 10.d0**(A_0+A_1)

!-- CLUDGE FOR NOW
      BDXKAPP_DT = ((10.d0**(P00))-(10.d0**(P01)))/
     $     (FreedmanTVAL(IT)-FreedmanTVAL(IT1))

      BDXKAPA_DT = ((10.d0**(A000))-(10.d0**(A001)))/
     $     (FreedmanTVAL(IT)-FreedmanTVAL(IT1))
      


      if(R00.eq.-7.d0) print *,'R00=-7.d0',TM1,PG1
      if(R10.eq.-7.d0) print *,'R10=-7.d0',TM1,PG1
      if(R11.eq.-7.d0) print *,'R11=-7.d0',TM1,PG1
      if(R01.eq.-7.d0) print *,'R01=-7.d0',TM1,PG1

      IF(BKAP_R.LE.0.d0.OR.BKAP_P.LE.0.d0.OR.BKAP_A.LE.0.d0) THEN
         print *,'neg/zero BKAP_R,BKAP_P,BKAP_A',BKAP_R,BKAP_P,BKAP_A
         PRINT *,'TM1,PG1=',TM1,PG1
         print *,'tval=',FreedmanTVAL(IT),FreedmanTVAL(IT1)
         print *,'pval=',FreedmanPVAL(IW),FreedmanPVAL(IW1)
         print *,'init Ws',WS1,WS,WS2,WS3
         print *,'ZXX=',R00,R10,R11,R01
         print *,'ZXX=',P00,P10,P11,P01
         print *,'ZXX=',A000,A010,A011,A001
         print *,'ZXX=',A100,A110,A111,A101
         print *,'final Ws',R_WS4,R_WS5
         print *,'final Ws',P_WS4,P_WS5
         print *,'final Ws',A_WS4,A_WS5
         print *, IT,IT1,IW,IW1

         write(proc_num,'(i3.3)') myid
         filename = 'SNAP_'//TRIM(ADJUSTL(proc_num))//'.'//
     $        '100'
         OPEN (12,FILE=filename,STATUS='UNKNOWN',
     $        FORM='UNFORMATTED')
         REWIND(12)
         CALL PRINTOUT(1)
         CLOSE(12)
         stop
      ENDIF
      RETURN
      END SUBROUTINE columnOPC

      SUBROUTINE READCOLUMNDEPTHOPCTABLE
!--- ikaptyp = 6 : Freedman(2007) opacities (updated 4/2008) normal table, column opacities
!--- ikaptyp = 7 : Freedman(2007) opacities (updated 4/2008) without TIVO, column opacities
      USE opc_init
      USE mpi_var_init
      USE input_init
      IMPLICIT NONE
      include "mpif.h"
      character text*20
      integer :: i,j,n_levels
      integer :: Tindex,Pindex,Cindex,N_entries
      double precision :: readT,readP,readC
      double precision :: readkapABS,Tabsorption_read
      double precision :: level_min,level_max
      double precision :: tmp_val
      if(ikaptyp.eq.6) then
         open(12,file="Freedman6117_column.dat",status='old')
         if(myid.eq.0) print *,'reading Freedman6117_column.dat'
      elseif(ikaptyp.eq.7) then
         open(12,file="Freedman6117_column_noTIVO.dat",status='old')
         if(myid.eq.0) print *,'reading Freedman6117_column_noTiVO.dat'
      else
         print *,'ikaptyp=',ikaptyp,'not avail in readcolumnopctable' 
         stop
      endif
      N_entries=5888
!-- Initially set everything to -7 to account for missing values
      ColumnXKAP_abs = -7.d0
!-- Skip header
      DO i=1,10
         READ(12,*) text
         if(myid.eq.0) then
            print *,text
         endif
      ENDDO
      READ(12,*) Tabsorption_read
      if(myid.eq.0) print *,'Absorption opacities calculated 
     %     with T_star=',Tabsorption_read

      READ(12,*) n_levels
      if(myid.eq.0) print *,'Absorption opacities calculated 
     %     with ',n_levels,' column-levels'

      READ(12,*) level_min,level_max
      if(myid.eq.0) print *,'levels run from',
     %     level_min,'to',level_max
      
!-- read in actual values
      DO i=1,N_entries
         READ(12,*) readT,readP,readC,readkapABS
         do j=1,42
            if(FreedmanTVAL(j).eq.readT) goto 100
         enddo
 100     CONTINUE
         Tindex = j
         do j=1,18
            if(FreedmanPVAL(j).eq.readP) goto 101
         enddo
 101     CONTINUE
         Pindex = j
         do j=1,8
            if(FreedmanCVAL(j).eq.readC) goto 102
         enddo
 102     CONTINUE
         Cindex = j
         ColumnXKAP_abs(Tindex,Pindex,Cindex)  = readkapABS
      enddo
      CLOSE(12)
!--- TEST PRINTOUT
c      if(myid.eq.0) then
c         print *,'Freedman Opacities:',FreedmanXKAP_plnk(1,2)
c         do Tindex=1,42
c            do Pindex=1,18
c               do Cindex=1,8
c                  write(*,'(5(1x,e12.4))') FreedmanTVAL(Tindex),
c     %                 FreedmanPVAL(Pindex),FreedmanCVAL(Cindex),
c     %                 ColumnXKAP_abs(Tindex,Pindex,Cindex)
c               enddo
c            enddo
c         enddo
c      endif

  
      do Tindex=1,42
         do Pindex=1,18
            tmp_val = ColumnXKAP_abs(Tindex,Pindex,8)
            do Cindex=7,1,-1
               if(ColumnXKAP_abs(Tindex,Pindex,Cindex).gt.tmp_val) then
                  ColumnXKAP_abs(Tindex,Pindex,Cindex)=tmp_val
               endif
               tmp_val = ColumnXKAP_abs(Tindex,Pindex,8)
            enddo
         enddo
      enddo


      RETURN
      END SUBROUTINE  READCOLUMNDEPTHOPCTABLE
