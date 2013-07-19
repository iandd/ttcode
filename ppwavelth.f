!----------
! IRAD=5 is radial flux only
! IRAD=6 is radial + explicit solver for horiz. components
!----------
      SUBROUTINE WAVELENGTH_RADIAL_HEATING
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      USE opc_init
      USE mpi_var_init
      USE mpi_grid_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k,bin,ip
      double precision :: dfdr(locNX,locNY,locNZ)
      double precision :: bottom_flux_sum,diffusivity_FAC
      double precision :: calc_extrapolated_local_BB
      double precision :: tau_switch_rad_routine

!-BOUNDARY CONDITIONS
      CALL BOUNDS

c      print *,'wave BF kappa',myid
!-CALCULATE THE ALL BINNED OPACITIES AND SOURCE FUNCTION XKAPB(b), xkapW(b), local_BB(b)
      CALL KAPPA
c      print *,'wave AF kappa',myid

!-CALCULATE ALL TAU[b] and DTAU(b). CALCULATED ON THE A-GRID
      tau_LW(:,:,:,:)  = 0.d0
      dtau_LW(:,:,:,:) = 0.d0
      DO K=1,locNZ
         DO J=1,locNY
             DO I=upperbnd(j,k),1,-1
               dtau_LW(i,j,k,:) = rh(i,j,k)*xkapLW(i,j,k,:)*dxa(i)
               tau_LW(i,j,k,:) = tau_LW(i+1,j,k,:)+dtau_LW(i,j,k,:)
            ENDDO
            dtau_LW(0,j,k,:) = dtau_LW(1,j,k,:)
            tau_LW(0,j,k,:) = tau_LW(1,j,k,:)+dtau_LW(0,j,k,:)
         ENDDO
      ENDDO
!-augment tau and dtau with a diffusivity-factor
      diffusivity_FAC=1.66
      tau_LW =  diffusivity_FAC*tau_LW
      
c-set where to make the switch from lambda-dep to diffusive
      tau_switch_rad_routine=2.5

!-CALCULATED F_UP
      flux_r_up(:,:,:,:) = 0.d0
      DO K=1,locNZ
         DO J=1,locNY
            DO bin=1,nwave_bins
               DO i=1,upperbnd(j,k)+1
                  if(dtau_LW(i-1,j,k,bin).ge.tau_switch_rad_routine.or.
     $                 i.eq.1) then !-optically thick region
                     flux_r_up(i,j,k,bin) = 
     $                    CALC_EXTRAPOLATED_LOCAL_BB(i,j,k,bin,.true.)
                  else !-optically thin region
!-the exponetially attenuated initial source
                     flux_r_up(i,j,k,bin)=local_BB(1,j,k,bin)* 
     $                    (exp(-(tau_LW(1,j,k,bin)-tau_LW(i,j,k,bin))))
!-integrate over the rest of the column below i and exponetially attenuated it
                     DO ip=1,i-1
                        flux_r_up(i,j,k,bin) = flux_r_up(i,j,k,bin) + 
     $                       local_BB(ip,j,k,bin)*
     $                 (exp(-(tau_LW(ip+1,j,k,bin)-tau_LW(i,j,k,bin))))*
     $                       (1.d0 - exp(-diffusivity_FAC*
     $                       dtau_LW(ip,j,k,bin)))
                     ENDDO
                  endif
               ENDDO
! constant luminosity above upperbnd
               do i=upperbnd(j,k)+2,locNX+1
                  flux_r_up(i,j,k,bin)=
     $                 flux_r_up(upperbnd(j,k)+1,j,k,bin)*
     $                 ((xxa(upperbnd(j,k)+1)/xxa(i))**2.0)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

!-CALCULATE F_DOWN
      flux_r_down(:,:,:,:)=0.d0
      DO K=1,locNZ
         DO J=1,locNY
            DO bin=1,nwave_bins
c               flux_r_down(upperbnd(j,k):locNX+1,j,k,bin)=
c     $              ((RSTAR*RSOL/semi)**(2.0))*bstar_bin(bin)
               DO i=upperbnd(j,k),1,-1
                  if(dtau_LW(i,j,k,bin).ge.tau_switch_rad_routine) then                      
                     flux_r_down(i,j,k,bin) = 
     $                    CALC_EXTRAPOLATED_LOCAL_BB(i,j,k,bin,.false.)
                  else

! I could put the exponentially attenuated stellar BB here
c                  flux_r_down(i,j,k,bin)=((RSTAR*RSOL/semi)**(2.0))*
c     $                 bstar_bin(bin)*
c     $           exp(-(tau_LW(i,j,k,bin)-tau_LW(upperbnd(j,k),j,k,bin)))
! second, integrate over the rest of the column above i and exponetially attenuated it
                     DO ip=upperbnd(j,k),i,-1
                        flux_r_down(i,j,k,bin) = flux_r_down(i,j,k,bin)+
     $                       local_BB(ip,j,k,bin)* 
     $                   (exp(-(tau_LW(i,j,k,bin)-tau_LW(ip,j,k,bin))))*
     $                       (1.d0-exp(-diffusivity_FAC*
     $                       dtau_LW(ip,j,k,bin)))
                     ENDDO
                  endif
               ENDDO
            ENDDO
         ENDDO
      ENDDO      
      
!-CALCULATE THE NET RADIAL FLUX (defined to be positive in the positive
!- radial dir) OVER ALL THE WAVELENGTH BINS
      flux_x = 0.d0
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX+1
               flux_x(i,j,k) = SUM(flux_r_up(i,j,k,:)) -
     $              SUM(flux_r_down(i,j,k,:))
            ENDDO
         ENDDO
      ENDDO

!-cALCULATE THE FLUX IN THE Y AND Z DIRCETIONS (THIS DOESN'T SEEM TO WORK)
c      ER=ARAD*(T**(4.0))
c      CALL PLANETSHIFT(6,-1)
cc!--   now correct the polar shifts if necessary
c      if(poles) then
c         if(proc_coords(3).eq.proc_dims(3)-1) then
c            CALL NORTH_POLE_SHIFT(6,-1)
c         endif
c         if(proc_coords(3).eq.0) then
c            CALL SOUTH_POLE_SHIFT(6,-1)
c         endif
c      endif
c      CALL CALCULATEFLUX        !-note: flux_x not calculated for irad=5
c      do k=1,locNZ
c         do j=1,locNY+1
c            FLUX_Y(upperbnd(j,k):locNX,j,k) = 0.d0
c         enddo
c      enddo
c      do k=1,locNZ+1
c         do j=1,locNY
c            FLUX_Z(upperbnd(j,k):locNX,j,k) = 0.d0
c         enddo
c      enddo
      FLUX_Y=0.d0
      FLUX_Z=0.d0

!-calculate the div(F) term using the fluxes
      CALL CALCULATEDIVF_WITHF(dfdr,1.d0)

!-PRINT OUT LOCAL_BB AT THE BOTTOM SUMED OVER ALL THE WAVELENGTH BINS
      if(myid.eq.0.and.III.eq.1) then
         write(*,'(A,1(1x,e11.5))') ' Flux from the interior=',
     $        flux_x(1,1,1)
      ENDIF
!-PRINT OUT LOCAL_BB AT THE BOTTOM SUMED OVER ALL THE WAVELENGTH BINS
      if(myid.eq.0.and.III.eq.1) then
         print *,''
         write(*,'(A,1(1x,e11.5))') ' Flux from the exterior=',
     $        flux_x(upperbnd(1,1),1,1)
         write(*,'(A,2(1x,e11.5))') ' sum(bstar_bin),sig*Tstar^4=',
     $        ((RSTAR*RSOL/semi)**(2.0))*sum(bstar_bin(:)),
     $        ((RSTAR*RSOL/semi)**(2.0))*sbconst*(5040.**4.0)

         write(*,'(A,1(1x,e11.5))') ' sum(fup(ub))=',
     $        sum(flux_r_up(upperbnd(1,1),1,1,:))
         write(*,'(A,1(1x,e11.5))') ' sum(fdown(ub))=',
     $        sum(flux_r_down(upperbnd(1,1),1,1,:))


         write(*,'(A,1(1x,e11.5))') ' sum(fup(1))=',
     $        sum(flux_r_up(1,1,1,:))
         write(*,'(A,1(1x,e11.5))') ' sum(fdown(1))=',
     $        sum(flux_r_down(1,1,1,:))

      endif


!-UPDATE T WITH d(Fr)/dr, AND STELLAR HEATING
      if(F_INCIDENT.gt.1) then  !-- implicit stellar heating
         CALL STELLARINPUT_HEATING
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,upperbnd(j,k)
                  T(I,J,K)=T(I,J,K) + DELT*(
     $                 -dfdr(i,j,k) + STELLARINPUT(i,j,k))/
     $                 (RH(I,J,K)*CV(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      else
!-UPDATE T WITH d(Fr)/dr
         DO K=1,locNZ
            DO J=1,locNY
               DO I=1,upperbnd(j,k)
                  T(I,J,K)=T(I,J,K) - DELT*dfdr(i,j,k)/
     $                 (RH(I,J,K)*CV(I,J,K))
               ENDDO
            ENDDO
         ENDDO
      endif


!-PRINT OUT FLUX INFORMATION
      if(III.eq.1) then 
         CALL PRINT_RADTRANSFER
      endif


!-CALL THE IMPLICIT SOLVER FOR Y and Z DIRECTIONS
      if(irad.eq.6) then 
         CALL RADDIF            !-difrx=0 for irad=6 in DIFF
      endif
!-BOUNDARY CONDITIONS
      CALL BOUNDS

      RETURN
      END SUBROUTINE WAVELENGTH_RADIAL_HEATING


      FUNCTION CALC_EXTRAPOLATED_LOCAL_BB(i_indx,j_indx,k_indx,bin_indx,
     $     Fupper)
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE rad_var_init
      USE grid_var_init
      USE opc_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i_indx,j_indx,k_indx,bin_indx
      double precision :: calc_extrapolated_local_BB
      double precision :: Tslope,path_length,Tlinear
      logical :: Fupper
      Tslope = (T(i_indx,j_indx,k_indx)-T(i_indx-1,j_indx,k_indx))/
     $     (xxb(i_indx)-xxb(i_indx-1))
      if(Fupper) then !-use kappa*rh from lower cell and extrapolate downward
         if(i_indx.eq.1) then
            path_length=1.d0/(rh(i_indx,j_indx,k_indx)*
     $           xkapLW(i_indx,j_indx,k_indx,bin_indx))
         else
            path_length=1.d0/(rh(i_indx-1,j_indx,k_indx)*
     $           xkapLW(i_indx-1,j_indx,k_indx,bin_indx))
         endif
         Tlinear = Tslope*(xxa(i_indx) - path_length - xxb(i_indx-1))+
     $        T(i_indx-1,j_indx,k_indx)
      else !-use kappa*rh from upper cell and extrapolate upward
         path_length=1.d0/(rh(i_indx,j_indx,k_indx)*
     $        xkapLW(i_indx,j_indx,k_indx,bin_indx))
         Tlinear = Tslope*(xxa(i_indx) + path_length - xxb(i_indx-1))+
     $        T(i_indx-1,j_indx,k_indx)
      endif
!-perhaps eventually I could pull out just the local_BB portion out of
!-burrowsOPC, it would probably be much faster
      call BurrowsOPC(Tlinear,RH(i_indx,j_indx,k_indx),
     $     xkapP(i_indx,j_indx,k_indx),xkapR(i_indx,j_indx,k_indx),
     $     xkapA(i_indx,j_indx,k_indx),xkapW(i_indx,j_indx,k_indx,:),
     $     xkapLW(i_indx,j_indx,k_indx,:),
     $     local_BB(i_indx,j_indx,k_indx,:),
     $     dxkapP_dT(i_indx,j_indx,k_indx),
     $     dxkapA_dT(i_indx,j_indx,k_indx),i_indx,j_indx,k_indx)
      calc_extrapolated_local_bb=local_BB(i_indx,j_indx,k_indx,bin_indx)
      RETURN
      END FUNCTION CALC_EXTRAPOLATED_LOCAL_BB


      SUBROUTINE PRINT_RADTRANSFER
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,j,k,bin
      character proc_num*3
      character filename*12
      character print_num*3

      write(print_num,'(i3.3)') III
      write(proc_num,'(i3.3)') myid
      
      filename = 'WLTH_'//TRIM(ADJUSTL(proc_num))//'.'//
     $     print_num
      OPEN (13,FILE=filename,STATUS='UNKNOWN',
     $     FORM='UNFORMATTED')
      WRITE (13) ((((flux_r_up(i,j,k,bin),I=1,locNX),J=1,locNY),
     %     K=1,locNZ),bin=1,nwave_bins)
      WRITE (13) ((((flux_r_down(i,j,k,bin),I=1,locNX),J=1,locNY),
     %     K=1,locNZ),bin=1,nwave_bins)
      CLOSE(13)
      RETURN
      END SUBROUTINE PRINT_RADTRANSFER





      SUBROUTINE GREY_TWOSTREAM
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      USE opc_init
      USE mpi_var_init
      USE mpi_grid_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k,ip
      double precision :: dfdr(locNX,locNY,locNZ)
      double precision :: diffusivity_FAC
      double precision :: tau_switch_rad_routine
      double precision :: dtau_grey(locNX,locNY,locNZ)
!-BOUNDARY CONDITIONS
      CALL BOUNDS
!-OPACITY
      CALL KAPPA
!-CALCULATE TAU and DTAU. CALCULATED ON THE A-GRID
c      tau_grey(:,:,:)  = 0.d0
c      dtau_grey(:,:,:) = 0.d0
      DO K=1,locNZ
         DO J=1,locNY
            tau_grey(locNX+1,j,k) = 0.d0
            DO I=locNX,1,-1
               dtau_grey(i,j,k) = rh(i,j,k)*xkapP(i,j,k)*dxa(i)
               tau_grey(i,j,k) = tau_grey(i+1,j,k)+dtau_grey(i,j,k)
            ENDDO
         ENDDO
      ENDDO
!-augment tau and dtau with a diffusivity-factor
      diffusivity_FAC=1.66
      tau_grey =  diffusivity_FAC*tau_grey

!-CALCULATED F_UP
c      flux_r_up_grey(:,:,:) = 0.d0
      DO K=1,locNZ
         DO J=1,locNY
            flux_r_up_grey(1,j,k)=SBCONST*(T(0,j,k)**4.d0)
            DO i=2,locNX+1
!-the exponetially attenuated initial source
               flux_r_up_grey(i,j,k)=SBCONST*(T(0,j,k)**4.d0)*
     $              (exp(-(tau_grey(1,j,k)-tau_grey(i,j,k))))
c!-integrate over the rest of the column below i and exponetially attenuated it
               DO ip=1,i-1
                  flux_r_up_grey(i,j,k) = flux_r_up_grey(i,j,k) + 
     $                 SBCONST*(T(ip,j,k)**4.0)*
     $                 (exp(-(tau_grey(ip+1,j,k)-tau_grey(i,j,k))))*
     $                 (1.d0 - exp(-diffusivity_FAC*
     $                 dtau_grey(ip,j,k)))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
!-CALCULATE F_DOWN
      flux_r_down_grey(:,:,:)=0.d0
      DO K=1,locNZ
         DO J=1,locNY
            flux_r_down_grey(locNX+1,j,k) = 0.d0
c            flux_r_down_grey(locNX,j,k) = 0.d0
            DO i=locNX,1,-1
c            DO i=locNX-1,1,-1
! integrate over the rest of the column above i and exponetially attenuated it
               DO ip=locNX,i,-1
c               DO ip=locNX-1,i,-1
                  flux_r_down_grey(i,j,k) = flux_r_down_grey(i,j,k)+
     $                 SBCONST*(T(ip,j,k)**4.d0)* 
     $                 (exp(-(tau_grey(i,j,k)-tau_grey(ip,j,k))))*
     $                 (1.d0-exp(-diffusivity_FAC*
     $                 dtau_grey(ip,j,k)))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
!-CALCULATE THE NET RADIAL FLUX (defined to be positive in the positive
!- radial dir) OVER ALL THE WAVELENGTH BINS
      flux_x = 0.d0
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX+1
               flux_x(i,j,k) = flux_r_up_grey(i,j,k) - 
     $              flux_r_down_grey(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      FLUX_Y=0.d0
      FLUX_Z=0.d0

!-calculate the div(F) term using the fluxes
      CALL CALCULATEDIVF_WITHF(dfdr,1.d0)

!-UPDATE T WITH d(Fr)/dr
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,upperbnd(j,k)
               T(I,J,K)=T(I,J,K) - DELT*(
     $              dfdr(i,j,k))/(RH(I,J,K)*CV(I,J,K))
            ENDDO
         ENDDO
      ENDDO


!-PRINT OUT FLUX INFORMATION

      IF ((IREA/10).ne.0 .and. 
     &     ((III-1)/(IREA/10))*(IREA/10).eq.(III-1)) THEN
c      if(III.eq.1) then          
         CALL PRINT_GREY_RADTRANSFER
      endif


!-BOUNDARY CONDITIONS
      CALL BOUNDS
      RETURN
      END SUBROUTINE GREY_TWOSTREAM




      SUBROUTINE PRINT_GREY_RADTRANSFER
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,j,k,bin
      character proc_num*3
      character filename*12
      character print_num*3

      write(print_num,'(i3.3)') num_print-1
      write(proc_num,'(i3.3)') myid
      
      filename = 'GREY_'//TRIM(ADJUSTL(proc_num))//'.'//
     $     print_num
      OPEN (13,FILE=filename,STATUS='UNKNOWN',
     $     FORM='UNFORMATTED')
      WRITE (13) (((tau_grey(i,j,k),I=1,locNX),J=1,locNY),
     %     K=1,locNZ)
      WRITE (13) (((flux_r_up_grey(i,j,k),I=1,locNX),J=1,locNY),
     %     K=1,locNZ)
      WRITE (13) (((flux_r_down_grey(i,j,k),I=1,locNX),J=1,locNY),
     %     K=1,locNZ)
      CLOSE(13)
      RETURN
      END SUBROUTINE PRINT_GREY_RADTRANSFER

