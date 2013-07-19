
      SUBROUTINE COUPLED_ENERGY_EQS
      USE input_init
      USE fluid_var_init
      USE force_var_init
      USE rad_var_init
      USE grid_var_init
      USE deltaER_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE global_constants
      USE hypre_var_init
      USE sor_var_init
      IMPLICIT NONE
      include "mpif.h"
      double precision :: locMAX_Q,locMAX_ER
      double precision :: global_deltaER,global_deltaQ
      double precision :: delta_max
      double precision :: delta_Q(locNX,locNY,locNZ)
      double precision :: Q(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: Q_old(-1:locNX+2,-1:locNY+2,-1:locNZ+2)
      double precision :: divF(locNX,locNY,locNZ)
      double precision :: divF_old(locNX,locNY,locNZ)
      double precision :: B_body_c(locNX,locNY,locNZ)
      double precision :: dFQ_dQ(locNX,locNY,locNZ)
      double precision :: dFR_dQ(locNX,locNY,locNZ)
      double precision :: theta_param,Q1,Q2,FluxR
      double precision :: Tincident,HEATTERM
      double precision :: bottomDQDR,Qupper,Qinitbottom,bottomDERDR
      double precision :: Fluxsum_top,kappaT_bottom
      double precision :: CALC_UNIVERSAL_BOTTOM_DERIV
      integer :: stat(MPI_STATUS_SIZE)

      integer :: i,j,k,indx,xmax_index,xmin_index
      double precision :: max_val,max_Qval,CALC_TINCIDENT
      integer :: max_iloc,max_jloc,max_kloc
      character non_converge_file*12
      character proc_num*3

!--- Set the convergence tolerance for Q
c       delta_max = 0.05d0       
      delta_max = 0.005d0       !-- (THIS WORKS WITH SOR, HYPRE worked for 0.05) !-req. for wasp12?
c      delta_max = 0.0005d0       !-- (THIS WORKS WITH SOR, HYPRE worked for 0.05)

!--- convergence parameter for updating Q and ER
c      conv_parm = 1.0
c      conv_parm = 0.5
c      conv_parm = 0.25
      conv_parm = 0.1

!--- zero all changes (needed for the points outside normal bounds)
      delta_ER  = 0.d0
      delta_Q   = 0.d0

      energy_ITER = 0 !--- define here for divf calculation

!--- Set type of update type - Careful here, time-centered may produce
!---- instabilities or not converge as well
c      theta_param = 0.5d0 ! time-centered (Crank-Nicolson) update
      theta_param = 1.d0 ! backwards Euler

!---  T (ie Q) boundary conditions
      CALL BOUNDS
!---  SET Q
      Q = RH*CV*T

!--- Calculate the stellar heating term
      if(F_INCIDENT.gt.1) then !-- implicit stellar heating
         CALL KAPPA
         CALL STELLARINPUT_HEATING
      endif

!--- CALCULATE THE FLUXES. REQUIRED IF USING E=2F/c 
      CALL CALCULATEFLUX
      CALL FLUX_BC

!-Q boundary conditions
!--upper Boundary
      first_Tincident = .true.      
      if(MPIleft.eq.MPI_PROC_NULL) then
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
!-integrate over the tau through the cell (this puts back the mu_star)
!-this should be valid for both the grey and wavelength dep. heating
c               Fluxsum_top = flux_x(upperbnd(j,k)-1,J,K)+
c     $              flux_y(upperbnd(j,k)-1,J,K)+
c     $              flux_z(upperbnd(j,k)-1,J,K)
c               Qupper = rh(upperbnd(j,k)-1,J,K)*cv(upperbnd(j,k)-1,J,K)*
cc     %              (( (abs(flux_x(upperbnd(j,k)-1,J,K))/(2.d0*SBCONST))
c     %              (((abs(Fluxsum_top)/(2.d0*SBCONST))
c     $              + (stellarinput(upperbnd(j,k)-1,j,k)/
c     $              (4.d0*SBCONST*rh(upperbnd(j,k)-1,j,k)*
c     $              xkapP(upperbnd(j,k)-1,j,k))) )**(0.25))
               Qupper = rh(upperbnd(j,k)-1,J,K)*cv(upperbnd(j,k)-1,J,K)*
     $              T(upperbnd(j,k),j,k)
               DO I=upperbnd(j,k),locNX+2
                  Q(i,j,k) = Qupper !--dQ/dr = 0
                  ER(I,J,K)=abs(2.d0*flux_x(i,j,k)/FD) !--Er = 2F/c
!--Er = 2F/c + 2*Finc: MARSHAK BOUNDARY CONDITIONS
c                  ER(I,J,K)=abs(2.d0*flux_x(i,j,k)/FD) + 
c     $                 (2.d0*SBCONST*(CALC_TINCIDENT(j,k)**4.0)/FD)
               ENDDO
            ENDDO
         ENDDO
      endif

!---- LEFT BOUNDARY (ie XMIN)
      if(MPIleft.eq.MPI_PROC_NULL) then
         if(BCXMIN.eq.8)        !-CALCULATE UNIVERSAL dQ/dr
     $        bottomDQDR = CALC_UNIVERSAL_BOTTOM_DERIV(.true.)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               if(BCXMIN.eq.5.or.BCXMIN.eq.7) then !-- fixed Q/ER derivatives
                  Q(0,j,K)   = Q(1,J,K) - 
     %                 RH(1,j,k)*cv(1,j,k)*D_T*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - 
     %                 RH(1,j,k)*cv(1,j,k)*D_T*(xxb(1)-xxb(-1))
                  ER(0,j,K)   = ER(1,J,K) -  
     %                 4.d0*ARAD*(T(1,j,k)**(3.d0))*D_T*(xxb(1)-xxb(0))
                  ER(-1,j,K)  = ER(1,J,K) -
     %                 4.d0*ARAD*(T(1,j,k)**(3.d0))*D_T*(xxb(1)-xxb(-1))
               elseif(BCXMIN.eq.6) then !-- fixed FLUX... fixed Q/ER values and derivative
                  Qinitbottom = RHinitbottom*cv(1,j,k)*Tinitbottom
                  bottomDQDR = -3.d0*Fbottom*xkapR(1,j,k)*
     $                 (RH(1,j,k)**5.d0)*(cv(1,j,k)**4.d0)/
     $                 (4.d0*ARAD*FD*(Qinitbottom**3.d0))+ 
     %                 (Qinitbottom*D_RH/rh(1,j,k))
                  Q(1,j,K)   = Qinitbottom
                  Q(2,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(2))
                  Q(0,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(-1))
!--   optically thick = ER=aT^4
                  ER(2,j,K)   = ARAD*((Q(2,j,k)/
     $                 (RH(2,j,k)*cv(2,j,k)))**(4.d0))
                  ER(1,j,K)   = ARAD*((Q(1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(0,j,k)*cv(0,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(-1,j,k)*cv(-1,j,k)))**(4.d0))
               elseif(BCXMIN.eq.8) then !--fixed flux, use universal dQ/dr (defined above), optically thick ER
                  Q(0,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(-1))
                  ER(1,j,K)   = ARAD*((Q(1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(0,j,k)*cv(0,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(-1,j,k)*cv(-1,j,k)))**(4.d0))
               elseif(BCXMIN.eq.9) then !--fixed flux, but free T and Q
                  if(ieta.eq.5) then !-fixed PR and viscosity at the bottom
                     bottomDQDR = -RH(1,j,k)*cv(1,j,k)*Fbottom/
     $                    (RH(1,j,k)*cv(1,j,k)*(XNUE/PR_NUM))
                  else
                     bottomDQDR = -RH(1,j,k)*cv(1,j,k)*
     $                    3.d0*Fbottom*xkapR(1,j,k)*RH(1,j,k)/
     $                    (4.d0*ARAD*FD*(T(1,j,k)**3.d0))
                  endif
                  Q(0,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(-1))
!-- optically thick = ER=aT^4
                  ER(1,j,K)   = ARAD*((Q(1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(0,j,k)*cv(0,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(-1,j,k)*cv(-1,j,k)))**(4.d0))
               elseif(BCXMIN.eq.11) then !--fixed values and derivatives set by ppinit
                  bottomDQDR = cv(1,j,k)* (TinitBottom*D_RH +
     $                 rhinitbottom*D_T)
                  Q(1,j,K)   = cv(1,j,k)*RHinitbottom*Tinitbottom
                  Q(2,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(2))
                  Q(0,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(-1))
!--   optically thick = ER=aT^4
                  ER(2,j,K)   = ARAD*((Q(2,j,k)/
     $                 (RH(2,j,k)*cv(2,j,k)))**(4.d0))
                  ER(1,j,K)   = ARAD*((Q(1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(0,j,k)*cv(0,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(-1,j,k)*cv(-1,j,k)))**(4.d0))
               elseif(BCXMIN.eq.12) then !--fixed ER and T
                  bottomDQDR = cv(1,j,k)*rh(1,j,k)*D_T
                  Q(0,j,K)   = Q(1,J,K)+bottomDQDR*(xxb(0)-xxb(1))
                  Q(-1,j,K)  = Q(1,J,K)+bottomDQDR*(xxb(-1)-xxb(1))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
c                  Q(0,j,K)   = RH(0,j,k)*cv(1,j,k)*T(0,J,K)
c                  Q(-1,j,K)  = RH(-1,j,k)*cv(1,j,k)*T(-1,J,K)
c                  ER(0,j,K)   = ARAD*(T(0,J,K)**(4.d0))
c                  ER(-1,j,K)   = ARAD*(T(-1,J,K)**(4.d0))
               else
                  print *,'bcxmin=',bcxmin,'not yet defined: PPENERGY'
               endif
            enddo
         enddo
      endif

!---- z boundaries
      if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
         DO J=-1,locNY+2
            DO i=-1,locNX+2
               Q(i,j,1)=Q(i,j,2)
               Q(i,j,0)=Q(i,j,1)
               Q(i,j,-1)=Q(i,j,1)
               ER(I,j,1)=ER(I,j,2)
               ER(I,j,0)=ER(I,j,1)
               ER(I,j,-1)=ER(I,j,1)
            ENDDO
         ENDDO
      endif
      if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL) then
         DO J=-1,locNY+2
            DO I=-1,locNX+2
               Q(I,j,locNZ)=Q(I,j,locNZ-1)
               Q(i,j,locNZ+1)=Q(i,j,locNZ)
               Q(i,j,locNZ+2)=Q(i,j,locNZ)
               ER(i,j,locNZ)=ER(i,j,locNZ-1)
               ER(i,j,locNZ+1)=ER(i,j,locNZ)
               ER(i,j,locNZ+2)=ER(i,j,locNZ)
            ENDDO
         ENDDO
      endif

!--- Call ER boundary passing routine  (for divF calculation)
      if(modtyp.eq.2) then
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
      elseif(modtyp.eq.6) then
         CALL SHIFTVAR_AXISYMMETRIC(6)
      else
         CALL NEWSHIFTVAR3D(6)
      endif
      
!--- save the original Q and ER
      Q_old = Q
      ER_old = ER

!--- Save the \divF at n for FQ calculation (and opacities)
c      CALL CALCULATEDIVF(divf_old,theta_param)
      CALL CALCULATEFLUX
      CALL FLUX_BC
      CALL CALCULATEDIVF_WITHF(divf_old,theta_param)

!----------------------------------------------------------
!--- Iterate until ITMAX or convergence
 4    CONTINUE
      energy_ITER = energy_ITER + 1
c      if(myid.eq.0) print *,'energyiter=',energy_ITER

!--- Call ER boundary passing routine  (for divF calculation)
      if(modtyp.eq.2) then
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
      elseif(modtyp.eq.6) then
         CALL SHIFTVAR_AXISYMMETRIC(6)
      else
         CALL NEWSHIFTVAR3D(6)
      endif

!--- Calculate the divF term for n+1 (and opacities)
c      CALL CALCULATEDIVF(divF,theta_param)
      CALL CALCULATEFLUX
      CALL FLUX_BC
      CALL CALCULATEDIVF_WITHF(divF,theta_param)

      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX
!--- Calculate the Blackbody flux B for term n+1 : b_body_c = 4*sigma*T^4/c = arad*T^4
               B_body_c(I,J,K) = arad*((Q(I,J,K)/
     1              (RH(I,J,K)*CV(I,J,K)))**(4.d0))
!--- Calculate the F parameters
               if(F_INCIDENT.gt.1) then !-- implicit stellar heating
c-old grey method (now its pre-computed)
c                  Tincident = CALC_TINCIDENT(j,k)
c                  STELLARINPUT(I,J,K) = SBCONST*
c     $                 ((Tincident**4.d0)/mu_star(j,k))*
c     &                 exp(-STELLARTAU(I,J,K)/mu_star(j,k))*
c     $                 dtaus(i,j,k)/dxb(i)
                  HEATTERM = STELLARINPUT(I,J,K)/(RH(I,J,K)*
     &                 xkapA(I,J,K))
                  FQ(i,j,k) = Q(I,J,K) - Q_old(I,J,K) + DELT*
     1                 ( (1-theta_param)*(divF_old(I,J,K) - divF(I,J,K))
     2                 + (RH(I,J,K)*XkapP(I,J,K)*FD*
     3                 (B_body_c(I,J,K)-ER(I,J,K))) - 
     4                 STELLARINPUT(I,J,K) )
               else
                  FQ(i,j,k) = Q(I,J,K) - Q_old(I,J,K) + DELT*
     1                 ( (1-theta_param)*(divF_old(I,J,K) - divF(I,J,K))
     2                 + RH(I,J,K)*XkapP(I,J,K)*FD*
     3                 (B_body_c(I,J,K)-ER(I,J,K)))
               endif
               FR(I,J,K)  = ER(I,J,K) - ER_old(I,J,K) + DELT*(
     1              divF(I,J,K) - RH(I,J,K)*XkapP(I,J,K)*FD*
     2              (B_body_c(I,J,K)-ER(I,J,K)) )
!--- Calculate the LOCALLY EVALUATED Jacobian componenets and eta
               if(F_INCIDENT.gt.1) then !-- implicit stellar heating
                  dFQ_dQ(I,J,K) = 1.d0 + (delt/cv(I,J,K))*(xkapP(I,J,K)*
     1                 4.d0*ARAD*FD*((Q(I,J,K)/(RH(I,J,K)*
     2                 CV(I,J,K)))**(3.d0)) + 
     3                 (FD*(B_body_c(I,J,K)-ER(I,J,K))*dxkapP_dT(I,J,K))
     4                 - HEATTERM*dxkapA_dT(i,j,k) ) 
               else !-- no stellar heating
                  dFQ_dQ(I,J,K) = 1.d0 + (delt/cv(I,J,K))*(xkapP(I,J,K)*
     1                 4.d0*ARAD*FD*((Q(I,J,K)/(RH(I,J,K)*
     2                 CV(I,J,K)))**(3.d0)) + 
     3                 FD*(B_body_c(I,J,K)-ER(I,J,K))*dxkapP_dT(I,J,K) )
               endif
               dFR_dQ(I,J,K) = -(delt/cv(I,J,K))*(xkapP(I,J,K)*
     1              4.d0*ARAD*FD*((Q(I,J,K)/(RH(I,J,K)*
     2              CV(I,J,K)))**(3.d0)) + 
     3              FD*(B_body_c(I,J,K)-ER(I,J,K))*dxkapP_dT(I,J,K) )
               if(dFQ_dQ(I,J,K).ne.0.d0) then
                  Jacob_eta(I,J,K) = -dFR_dQ(I,J,K)/dFQ_dQ(I,J,K)
               else
                  print *,'JACOB_ETA WOULD BE INFINITY:',
     $                 dFR_dQ(I,J,K),dFQ_dQ(I,J,K)
                  Jacob_eta(I,J,K) = 0.d0                  
               endif
!--- CALCULATE THE RHS for the \delta_ER Equation (cannot do as matrix)
               RHS_ER(I,J,K) = -FR(I,J,K) - Jacob_eta(I,J,K)*FQ(I,J,K)
            enddo
         enddo
      enddo

!--- RE-FILL E7 AND E8 NOW THAT JACOB_ETA IS KNOWN
      CALL FILL_EVALUES(theta_param)


!--- SOLVE FOR \delta_ER
      IF (IRAD.EQ.2) THEN !---  USE HYPRE
c         if(myid.eq.0) print *,'starting hypre',III
         CALL SOLVE_WITH_HYPRE
         if(modtyp.eq.2) then
            CALL PLANETSHIFT(8,-1)
!--   now correct the polar shifts if necessary
            if(poles) then
               if(proc_coords(3).eq.proc_dims(3)-1) then
                  CALL NORTH_POLE_SHIFT(8,-1)
               endif
               if(proc_coords(3).eq.0) then
                  CALL SOUTH_POLE_SHIFT(8,-1)
               endif
            endif
         elseif(modtyp.eq.6) then
            CALL SHIFTVAR_AXISYMMETRIC(8)
         else
            CALL NEWSHIFTVAR3D(8)
         endif
      ELSEIF(IRAD.eq.3) THEN !---  USE SOR
         if(energy_iter.eq.1.and.(mod(III,10000).eq.0.or.III.eq.1)) then
            call adjust_sor_param
         endif
         CALL SOLVE_WITH_SOR(.true.)
      ENDIF

!--- Calculate Algebreic \deltaQ, and find max of \delta's
      locMAX_Q  = 0.d0
      locMAX_ER = 0.d0
      max_val = 0.d0
      max_Qval = 0.d0
      do k=1,locNZ
         do j=1,locNY
            if(MPIright.eq.MPI_PROC_NULL) then
               xmax_index = upperbnd(J,K)-1 !-- don't include upperbnd (cause its set with BC)
            else
               xmax_index = locNX
            endif
            if(BCXMIN.eq.6.or.BCXMIN.eq.11) then 
               xmin_index = 3 !-- bottom (i=1,2) is also set below (Apr 2010)
            elseif(BCXMIN.eq.12) then 
               xmin_index = 1
            else
               xmin_index = 2 !-- bottom (i=1) is also set below (May 6,2009)
            endif
            DO I=xmin_index,xmax_index
               Q1 = -delt*(1.d0-theta_param)
               Q2 = delt*RH(i,j,k)*xkapP(i,j,k)*FD
               delta_Q(i,j,k)  = (1.d0/dFQ_dQ(i,j,k))*(
     $              Q1*EVAL_1_6(i,j,k)*delta_ER(i,j,k) -
     $              Q1*EVAL1(i,j,k)*delta_ER(i+1,j,k) -
     $              Q1*EVAL2(i,j,k)*delta_ER(i-1,j,k) -
     $              Q1*EVAL3(i,j,k)*delta_ER(i,j+1,k) -
     $              Q1*EVAL4(i,j,k)*delta_ER(i,j-1,k) -
     $              Q1*EVAL5(i,j,k)*delta_ER(i,j,k+1) -
     $              Q1*EVAL6(i,j,k)*delta_ER(i,j,k-1) +
     $              Q2*delta_ER(i,j,k) - FQ(i,j,k) )

!---------- TEST TEST TEST---------- TEsT TEST TEST---------- TEST TEST TEST
               if( (Q(i,j,k)+conv_parm*delta_Q(i,j,k)) .le.0.d0) then 
                  PRINT *,'ERROR IN DELTA_Q: III, energy_iter=',
     $                 III,energy_iter
                  print *,'non-positive value of Q',
     $                 (Q(i,j,k)+conv_parm*delta_Q(i,j,k))
                  print *,'III,num_iter,myid',III,num_iter,myid
                  print *,'Q/rh',Q(i,j,k),RH(i,j,k),myid
                  print *,'(i,j,k,myid) of error',i,j,k,myid
                  print *,'(upperbnd,myid) of error',upperbnd(j,k)
     $                 ,myid
                  
                  write(*,'(A,11(1x,e12.6))') 
     $                 'deltaQ:',delta_Q(i,j,k), (1.d0/dFQ_dQ(i,j,k)),
     $                 Q1*EVAL_1_6(i,j,k)*delta_ER(i,j,k),
     $                 -Q1*EVAL1(i,j,k)*delta_ER(i+1,j,k),
     $                 -Q1*EVAL2(i,j,k)*delta_ER(i-1,j,k),
     $                 -Q1*EVAL3(i,j,k)*delta_ER(i,j+1,k),
     $                 -Q1*EVAL4(i,j,k)*delta_ER(i,j-1,k),
     $                 -Q1*EVAL5(i,j,k)*delta_ER(i,j,k+1),
     $                 -Q1*EVAL6(i,j,k)*delta_ER(i,j,k-1), 
     $                 Q2*delta_ER(i,j,k),-FQ(i,j,k)
                  write(*,'(A,I8,10(1x,e12.6))') 
     $                 'deltaQ(i,j,k) factors: ',myid,
     $                 Q2*delta_ER(i,j,k),FQ(i,j,k),
     $                 delt,RH(i,j,k),xkapP(i,j,k),
     $                 Q(I,J,K)-Q_old(I,J,K),
     3                 B_body_c(I,J,K),ER(I,J,K),
     4                 STELLARINPUT(I,J,K),T(i,j,k)
                  write(*,'(A,7(1x,e12.6))') 
     $                 'ER:',ER(i,j,k),ER(i+1,j,k),ER(i-1,j,k),
     $                 ER(i,j+1,k),ER(i,j-1,k),ER(i,j,k+1),
     $                 ER(i,j,k-1)
                  print *,'Q/Qupper',Q(i,j,k),Qupper,
     &                 upperbnd(j,k)
                  print *,'Tupper=Q/(rh*cv)',Qupper/(
     $                 rh(upperbnd(j,k)-1,J,K)*cv(upperbnd(j,k)-1,J,K))
                  write(*,'(A,7(1x,e12.6))') 'Qupper factors:',
     $                 Qupper,rh(upperbnd(j,k)-1,J,K), 
     $                 cv(upperbnd(j,k)-1,J,K),
     %                 flux_x(upperbnd(j,k)-1,J,K),
     $                 stellarinput(upperbnd(j,k)-1,j,k),
     $                 xkapP(upperbnd(j,k)-1,j,k),mu_star(j,k)
                  write(*,'(A,9(1x,e12.6))') 
     1                 'evals:',EVAL7(i,j,k),
     1                 EVAL8(i,j,k)*EVAL_1_6(i,j,k),
     $                 -EVAL8(i,j,k)*EVAL2(i,j,k),
     1                 -EVAL8(i,j,k)*EVAL1(i,j,k),
     1                 -EVAL8(i,j,k)*EVAL4(i,j,k),
     1                 -EVAL8(i,j,k)*EVAL3(i,j,k),
     1                 -EVAL8(i,j,k)*EVAL6(i,j,k),
     1                 -EVAL8(i,j,k)*EVAL5(i,j,k),
     1                 EVAL8(i,j,k)                 
                  write(*,'(A,7(1x,e12.6))') 
     1                 'FX:',flux_x(i,j,k),
     1                 flux_x(i+1,j,k),
     1                 flux_x(i-1,j,k),
     1                 flux_x(i,j+1,k),
     1                 flux_x(i,j-1,k),
     1                 flux_x(i,j,k+1),
     1                 flux_x(i,j,k-1)
                  write(*,'(A,7(1x,e12.6))') 
     1                 'FY:',flux_y(i,j,k),
     1                 flux_y(i+1,j,k),
     1                 flux_y(i-1,j,k),
     1                 flux_y(i,j+1,k),
     1                 flux_y(i,j-1,k),
     1                 flux_y(i,j,k+1),
     1                 flux_y(i,j,k-1)
                  write(*,'(A,7(1x,e12.6))') 
     1                 'FZ:',flux_z(i,j,k),
     1                 flux_z(i+1,j,k),
     1                 flux_z(i-1,j,k),
     1                 flux_z(i,j+1,k),
     1                 flux_z(i,j-1,k),
     1                 flux_z(i,j,k+1),
     1                 flux_z(i,j,k-1)
                  write(*,'(A,7(1x,e12.6))') 
     1                 'deltaER:',delta_ER(i,j,k),
     1                 delta_ER(i+1,j,k),
     1                 delta_ER(i-1,j,k),
     1                 delta_ER(i,j+1,k),
     1                 delta_ER(i,j-1,k),
     1                 delta_ER(i,j,k+1),
     1                 delta_ER(i,j,k-1)
                  write(*,'(A,7(1x,e12.6))') 
     1                 'RHS:',RHS_ER(i,j,k),
     1                 RHS_ER(i+1,j,k),
     1                 RHS_ER(i-1,j,k),
     1                 RHS_ER(i,j+1,k),
     1                 RHS_ER(i,j-1,k),
     1                 RHS_ER(i,j,k+1),
     1                 RHS_ER(i,j,k-1)
                  write(*,'(A,7(1x,e12.6))') 
     1                 'FQ:',FQ(i,j,k),
     1                 FQ(i+1,j,k),
     1                 FQ(i-1,j,k),
     1                 FQ(i,j+1,k),
     1                 FQ(i,j-1,k),
     1                 FQ(i,j,k+1),
     1                 FQ(i,j,k-1)
                  write(*,'(A,7(1x,e12.6))') 
     1                 'FR:',FR(i,j,k),
     1                 FR(i+1,j,k),
     1                 FR(i-1,j,k),
     1                 FR(i,j+1,k),
     1                 FR(i,j-1,k),
     1                 FR(i,j,k+1),
     1                 FR(i,j,k-1)
                  call clean_stop
                  
               endif




!--- fraction change
               locMAX_Q = MAX(locMAX_Q,abs(delta_Q(i,j,k)/Q(i,j,k)))
               locMAX_ER= MAX(locMAX_ER,abs(delta_ER(i,j,k)/ER(i,j,k)))
               
               if(energy_ITER.EQ.ITMAX.and.
     $              abs(locMAX_Q).gt.abs(max_val)) then
                  max_val = delta_Q(i,j,k)/Q(i,j,k) 
                  max_Qval = Q(i,j,k)
                  max_iloc = i
                  max_jloc = j
                  max_kloc = k
               endif

            enddo
         enddo
      enddo
      
c      if(myid.eq.0) then
c         write(*,'(I8,I8,2(1x,e12.6),4(1x,I6))') 
c     %        num_iter,energy_ITER,max_val,max_Qval,max_iloc,max_jloc,
c     %        max_kloc,upperbnd(max_jloc,max_kloc)
c      endif


!--- Update ER and Q
      DO K=1,locNZ
         DO J=1,locNY
            if(MPIright.eq.MPI_PROC_NULL) then
               xmax_index = upperbnd(J,K)-1 !-- don't include upperbnd (cause its set with BC)
            else
               xmax_index = locNX
            endif
            if(BCXMIN.eq.6.or.BCXMIN.eq.11) then 
               xmin_index = 3 !-- bottom (i=1,2) is also set below (Apr 2010)
            elseif(BCXMIN.eq.12) then 
               xmin_index = 1
            else
               xmin_index = 2 !-- bottom (i=1) is also set below (May 6,2009)
            endif
            DO I=xmin_index,xmax_index
               ER(i,j,k) = ER(i,j,k) + conv_parm*delta_ER(i,j,k)
               Q(i,j,k)  = Q(i,j,k) + conv_parm*delta_Q(i,j,k)
!--- LAPLACIAN TEST CASE
c               ER(i,j,k) = ER_OLD(i,j,k)
c               Q(i,j,k)  = Q_OLD(i,j,k)
            ENDDO
         ENDDO
      ENDDO

!--- BOUNDARY CONDITION ON Q AND ER
!---- LEFT BOUNDARY (ie XMIN)
      first_Tincident = .true.
      if(MPIleft.eq.MPI_PROC_NULL) then
         if(BCXMIN.eq.8)        !-CALCULATE UNIVERSAL dQ/dr
     $        bottomDQDR=CALC_UNIVERSAL_BOTTOM_DERIV(.true.)
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
               if(BCXMIN.eq.5.or.BCXMIN.eq.7) then !-- fixed Q/ER derivatives
                  Q(0,j,K)   = Q(1,J,K) - 
     %                 RH(1,j,k)*cv(1,j,k)*D_T*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - 
     %                 RH(1,j,k)*cv(1,j,k)*D_T*(xxb(1)-xxb(-1))
                  ER(0,j,K)   = ER(1,J,K) -  
     %                 4.d0*ARAD*(T(1,j,k)**(3.d0))*D_T*(xxb(1)-xxb(0))
                  ER(-1,j,K)  = ER(1,J,K) -
     %                 4.d0*ARAD*(T(1,j,k)**(3.d0))*D_T*(xxb(1)-xxb(-1))
               elseif(BCXMIN.eq.6) then !-- fixed FLUX... fixed Q/ER values and derivative
                  Qinitbottom = RHinitbottom*cv(1,j,k)*Tinitbottom
                  bottomDQDR = -3.d0*Fbottom*xkapR(1,j,k)*
     $                 (RH(1,j,k)**5.d0)*(cv(1,j,k)**4.d0)/
     $                 (4.d0*ARAD*FD*(Qinitbottom**3.d0))+ 
     %                 (Qinitbottom*D_RH/rh(1,j,k))
                  Q(1,j,K)   = Qinitbottom
                  Q(2,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(2))
                  Q(0,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(-1))
!--   optically thick = ER=aT^4
                  ER(2,j,K)   = ARAD*((Q(2,j,k)/
     $                 (RH(2,j,k)*cv(2,j,k)))**(4.d0))
                  ER(1,j,K)   = ARAD*((Q(1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(0,j,k)*cv(0,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(-1,j,k)*cv(-1,j,k)))**(4.d0))
               elseif(BCXMIN.eq.8) then !--fixed flux, use universal dQ/dr (defined above), optically thick ER
                  Q(0,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(-1))
                  ER(1,j,K)   = ARAD*((Q(1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(0,j,k)*cv(0,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(-1,j,k)*cv(-1,j,k)))**(4.d0))
               elseif(BCXMIN.eq.9) then !--fixed flux, but free T and Q
                  if(ieta.eq.5) then !-fixed PR and viscosity at the bottom
                     bottomDQDR = -RH(1,j,k)*cv(1,j,k)*Fbottom/
     $                    (RH(1,j,k)*cv(1,j,k)*(XNUE/PR_NUM))
                  else
                     bottomDQDR = -RH(1,j,k)*cv(1,j,k)*
     $                    3.d0*Fbottom*xkapR(1,j,k)*RH(1,j,k)/
     $                    (4.d0*ARAD*FD*(T(1,j,k)**3.d0))
                  endif
                  Q(0,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(-1))
!-- optically thick = ER=aT^4
                  ER(1,j,K)   = ARAD*((Q(1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(0,j,k)*cv(0,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(-1,j,k)*cv(-1,j,k)))**(4.d0))
               elseif(BCXMIN.eq.11) then !--fixed values and derivatives set by ppinit
                  bottomDQDR = cv(1,j,k)* (TinitBottom*D_RH +
     $                 rhinitbottom*D_T)
                  Q(1,j,K)   = cv(1,j,k)*RHinitbottom*Tinitbottom
                  Q(2,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(2))
                  Q(0,j,K)   = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(0))
                  Q(-1,j,K)  = Q(1,J,K) - bottomDQDR*(xxb(1)-xxb(-1))
!--   optically thick = ER=aT^4
                  ER(2,j,K)   = ARAD*((Q(2,j,k)/
     $                 (RH(2,j,k)*cv(2,j,k)))**(4.d0))
                  ER(1,j,K)   = ARAD*((Q(1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(0,j,k)*cv(0,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(-1,j,k)*cv(-1,j,k)))**(4.d0))
               elseif(BCXMIN.eq.12) then !--fixed flux, T, and rh
                  bottomDQDR = cv(1,j,k)*rh(1,j,k)*D_T
                  Q(0,j,K)   = Q(1,J,K)+bottomDQDR*(xxb(0)-xxb(1))
                  Q(-1,j,K)  = Q(1,J,K)+bottomDQDR*(xxb(-1)-xxb(1))
                  ER(0,j,K)   = ARAD*((Q(0,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
                  ER(-1,j,K)  = ARAD*((Q(-1,j,k)/
     $                 (RH(1,j,k)*cv(1,j,k)))**(4.d0))
c                  bottomDERDR = -Fbottom*
c     $                 3.d0*xkapRinitbottom*RHinitbottom/FD
c                  ER(0,j,K)   = ER(1,J,K)+bottomDERDR*(xxb(0)-xxb(1))
c                  ER(-1,j,K)  = ER(1,J,K)+bottomDERDR*(xxb(-1)-xxb(1))
               else
                  print *,'bcxmin=',bcxmin,'not yet defined: PPENERGY'
               endif
            enddo
         enddo
      endif

!-TRY ADDING THIS CALC FLUX/BC HERE SO NEW Q/ER BC IS CORRECT
c      CALL CALCULATEFLUX
c      CALL FLUX_BC

!-- UPPER BOUNDARY----------------------------------
      if(MPIleft.eq.MPI_PROC_NULL) then
         DO K=-1,locNZ+2
            DO J=-1,locNY+2
c- Q = T*(cv*rh) = (ER/a+gamma^4*Tirr^4*e(-tau_star/mu))^(0.25) *(cv*rh) with E=2F/c
c               Fluxsum_top = flux_x(upperbnd(j,k)-1,J,K)+
c     $              flux_y(upperbnd(j,k)-1,J,K)+
c     $              flux_z(upperbnd(j,k)-1,J,K)
c               Qupper = rh(upperbnd(j,k)-1,J,K)*cv(upperbnd(j,k)-1,J,K)*
cc     %              (( (abs(flux_x(upperbnd(j,k)-1,J,K))/(2.d0*SBCONST))
c     %              (((abs(Fluxsum_top)/(2.d0*SBCONST))
c     $              + (stellarinput(upperbnd(j,k)-1,j,k)/
c     $              (4.d0*SBCONST*rh(upperbnd(j,k)-1,j,k)*
c     $              xkapP(upperbnd(j,k)-1,j,k))) )**(0.25))
               Qupper = rh(upperbnd(j,k)-1,J,K)*cv(upperbnd(j,k)-1,J,K)*
     $              T(upperbnd(j,k),j,k)
               DO I=upperbnd(j,k),locNX+2
                  Q(i,j,k) = Qupper !-isothermal
                  ER(I,J,K)=abs(2.d0*flux_x(i,j,k)/FD) !--   Er = 2F/c
!--Er = 2F/c + 2*Finc: MARSHAK BOUNDARY CONDITIONS
c                  ER(I,J,K)=abs(2.d0*flux_x(i,j,k)/FD) + 
c     $                 (2.d0*SBCONST*(CALC_TINCIDENT(j,k)**4.d0)/FD)
               ENDDO
            ENDDO
         ENDDO
      endif


!---- z boundaries
      if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL) then
         DO J=-1,locNY+2
            DO i=-1,locNX+2
!-- ADDED 6/3/08
               Q(i,j,1)=Q(i,j,2)
               ER(i,j,1)=ER(i,j,2)
!-- zero derivative at the edges
               Q(i,j,0)=Q(i,j,1)
               Q(i,j,-1)=Q(i,j,1)
               ER(i,j,0)=ER(i,j,1)
               ER(i,j,-1)=ER(i,j,1)
            ENDDO
         ENDDO
      endif
      if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL) then
         DO J=-1,locNY+2
            DO I=-1,locNX+2
!-- ADDED 6/3/08
               Q(I,j,locNZ)=Q(I,j,locNZ-1)
               ER(I,j,locNZ)=ER(I,j,locNZ-1)
!-- zero derivative at the edges
               Q(i,j,locNZ+1)=Q(i,j,locNZ)
               Q(i,j,locNZ+2)=Q(i,j,locNZ)
               ER(i,j,locNZ+1)=ER(i,j,locNZ)
               ER(i,j,locNZ+2)=ER(i,j,locNZ)
            ENDDO
         ENDDO
      endif

!--- Calculate the largest change in Q, save as global_deltaQ
      call MPI_ALLREDUCE(locMAX_Q,global_deltaQ,1,
     %     MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(locMAX_ER,global_deltaER,1,
     %     MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

c      if(myid.eq.0) then
c         write(*,'(A,3(I8,1x),3(1x,e12.6))') 
c     %     'matrix test1_2',energy_iter,num_iter,sor_iter,global_deltaQ,
c     $     global_deltaER,delta_max
c         print *,'----------------------'
c         print *,''
c      endif

!--- Check if maximum iterations are exceeded
      IF (energy_ITER.EQ.ITMAX) THEN
         if(myid.eq.0) then
            print *,''
            print *,'Max iteration reached in coupled_energy_eqs'
            write(*,'(I8,1x,3(1x,e12.6))') 
     %           energy_iter,global_deltaQ,global_deltaER,delta_max

         endif
         write(*,'(A,I8,1x,e12.6,4(1x,I8),1x,e12.6)') 
     %        'ind. proc. values:',myid,max_val,max_iloc,max_jloc,
     %        max_kloc,upperbnd(max_jloc,max_kloc),
     %        Q(max_iloc,max_jloc,max_kloc)/
     %        (rh(max_iloc,max_jloc,max_kloc)*
     %        cv(max_iloc,max_jloc,max_kloc))

c- writeout a snapshot of the simulation
c         write(proc_num,'(i3.3)') myid
c         non_converge_file = 'NCVG_'//TRIM(ADJUSTL(proc_num))//'.005'
c         OPEN (12,FILE=non_converge_file,STATUS='UNKNOWN',
c     $        FORM='UNFORMATTED')
c         REWIND(12)
c         CALL PRINTOUT(1)
c         CLOSE(12)        
         
         goto 500
      end if
      
!--- Check Convergence on Q
      IF (global_deltaQ.GT.delta_max) GOTO 4
c      if(myid.eq.0) then
c         write(*,'(A,1x,I8,1x,3(1x,e12.6))') 
c     %        'CONVERGED PPENERGY:',energy_iter,global_deltaQ,
c     %        global_deltaER,delta_max
c      endif

 500  CONTINUE      
!--- PUT TOTAL CHANGE IN DELTA_ER
      delta_ER = ER - ER_OLD
!--- Fill T-array
      T = Q/(RH*CV)
      RETURN
      END SUBROUTINE COUPLED_ENERGY_EQS

      SUBROUTINE CALCULATEDIVF(DIVF_IN,theta_param)
!--- This routine is only checked for spherical coordiates
!--- Check geomety terms in gradient and divergence for ncosys.ne.2
      USE fluid_var_init
      USE mpi_var_init
      USE deltaER_var_init
      IMPLICIT NONE
      integer :: i,j,k
      double precision :: DIVF_IN(locNX,locNY,locNZ)
      double precision :: theta_param
!--- Calculate the FLD parameter difrx = (1/DX_b)*lambdax*c/(rho*(kappa+sig))
!--- This also calculates opacites through DIFF and FLIM calls
      CALL DIFF
!---  Fill THE EVAL arrays for the FLUX equation AND the \delta_ER Equation
      CALL FILL_EVALUES(theta_param)
      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX
               DIVF_IN(i,j,k) = -EVAL_1_6(i,j,k)*ER(i,j,k) + 
     $              EVAL1(i,j,k)*ER(i+1,j,k) + 
     $              EVAL2(i,j,k)*ER(i-1,j,k) + 
     $              EVAL3(i,j,k)*ER(i,j+1,k) + 
     $              EVAL4(i,j,k)*ER(i,j-1,k) + 
     $              EVAL5(i,j,k)*ER(i,j,k+1) + 
     $              EVAL6(i,j,k)*ER(i,j,k-1)
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE CALCULATEDIVF


      SUBROUTINE CALCULATEDIVF_WITHF(DIVF_IN,theta_param)
!--- This routine is only checked for spherical coordiates
      USE input_init
      USE fluid_var_init
      USE deltaER_var_init
      USE grid_var_init
      USE rad_var_init
      USE mpi_var_init
      IMPLICIT NONE
      integer i,j,k
      double precision :: DIVF_IN(locNX,locNY,locNZ)
      double precision :: theta_param
      SELECT CASE(NCOSYS)
      CASE (0)
         do k=1,locNZ
            do j=1,locNY
               do i=1,locNX
                  DIVF_IN(i,j,k) = 
     %                 ((FLUX_X(i+1,j,k)-FLUX_X(i,j,k))/dxa(i))+
     %                 ((FLUX_Y(i,j+1,k)-FLUX_Y(i,j,k))/dya(j))+
     %                 ((FLUX_Z(i,j,k+1)-FLUX_Z(i,j,k))/dza(k))
               enddo
            enddo
         enddo
      CASE (1)
         print *,'Check geomety terms in gradient and divergence'
         print *,'for ncosys=1:CALCULATEDIVF_WITH F'
         call clean_stop
      CASE (2)
         do k=1,locNZ
            do j=1,locNY
               do i=1,locNX
                  DIVF_IN(i,j,k) = 
     %                 ((surxa(i+1)/(surxb(i)*dxa(i)))*FLUX_X(i+1,j,k)-
     %                 (surxa(i)/(surxb(i)*dxa(i)))*FLUX_X(i,j,k)) +
     %                 ((1.d0/(xxb(i)*surzb(k)*dya(j)))*FLUX_Y(i,j+1,k)-
     %                 (1.d0/(xxb(i)*surzb(k)*dya(j)))*FLUX_Y(i,j,k)) +
     %                 ((surza(k+1)/(xxb(i)*surzb(k)*dza(k)))*
     %                 FLUX_Z(i,j,k+1)-
     %                (surza(k)/(xxb(i)*surzb(k)*dza(k)))*FLUX_Z(i,j,k))
               enddo
            enddo
         enddo
      CASE (3:)
         print *,'NCOSYS not available in CALCULATEDIVF_WITH F',NCOSYS
         call clean_stop
      CASE DEFAULT
         print *,'NCOSYS not available in CALCULATEDIVF_WITH F',NCOSYS
         call clean_stop
      END SELECT
!---  Fill THE EVAL arrays for the \delta_ER Equation
      CALL FILL_EVALUES(theta_param)
      RETURN
      END SUBROUTINE CALCULATEDIVF_WITHF

      SUBROUTINE CALCULATEFLUX
!--- This routine is only checked for spherical coordiates
!--- Check geomety terms in gradient and divergence for ncosys.ne.2
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE rad_var_init
      USE deltaER_var_init
      IMPLICIT NONE
      integer :: i,j,k
!--- Calculate the FLD parameter difrx = (1/DX_b)*lambdax*c/(rho*(kappa+sig))
!--- This also calculates opacites through DIFF and FLIM calls
      CALL DIFF
!--- flux in the x-direction
      if(irad.ne.5.and.irad.ne.6) then !-radial flux calculated in ppwavelth for 5 and 6
         do k=1,locNZ
            do j=1,locNY
               do i=1,locNX+1
                  FLUX_X(i,j,k) = -difrx(i,j,k)*(ER(i,j,k)-ER(i-1,j,k))
               enddo
            enddo
         enddo
      endif
!--- flux in the y-direction
      do k=1,locNZ
         do j=1,locNY+1
            do i=1,locNX
               FLUX_Y(i,j,k) = -difry(i,j,k)*(ER(i,j,k)-ER(i,j-1,k))
            enddo
         enddo
      enddo
!--- flux in the z-direction
      do k=1,locNZ+1
         do j=1,locNY
            do i=1,locNX
               FLUX_Z(i,j,k) = -difrz(i,j,k)*(ER(i,j,k)-ER(i,j,k-1))
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE CALCULATEFLUX

      SUBROUTINE FLUX_BC
!--- This routine is only checked for spherical coordiates
!--- Check geomety terms in gradient and divergence for ncosys.ne.2
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE deltaER_var_init
      USE global_constants
      USE grid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k 
      double precision :: CALC_TINCIDENT

!-required if using F = F*(r(ub-1)/r(ub))^2 + Fstar*dtau
c      CALL CALCSTELLARTAU
      do k=1,locNZ
         do j=1,locNY
            if(MPIright.eq.MPI_PROC_NULL) then
               do i=upperbnd(j,k),locNX+1
!-- F=cE
c                  FLUX_X(i,j,k) = FD*ER(upperbnd(j,k)-1,j,k)
!-- F=cE/2
                  FLUX_X(i,j,k) = FD*ER(upperbnd(j,k)-1,j,k)/2.d0
!- F = cE/2 - 2*Finc: MARSHAK BOUNDARY CONDITIONS
c                  FLUX_X(I,J,K)= (FD*ER(upperbnd(j,k)-1,j,k)/2.d0) -
c     $                 (2.d0*SBCONST*(CALC_TINCIDENT(j,k)**4.d0)/FD)
!-- dF/dr=0
c                  FLUX_X(i,j,k) = FLUX_X(upperbnd(j,k)-1,j,k)
!-- F*r^2_ub = F*r^2_(ub-1)... ie lum = constant
c                  FLUX_X(i,j,k) =((xxa(upperbnd(j,k)-1)/xxa(i))**2)
c     $                 *FLUX_X(upperbnd(j,k)-1,j,k)
!-- F = F*(r(ub-1)/r(ub))^2 + Fstar*dtau
c                  FLUX_X(i,j,k)=(((xxa(upperbnd(j,k)-1)**2)/(xxa(i)**2))
c     $                 *FLUX_X(upperbnd(j,k)-1,j,k) +
c     $                 SBCONST*(CALC_TINCIDENT(j,k)**4)*
c     $                 dtaus(upperbnd(j,k)-1,j,k))
!-- F=input stellar flux + interior flux
c                  FLUX_X(i,j,k)=SBCONST*(CALC_TINCIDENT(j,k)**4) + 
c     $                 FBottom
!--y and z fluxes
                  if(i.ne.locNX+1) then
                     FLUX_Y(i,j,k) = 0.d0
                     FLUX_Z(i,j,k) = 0.d0
                  endif
               enddo
            endif
            if(MPIleft.eq.MPI_PROC_NULL) then
c               FLUX_X(2,j,k) = Fbottom
c               FLUX_X(1,j,k) = Fbottom*((xxa(2)/xxa(1))**2.d0)
c               FLUX_Y(1,j,k) = 0.d0
c               FLUX_Z(1,j,k) = 0.d0
               FLUX_X(1,j,k) = Fbottom
c               FLUX_X(1,j,k) = 0.d0
            endif
         enddo
      enddo
      RETURN
      END SUBROUTINE FLUX_BC

      SUBROUTINE FILL_EVALUES(theta_param)
!--- This routine fills the the E-constants for the matrix inversion
!---  for solving for \delta_ER
      USE grid_var_init
      USE fluid_var_init
      USE rad_var_init
      USE deltaER_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      double precision :: theta_param
      integer :: i,j,k
      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX
               EVAL1(i,j,k) = -surxa(i+1)*difrx(i+1,j,k)/
     %              (surxb(i)*dxa(i))
               EVAL2(i,j,k) = -surxa(i)*difrx(i,j,k)/
     %              (surxb(i)*dxa(i))
               EVAL3(i,j,k) = -difry(i,j+1,k)/
     %              (xxb(i)*surzb(k)*dya(j))
               EVAL4(i,j,k) = -difry(i,j,k)/
     %              (xxb(i)*surzb(k)*dya(j))
               EVAL5(i,j,k) = -difrz(i,j,k+1)*surza(k+1)/
     %              (xxb(i)*surzb(k)*dza(k))
               EVAL6(i,j,k) = -difrz(i,j,k)*surza(k)/
     %              (xxb(i)*surzb(k)*dza(k))              
               EVAL7(i,j,k) = 1.d0 + (1.d0-Jacob_eta(i,j,k))*delt*
     %              rh(i,j,k)*xkapP(i,j,k)*FD
               EVAL8(i,j,k) = (Jacob_eta(i,j,k)*
     %              (1.d0-theta_param)-1.d0)*delt
               EVAL_1_6(i,j,k) = EVAL1(i,j,k)+EVAL2(i,j,k)+EVAL3(i,j,k)+
     %              EVAL4(i,j,k)+EVAL5(i,j,k)+EVAL6(i,j,k)
c               EVAL1(i,j,k) = 1.d0
c               EVAL2(i,j,k) = 1.d0
c               EVAL3(i,j,k) = 1.d0
c               EVAL4(i,j,k) = 1.d0
c               EVAL5(i,j,k) = 1.d0
c               EVAL6(i,j,k) = 1.d0
c               EVAL7(i,j,k) = 6.d0
c               EVAL8(i,j,k) = 1.d0
c               EVAL_1_6(i,j,k) = 0.d0
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE FILL_EVALUES



