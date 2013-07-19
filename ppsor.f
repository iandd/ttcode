
      SUBROUTINE SOLVE_WITH_SOR(print_max_values)
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE deltaER_var_init
      USE sor_var_init
      USE rad_var_init

      USE grid_var_init

      IMPLICIT NONE
      include "mpif.h"
      integer i,j,k,IMIN
      INTEGER,PARAMETER :: ITMIN=1,num_indx=3
      double precision :: resid,RNORM,RHSNORM
      double precision :: max_resid
      integer :: max_loc(3),zmax_sor_index
      logical :: print_max_values
!--- SETUP PROBLEM WHEN SOLVING FOR delta_ER ... (ie irad=3)
      if(irad.eq.3) then
!--- Fill the matrix coefficients
         CALL FILL_SOR_MATRIX
!---  Fill the RHS vector
         CALL B_VECTOR_BOUNDARIES
!---  Initilize delta_ER to small value
c         delta_ER = 0.001*ER
         delta_ER = 0.d0
      endif

**  RHS - Norm
      RHSNORM=SUM(ABS(RHS_ER))
!-- Find the global_RHSNORM
      call MPI_ALLREDUCE(RHSNORM,global_RHSNORM,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

!-- Set the maximum z-index for the SOR solver      
      if(poles.and.proc_coords(3).eq.proc_dims(3)-1) then
         zmax_sor_index = locNZ-1
      else
         zmax_sor_index = locNZ
      endif


      SOR_ITER = 0
 4    CONTINUE
      SOR_ITER = SOR_ITER + 1

c      if(myid.eq.0) print *,'soriter=',sor_ITER

!-- Zero the convergence parameter each full iteration
      rnorm = 0.d0

!-- Zero the max resid for the final integration debug
      max_resid = 0.0

**  ODD - EVEN   *************
      DO K=1,zmax_sor_index
         DO J=1,locNY
            IMIN=2-MOD(J+K,2)
            DO I=IMIN,locNX,2
c               DO I=IMIN,upperbnd(j,k)-1,2
               resid = SOR_A(i,j,k)*delta_ER(i,j,k) + 
     1              SOR_B(i,j,k)*delta_ER(i-1,j,k) + 
     2              SOR_C(i,j,k)*delta_ER(i+1,j,k) + 
     3              SOR_D(i,j,k)*delta_ER(i,j-1,k) + 
     4              SOR_E(i,j,k)*delta_ER(i,j+1,k) + 
     5              SOR_F(i,j,k)*delta_ER(i,j,k-1) + 
     6              SOR_G(i,j,k)*delta_ER(i,j,k+1) - 
     7              RHS_ER(i,j,k)

               RNORM=RNORM+abs(resid)
               delta_ER(I,J,K) = delta_ER(I,J,K) - 
     %              SORPARAM*resid/SOR_A(I,J,K)

               
               IF (SOR_ITER.EQ.ITMAX) THEN
                  if(abs(resid).gt.max_resid) then
                     max_resid = abs(resid)
                     max_loc(1) = i
                     max_loc(2) = j
                     max_loc(3) = k
                  endif
               endif
            ENDDO
         ENDDO
      ENDDO
!--- Call delta_ER boundary passing routine
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
      
      DO K=1,zmax_sor_index
         DO J=1,locNY
            IMIN=2-MOD(J+K+1,2)
c            DO I=IMIN,upperbnd(j,k)-1,2
            DO I=IMIN,locNX,2
               resid = SOR_A(i,j,k)*delta_ER(i,j,k) + 
     1              SOR_B(i,j,k)*delta_ER(i-1,j,k) + 
     2              SOR_C(i,j,k)*delta_ER(i+1,j,k) + 
     3              SOR_D(i,j,k)*delta_ER(i,j-1,k) + 
     4              SOR_E(i,j,k)*delta_ER(i,j+1,k) + 
     5              SOR_F(i,j,k)*delta_ER(i,j,k-1) + 
     6              SOR_G(i,j,k)*delta_ER(i,j,k+1) - 
     7              RHS_ER(i,j,k)

               RNORM=RNORM+abs(resid)
               delta_ER(I,J,K) = delta_ER(I,J,K) - 
     %              SORPARAM*resid/SOR_A(I,J,K)

               
               IF (SOR_ITER.EQ.ITMAX) THEN
                  if(abs(resid).gt.max_resid) then
                     max_resid = abs(resid)
                     max_loc(1) = i
                     max_loc(2) = j
                     max_loc(3) = k
                  endif
               endif
            ENDDO
         ENDDO
      ENDDO
!---  Call delta_ER boundary passing routine
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

!--- find the global_RNORM
      call MPI_ALLREDUCE(RNORM,global_RNORM,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

!--- Check if maximum iterations are exceeded
!- print, except for when called from adjust_sor
      IF (SOR_ITER.EQ.ITMAX) THEN
         if(myid.eq.0.and.print_max_values) then
            write(*,'(A,2(1x,I8),4(1x,e12.6))') 
     $           'MAX COUPLED SOR ITERATION',
     $           num_iter,energy_iter,global_RNORM/global_RHSNORM,
     %           global_RNORM,global_RHSNORM,EPSMAT
         endif
         if(print_max_values) then
            write(*,'(A,I8,1x,e12.6,1x,e12.6,1x,I8,1x,I8,1x,I8,1x,I8)') 
     %           '     ',myid,RNORM/RHSNORM,max_resid,max_loc(1),
     %           max_loc(2),max_loc(3),upperbnd(max_loc(2),max_loc(3))
         endif
         goto 500
      endif
      
      IF (global_RNORM/global_RHSNORM.GT.EPSMAT.OR.
     $     SOR_ITER.LT.ITMIN) GOTO 4

c      if(myid.eq.0) then
c         write(*,'(A,2(I8,1x),e12.6)') 
c     %        'condition1, SOR',SOR_iter,energy_iter,
c     %        global_RNORM/global_RHSNORM
c      endif

 500  CONTINUE
      RETURN
      END SUBROUTINE SOLVE_WITH_SOR
      
      SUBROUTINE FILL_SOR_MATRIX
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE deltaER_var_init
      USE sor_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer i,j,k
      do k = 1,locNZ
         do j = 1,locNY
            do i = 1,locNX
               SOR_A(i,j,k)   = EVAL7(i,j,k)+EVAL8(i,j,k)*
     $              EVAL_1_6(i,j,k)
               SOR_B(i,j,k) = -EVAL8(i,j,k)*EVAL2(i,j,k)
               SOR_C(i,j,k) = -EVAL8(i,j,k)*EVAL1(i,j,k)
               SOR_D(i,j,k) = -EVAL8(i,j,k)*EVAL4(i,j,k)
               SOR_E(i,j,k) = -EVAL8(i,j,k)*EVAL3(i,j,k)
               SOR_F(i,j,k) = -EVAL8(i,j,k)*EVAL6(i,j,k)
               SOR_G(i,j,k) = -EVAL8(i,j,k)*EVAL5(i,j,k)
!-- values beyond upperbnd
               if(MPIright.eq.MPI_PROC_NULL.and.i.ge.upperbnd(j,k)) then !--------right boundary
                  SOR_A(i,j,k) = 1.d0  
                  SOR_B(i,j,k) = 0.d0
                  SOR_C(i,j,k) = 0.d0
                  SOR_D(i,j,k) = 0.d0
                  SOR_E(i,j,k) = 0.d0
                  SOR_F(i,j,k) = 0.d0
                  SOR_G(i,j,k) = 0.d0
               endif
!-- neighbor values beyond upperbnd
               if(MPIright.eq.MPI_PROC_NULL.and.
     %              i.eq.upperbnd(j,k)-1) then !--------right-neighbor boundary
                  SOR_C(i,j,k) = 0.d0
               endif

!---- BOUNDARIES
c               if( (MPIleft.eq.MPI_PROC_NULL.and.i.eq.1) .or. !---------left boundary
               if(
     1              (MPIright.eq.MPI_PROC_NULL.and.i.eq.locNX) .or. !----right boundary
     2              (MPIlower.eq.MPI_PROC_NULL.and.j.eq.1) .or.!--------lower boundary
     3              (MPIupper.eq.MPI_PROC_NULL.and.j.eq.locNY) .or. !---upper boundary
     4              (.not.poles.and.MPIbelow.eq.MPI_PROC_NULL.and.
     5              k.eq.1) .or. !-------below boundary
     6              (.not.poles.and.MPIabove.eq.MPI_PROC_NULL.and.
     7              k.eq.locNZ) ) then !-above boundary
                  SOR_A(i,j,k) = 1.d0  
                  SOR_B(i,j,k) = 0.d0
                  SOR_C(i,j,k) = 0.d0
                  SOR_D(i,j,k) = 0.d0
                  SOR_E(i,j,k) = 0.d0
                  SOR_F(i,j,k) = 0.d0
                  SOR_G(i,j,k) = 0.d0
               endif
!---- NEIGHBOR BOUNDARY VALUES
!--------left-neighbor boundary
c               if(MPIleft.eq.MPI_PROC_NULL.and.i.eq.2) then
c                  SOR_B(i,j,k) = 0.d0
c               endif
!--------right-neighbor boundary
               if(MPIright.eq.MPI_PROC_NULL.and.i.eq.locNX-1) then
                  SOR_C(i,j,k) = 0.d0
               endif
!--------lower-neighbor boundary
               if(MPIlower.eq.MPI_PROC_NULL.and.j.eq.2) then
                  SOR_D(i,j,k) = 0.d0
                  print *,'IN BAD SPOT NEIGHBOR'
               endif
!--------upper-neighbor boundary
               if(MPIupper.eq.MPI_PROC_NULL.and.j.eq.locNY-1) then
                  SOR_E(i,j,k) = 0.d0
                  print *,'IN BAD SPOT NEIGHBOR'
               endif
!--------below-neighbor boundary
               if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL.and.
     $              k.eq.2) then
                  SOR_F(i,j,k) = 0.d0
               endif
!--------above-neighbor boundary
               if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL.and.
     $              k.eq.locNZ-1) then
                  SOR_G(i,j,k) = 0.d0
               endif
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE FILL_SOR_MATRIX


      SUBROUTINE B_VECTOR_BOUNDARIES
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE deltaER_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer i,j,k
!--- FILL THE VECTORS
      do k = 1,locNZ
         do j = 1,locNY
            do i = 1,locNX
!-- values beyond upperbnd
               if(MPIright.eq.MPI_PROC_NULL.and.i.ge.upperbnd(j,k)) then !--------right boundary
                  RHS_ER(i,j,k) = 0.d0
               endif
!--------left boundary
c               if(MPIleft.eq.MPI_PROC_NULL.and.i.eq.1) then
c                  RHS_ER(i,j,k) = 0.d0
c               endif
!--------right boundary
               if(MPIright.eq.MPI_PROC_NULL.and.i.eq.locNX) then
                  RHS_ER(i,j,k) = 0.d0
               endif
!--------lower boundary
              if(MPIlower.eq.MPI_PROC_NULL.and.j.eq.1) then
                  RHS_ER(i,j,k) = 0.d0
                  print *,'IN BAD SPOT'
               endif
!--------upper boundary
               if(MPIupper.eq.MPI_PROC_NULL.and.j.eq.locNY) then
                  RHS_ER(i,j,k) = 0.d0
                  print *,'IN BAD SPOT'
               endif
!--------below boundary
               if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL.and.
     $              k.eq.1) then
                  RHS_ER(i,j,k) = 0.d0
               endif
!--------above boundary
               if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL.and.
     $              k.eq.locNZ) then
                  RHS_ER(i,j,k) = 0.d0
               endif
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE B_VECTOR_BOUNDARIES


      SUBROUTINE ADJUST_SOR_PARAM
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE sor_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,save_soriter(50)
      double precision :: save_sorparam(50)
      SORPARAM = 1.01
      do i=1,50
         CALL SOLVE_WITH_SOR(.false.)
         save_soriter(i) = SOR_ITER
         save_sorparam(i) = SORPARAM
         SORPARAM=SORPARAM+0.02
      enddo
      SORPARAM=save_sorparam(minloc(save_soriter,1))
      if(myid.eq.0) then !.and.mod(III,iprstep)) then
         write(*,'(A,I8,1x,e12.6,1x,I8)') 'SOR PARAM=',III,SORPARAM,
     $        save_soriter(minloc(save_soriter,1))
      endif
      RETURN
      END SUBROUTINE ADJUST_SOR_PARAM
