      SUBROUTINE DIFFUSION(scalar_num)
**     **Implicit** Diffusion (SOR)
      USE dif_var_init
      USE scalar_var_init
      USE fluid_var_init
      USE global_constants
      USE input_init
      USE sor_var_init
      USE mpi_var_init
      USE grid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,K,IMIN,scalar_num
      INTEGER,PARAMETER :: ITMIN=5,num_indx=3
      INTEGER,PARAMETER :: SORITMAX=500 
      double precision :: resid,RNORM,RHSNORM
      double precision :: max_resid
      integer :: max_loc(3)
      double precision :: global_RNORM,global_RHSNORM
!-- Set pointer to the scalar
      DIFF_PTR => pass_sc

!-- Set rhs=old value
      do k = 1,locNZ
         do j = 1,locNY
            do i = 1,locNX
               RHS_DIF(i,j,k) = DIFF_PTR(i,j,k)
            enddo
         enddo
      enddo

!-- Calculate diffusion coeff. and SOR matrix elements
      CALL DIFFUSION_MAT_ELEM !- set diffusion boundary conditions here..

**  RHS - Norm
      RHSNORM=SUM(ABS(RHS_DIF))
!-- Find the global_RHSNORM
      call MPI_ALLREDUCE(RHSNORM,global_RHSNORM,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

!-- Zero the max resid for the final integration debug
      max_resid = 0.0

      SOR_ITER = 0
 4    CONTINUE
      SOR_ITER = SOR_ITER + 1

!-- Zero the convergence parameter each full iteration
      rnorm = 0.d0

**  ODD - EVEN   *************
      DO K=1,locNZ
         DO J=1,locNY
            IMIN=2-MOD(J+K,2)
               DO I=IMIN,locNX,2 
                  resid = SOR_A(i,j,k)*DIFF_PTR(i,j,k) + 
     1                 SOR_B(i,j,k)*DIFF_PTR(i-1,j,k) + 
     2                 SOR_C(i,j,k)*DIFF_PTR(i+1,j,k) + 
     3                 SOR_D(i,j,k)*DIFF_PTR(i,j-1,k) + 
     4                 SOR_E(i,j,k)*DIFF_PTR(i,j+1,k) + 
     5                 SOR_F(i,j,k)*DIFF_PTR(i,j,k-1) + 
     6                 SOR_G(i,j,k)*DIFF_PTR(i,j,k+1) - 
     7                 RHS_DIF(i,j,k)
                  RNORM=RNORM+abs(resid)
                  DIFF_PTR(I,J,K) = DIFF_PTR(I,J,K) - 
     %                 SORPARAM*resid/SOR_A(I,J,K)
                  IF (SOR_ITER.EQ.SORITMAX) THEN
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
!--- Call diffusion boundary passing routine
      if(modtyp.eq.2) then
         CALL PLANETSHIFT(10,-1)
!--   now correct the polar shifts if necessary
         if(poles) then
            if(proc_coords(3).eq.proc_dims(3)-1) then
               CALL NORTH_POLE_SHIFT(10,-1)
            endif
            if(proc_coords(3).eq.0) then
               CALL SOUTH_POLE_SHIFT(10,-1)
            endif
         endif
      else
         CALL NEWSHIFTVAR3D(10)
      endif

      DO K=1,locNZ
         DO J=1,locNY
            IMIN=2-MOD(J+K+1,2)
            DO I=IMIN,locNX,2
                  resid = SOR_A(i,j,k)*DIFF_PTR(i,j,k) + 
     1                 SOR_B(i,j,k)*DIFF_PTR(i-1,j,k) + 
     2                 SOR_C(i,j,k)*DIFF_PTR(i+1,j,k) + 
     3                 SOR_D(i,j,k)*DIFF_PTR(i,j-1,k) + 
     4                 SOR_E(i,j,k)*DIFF_PTR(i,j+1,k) + 
     5                 SOR_F(i,j,k)*DIFF_PTR(i,j,k-1) + 
     6                 SOR_G(i,j,k)*DIFF_PTR(i,j,k+1) - 
     7                 RHS_DIF(i,j,k)
                  RNORM=RNORM+abs(resid)
                  DIFF_PTR(I,J,K) = DIFF_PTR(I,J,K) - 
     %                 SORPARAM*resid/SOR_A(I,J,K)
                  
                  IF (SOR_ITER.EQ.SORITMAX) THEN
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
!--- Call diffusion boundary passing routine
      if(modtyp.eq.2) then
         CALL PLANETSHIFT(10,-1)
!--   now correct the polar shifts if necessary
         if(poles) then
            if(proc_coords(3).eq.proc_dims(3)-1) then
               CALL NORTH_POLE_SHIFT(10,-1)
            endif
            if(proc_coords(3).eq.0) then
               CALL SOUTH_POLE_SHIFT(10,-1)
            endif
         endif
      else
         CALL NEWSHIFTVAR3D(10)
      endif

!--- find the global_RNORM
      call MPI_ALLREDUCE(RNORM,global_RNORM,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

!--- Check if maximum iterations are exceeded
      IF (SOR_ITER.EQ.SORITMAX) THEN
         if(myid.eq.0) then
            print *,''
            print *,'Max iteration reached: DIF SOR matrix solver'
            print *,sor_iter,num_iter,global_RNORM/global_RHSNORM,
     %           global_RNORM,global_RHSNORM,EPSMAT
         endif
         write(*,'(A,I8,1x,e12.6,1x,e12.6,1x,I8,1x,I8,1x,I8,1x,I8)') 
     %        '     ',myid,RNORM,max_resid,max_loc(1),
     %        max_loc(2),max_loc(3),upperbnd(max_loc(2),max_loc(3))
         goto 500
      end if

      IF (global_RNORM/global_RHSNORM.GT.EPSMAT.OR.
     $     SOR_ITER.LT.ITMIN) GOTO 4

c      if(myid.eq.0) then
c         write(*,'(A,I8,2(1x,e12.6))') 
c     %        'converged, DIFF sor',SOR_iter,
c     %        global_RNORM/global_RHSNORM,EPSMAT
c      endif

 500  CONTINUE
      NULLIFY(DIFF_PTR)
      RETURN
      END SUBROUTINE DIFFUSION

      SUBROUTINE DIFFUSION_MAT_ELEM
**    Sets Matrixelements for Diffusion calculation
      USE input_init
      USE grid_var_init
      USE dif_var_init
      USE fluid_var_init
      USE sor_var_init
      USE mpi_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      INTEGER :: I,J,K
**  Calculate diffusion coefficients 
      CALL DIFFUSION_COEF
**  FILL THE DVALUES
      CALL FILL_DVALUES
**  NOW FILL THE SOR CONSTANTS
      SOR_A = 1.d0 + DVAL_1_6
      SOR_B = -DVAL2
      SOR_C = -DVAL1
      SOR_D = -DVAL4
      SOR_E = -DVAL3
      SOR_F = -DVAL6
      SOR_G = -DVAL5


**now set boundaries

!----------XMIN------------------------------------------------------------
!-- Fixed temperature derivative at the bottom (i=1). Corresponds to BOUNDXMIN(0)
c      if(MPIleft.eq.MPI_PROC_NULL) then !---------left boundary
c         SOR_A(1,1:locNY,1:locNZ) = 1.d0 + DVAL1(1,1:locNY,1:locNZ)
c         SOR_B(1,1:locNY,1:locNZ) = 0.d0
c         SOR_D(1,1:locNY,1:locNZ) = 0.d0
c         SOR_E(1,1:locNY,1:locNZ) = 0.d0
c         SOR_F(1,1:locNY,1:locNZ) = 0.d0
c         SOR_G(1,1:locNY,1:locNZ) = 0.d0
c         RHS_DIF(1,1:locNY,1:locNZ) = RHS_DIF(1,1:locNY,1:locNZ) + 
c     $        D_T*(xxb(0)-xxb(1))*DVAL2(1,1:locNY,1:locNZ)
c      endif

!----------XMAX------------------------------------------------------------
!-- Zero temperature gradient at the top (i=locNX) Corresponds to BOUNDXMAX(8)
      if(MPIright.eq.MPI_PROC_NULL) then 
         SOR_A(locNX,1:locNY,1:locNZ)=1.d0+DVAL2(locNX,1:locNY,1:locNZ)
         SOR_B(locNX,1:locNY,1:locNZ)= -DVAL2(locNX,1:locNY,1:locNZ)
         SOR_C(locNX,1:locNY,1:locNZ)= 0.d0
         SOR_D(locNX,1:locNY,1:locNZ)= 0.d0
         SOR_E(locNX,1:locNY,1:locNZ)= 0.d0
         SOR_F(locNX,1:locNY,1:locNZ)= 0.d0
         SOR_G(locNX,1:locNY,1:locNZ)= 0.d0
      endif

      RETURN
      END SUBROUTINE DIFFUSION_MAT_ELEM


      SUBROUTINE DIFFUSION_COEF
**    calculate radiative diffusion coefficients
      USE input_init
      USE fluid_var_init
      USE dif_var_init
      USE grid_var_init
      USE global_constants
      USE mpi_var_init
      USE global_constants
      IMPLICIT NONE
      include "mpif.h"
      INTEGER :: I,J,k
      double precision :: diffusion_const
!-- set diffusion constant
c      diffusion_const = (DXA(1)**(2.0))/(0.01*DAY)
      diffusion_const = xnue/Pr_num
      if(myid.eq.0.and.III.eq.1) then
          print *,'Prandtl Number=',Pr_num
         print *,'diffusion constant=',diffusion_const
         print *,'Fbot_T=-k*D_T =',-diffusion_const*D_T
         print *,'Fbot_epsilon=-rh*cv*k*D_T =',-rh(1,1,1)*cv(1,1,1)*
     $        diffusion_const*D_T
      endif
** difrx_DIF = DIFF_COEFF*(1/DX_b) - defined on iA
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX+1
               difrx_DIF(I,J,K)= diffusion_const/DXB(I)
            ENDDO
         ENDDO
      ENDDO
** difry = diff_coeff*(1/DY*r*cos(theta))- defined on jA
      DO K=1,locNZ
         DO J=1,locNY+1
            DO I=1,locNX
               difry_DIF(I,J,K)= diffusion_const/ 
     &              (DYB(J)*geoxg(i)*geozg(k))
            ENDDO
         ENDDO
      ENDDO
** difrz = diff_coeff*(1/DZ*r) - defined on kA
      DO K=1,locNZ+1
         DO J=1,locNY
            DO I=1,locNX
               difrz_DIF(I,J,K)= diffusion_const/(DZB(K)*geoxh(i))
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DIFFUSION_COEF
 

      SUBROUTINE FILL_DVALUES
!--- This routine fills the the D-constants for the matrix inversion
!---  for solving diffusion equation (includes dt)
      USE grid_var_init
      USE fluid_var_init
      USE dif_var_init
      IMPLICIT NONE
      double precision :: theta_param
      integer :: i,j,k
      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX 
               DVAL1(i,j,k) = delt*surxa(i+1)*difrx_DIF(i+1,j,k)/
     %              (surxb(i)*dxa(i))
               DVAL2(i,j,k) = delt*surxa(i)*difrx_DIF(i,j,k)/
     %              (surxb(i)*dxa(i))
               DVAL3(i,j,k) = delt*difry_DIF(i,j+1,k)/
     %              (geoxg(i)*surzb(k)*dya(j))
               DVAL4(i,j,k) = delt*difry_DIF(i,j,k)/
     %              (geoxg(i)*surzb(k)*dya(j))
               DVAL5(i,j,k) = delt*difrz_DIF(i,j,k+1)*surza(k+1)/
     %              (geoxg(i)*surzb(k)*dza(k))
               DVAL6(i,j,k) = delt*difrz_DIF(i,j,k)*surza(k)/
     %              (geoxg(i)*surzb(k)*dza(k))              
               DVAL_1_6(i,j,k) = DVAL1(i,j,k)+DVAL2(i,j,k)+DVAL3(i,j,k)+
     %              DVAL4(i,j,k)+DVAL5(i,j,k)+DVAL6(i,j,k)
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE FILL_DVALUES

