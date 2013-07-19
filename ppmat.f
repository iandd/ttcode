      SUBROUTINE SOLVMAT
*     SOR-Method for solving the Matrixequations
*     SORPARAM is the overrelaxation parameter set in ipp
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      IMPLICIT NONE
      include "mpif.h"
      INTEGER :: I,J,K,ITER,IMIN
      INTEGER,PARAMETER :: ITMIN=5 
      DOUBLE PRECISION :: XS(locNX,locNY,locNZ)
      DOUBLE PRECISION :: OM,RHSNORM,RNORM
      DOUBLE PRECISION :: test1,test2,test3,test4,test5


      if(mod(III,iprstep).eq.0.or.III.eq.1) then
         call adjust_sor_param
      endif


**  Fixed Omega (OM used in iteration)
      OM=SORPARAM
**  RHS - Norm
      RHSNORM=0.d0
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               RHSNORM=RHSNORM+RHS(I,J,K)*RHS(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      RHSNORM=SQRT(RHSNORM)
!-- Find the global_RHSNORM
      call MPI_ALLREDUCE(RHSNORM,global_RHSNORM,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)


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
      elseif(modtyp.eq.6) then
         CALL SHIFTVAR_AXISYMMETRIC(6)
      else  !-- shifting routine for everything else
         IF(locNZ.eq.1) THEN
            CALL SHIFTVAR2D(7)
         ENDIF 
         IF(locNZ.ne.1) THEN
c            CALL SHIFTVAR3D(6)
            CALL NEWSHIFTVAR3D(6)
         ENDIF
      endif

!-------------------
!--   START OF THE SOLVE-MATRIX ITERATIONS
!-------------------
      ITER = 0
 4    CONTINUE
      ITER = ITER + 1

**  ODD - EVEN   *************
      DO K=1,locNZ
         DO J=1,locNY
            IMIN=2-MOD(J+K,2)
            DO I=IMIN,locNX,2
               XS(I,J,K) = B(I,J,K)*ER(I,J,K) -RHS(I,J,K)
     &              +AX(I,J,K)*ER(I-1,J,K)+CX(I,J,K)*ER(I+1,J,K)
     &              +AY(I,J,K)*ER(I,J-1,K)+CY(I,J,K)*ER(I,J+1,K)
     &              +AZ(I,J,K)*ER(I,J,K-1)+CZ(I,J,K)*ER(I,J,K+1)
               ER(I,J,K) = ER(I,J,K) - OM * XS(I,J,K)/B(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      
      if(irad.ne.5.and.irad.ne.6) CALL RADIATIVEBOUND
      

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
      elseif(modtyp.eq.6) then
         CALL SHIFTVAR_AXISYMMETRIC(6)
      else  !-- shifting routine for everything else
         IF(locNZ.eq.1) THEN
            CALL SHIFTVAR2D(7)
         ENDIF 
         IF(locNZ.ne.1) THEN
            CALL NEWSHIFTVAR3D(6)
c            CALL SHIFTVAR3D(6)
         ENDIF
      endif

      DO K=1,locNZ
         DO J=1,locNY
            IMIN=2-MOD(J+K+1,2)
            DO I=IMIN,locNX,2
               XS(I,J,K) = B(I,J,K)*ER(I,J,K) -RHS(I,J,K)
     &              +AX(I,J,K)*ER(I-1,J,K)+CX(I,J,K)*ER(I+1,J,K)
     &              +AY(I,J,K)*ER(I,J-1,K)+CY(I,J,K)*ER(I,J+1,K)
     &              +AZ(I,J,K)*ER(I,J,K-1)+CZ(I,J,K)*ER(I,J,K+1)
               ER(I,J,K) = ER(I,J,K) - OM * XS(I,J,K)/B(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      if(irad.ne.5.and.irad.ne.6) CALL RADIATIVEBOUND

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
      elseif(modtyp.eq.6) then
         CALL SHIFTVAR_AXISYMMETRIC(6)
      else  !-- shifting routine for everything else
         IF(locNZ.eq.1) THEN
            CALL SHIFTVAR2D(7)
         ENDIF 
         IF(locNZ.ne.1) THEN
            CALL NEWSHIFTVAR3D(6)
c            CALL SHIFTVAR3D(6)
         ENDIF
      endif

**  Residual Norm
      RNORM=0.d0
      DO K=1,locNZ
         DO J=1,locNY
            DO I=1,locNX
               RNORM=RNORM+XS(I,J,K)*XS(I,J,K) 
            ENDDO
         ENDDO
      ENDDO
!--- find the global_RNORM
      call MPI_ALLREDUCE(RNORM,global_RNORM,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      global_RNORM=SQRT(global_RNORM)

!--- Check if maximum iterations are exceeded
      IF (ITER.EQ.ITMAX) THEN
         if(myid.eq.0) then
            write(*,'(A,2(1x,I8),5(1x,e12.6))') 
     $           'Max iteration reached in SOR matrix solver',
     $           III,iter,global_RNORM/global_RHSNORM,
     $           global_RNORM,global_RHSNORM,EPSMAT,OM
         endif
c         write(*,'(A,1(I5),6(e19.8))') 'ind:',myid,
c     %        SQRT(RNORM)/sqrt(RHSNORM),
c     %        test1,test2,test3,test4,test5
         goto 500
      end if
      
      IF (global_RNORM/global_RHSNORM.GT.EPSMAT.OR.ITER.LT.ITMIN) GOTO 4
c      if(myid.eq.0) then
c         print *,'condition1',iter,global_RNORM/global_RHSNORM
c      endif

 500  CONTINUE

      RETURN
      END SUBROUTINE SOLVMAT

