
      SUBROUTINE GASCOOL
! Artifically cools the gas by relaxing the temperature back to a
!  pre-defined value defined by Trelax
      USE input_init
      USE fluid_var_init
      USE relax_var_init
      IMPLICIT NONE
      integer :: i,j,k
      double precision ::epsilon
      epsilon = 0.01
      DO k=1,locNZ
         DO j=1,locNY
            DO i=1,locNX
c               T(i,j,k) = Trelax(i,j,k) + (T(i,j,k)-Trelax(i,j,k))*
c     $              exp(-1.d0*III/2000.d0)
c               ER(i,j,k) = ERrelax(i,j,k) + (ER(i,j,k)-ERrelax(i,j,k))*
c     $              exp(-1.d0*III/2000.d0)
               T(i,j,k) = T(i,j,k) - epsilon*(T(i,j,k)-Trelax(i,j,k))*
     $              exp(-1.d0*III/2000.d0)
               ER(i,j,k) = ER(i,j,k)-epsilon*(ER(i,j,k)-ERrelax(i,j,k))*
     $              exp(-1.d0*III/2000.d0)

            enddo
         enddo
      enddo
      CALL BOUNDS
      RETURN
      END SUBROUTINE GASCOOL

      SUBROUTINE READRELAX
      USE input_init
      USE global_var_init
      USE fluid_var_init
      USE mpi_var_init
      USE relax_var_init
      IMPLICIT NONE
      integer :: i,j,k,loci,locj,lock
      integer :: dummy_global_NX,dummy_global_NY,dummy_global_NZ
      integer :: dummy_locNX,dummy_locNY,dummy_locNZ,dummy_numprocs
      double precision :: dummy_XXB(global_NX),dummy_XYB(global_NY)
      double precision :: dummy_XZB(global_NY)
      double precision :: dummy(global_NX,global_NY,global_NZ)
      double precision :: global_Trelax(global_NX,global_NY,global_NZ)
      double precision :: global_ERrelax(global_NX,global_NY,global_NZ)
!-- read in the T and ER from a file
      if(myid.eq.0) then
         print *,'READING RELAXATION FILE'
      endif
      OPEN (12,FILE='RELAX',STATUS='OLD',FORM='UNFORMATTED')
      rewind(12)
      READ (12) dummy_global_NX
      READ (12) dummy_global_NY
      READ (12) dummy_global_NZ
      READ (12) dummy_locNX,dummy_locNY,dummy_locNZ,dummy_numprocs
      READ (12) ZEIT,DELT,NUM_ITER
      READ (12) (dummy_XXB(I),I=1,global_NX)
      READ (12) (dummy_XYB(J),J=1,global_NY)
      READ (12) (dummy_XZB(J),K=1,global_NZ)
      READ (12) (((dummy(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      READ (12) (((dummy(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      READ (12) (((dummy(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      READ (12) (((dummy(I,J,K),I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      READ (12) (((global_Trelax(I,J,K) ,I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      READ (12) (((global_ERrelax(I,J,K) ,I=1,global_NX),J=1,global_NY),
     %     K=1,global_NZ)
      CLOSE(12)
!-- Initilize local copies of the variables
      do k=1,locNZ
         lock = k+locNZ*proc_coords(3)
         if(lock.gt.global_NZ.or.lock.lt.1) then
            print *,'ERROR IN INIT-k'
            stop
         endif
         do j=1,locNY
            locj = j+locNY*proc_coords(2)
            if(locj.gt.global_NY.or.locj.lt.1) then
               print *,'ERROR IN INIT-j'
               stop
            endif
            do i=1,locNX
               loci = i+locNX*proc_coords(1)
               if(loci.gt.global_NX.or.loci.lt.1) then
                  print *,'ERROR IN INIT-i'
                  stop
               endif
               Trelax(i,j,k)  = global_Trelax(loci,locj,lock)
               ERrelax(i,j,k) = global_ERrelax(loci,locj,lock)
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE READRELAX
