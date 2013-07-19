      SUBROUTINE SOLVE_WITH_HYPRE
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      USE deltaER_var_init
      USE mpi_var_init
      USE global_constants
      USE hypre_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer i,j,k,indx
!----  CREATE THE HYPRE MATRIX
      CALL CREATEHYPREMATRIX
!----  FILL THE MATRIX (AND SET BOUNDARIES)
      CALL FILLHYPREMATRIX
!----  ASSEMBLE THE MATRIX
      call bHYPRE_SStructMatrix_Assemble_f( H_matrix, ierr, except)
!----  CREATE THE X AND B HYPRE VECTORS
      CALL CREATEHYPREVECTORS
!----  FILL THE VECTORS (AND SET BOUNDARIES)
      CALL FILLHYPREVECTORS
!---  ASSEMBLE VECTORS
      call bHYPRE_SStructVector_Assemble_f( b_vec, ierr, except )
      call bHYPRE_SStructVector_Assemble_f( x_vec, ierr, except )

!--- Cast objects from SSTRUC class to STRUC class so we can use
!---  the STRUC solvers
      call bHYPRE_SStructMatrix_GetObject_f( H_matrix, dummy, ierr,
     %     except )
      call bHYPRE_StructMatrix__cast_f( dummy, sA, except )
      call sidl_BaseInterface_deleteRef_f( dummy, except )
      call bHYPRE_SStructVector_GetObject_f(b_vec, dummy,ierr,except)
      call bHYPRE_Vector__cast_f( dummy, vb, except )
      call sidl_BaseInterface_deleteRef_f( dummy, except )
      call bHYPRE_SStructVector_GetObject_f(x_vec, dummy,ierr,except)
      call bHYPRE_Vector__cast_f( dummy, vx, except )
      call sidl_BaseInterface_deleteRef_f( dummy, except )

!--- Cast and create an empty PCG Struct solver 
      call bHYPRE_Operator__cast_f( sA, opA, except )
      call bHYPRE_PCG_Create_f( bHYPRE_mpicomm, opA, PCGsolver,
     1     except )

!--- Set PCG parameters
      tol = 1.0e-6
      call bHYPRE_PCG_SetDoubleParameter_f( PCGsolver, "Tolerance",
     1     tol, ierr, except )
!-- 2: prints info on hypre iterations
c      call bHYPRE_PCG_SetIntParameter_f( PCGsolver, "PrintLevel",
c     1     2, ierr, except )
      call bHYPRE_PCG_SetIntParameter_f( PCGsolver, "PrintLevel",
     1     1, ierr, except )
      call bHYPRE_PCG_SetIntParameter_f( PCGsolver, "MaxIter",
     1     50, ierr, except )
!---  Create the Struct SMG solver for use as a preconditioner
      call bHYPRE_StructSMG_Create_f( bHYPRE_mpicomm, sA, SMGprecond,
     1     except )

!--- Set SMG parameters
      call bHYPRE_StructSMG_SetIntParameter_f( SMGprecond,
     1     "MaxIter", 1, ierr, except )
      tol = 0.0
      call bHYPRE_StructSMG_SetDoubleParameter_f( SMGprecond,
     1     "Tolerance", tol, ierr, except )
      call bHYPRE_StructSMG_SetIntParameter_f( SMGprecond,
     1     "ZeroGuess", 1, ierr, except )
      call bHYPRE_StructSMG_SetIntParameter_f( SMGprecond,
     1     "NumPreRelax", 1, ierr, except )
      call bHYPRE_StructSMG_SetIntParameter_f( SMGprecond,
     1     "NumPostRelax", 1, ierr, except )

!---  Cast preconditioner
      call bHYPRE_Solver__cast_f( SMGprecond, precond, except )
      call bHYPRE_PCG_SetPreconditioner_f( PCGsolver, precond,
     1     ierr, except )

!---  Setup and Solve
      call bHYPRE_PCG_Setup_f( PCGsolver, vb, vx, ierr, except )
      call bHYPRE_PCG_Apply_f( PCGsolver, vb, vx, ierr, except )

!---  Extract solution (note use SSTRUC pointer, not the cast pointer)
      nvalues = locNX*locNY*locNZ
      call bHYPRE_SStructVector_GetBoxValues_f(x_vec, hypre_part,
     1     ilower, iupper, hypre_ndim, hypre_var, answer, nvalues,
     2     ierr, except)
      indx = 1
      do k = 1,locNZ
         do j = 1,locNY
            do i = 1,locNX
               delta_ER(i,j,k) = answer(indx)
               indx = indx + 1
            enddo
         enddo
      enddo

!----- Free Preconditioner and Solver Memory
      call bHYPRE_Operator_deleteRef_f( opA, except )
      call bHYPRE_Vector_deleteRef_f( vx, except )
      call bHYPRE_Vector_deleteRef_f( vb, except )
      call bHYPRE_StructMatrix_deleteRef_f( sA, except )
      call bHYPRE_Solver_deleteRef_f( precond, except )
      call bHYPRE_StructSMG_deleteRef_f( SMGprecond, except )
      call bHYPRE_PCG_deleteRef_f( PCGsolver, except )
!----- Free Vector and Matrix Memory
      call bHYPRE_SStructVector_deleteRef_f( x_vec, except )
      call bHYPRE_SStructVector_deleteRef_f( b_vec, except )
      call bHYPRE_SStructMatrix_deleteRef_f( H_matrix, except )

      RETURN
      END SUBROUTINE SOLVE_WITH_HYPRE

      SUBROUTINE FILLHYPREVECTORS
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE hypre_var_init
      USE deltaER_var_init
      IMPLICIT NONE
      include "mpif.h"
      include 'bHYPRE_SStructVariable.inc'
      integer indx,i,j,k
!--- Initilize the hypre vectors
      call bHYPRE_SStructVector_Initialize_f( b_vec, ierr, except )
      call bHYPRE_SStructVector_Initialize_f( x_vec, ierr, except )

!--- FILL THE VECTORS
      nvalues=locNX*locNY*locNZ
      indx = 1
      do k = 1,locNZ
         do j = 1,locNY
            do i = 1,locNX
               x_values(indx)   = 1.d0
               b_values(indx)   = RHS_ER(i,j,k)
c--- LAPLACIAN TEST CASE
c               b_values(indx)   = 1.d0
!-- values beyond upperbnd
               if(MPIright.eq.MPI_PROC_NULL.and.i.ge.upperbnd(j,k)) then !--------right boundary
c               if(MPIright.eq.MPI_PROC_NULL.and.i.ge.
c     $              upperbnd(j,k)-10) then !--------right boundary
                  b_values(indx)   = 0.d0
               endif
!-- neighbor values beyond upperbnd
c               if(MPIright.eq.MPI_PROC_NULL.and.
c     %              i.eq.upperbnd(j,k)-1) then !--------right-neighbor boundary
c                  b_values(indx)=b_values(indx)+0.d0
c               endif

!--------left boundary
               if(MPIleft.eq.MPI_PROC_NULL.and.i.eq.1) then
                  b_values(indx)=0.d0
               endif
!--------left-neighbor boundary
c               if(MPIleft.eq.MPI_PROC_NULL.and.i.eq.2) then
c                  b_values(indx)=b_values(indx)+0.d0
c               endif
!--------right boundary
               if(MPIright.eq.MPI_PROC_NULL.and.i.eq.locNX) then
                  b_values(indx)=0.d0
               endif
!--------right-neighbor boundary
c               if(MPIright.eq.MPI_PROC_NULL.and.i.eq.locNX-1) then
c                  b_values(indx)=b_values(indx)+0.d0
c               endif
!--------lower boundary
               if(MPIlower.eq.MPI_PROC_NULL.and.j.eq.1) then
                  b_values(indx)=0.d0
                  print *,'IN BAD SPOT'
               endif
!--------lower-neighbor boundary
c               if(MPIlower.eq.MPI_PROC_NULL.and.j.eq.2) then
c                  b_values(indx)=b_values(indx)+0.d0
c                  print *,'IN BAD SPOT'
c               endif
!--------upper boundary
               if(MPIupper.eq.MPI_PROC_NULL.and.j.eq.locNY) then
                  b_values(indx)=0.d0
                  print *,'IN BAD SPOT'
               endif
!--------upper-neighbor boundary
c               if(MPIupper.eq.MPI_PROC_NULL.and.j.eq.locNY-1) then
c                  b_values(indx)=b_values(indx)+0.d0
c                  print *,'IN BAD SPOT'
c               endif
!--------below boundary
               if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL.and.
     $              k.eq.1) then
                  b_values(indx)=0.d0
               endif
!--------below-neighbor boundary
c               if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL.and.
c     $              k.eq.2) then
c                  b_values(indx)=b_values(indx)+0.d0
c               endif
!--------above boundary
               if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL.and.
     $              k.eq.locNZ) then
                  b_values(indx)=0.d0
               endif
!--------above-neighbor boundary
c               if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL.and.
c     $              k.eq.locNZ-1) then
c                  b_values(indx)=b_values(indx)+0.d0
c               endif
!-- increment indx
               indx = indx + 1
            enddo
         enddo
      enddo
!--- fill the hypre vectors
      call bHYPRE_SStructVector_SetBoxValues_f( b_vec, hypre_part,
     1     ilower, iupper, hypre_ndim, hypre_var, b_values, nvalues,
     2     ierr, except )
      call bHYPRE_SStructVector_SetBoxValues_f( x_vec, hypre_part,
     1     ilower, iupper, hypre_ndim, hypre_var, x_values, nvalues,
     2     ierr, except )
      RETURN
      END SUBROUTINE FILLHYPREVECTORS

      SUBROUTINE FILLHYPREMATRIX
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE hypre_var_init
      USE deltaER_var_init
      IMPLICIT NONE
      include "mpif.h"
      include 'bHYPRE_SStructVariable.inc'
      integer indx,i,j,k
!--- Initilize the hypre vectors
      call bHYPRE_SStructMatrix_Initialize_f( H_matrix, ierr, except )
      nentries=7
!--- Assign values to the matrix but don't Assemble until boundary 
!      conditions are set
      nvalues=locNX*locNY*locNZ*nentries
      indx = 1
      do k = 1,locNZ
         do j = 1,locNY
            do i = 1,locNX
               H_values(indx)   = EVAL7(i,j,k)+EVAL8(i,j,k)*
     $              EVAL_1_6(i,j,k)
               H_values(indx+1) = -EVAL8(i,j,k)*EVAL2(i,j,k)
               H_values(indx+2) = -EVAL8(i,j,k)*EVAL1(i,j,k)
               H_values(indx+3) = -EVAL8(i,j,k)*EVAL4(i,j,k)
               H_values(indx+4) = -EVAL8(i,j,k)*EVAL3(i,j,k)
               H_values(indx+5) = -EVAL8(i,j,k)*EVAL6(i,j,k)
               H_values(indx+6) = -EVAL8(i,j,k)*EVAL5(i,j,k)

c--- LAPLACIAN TEST CASE
c               H_values(indx)   = 6.0
c               H_values(indx+1) = -1.0
c               H_values(indx+2) = -1.0
c               H_values(indx+3) = -1.0
c               H_values(indx+4) = -1.0
c               H_values(indx+5) = -1.0
c               H_values(indx+6) = -1.0

!-- values at and beyond upperbnd
               if(MPIright.eq.MPI_PROC_NULL.and.i.ge.upperbnd(j,k)) then !--------right boundary
c               if(MPIright.eq.MPI_PROC_NULL.and.i.ge.
c     $              upperbnd(j,k)-10) then !--------right boundary
                  H_values(indx) = 1.0
                  H_values(indx+1) = 0.0
                  H_values(indx+2) = 0.0
                  H_values(indx+3) = 0.0
                  H_values(indx+4) = 0.0
                  H_values(indx+5) = 0.0
                  H_values(indx+6) = 0.0
               endif
!-- neighbor values beyond upperbnd
               if(MPIright.eq.MPI_PROC_NULL.and.
     %              i.eq.upperbnd(j,k)-1) then !--------right-neighbor boundary
c               if(MPIright.eq.MPI_PROC_NULL.and.
c     %              i.eq.upperbnd(j,k)-11) then !--------right-neighbor boundary
                  H_values(indx+2) = 0.0
               endif

!---- BOUNDARIES
               if( (MPIleft.eq.MPI_PROC_NULL.and.i.eq.1) .or. !---------left boundary
     1              (MPIright.eq.MPI_PROC_NULL.and.i.eq.locNX) .or.!----right boundary
     2              (MPIlower.eq.MPI_PROC_NULL.and.j.eq.1) .or.!--------lower boundary
     3              (MPIupper.eq.MPI_PROC_NULL.and.j.eq.locNY) .or. !---upper boundary
     4              (.not.poles.and.MPIbelow.eq.MPI_PROC_NULL.and.
     5              k.eq.1) .or. !-------below boundary
     6              (.not.poles.and.MPIabove.eq.MPI_PROC_NULL.and.
     7              k.eq.locNZ) ) then !-above boundary
                  H_values(indx) = 1.0
                  H_values(indx+1) = 0.0
                  H_values(indx+2) = 0.0
                  H_values(indx+3) = 0.0
                  H_values(indx+4) = 0.0
                  H_values(indx+5) = 0.0
                  H_values(indx+6) = 0.0
               endif
!---- NEIGHBOR BOUNDARY VALUES
!--------left-neighbor boundary
               if(MPIleft.eq.MPI_PROC_NULL.and.i.eq.2) then
                  H_values(indx+1) = 0.0
               endif
!--------right-neighbor boundary
               if(MPIright.eq.MPI_PROC_NULL.and.i.eq.locNX-1) then
                  H_values(indx+2) = 0.0
               endif
!--------lower-neighbor boundary
               if(MPIlower.eq.MPI_PROC_NULL.and.j.eq.2) then
                  H_values(indx+3) = 0.0
                  print *,'IN BAD SPOT NEIGHBOR'
               endif
!--------upper-neighbor boundary
               if(MPIupper.eq.MPI_PROC_NULL.and.j.eq.locNY-1) then
                  H_values(indx+4) = 0.0
                  print *,'IN BAD SPOT NEIGHBOR'
               endif
!--------below-neighbor boundary
               if(.not.poles.and.MPIbelow.eq.MPI_PROC_NULL.and.
     $              k.eq.2) then
                  H_values(indx+5) = 0.0
               endif
!--------above-neighbor boundary
               if(.not.poles.and.MPIabove.eq.MPI_PROC_NULL.and.
     $              k.eq.locNZ-1) then
                  H_values(indx+6) = 0.0
               endif
               indx = indx + 7
            enddo
         enddo
      enddo
      call bHYPRE_SStructMatrix_SetBoxValues_f( H_matrix, hypre_part,
     1     ilower, iupper, hypre_ndim, hypre_var, nentries, 
     2     stencil_indices, H_values, nvalues, ierr, except )
      return
      end SUBROUTINE FILLHYPREMATRIX

      SUBROUTINE CREATEHYPREVECTORS
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE hypre_var_init
      IMPLICIT NONE
      include 'bHYPRE_SStructVariable.inc'
!--- SETUP THE X AND B HYPRE VECTORS
      call bHYPRE_SStructVector_Create_f( bHYPRE_mpicomm, grid, b_vec,
     1     except )
      call bHYPRE_SStructVector_Create_f( bHYPRE_mpicomm, grid, x_vec,
     1     except )
!--- Set the object type for the vectors to be the struct type
      object_type = HYPRE_STRUCT;
      call bHYPRE_SStructVector_SetObjectType_f( b_vec, object_type,
     1     ierr, except)
      call bHYPRE_SStructVector_SetObjectType_f( x_vec, object_type,
     1     ierr, except)
      RETURN
      END SUBROUTINE CREATEHYPREVECTORS

      SUBROUTINE CREATEHYPREMATRIX
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE hypre_var_init
      IMPLICIT NONE
      include 'bHYPRE_SStructVariable.inc'
!--- Create the empty SSTRUC matrix object
      call bHYPRE_SStructMatrix_Create_f( bHYPRE_mpicomm, graph,
     1     H_matrix, except )
!--- Set the object type (=HYPRE_STRUCT, ie a purely structured example)
      object_type = HYPRE_STRUCT
      call bHYPRE_SStructMatrix_SetObjectType_f( H_matrix, object_type,
     1     ierr, except )
      RETURN
      END SUBROUTINE CREATEHYPREMATRIX

      SUBROUTINE SETUPHYPREGRID
!-routine to setup a standard hypre grid without the special condition
! for the poles
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE hypre_var_init
      IMPLICIT NONE
      include "mpif.h"
      include 'bHYPRE_SStructVariable.inc'
      if(myid.eq.0) then 
         print *,'Setting up HYPRE grid WITHOUT polar conditions'
      endif
      hypre_ndim = 3
      hypre_nparts = 1
      hypre_part = 0
      hypre_var = 0
      call bHYPRE_SStructGrid_Create_f(bHYPRE_mpicomm, hypre_ndim,
     1     hypre_nparts, grid, except)
!--- Set the extents of the grid using a single global grid over all
!     the processors. In this style setup this step is what tells hypre
!     which processors are its neighboors
      ilower(1) =proc_coords(1)*locNX+1
      ilower(2) =proc_coords(2)*locNY+1
      ilower(3) =proc_coords(3)*locNZ+1
      iupper(1) =(proc_coords(1)+1)*locNX
      iupper(2) =(proc_coords(2)+1)*locNY
      iupper(3) =(proc_coords(3)+1)*locNZ
      call bHYPRE_SStructGrid_SetExtents_f( grid, hypre_part,
     1     ilower, iupper, hypre_ndim, ierr, except )
c         write(*,'(I8,A,1x,3(1x,I8),A,A,3(1x,I8),A)')
c     1     myid,'(',ilower(1),ilower(2),ilower(3),')',
c     $     '(',iupper(1),iupper(2),iupper(3),')'
!--- Set the variable(=0) type(=CELL) and number(=1) on each part
      vartypes(1) = CELL
      call bHYPRE_SStructGrid_SetVariable_f( grid, 0, hypre_var, 1,
     1     vartypes(1), ierr, except )
      call bHYPRE_SStructGrid_Assemble_f( grid, ierr, except )
      RETURN
      END SUBROUTINE SETUPHYPREGRID


      SUBROUTINE SETUPHYPREGRID_PARTS
!-- Setup hypre using numproc 'parts'. You must first call the normal setuphypregrid,
!     then this routine adds information about the polar processors neighboring's cells
!     setting them with SetNeighborBox_f.
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE hypre_var_init
      IMPLICIT NONE
      include "mpif.h"
      include 'bHYPRE_SStructVariable.inc'
      integer*4 ilower_my_indx(3),iupper_my_indx(3)
      integer*4 ilower_neigh_indx(3),iupper_neigh_indx(3)
      integer*4 hypre_index_map(3)
      integer neigh_proc_yindx
      if(myid.eq.0) then 
         print *,'Setting up HYPRE grid WITH polar conditions'
      endif
!- set parameters for the hypre grid, mimicing the mpi grid
      hypre_ndim = 3
      hypre_nparts = 1
      hypre_part = 0
      hypre_var = 0

!--- Check that there are no processors to left or right      
      if(MPIleft.ne.MPI_PROC_NULL.or.
     $     MPIright.ne.MPI_PROC_NULL) then !--left/right x-edges
         print *,'SETUPHYPREGRID_PARTS is not setup for processors'
         print *,'  in X-DIRECTION'
         call clean_stop
      endif

!--- Setup the dimension index maping between the part. This is the same for all
      hypre_index_map(1)=0
      hypre_index_map(2)=1
      hypre_index_map(3)=2

!--- Setup spatial relation between processors in the +Z (above) direction
!     for the processors at the north pole
      if(poles.and.proc_coords(3).eq.proc_dims(3)-1) then
!- in local parts space
         ilower_my_indx(1) = proc_coords(1)*locNX+1
         ilower_my_indx(2) = proc_coords(2)*locNY+1
         ilower_my_indx(3) = (proc_coords(3)+1)*locNZ+1
         iupper_my_indx(1) = (proc_coords(1)+1)*locNX
         iupper_my_indx(2) = (proc_coords(2)+1)*locNY
         iupper_my_indx(3) = (proc_coords(3)+1)*locNZ+1
!- in neighboring parts space
         neigh_proc_yindx = mod(proc_coords(2)+
     $        proc_dims(2)/2,proc_dims(2))
c         write(*,'(A,1x,3(1x,I8))') 'NP neigh test: ',myid,
c     $        proc_coords(2),neigh_proc_yindx
         ilower_neigh_indx(1) =  proc_coords(1)*locNX+1
         ilower_neigh_indx(2) =  neigh_proc_yindx*locNY+1
         ilower_neigh_indx(3) =  (proc_coords(3)+1)*locNZ-2
         iupper_neigh_indx(1) =  (proc_coords(1)+1)*locNX
         iupper_neigh_indx(2) =  (neigh_proc_yindx+1)*locNY
         iupper_neigh_indx(3) =  (proc_coords(3)+1)*locNZ-2
c      write(*,'(I8,A,1x,3(1x,I8),A,A,3(1x,I8),A)')
c     1     myid,'north pole mine(',ilower_my_indx(1),
c     1        ilower_my_indx(2)
c     2        ,ilower_my_indx(3),')',
c     $     '(',iupper_my_indx(1),iupper_my_indx(2),
c     4        iupper_my_indx(3),')'
c      write(*,'(I8,A,1x,3(1x,I8),A,A,3(1x,I8),A)')
c     1     myid,'north pole neigh(',ilower_neigh_indx(1),
c     1        ilower_neigh_indx(2)
c     2        ,ilower_neigh_indx(3),')',
c     $     '(',iupper_neigh_indx(1),iupper_neigh_indx(2),
c     4        iupper_neigh_indx(3),')'
c         print *,'setting strucgrid,NP',myid
         call bHYPRE_SStructGrid_SetNeighborBox_f(grid,hypre_part,
     $        ilower_my_indx, iupper_my_indx, MPIabove_pole,
     $        ilower_neigh_indx,iupper_neigh_indx,
     $        hypre_index_map, ierr, except )
c         print *,'done setneigh NP',myid
      endif


!--- Setup spatial relation between processors in the -Z (below) direction
!     for the processors at the south pole
      if(poles.and.proc_coords(3).eq.0) then
!- in local parts space
         ilower_my_indx(1) = proc_coords(1)*locNX+1
         ilower_my_indx(2) = proc_coords(2)*locNY+1
         ilower_my_indx(3) = (proc_coords(3))*locNZ ! I could put 0 here since it always is
         iupper_my_indx(1) = (proc_coords(1)+1)*locNX
         iupper_my_indx(2) = (proc_coords(2)+1)*locNY
         iupper_my_indx(3) = (proc_coords(3))*locNZ ! I could put 0 here since it always is
!- in neighboring parts space
         neigh_proc_yindx = mod(proc_coords(2)+
     $        proc_dims(2)/2,proc_dims(2))  
c         write(*,'(A,1x,3(1x,I8))') 'SP neigh test: ',myid,
c     $        proc_coords(2),neigh_proc_yindx
         ilower_neigh_indx(1) =  proc_coords(1)*locNX+1
         ilower_neigh_indx(2) =  neigh_proc_yindx*locNY+1
         ilower_neigh_indx(3) =  proc_coords(3)*locNZ+1 ! I could put 1 here since it always is
         iupper_neigh_indx(1) =  (proc_coords(1)+1)*locNX
         iupper_neigh_indx(2) =  (neigh_proc_yindx+1)*locNY
         iupper_neigh_indx(3) =  proc_coords(3)*locNZ+1 ! I could put 1 here since it always is
c         write(*,'(I8,A,1x,3(1x,I8),A,A,3(1x,I8),A)')
c     1        myid,'south pole mine(',ilower_my_indx(1),
c     1        ilower_my_indx(2)
c     2        ,ilower_my_indx(3),')',
c     $        '(',iupper_my_indx(1),iupper_my_indx(2),
c     4        iupper_my_indx(3),')'
c         write(*,'(I8,A,1x,3(1x,I8),A,A,3(1x,I8),A)')
c     1        myid,'south pole neigh(',ilower_neigh_indx(1),
c     1        ilower_neigh_indx(2)
c     2        ,ilower_neigh_indx(3),')',
c     $        '(',iupper_neigh_indx(1),iupper_neigh_indx(2),
c     4        iupper_neigh_indx(3),')'
         call bHYPRE_SStructGrid_SetNeighborBox_f(grid,hypre_part,
     $        ilower_my_indx, iupper_my_indx, MPIabove_pole,
     $        ilower_neigh_indx,iupper_neigh_indx,
     $        hypre_index_map, ierr, except )
c         print *,'done setneigh SP',myid

      endif

      RETURN
      END SUBROUTINE SETUPHYPREGRID_PARTS


      SUBROUTINE SETUPHYPRESTENCIL
      USE input_init
      USE fluid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE hypre_var_init
      IMPLICIT NONE
      include 'bHYPRE_SStructVariable.inc' 
      integer :: j
      call bHYPRE_SStructStencil_Create_f(hypre_ndim,7,stencil,except)
!--- Assign stencil offsets
!---   (i,j,k)
      offsets(1,1) = 0
      offsets(2,1) = 0
      offsets(3,1) = 0
!---   (i-1,j,k)
      offsets(1,2) = -1
      offsets(2,2) = 0
      offsets(3,2) = 0
!---   (i+1,j,k)
      offsets(1,3) = 1
      offsets(2,3) = 0
      offsets(3,3) = 0
!---   (i,j-1,k)
      offsets(1,4) = 0
      offsets(2,4) = -1
      offsets(3,4) = 0
!---   (i,j+1,k)
      offsets(1,5) = 0
      offsets(2,5) = 1
      offsets(3,5) = 0
!---   (i,j,k-1)
      offsets(1,6) = 0
      offsets(2,6) = 0
      offsets(3,6) = -1
!---   (i,j,k+1)
      offsets(1,7) = 0
      offsets(2,7) = 0
      offsets(3,7) = 1
!--- Fill stencil and assign numerical values to the offsets for var=0
      do sentry = 1, 7
         call bHYPRE_SStructStencil_SetEntry_f( stencil, sentry-1,
     1        offsets(1,sentry), hypre_ndim, hypre_var, ierr, except )
      enddo
!--- Create the graph object
      call bHYPRE_SStructGraph_Create_f( bHYPRE_mpicomm, grid, graph,
     1     except )
!--- Connect the graph and stencil objects
      call bHYPRE_SStructGraph_SetStencil_f( graph, hypre_part, 
     1     hypre_var, stencil, ierr, except )

!--  For simulations that include the polar region, run the 
!     setNeighbor_box routines
      if(poles) then
         CALL SETUPHYPREGRID_PARTS
      endif

!--- Assemble the final graph
      call bHYPRE_SStructGraph_Assemble_f( graph, ierr, except )
!---  Set stencil indices
      nentries=7
      do j = 1, nentries
         stencil_indices(j) = j-1
      enddo
      RETURN
      END SUBROUTINE SETUPHYPRESTENCIL

