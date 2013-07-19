
      subroutine laplacetest
      USE input_init
      USE fluid_var_init
      USE rad_var_init
      USE grid_var_init
      USE deltaER_var_init
      USE mpi_var_init
      USE global_constants
      USE hypre_var_init
      IMPLICIT NONE
      integer i,j,k,indx

!----  CREATE THE HYPRE MATRIX
!--- Create the empty SSTRUC matrix object
      call bHYPRE_SStructMatrix_Create_f( bHYPRE_mpicomm, graph,
     1     H_matrix, except )
!--- Set the object type (=HYPRE_STRUCT, ie a purely structured example)
      object_type = HYPRE_STRUCT
      call bHYPRE_SStructMatrix_SetObjectType_f( H_matrix, object_type,
     1     ierr, except )

!----  CREATE THE X AND B HYPRE VECTORS
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


!----  FILL THE MATRIX (AND SET BOUNDARIES)
      CALL FILLHYPREMATRIX
!----  FILL THE VECTORS (AND SET BOUNDARIES)
      CALL FILLHYPREVECTORS

!---  Assemble vectors and matrix
      call bHYPRE_SStructVector_Assemble_f( b_vec, ierr, except )
      call bHYPRE_SStructVector_Assemble_f( x_vec, ierr, except )
      call bHYPRE_SStructMatrix_Assemble_f( H_matrix, ierr, except)

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
      return
      end
