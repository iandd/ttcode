
      SUBROUTINE SHIFTVAR2D(f)
!--- SHIFT variable in 2d
!    f=1-5 => T,RH,V,G,H
!    f=6   => ER
!    f=7   => ER w/out  -1 and N+2
!-------------------
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      USE dif_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: f,k,stat(MPI_STATUS_SIZE)
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: fluid_point
      if(f.eq.1) fluid_point => T
      if(f.eq.2) fluid_point => RH
      if(f.eq.3) fluid_point => V
      if(f.eq.4) fluid_point => G
      if(f.eq.5) fluid_point => H
      if(f.eq.6) fluid_point => ER
      if(f.eq.7) fluid_point => ER
      if(f.ge.8) then
         print *,'not avalible:SHIFTVAR2D'
         stop
      endif
!---- X-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(locNX,1:locNY,1:locNZ),
     %     locNY*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %     fluid_point(0,1:locNY,1:locNZ),locNY*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
      if(f.ne.3.and.f.ne.7) THEN
!V(-1,j,k) doesn't exist: don't pass for SOR
         call MPI_SENDRECV(fluid_point(locNX-1,1:locNY,1:locNZ),
     %        locNY*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %        fluid_point(-1,1:locNY,1:locNZ),locNY*locNZ,
     %        MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
      endif
!---- X-COORD, - DIRECTION 
      call MPI_SENDRECV(fluid_point(1,1:locNY,1:locNZ),
     %     locNY*locNZ,MPI_DOUBLE_PRECISION,MPIleft,1,
     %     fluid_point(locNX+1,1:locNY,1:locNZ),locNY*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIright,1,COMM_CART,stat,ierr)
      if(f.ne.7) then !Don't pass for SOR
         call MPI_SENDRECV(fluid_point(2,1:locNY,1:locNZ),
     %        locNY*locNZ,MPI_DOUBLE_PRECISION,MPIleft,1,
     %        fluid_point(locNX+2,1:locNY,1:locNZ),locNY*locNZ,
     %        MPI_DOUBLE_PRECISION,MPIright,1,COMM_CART,stat,ierr)
      endif
!---- Y-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,locNY,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     fluid_point(1:locNX,0,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.4.and.f.ne.7) then 
!G(i,-1,k) doesn't exist: don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,locNY-1,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %        fluid_point(1:locNX,-1,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %        COMM_CART,stat,ierr)
      endif
!---- Y-COORD, - DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,1,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     fluid_point(1:locNX,locNY+1,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.7) then !Don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,2,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %        fluid_point(1:locNX,locNY+2,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %        COMM_CART,stat,ierr)
      endif
!----- CORNER SHIFT FOR 2D PROCESSORS
      k=1
!----- lower left and upper right corner
      call MPI_SENDRECV(fluid_point(1,1,k),
     %     1,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %     fluid_point(locNX+1,locNY+1,k),
     %     1,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(locNX,locNY,k),
     %     1,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %     fluid_point(0,0,k),
     %     1,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %     COMM_CART,stat,ierr)
!----- upper left and lower right corner
      call MPI_SENDRECV(fluid_point(locNX,1,k),
     %     1,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %     fluid_point(0,locNY+1,k),
     %     1,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1,locNY,k),
     %     1,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %     fluid_point(locNX+1,0,k),
     %     1,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %     COMM_CART,stat,ierr)
!-----DONE SHIFTING
!--- RELEASE THE POINTER 
      NULLIFY(fluid_point)
      RETURN
      END SUBROUTINE SHIFTVAR2D


      SUBROUTINE SETUP2DMPIGRID
      USE input_init
      USE mpi_grid_init
      USE mpi_var_init
      USE fluid_var_init
      IMPLICIT NONE
      integer ::shift_coord,shift_dir
      include "mpif.h"
!                      (upper)
!          !-----!-----!------!------!
!          !     !     !      !      !
!          !  2  !  5  !  8   !  11  !
!          !-----!-----!------!------!
!       y  !     !     !      !      !
!(left) ^  !  1  !  4  !  7   !  10  ! (right)
!       !  !-----!-----!------!------!
!       !  !     !     !      !      !
!       !  !  0  !  3  !  6   !  9   !
!       !  !-----!-----!------!------!
!       !------->x    (lower)
!
      IF(locNZ.ne.1) then
         print *,'not yet defined in the z-direction:newshift'
         stop
      ENDIF
      IF(myid.eq.1) then
         print *,'Setting up 2D processor grid'
      ENDIF
!---- X-COORD, + DIRECTION
      shift_coord = 0
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIleft  = source
      MPIright  = dest
!---- Y-COORD, + DIRECTION
      shift_coord = 1
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIlower  = source
      MPIupper = dest

      if(MODTYP.eq.3.and.MODVER.eq.0) then !-PARTIAL DISK
         MPIupperright = MPIright+1
         MPIlowerright = MPIright-1
         MPIupperleft  = MPIleft+1
         MPIlowerleft  = MPIleft-1
!-- y boundary processors
         if(proc_coords(2).eq.0) then
            MPIlowerright = MPI_PROC_NULL
            MPIlowerleft  = MPI_PROC_NULL
         endif
         if(proc_coords(2).eq.proc_dims(2)-1) then
            MPIupperright = MPI_PROC_NULL
            MPIupperleft  = MPI_PROC_NULL
         endif
!-- all the variables that contain z-information (ie. above and below)
         MPIabove           = MPI_PROC_NULL
         MPIbelow           = MPI_PROC_NULL
         MPIupperabove      = MPI_PROC_NULL
         MPIupperbelow      = MPI_PROC_NULL
         MPIlowerabove      = MPI_PROC_NULL
         MPIlowerbelow      = MPI_PROC_NULL
         MPIaboveright      = MPI_PROC_NULL
         MPIbelowright      = MPI_PROC_NULL
         MPIaboveleft       = MPI_PROC_NULL
         MPIbelowleft       = MPI_PROC_NULL
         MPIupperaboveright = MPI_PROC_NULL
         MPIupperbelowright = MPI_PROC_NULL
         MPIloweraboveright = MPI_PROC_NULL
         MPIlowerbelowright = MPI_PROC_NULL
         MPIupperaboveleft  = MPI_PROC_NULL
         MPIupperbelowleft  = MPI_PROC_NULL
         MPIloweraboveleft  = MPI_PROC_NULL
         MPIlowerbelowleft  = MPI_PROC_NULL
      ELSEIF((MODTYP.eq.3.and.MODVER.eq.1).or.MODTYP.eq.4.or.
     %        MODTYP.eq.5) then !-FULL DISK
         MPIupperright = MPIupper + proc_dims(2)
         MPIlowerright = MPIlower + proc_dims(2)
         MPIupperleft  = MPIupper - proc_dims(2)
         MPIlowerleft  = MPIlower - proc_dims(2)
      else
         print *,'MODVER/MODTYP=',MODVER/MODTYP,
     %        ' is not defined:SETUP2DMPIGRID'
         stop
      endif

!-- x boundary processors
      if(proc_coords(1).eq.0) then
         MPIupperleft = MPI_PROC_NULL
         MPIlowerleft = MPI_PROC_NULL
      endif
      if(proc_coords(1).eq.proc_dims(1)-1) then
         MPIupperright = MPI_PROC_NULL
         MPIlowerright = MPI_PROC_NULL
      endif
c      print *,myid,MPIlower,MPIupper
c      print *,myid,MPIleft,MPIright
c      print *,myid,MPIupperright,MPIupperleft
c      print *,myid,MPIlowerright,MPIlowerleft
c      print *,myid,
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      stop
      RETURN
      END SUBROUTINE SETUP2DMPIGRID
      

      SUBROUTINE SHIFTVAR3D(f)
!--- SHIFT variable in 3D
!    f=1-5 => T,RH,V,G,H
!    f=6   => ER
!    f=7   => ER w/out  -1 and N+2
!-------------------
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: f,k,stat(MPI_STATUS_SIZE)
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: fluid_point
      if(f.eq.1) fluid_point => T
      if(f.eq.2) fluid_point => RH
      if(f.eq.3) fluid_point => V
      if(f.eq.4) fluid_point => G
      if(f.eq.5) fluid_point => H
      if(f.eq.6) fluid_point => ER
      if(f.eq.7) fluid_point => ER !(don't pass for:i=-1,locNX+2,j=-1,locNY+2)
      if(f.ge.8) then
         print *,'not avalible:SHIFTVAR3D'
         stop
      endif
!---- X-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(locNX,1:locNY,1:locNZ),
     %     locNY*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %     fluid_point(0,1:locNY,1:locNZ),locNY*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
      if(f.ne.3.and.f.ne.7) THEN
!V(-1,j,k) doesn't exist: don't pass for SOR
         call MPI_SENDRECV(fluid_point(locNX-1,1:locNY,1:locNZ),
     %        locNY*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %        fluid_point(-1,1:locNY,1:locNZ),locNY*locNZ,
     %        MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
      endif
!---- X-COORD, - DIRECTION 
      call MPI_SENDRECV(fluid_point(1,1:locNY,1:locNZ),
     %     locNY*locNZ,MPI_DOUBLE_PRECISION,MPIleft,1,
     %     fluid_point(locNX+1,1:locNY,1:locNZ),locNY*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIright,1,COMM_CART,stat,ierr)
      if(f.ne.7) then !Don't pass for SOR
         call MPI_SENDRECV(fluid_point(2,1:locNY,1:locNZ),
     %        locNY*locNZ,MPI_DOUBLE_PRECISION,MPIleft,1,
     %        fluid_point(locNX+2,1:locNY,1:locNZ),locNY*locNZ,
     %        MPI_DOUBLE_PRECISION,MPIright,1,COMM_CART,stat,ierr)
      endif
!---- Y-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,locNY,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     fluid_point(1:locNX,0,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.4.and.f.ne.7) then 
!G(i,-1,k) doesn't exist: don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,locNY-1,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %        fluid_point(1:locNX,-1,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %        COMM_CART,stat,ierr)
      endif
!---- Y-COORD, - DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,1,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     fluid_point(1:locNX,locNY+1,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.7) then !Don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,2,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %        fluid_point(1:locNX,locNY+2,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %        COMM_CART,stat,ierr)
      endif
!---- Z-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,1:locNY,locNZ),
     %     locNX*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     fluid_point(1:locNX,1:locNY,0),
     %     locNX*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.5.and.f.ne.7) then !H(i,j,-1) doesn't exist, don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,1:locNY,locNZ-1),
     %        locNX*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %        fluid_point(1:locNX,1:locNY,-1),
     %        locNX*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %        COMM_CART,stat,ierr)
      endif
!---- Z-COORD, - DIRECTION (JUST REVERSE SOURCE AND DEST)
      call MPI_SENDRECV(fluid_point(1:locNX,1:locNY,1),
     %     locNX*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     fluid_point(1:locNX,1:locNY,locNZ+1),
     %     locNX*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.7) then !Don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,1:locNY,2),
     %        locNX*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %        fluid_point(1:locNX,1:locNY,locNZ+2),
     %        locNX*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %        COMM_CART,stat,ierr)
      endif
!
!----- EDGE SHIFTING
!
!--------y-z plane
!--------- upperabove and lowerbelow
      call MPI_SENDRECV(fluid_point(1:locNX,locNY,locNZ),
     %     locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %     fluid_point(1:locNX,0,0),
     %     locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1:locNX,1,1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %     fluid_point(1:locNX,locNY+1,locNZ+1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %     COMM_CART,stat,ierr)
!--------- lowerabove and upperbelow
      call MPI_SENDRECV(fluid_point(1:locNX,1,locNZ),
     %     locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %     fluid_point(1:locNX,locNY+1,0),
     %     locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1:locNX,locNY,1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %     fluid_point(1:locNX,0,locNZ+1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %     COMM_CART,stat,ierr)
!--------x-y plane
!--------- upperright and lowerleft
      call MPI_SENDRECV(fluid_point(1,1,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %     fluid_point(locNX+1,locNY+1,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(locNX,locNY,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %     fluid_point(0,0,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %     COMM_CART,stat,ierr)
!--------- upperleft and lowerright
      call MPI_SENDRECV(fluid_point(1,locNY,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %     fluid_point(locNX+1,0,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(locNX,1,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %     fluid_point(0,locNY+1,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %     COMM_CART,stat,ierr)
!--------x-z plane
!--------- belowleft and aboveright
      call MPI_SENDRECV(fluid_point(1,1:locNY,1),
     %     locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %     fluid_point(locNX+1,1:locNY,locNZ+1),
     %     locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(locNX,1:locNY,locNZ),
     %     locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %     fluid_point(0,1:locNY,0),
     %     locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %     COMM_CART,stat,ierr)
!--------- belowright and aboveleft
      call MPI_SENDRECV(fluid_point(locNX,1:locNY,1),
     %     locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %     fluid_point(0,1:locNY,locNZ+1),
     %     locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1,1:locNY,locNZ),
     %     locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %     fluid_point(locNX+1,1:locNY,0),
     %     locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %     COMM_CART,stat,ierr)
!
!----- CORNER SHIFTING
!
!--------- upperaboveright and lowerbelowleft
      call MPI_SENDRECV(fluid_point(locNX,locNY,locNZ),
     %     1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %     fluid_point(0,0,0),
     %     1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1,1,1),
     %     1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %     fluid_point(locNX+1,locNY+1,locNZ+1),
     %     1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %     COMM_CART,stat,ierr)
!--------- loweraboveright and upperbelowleft
      call MPI_SENDRECV(fluid_point(locNX,1,locNZ),
     %     1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %     fluid_point(0,locNY+1,0),
     %     1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1,locNY,1),
     %     1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %     fluid_point(locNX+1,0,locNZ+1),
     %     1,MPI_DOUBLE_PRECISION,MPIloweraboveright ,1,
     %     COMM_CART,stat,ierr)
!--------- upperbelowright and loweraboveleft
      call MPI_SENDRECV(fluid_point(locNX,locNY,1),
     %     1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %     fluid_point(0,0,locNZ+1),
     %     1,MPI_DOUBLE_PRECISION,MPIloweraboveleft,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1,1,locNZ),
     %     1,MPI_DOUBLE_PRECISION,MPIloweraboveleft ,1,
     %     fluid_point(locNX+1,locNY+1,0),
     %     1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %     COMM_CART,stat,ierr)
!--------- lowerbelowright and upperaboveleft
      call MPI_SENDRECV(fluid_point(locNX,1,1),
     %     1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %     fluid_point(0,locNY+1,locNZ+1),
     %     1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1,locNY,locNZ),
     %     1,MPI_DOUBLE_PRECISION,MPIupperaboveleft ,1,
     %     fluid_point(locNX+1,0,0),
     %     1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %     COMM_CART,stat,ierr)
!-----DONE SHIFTING
!--- RELEASE THE POINTER 
      NULLIFY(fluid_point)
      RETURN
      END SUBROUTINE SHIFTVAR3D


      SUBROUTINE SHIFTVAR_AXISYMMETRIC(f)
!--- SHIFT variable axisymmetricly
!    f=1-5 => T,RH,V,G,H
!    f=6   => ER
!    f=7   => ER w/out  -1 and N+2
!-------------------
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_grid_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: f,k,stat(MPI_STATUS_SIZE)
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: fluid_point
      if(f.eq.1) fluid_point => T
      if(f.eq.2) fluid_point => RH
      if(f.eq.3) fluid_point => V
      if(f.eq.4) fluid_point => G
      if(f.eq.5) fluid_point => H
      if(f.eq.6) fluid_point => ER
      if(f.eq.7) then
         print *, 'f=7 (for ER) should never use axisymmetric shift'
         stop
      endif
      if(f.eq.8) fluid_point => delta_ER
      if(f.ge.9) then
         print *,'not avalible:SHIFTVAR_AXISYMMETRIC',f
         stop
      endif

!---- X-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(locNX,1,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %     fluid_point(0,1,1:locNZ),locNZ,
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
      if(f.ne.3) THEN
!-- V(-1,j,k) doesn't exist
         call MPI_SENDRECV(fluid_point(locNX-1,1,1:locNZ),
     %        locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %        fluid_point(-1,1,1:locNZ),locNZ,
     %        MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
      endif
!---- X-COORD, - DIRECTION 
      call MPI_SENDRECV(fluid_point(1,1,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIleft,1,
     %     fluid_point(locNX+1,1,1:locNZ),locNZ,
     %     MPI_DOUBLE_PRECISION,MPIright,1,COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(2,1,1:locNZ),
     %     locNZ,MPI_DOUBLE_PRECISION,MPIleft,1,
     %     fluid_point(locNX+2,1,1:locNZ),locNZ,
     %     MPI_DOUBLE_PRECISION,MPIright,1,COMM_CART,stat,ierr)


!---- DON'T DO Y-COORD, + DIRECTION
!---- DON'T DO Y-COORD, - DIRECTION

!---- Z-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,1,locNZ),
     %     locNX,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     fluid_point(1:locNX,1,0),
     %     locNX,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.5) then !H(i,j,-1) doesn't exist
         call MPI_SENDRECV(fluid_point(1:locNX,1,locNZ-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIabove,1,
     %        fluid_point(1:locNX,1,-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %        COMM_CART,stat,ierr)
      endif
!---- Z-COORD, - DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,1,1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     fluid_point(1:locNX,1,locNZ+1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1:locNX,1,2),
     %     locNX,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     fluid_point(1:locNX,1,locNZ+2),
     %     locNX,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     COMM_CART,stat,ierr)

!
!----- EDGE SHIFTING
!
!--------don't do y-z plane
!--------don't do x-y plane
!
!--------x-z plane
!--------- to belowleft
      call MPI_SENDRECV(fluid_point(1,1,1),
     %     1,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %     fluid_point(locNX+1,1,locNZ+1),
     %     1,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(2,1,2),
     %     1,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %     fluid_point(locNX+2,1,locNZ+2),
     %     1,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(1,1,2),
     %     1,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %     fluid_point(locNX+1,1,locNZ+2),
     %     1,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(2,1,1),
     %     1,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %     fluid_point(locNX+2,1,locNZ+1),
     %     1,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %     COMM_CART,stat,ierr)
!-- to aboveright
      call MPI_SENDRECV(fluid_point(locNX,1,locNZ),
     %     1,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %     fluid_point(0,1,0),
     %     1,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.3.and.f.ne.5) then !-- no V(-1,1,-1) or H(-1,1,-1)
         call MPI_SENDRECV(fluid_point(locNX-1,1,locNZ-1),
     %        1,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %        fluid_point(-1,1,-1),
     %        1,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %        COMM_CART,stat,ierr)
      endif
      if(f.ne.3) then  !-- no V(-1,1,0)
         call MPI_SENDRECV(fluid_point(locNX-1,1,locNZ),
     %        1,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %        fluid_point(-1,1,0),
     %        1,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %        COMM_CART,stat,ierr)
      endif
      if(f.ne.5) then  !-- no H(0,1,-1)
         call MPI_SENDRECV(fluid_point(locNX,1,locNZ-1),
     %        1,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %        fluid_point(0,1,-1),
     %        1,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %        COMM_CART,stat,ierr)
      endif
!--------- to aboveleft
      call MPI_SENDRECV(fluid_point(1,1,locNZ),
     %     1,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %     fluid_point(locNX+1,1,0),
     %     1,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.5) then !-- no H(*,1,-1)
         call MPI_SENDRECV(fluid_point(2,1,locNZ-1),
     %        1,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %        fluid_point(locNX+2,1,-1),
     %        1,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %        COMM_CART,stat,ierr)
         call MPI_SENDRECV(fluid_point(1,1,locNZ-1),
     %        1,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %        fluid_point(locNX+1,1,-1),
     %        1,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %        COMM_CART,stat,ierr)
      endif
      call MPI_SENDRECV(fluid_point(2,1,locNZ),
     %     1,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %     fluid_point(locNX+2,1,0),
     %     1,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %     COMM_CART,stat,ierr)

!---------to belowright
      call MPI_SENDRECV(fluid_point(locNX,1,1),
     %     1,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %     fluid_point(0,1,locNZ+1),
     %     1,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.3) then !-- no V(-1,*,*) 
         call MPI_SENDRECV(fluid_point(locNX-1,1,2),
     %        1,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %        fluid_point(-1,1,locNZ+2),
     %        1,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %        COMM_CART,stat,ierr)
         call MPI_SENDRECV(fluid_point(locNX-1,1,1),
     %        1,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %        fluid_point(-1,1,locNZ+1),
     %        1,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %        COMM_CART,stat,ierr)
      endif
      call MPI_SENDRECV(fluid_point(locNX,1,2),
     %     1,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %     fluid_point(0,1,locNZ+2),
     %     1,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %     COMM_CART,stat,ierr)
!-----DONE SHIFTING
!--- RELEASE THE POINTER 
      NULLIFY(fluid_point)
      RETURN
      END SUBROUTINE SHIFTVAR_AXISYMMETRIC


      SUBROUTINE SETUP3DMPIGRID
      USE input_init
      USE mpi_grid_init
      USE mpi_var_init
      USE fluid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer ::shift_coord,shift_dir
      integer :: stat(MPI_STATUS_SIZE)
!------------------------------------
!     LEVEL 1 (X=0 : left)
!                    (upper)
!          !-----!-----!------!
!          !     !     !      !
!          !  6  !  7  !  8   !
!          !-----!-----!------!
!       Y  !     !     !      !
!(below)^  !  3  !  4  !  5   ! (above)
!       !  !-----!-----!------!
!       !  !     !     !      !
!       !  !  0  !  1  !  2   !
!       !  !-----!-----!------!
!       !------->Z   (lower)
!
!     LEVEL 2 (X=1 : middle)
!                    (upper)
!          !-----!-----!------!
!          !     !     !      !
!          !  15 !  16 !  17  !
!          !-----!-----!------!
!       Y  !     !     !      !
!(below)^  !  12 !  13 !  14  ! (above)
!       !  !-----!-----!------!
!       !  !     !     !      !
!       !  !  9  !  10 ! 11   !
!       !  !-----!-----!------!
!       !------->Z   (lower)
!
!     LEVEL 3 (X=2 : right)
!                    (upper)
!          !-----!-----!------!
!          !     !     !      !
!          !  24 !  25 !  26  !
!          !-----!-----!------!
!       Y  !     !     !      !
!(below)^  !  21 !  22 !  23  ! (above)
!       !  !-----!-----!------!
!       !  !     !     !      !
!       !  !  18 !  19 !  20  !
!       !  !-----!-----!------!
!       !------->Z   (lower)
!

      IF(myid.eq.0) then
         print *,'Setting up 3D processor grid'
      ENDIF

!---- X-COORD, + DIRECTION
      shift_coord = 0
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIleft  = source
      MPIright  = dest
!---- Y-COORD, + DIRECTION
      shift_coord = 1
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIlower  = source
      MPIupper = dest
!---- Z-COORD, + DIRECTION
      shift_coord = 2
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIbelow = source
      MPIabove = dest

!-- y-z plane
      MPIupperabove = MPIupper+1
      MPIupperbelow = MPIupper-1
      MPIlowerabove = MPIlower+1
      MPIlowerbelow = MPIlower-1
!-- y-x plane
      MPIlowerright = MPIlower+proc_dims(2)*proc_dims(3)
      MPIlowerleft  = MPIlower-proc_dims(2)*proc_dims(3)
      MPIupperright = MPIupper+proc_dims(2)*proc_dims(3)
      MPIupperleft  = MPIupper-proc_dims(2)*proc_dims(3)
!-- z-x plane
      MPIaboveright = MPIright+1
      MPIbelowright = MPIright-1
      MPIaboveleft = MPIleft+1
      MPIbelowleft = MPIleft-1
!-- right corners
      MPIupperaboveright = MPIupperright+1
      MPIupperbelowright = MPIupperright-1
      MPIloweraboveright = MPIlowerright+1
      MPIlowerbelowright = MPIlowerright-1
!-- left corners
      MPIupperaboveleft = MPIupperleft+1
      MPIupperbelowleft = MPIupperleft-1
      MPIloweraboveleft = MPIlowerleft+1
      MPIlowerbelowleft = MPIlowerleft-1
!-- zero the pole designations as a start
      MPIbelow_pole = MPI_PROC_NULL
      MPIabove_pole = MPI_PROC_NULL
!-- x-boundaries
      if(.not.periods(1).and.proc_coords(1).eq.0) then
         MPIleft  = MPI_PROC_NULL
         MPIlowerleft = MPI_PROC_NULL
         MPIupperleft = MPI_PROC_NULL
         MPIaboveleft = MPI_PROC_NULL
         MPIbelowleft = MPI_PROC_NULL
         MPIupperaboveleft = MPI_PROC_NULL
         MPIupperbelowleft = MPI_PROC_NULL
         MPIloweraboveleft = MPI_PROC_NULL
         MPIlowerbelowleft = MPI_PROC_NULL
      endif
      if(.not.periods(1).and.proc_coords(1).eq.proc_dims(1)-1) then
         MPIright = MPI_PROC_NULL
         MPIlowerright = MPI_PROC_NULL
         MPIupperright = MPI_PROC_NULL
         MPIaboveright = MPI_PROC_NULL
         MPIbelowright = MPI_PROC_NULL
         MPIupperaboveright = MPI_PROC_NULL
         MPIupperbelowright = MPI_PROC_NULL
         MPIloweraboveright = MPI_PROC_NULL
         MPIlowerbelowright = MPI_PROC_NULL
      endif

!-- y-boundaries
      if(.not.periods(2).and.proc_coords(2).eq.0) then
         MPIlower = MPI_PROC_NULL
         MPIlowerright = MPI_PROC_NULL
         MPIlowerleft = MPI_PROC_NULL
         MPIlowerabove = MPI_PROC_NULL
         MPIlowerbelow = MPI_PROC_NULL
         MPIloweraboveright = MPI_PROC_NULL
         MPIlowerbelowright = MPI_PROC_NULL
         MPIloweraboveleft = MPI_PROC_NULL
         MPIlowerbelowleft = MPI_PROC_NULL
      endif
      if(.not.periods(2).and.proc_coords(2).eq.proc_dims(2)-1) then
         MPIupper = MPI_PROC_NULL
         MPIupperright = MPI_PROC_NULL
         MPIupperleft = MPI_PROC_NULL
         MPIupperabove = MPI_PROC_NULL
         MPIupperbelow = MPI_PROC_NULL
         MPIupperaboveright = MPI_PROC_NULL
         MPIupperbelowright = MPI_PROC_NULL
         MPIupperaboveleft = MPI_PROC_NULL
         MPIupperbelowleft = MPI_PROC_NULL
      endif
!-- z-boundaries
!--   South pole, solid boundary (This is implemented even 
!      for simulations with the poles so the sendrecive calls
!      work for all the other guys 
      if(.not.periods(3).and.proc_coords(3).eq.0) then
         MPIbelow = MPI_PROC_NULL
         MPIupperbelow = MPI_PROC_NULL
         MPIlowerbelow = MPI_PROC_NULL
         MPIbelowright = MPI_PROC_NULL
         MPIbelowleft = MPI_PROC_NULL
         MPIupperbelowright = MPI_PROC_NULL
         MPIlowerbelowright = MPI_PROC_NULL
         MPIupperbelowleft = MPI_PROC_NULL
         MPIlowerbelowleft = MPI_PROC_NULL
      endif
!--   North pole, solid boundary (This is implemented even 
!      for simulations with the poles so the sendrecive calls
!      work for all the other guys 
      if(.not.periods(3).and.proc_coords(3).eq.proc_dims(3)-1) then
         MPIabove = MPI_PROC_NULL
         MPIupperabove = MPI_PROC_NULL
         MPIlowerabove = MPI_PROC_NULL
         MPIaboveright = MPI_PROC_NULL
         MPIaboveleft = MPI_PROC_NULL
         MPIupperaboveright = MPI_PROC_NULL
         MPIloweraboveright = MPI_PROC_NULL
         MPIupperaboveleft = MPI_PROC_NULL
         MPIloweraboveleft = MPI_PROC_NULL
      endif
!--   South pole, include pole
      if(poles.and.proc_coords(3).eq.0) then
         if(myid.lt.numprocs/2) then
            MPIbelow_pole = myid + (numprocs/2)
         else
            MPIbelow_pole = myid - (numprocs/2)
         endif
      endif
!--   North pole, include pole
      if(poles.and.proc_coords(3).eq.proc_dims(3)-1) then
         if(myid.lt.numprocs/2) then
            MPIabove_pole = myid + (numprocs/2)
         else
            MPIabove_pole = myid - (numprocs/2)
         endif
      endif
!-----------------
!
!--- For debugging purposes ----
c      print *,'1',myid,MPIlowerbelow,MPIbelow,MPIupperbelow
c      print *,'2',myid,MPIupper,MPIupperabove,MPIabove
c      print *,'3',myid,MPIlowerabove,MPIlower,MPIright
c      if(myid.eq.12) then
c         print *,'0',myid,proc_coords(1),proc_coords(2),proc_coords(3)
c         print *,'1',MPIlower,MPIupper
c         print *,'2',MPIleft,MPIright
c         print *,'3',MPIbelow,MPIabove
c         print *,'4',MPIupperabove,MPIupperbelow
c         print *,'5',MPIlowerabove,MPIlowerbelow
c         print *,'6',MPIaboveright,MPIbelowright
c         print *,'7',MPIaboveleft,MPIbelowleft 
c         print *,'8',MPIupperright,MPIupperleft
c         print *,'9',MPIlowerright,MPIlowerleft
c         print *,'10',MPIupperaboveright,MPIupperbelowright
c         print *,'11',MPIloweraboveright,MPIlowerbelowright
c         print *,'12',MPIupperaboveleft,MPIupperbelowleft
c         print *,'13',MPIloweraboveleft,MPIlowerbelowleft
c      endif
c      print *,myid,MPIright,MPIleft
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      call clean_stop
      RETURN
      END SUBROUTINE SETUP3DMPIGRID


      SUBROUTINE SETUP_AXISYMMETRIC
      USE input_init
      USE mpi_grid_init
      USE mpi_var_init
      USE fluid_var_init
      IMPLICIT NONE
      integer ::shift_coord,shift_dir
      include "mpif.h"
!                    (above)
!          !-----!-----!------!
!          !     !     !      !
!          !  2  !  5  !  8   !
!          !-----!-----!------!
!       Z  !     !     !      !
!(left) ^  !  1  !  4  !  7   ! (right)
!       !  !-----!-----!------!
!       !  !     !     !      !
!       !  !  0  !  3  !  6   !
!       !  !-----!-----!------!
!       !------->X   (below)
!
!--   Anything with upper or lower should be = MPI_PROC_NULL
!
      IF(myid.eq.1) then
         print *,'Setting up AXISYMETRIC processor grid'
      ENDIF

!---- X-COORD, + DIRECTION
      shift_coord = 0
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIleft  = source
      MPIright  = dest
!---- Y-COORD, + DIRECTION
      MPIlower = MPI_PROC_NULL
      MPIupper = MPI_PROC_NULL
!---- Z-COORD, + DIRECTION
      shift_coord = 2
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIbelow = source
      MPIabove = dest

!-- y-z plane
      MPIupperabove = MPI_PROC_NULL
      MPIupperbelow = MPI_PROC_NULL
      MPIlowerabove = MPI_PROC_NULL
      MPIlowerbelow = MPI_PROC_NULL
!-- y-x plane
      MPIlowerright = MPI_PROC_NULL
      MPIlowerleft  = MPI_PROC_NULL
      MPIupperright = MPI_PROC_NULL
      MPIupperleft  = MPI_PROC_NULL
!-- z-x plane
      MPIaboveright = MPIright+1
      MPIbelowright = MPIright-1
      MPIaboveleft  = MPIleft+1
      MPIbelowleft  = MPIleft-1
!-- right corners
      MPIupperaboveright = MPI_PROC_NULL
      MPIupperbelowright = MPI_PROC_NULL
      MPIloweraboveright = MPI_PROC_NULL
      MPIlowerbelowright = MPI_PROC_NULL
!-- left corners
      MPIupperaboveleft = MPI_PROC_NULL
      MPIupperbelowleft = MPI_PROC_NULL
      MPIloweraboveleft = MPI_PROC_NULL
      MPIlowerbelowleft = MPI_PROC_NULL

!-- x-boundaries
      if(proc_coords(1).eq.0) then
         MPIleft  = MPI_PROC_NULL
         MPIlowerleft = MPI_PROC_NULL
         MPIupperleft = MPI_PROC_NULL
         MPIaboveleft = MPI_PROC_NULL
         MPIbelowleft = MPI_PROC_NULL
         MPIupperaboveleft = MPI_PROC_NULL
         MPIupperbelowleft = MPI_PROC_NULL
         MPIloweraboveleft = MPI_PROC_NULL
         MPIlowerbelowleft = MPI_PROC_NULL
      endif
      if(proc_coords(1).eq.proc_dims(1)-1) then
         MPIright = MPI_PROC_NULL
         MPIlowerright = MPI_PROC_NULL
         MPIupperright = MPI_PROC_NULL
         MPIaboveright = MPI_PROC_NULL
         MPIbelowright = MPI_PROC_NULL
         MPIupperaboveright = MPI_PROC_NULL
         MPIupperbelowright = MPI_PROC_NULL
         MPIloweraboveright = MPI_PROC_NULL
         MPIlowerbelowright = MPI_PROC_NULL
      endif

!-- z-boundaries
      if(proc_coords(3).eq.0) then
         MPIbelow = MPI_PROC_NULL
         MPIupperbelow = MPI_PROC_NULL
         MPIlowerbelow = MPI_PROC_NULL
         MPIbelowright = MPI_PROC_NULL
         MPIbelowleft = MPI_PROC_NULL
         MPIupperbelowright = MPI_PROC_NULL
         MPIlowerbelowright = MPI_PROC_NULL
         MPIupperbelowleft = MPI_PROC_NULL
         MPIlowerbelowleft = MPI_PROC_NULL
      endif
      if(proc_coords(3).eq.proc_dims(3)-1) then
         MPIabove = MPI_PROC_NULL
         MPIupperabove = MPI_PROC_NULL
         MPIlowerabove = MPI_PROC_NULL
         MPIaboveright = MPI_PROC_NULL
         MPIaboveleft = MPI_PROC_NULL
         MPIupperaboveright = MPI_PROC_NULL
         MPIloweraboveright = MPI_PROC_NULL
         MPIupperaboveleft = MPI_PROC_NULL
         MPIloweraboveleft = MPI_PROC_NULL
      endif
!--- For debugging purposes ----
c      print *,'1',myid,MPIlowerbelow,MPIbelow,MPIupperbelow
c      print *,'2',myid,MPIupper,MPIupperabove,MPIabove
c      print *,'3',myid,MPIlowerabove,MPIlower,MPIright
c      if(myid.eq.13) then
c         print *,'0',myid,proc_coords(1),proc_coords(2),proc_coords(3)
c         print *,'1',MPIlower,MPIupper
c         print *,'2',MPIleft,MPIright
c         print *,'3',MPIbelow,MPIabove
c         print *,'4',MPIupperabove,MPIupperbelow
c         print *,'5',MPIlowerabove,MPIlowerbelow
c         print *,'6',MPIaboveright,MPIbelowright
c         print *,'7',MPIaboveleft,MPIbelowleft 
c         print *,'8',MPIupperright,MPIupperleft
c         print *,'9',MPIlowerright,MPIlowerleft
c         print *,'10',MPIupperaboveright,MPIupperbelowright
c         print *,'11',MPIloweraboveright,MPIlowerbelowright
c         print *,'12',MPIupperaboveleft,MPIupperbelowleft
c         print *,'13',MPIloweraboveleft,MPIlowerbelowleft
c      endif
c      print *,myid,MPIright,MPIleft
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      call clean_stop
      RETURN
      END SUBROUTINE SETUP_AXISYMMETRIC

      SUBROUTINE PLANETSHIFT(f,scalar_num)
!--- SHIFT variable for standard planet case
!--- This routine is for use with hot-jupiter (modtyp=2) simulations.
!--- They are different due to the mpi setup. All mpi designations that
!--- include right and left are =MPI_PROC_NULL. This routine just passes 
!--- in the y and z directions. Consider using this for the 2d shift too...
!    f=1-5 => T,RH,V,G,H
!    f=6   => ER
!    f=7   => not availble
!    f=8   => delta_ER
!    f=9   => pass_sc
!    f=10   => DIFF_PTR
!    f=11   => GD
!-------------------
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_grid_init
      USE dif_var_init
      USE scalar_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: f,k,stat(MPI_STATUS_SIZE),scalar_num
      integer :: low_indx,add_indx
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: fluid_point
      if(f.eq.1) fluid_point => T
      if(f.eq.2) fluid_point => RH
      if(f.eq.3) fluid_point => V
      if(f.eq.4) fluid_point => G
      if(f.eq.5) fluid_point => H
      if(f.eq.6) fluid_point => ER
      if(f.eq.7) then
         print *, 'f=7 (for ER) should never use planetshift'
         call clean_stop
      endif
      if(f.eq.8) fluid_point => delta_ER
c--now this is not used.. cause it won't let me point a 3d pointer array to a 4d array
c  instead, I set pass_sc(:,:,:) = scalar(n,:,:,:) in ppscalar.f 
c      if(f.eq.9) fluid_point => scalar(scalar_num,:,:,:)
      if(f.eq.9) fluid_point => pass_sc
      if(f.eq.10) fluid_point => DIFF_PTR
      if(f.eq.11) fluid_point => GD
      if(f.ge.12) then          !-note 11 is used in north/south_pole_shift
         print *,'not avalible:planetshift'
         call clean_stop
      endif
!--set array limits for radial velocity
!   (V is special because of 2d domain decomp.)
      if(f.eq.3) then !- for V(0:locNX+2,-1:locNY+2,-1:locNZ+2)
         low_indx = 0
         add_indx = 3
      else !- for everything else
         low_indx = -1
         add_indx = 4
      endif
!--- from lowerbelow
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,locNY,locNZ),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %     fluid_point(low_indx:locNX+2,0,0),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.4.and.f.ne.5.and.f.ne.11) then !- G/GD(*,-1,-1) and H(*,-1,-1) don't exist
         call MPI_SENDRECV(
     %        fluid_point(low_indx:locNX+2,locNY-1,locNZ-1)
     %        ,locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %        fluid_point(low_indx:locNX+2,-1,-1),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %        COMM_CART,stat,ierr)
      endif
      if(f.ne.5) then           !- H(*,0,-1) doesn't exist
         call MPI_SENDRECV(
     %        fluid_point(low_indx:locNX+2,locNY,locNZ-1),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %        fluid_point(low_indx:locNX+2,0,-1),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %        COMM_CART,stat,ierr)
      endif
      if(f.ne.4.and.f.ne.11) then !- G/GD(*,-1,0) doesn't exist
         call MPI_SENDRECV(
     %        fluid_point(low_indx:locNX+2,locNY-1,locNZ),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %        fluid_point(low_indx:locNX+2,-1,0),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %        COMM_CART,stat,ierr)      
      endif
!--- from lower
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,locNY,1:locNZ),
     %     (locNX+add_indx)*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     fluid_point(low_indx:locNX+2,0,1:locNZ),
     %     (locNX+add_indx)*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.4.and.f.ne.11) then !- G/GD(*,-1,1:locNZ) doesn't exist
         call MPI_SENDRECV(fluid_point(low_indx:locNX+2,locNY-1,1:locNZ)
     %        ,(locNX+add_indx)*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %        fluid_point(low_indx:locNX+2,-1,1:locNZ),
     %        (locNX+add_indx)*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %        COMM_CART,stat,ierr)
      endif
!--- from lowerabove
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,locNY,1),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %     fluid_point(low_indx:locNX+2,0,locNZ+1),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.4.and.f.ne.11) then !- G/GD(*,-1,locNZ+1 and locNZ+2) don't exist
         call MPI_SENDRECV(fluid_point(low_indx:locNX+2,locNY-1,2),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        fluid_point(low_indx:locNX+2,-1,locNZ+2),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        COMM_CART,stat,ierr)
         call MPI_SENDRECV(fluid_point(low_indx:locNX+2,locNY-1,1),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        fluid_point(low_indx:locNX+2,-1,locNZ+1),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        COMM_CART,stat,ierr)
      endif
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,locNY,2),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %     fluid_point(low_indx:locNX+2,0,locNZ+2),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %     COMM_CART,stat,ierr)
!--- from above
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,1:locNY,1),
     %     (locNX+add_indx)*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     fluid_point(low_indx:locNX+2,1:locNY,locNZ+1),
     %     (locNX+add_indx)*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,1:locNY,2),
     %     (locNX+add_indx)*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     fluid_point(low_indx:locNX+2,1:locNY,locNZ+2),
     %     (locNX+add_indx)*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     COMM_CART,stat,ierr)
!--- from upperabove
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,1,1),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %     fluid_point(low_indx:locNX+2,locNY+1,locNZ+1),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,2,2),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %     fluid_point(low_indx:locNX+2,locNY+2,locNZ+2),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,2,1),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %     fluid_point(low_indx:locNX+2,locNY+2,locNZ+1),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,1,2),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %     fluid_point(low_indx:locNX+2,locNY+1,locNZ+2),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %     COMM_CART,stat,ierr)
!--- from upper
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,1,1:locNZ),
     %     (locNX+add_indx)*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     fluid_point(low_indx:locNX+2,locNY+1,1:locNZ),
     %     (locNX+add_indx)*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     COMM_CART,stat,ierr)
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,2,1:locNZ),
     %     (locNX+add_indx)*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     fluid_point(low_indx:locNX+2,locNY+2,1:locNZ),
     %     (locNX+add_indx)*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     COMM_CART,stat,ierr)
!--- from upperbelow
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,1,locNZ),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %     fluid_point(low_indx:locNX+2,locNY+1,0),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.5) then !- H(*,locNY+1/locNY+2,-1) don't exist
         call MPI_SENDRECV(fluid_point(low_indx:locNX+2,2,locNZ-1),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        fluid_point(low_indx:locNX+2,locNY+2,-1),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        COMM_CART,stat,ierr)
         call MPI_SENDRECV(fluid_point(low_indx:locNX+2,1,locNZ-1),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        fluid_point(low_indx:locNX+2,locNY+1,-1),
     %        locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        COMM_CART,stat,ierr)
      endif
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,2,locNZ),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %     fluid_point(low_indx:locNX+2,locNY+2,0),
     %     locNX+add_indx,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %     COMM_CART,stat,ierr)
!--- from below
      call MPI_SENDRECV(fluid_point(low_indx:locNX+2,1:locNY,locNZ),
     %     (locNX+add_indx)*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     fluid_point(low_indx:locNX+2,1:locNY,0),
     %     (locNX+add_indx)*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.5) then !- H(*,1:locNY,-1) doesn't exist
         call MPI_SENDRECV(fluid_point(low_indx:locNX+2,1:locNY,locNZ-1)
     %        ,(locNX+add_indx)*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %        fluid_point(low_indx:locNX+2,1:locNY,-1),
     %        (locNX+add_indx)*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %        COMM_CART,stat,ierr)
      endif
!-----DONE SHIFTING
!--- RELEASE THE POINTER 
      NULLIFY(fluid_point)
      RETURN
      END SUBROUTINE PLANETSHIFT




      SUBROUTINE SETUPNEW3DMPIGRID
      USE input_init
      USE mpi_grid_init
      USE mpi_var_init
      USE fluid_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: shift_coord,shift_dir
      integer :: stat(MPI_STATUS_SIZE)
!
!     LEVEL 1 (X=0 : left)
!                    (upper)
!          !-----!-----!------!
!          !     !     !      !
!          !  6  !  7  !  8   !
!          !-----!-----!------!
!       Y  !     !     !      !
!(below)^  !  3  !  4  !  5   ! (above)
!       !  !-----!-----!------!
!       !  !     !     !      !
!       !  !  0  !  1  !  2   !
!       !  !-----!-----!------!
!       !------->Z   (lower)
!
!     LEVEL 2 (X=1 : middle)
!                    (upper)
!          !-----!-----!------!
!          !     !     !      !
!          !  15 !  16 !  17  !
!          !-----!-----!------!
!       Y  !     !     !      !
!(below)^  !  12 !  13 !  14  ! (above)
!       !  !-----!-----!------!
!       !  !     !     !      !
!       !  !  9  !  10 ! 11   !
!       !  !-----!-----!------!
!       !------->Z   (lower)
!
!     LEVEL 3 (X=2 : right)
!                    (upper)
!          !-----!-----!------!
!          !     !     !      !
!          !  24 !  25 !  26  !
!          !-----!-----!------!
!       Y  !     !     !      !
!(below)^  !  21 !  22 !  23  ! (above)
!       !  !-----!-----!------!
!       !  !     !     !      !
!       !  !  18 !  19 !  20  !
!       !  !-----!-----!------!
!       !------->Z   (lower)
!

      IF(myid.eq.1) then
         print *,'Setting up NEW 3D processor grid'
         
      ENDIF

!---- X-COORD, + DIRECTION
      shift_coord = 0
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIleft  = source
      MPIright  = dest
!---- Y-COORD, + DIRECTION
      shift_coord = 1
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIlower  = source
      MPIupper = dest
!---- Z-COORD, + DIRECTION
      shift_coord = 2
      shift_dir = 1
      call MPI_CART_SHIFT(COMM_CART,shift_coord,shift_dir,source,
     %     dest,ierr)
      MPIbelow = source
      MPIabove = dest

!--now pass this info around to fill other variables.
!--FIRST set all to MPI_PROC_NULL, for those that are never re-set
      MPIupperabove      = MPI_PROC_NULL
      MPIupperbelow      = MPI_PROC_NULL
      MPIlowerabove      = MPI_PROC_NULL
      MPIlowerbelow      = MPI_PROC_NULL
      MPIaboveright      = MPI_PROC_NULL
      MPIbelowright      = MPI_PROC_NULL
      MPIaboveleft       = MPI_PROC_NULL
      MPIbelowleft       = MPI_PROC_NULL
      MPIupperright      = MPI_PROC_NULL
      MPIupperleft       = MPI_PROC_NULL
      MPIlowerright      = MPI_PROC_NULL
      MPIlowerleft       = MPI_PROC_NULL
      MPIupperaboveright = MPI_PROC_NULL
      MPIupperbelowright = MPI_PROC_NULL
      MPIloweraboveright = MPI_PROC_NULL
      MPIlowerbelowright = MPI_PROC_NULL
      MPIupperaboveleft  = MPI_PROC_NULL
      MPIupperbelowleft  = MPI_PROC_NULL
      MPIloweraboveleft  = MPI_PROC_NULL
      MPIlowerbelowleft  = MPI_PROC_NULL

!--within your x-plane
!--set MPIlowerbelow
      call MPI_SENDRECV(MPIbelow,1,MPI_INTEGER,MPIupper,1,
     %     MPIlowerbelow,1,MPI_INTEGER,MPIlower,1,COMM_CART,stat,ierr)
!--set MPIlowerabove
      call MPI_SENDRECV(MPIabove,1,MPI_INTEGER,MPIupper,1,
     %     MPIlowerabove,1,MPI_INTEGER,MPIlower,1,COMM_CART,stat,ierr)
!--set MPIupperbelow
      call MPI_SENDRECV(MPIbelow,1,MPI_INTEGER,MPIlower,1,
     %     MPIupperbelow,1,MPI_INTEGER,MPIupper,1,COMM_CART,stat,ierr)
!--set MPIupperabove
      call MPI_SENDRECV(MPIabove,1,MPI_INTEGER,MPIlower,1,
     %     MPIupperabove,1,MPI_INTEGER,MPIupper,1,COMM_CART,stat,ierr)
!---from the plane to the right (increasing proc_coord(1))
!--set MPIaboveright
      call MPI_SENDRECV(MPIabove,1,MPI_INTEGER,MPIleft,1,
     %     MPIaboveright,1,MPI_INTEGER,MPIright,1,COMM_CART,stat,ierr)
!--set MPIbelowright
      call MPI_SENDRECV(MPIbelow,1,MPI_INTEGER,MPIleft,1,
     %     MPIbelowright,1,MPI_INTEGER,MPIright,1,COMM_CART,stat,ierr)
!--set MPIupperright
      call MPI_SENDRECV(MPIupper,1,MPI_INTEGER,MPIleft,1,
     %     MPIupperright,1,MPI_INTEGER,MPIright,1,COMM_CART,stat,ierr)
!--set MPIlowerright
      call MPI_SENDRECV(MPIlower,1,MPI_INTEGER,MPIleft,1,
     %     MPIlowerright,1,MPI_INTEGER,MPIright,1,COMM_CART,stat,ierr)
!--set MPIupperaboveright
      call MPI_SENDRECV(MPIupperabove,1,MPI_INTEGER,MPIleft,1,
     %     MPIupperaboveright,1,MPI_INTEGER,MPIright,1,COMM_CART,
     %     stat,ierr)
!--set MPIupperbelowright
      call MPI_SENDRECV(MPIupperbelow,1,MPI_INTEGER,MPIleft,1,
     %     MPIupperbelowright,1,MPI_INTEGER,MPIright,1,COMM_CART,
     %     stat,ierr)
!--set MPIloweraboveright
      call MPI_SENDRECV(MPIlowerabove,1,MPI_INTEGER,MPIleft,1,
     %     MPIloweraboveright,1,MPI_INTEGER,MPIright,1,COMM_CART,
     %     stat,ierr)
!--set MPIlowerbelowright
      call MPI_SENDRECV(MPIlowerbelow,1,MPI_INTEGER,MPIleft,1,
     %     MPIlowerbelowright,1,MPI_INTEGER,MPIright,1,COMM_CART,
     %     stat,ierr)
!---from the plane to the left (decreasing proc_coord(1))
!--set MPIaboveleft
      call MPI_SENDRECV(MPIabove,1,MPI_INTEGER,MPIright,1,
     %     MPIaboveleft,1,MPI_INTEGER,MPIleft,1,COMM_CART,stat,ierr)
!--set MPIbelowleft
      call MPI_SENDRECV(MPIbelow,1,MPI_INTEGER,MPIright,1,
     %     MPIbelowleft,1,MPI_INTEGER,MPIleft,1,COMM_CART,stat,ierr)
!--set MPIupperleft
      call MPI_SENDRECV(MPIupper,1,MPI_INTEGER,MPIright,1,
     %     MPIupperleft,1,MPI_INTEGER,MPIleft,1,COMM_CART,stat,ierr)
!--set MPIlowerleft
      call MPI_SENDRECV(MPIlower,1,MPI_INTEGER,MPIright,1,
     %     MPIlowerleft,1,MPI_INTEGER,MPIleft,1,COMM_CART,stat,ierr)
!--set MPIupperaboveleft
      call MPI_SENDRECV(MPIupperabove,1,MPI_INTEGER,MPIright,1,
     %     MPIupperaboveleft,1,MPI_INTEGER,MPIleft,1,COMM_CART,
     %     stat,ierr)
!--set MPIupperbelowleft
      call MPI_SENDRECV(MPIupperbelow,1,MPI_INTEGER,MPIright,1,
     %     MPIupperbelowleft,1,MPI_INTEGER,MPIleft,1,COMM_CART,
     %     stat,ierr)
!--set MPIloweraboveleft
      call MPI_SENDRECV(MPIlowerabove,1,MPI_INTEGER,MPIright,1,
     %     MPIloweraboveleft,1,MPI_INTEGER,MPIleft,1,COMM_CART,
     %     stat,ierr)
!--set MPIlowerbelowleft
      call MPI_SENDRECV(MPIlowerbelow,1,MPI_INTEGER,MPIright,1,
     %     MPIlowerbelowleft,1,MPI_INTEGER,MPIleft,1,COMM_CART,
     %     stat,ierr)
c      if(myid.eq.8) then
c         print *,'0',myid,proc_coords(1),proc_coords(2),proc_coords(3)
c         print *,'1',MPIlower,MPIupper
c         print *,'2',MPIleft,MPIright
c         print *,'3',MPIbelow,MPIabove
c         print *,'4',MPIupperabove,MPIupperbelow
c         print *,'5',MPIlowerabove,MPIlowerbelow
c         print *,'6',MPIaboveright,MPIbelowright
c         print *,'7',MPIaboveleft,MPIbelowleft 
c         print *,'8',MPIupperright,MPIupperleft
c         print *,'9',MPIlowerright,MPIlowerleft
c         print *,'10',MPIupperaboveright,MPIupperbelowright
c         print *,'11',MPIloweraboveright,MPIlowerbelowright
c         print *,'12',MPIupperaboveleft,MPIupperbelowleft
c         print *,'13',MPIloweraboveleft,MPIlowerbelowleft
c      endif
      RETURN
      END SUBROUTINE SETUPNEW3DMPIGRID


      SUBROUTINE NEWSHIFTVAR3D(f)
!    f=1-5 => T,RH,V,G,H
!    f=6   => ER
!    f=7   => ER w/out  -1 and N+2
!-------------------
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_grid_init
      USE dif_var_init
      USE scalar_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: f,k,stat(MPI_STATUS_SIZE)
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: fluid_point
      if(f.eq.1) fluid_point => T
      if(f.eq.2) fluid_point => RH
      if(f.eq.3) fluid_point => V
      if(f.eq.4) fluid_point => G
      if(f.eq.5) fluid_point => H
      if(f.eq.6) fluid_point => ER
      if(f.eq.7) fluid_point => ER !(don't pass for:i=-1,locNX+2,j=-1,locNY+2)
      if(f.eq.8) fluid_point => delta_ER
      if(f.eq.9) fluid_point => pass_sc
      if(f.eq.10) fluid_point => DIFF_PTR
      if(f.ge.11) then
         print *,'not yet avalible:NEWSHIFTVAR3D'
         call clean_stop
      endif
!---- X-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(locNX,1:locNY,1:locNZ),
     %     locNY*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %     fluid_point(0,1:locNY,1:locNZ),locNY*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
      if(f.ne.3.and.f.ne.7) THEN
!V(-1,j,k) doesn't exist: don't pass for SOR
         call MPI_SENDRECV(fluid_point(locNX-1,1:locNY,1:locNZ),
     %        locNY*locNZ,MPI_DOUBLE_PRECISION,MPIright,1,
     %        fluid_point(-1,1:locNY,1:locNZ),locNY*locNZ,
     %        MPI_DOUBLE_PRECISION,MPIleft,1,COMM_CART,stat,ierr)
      endif
!---- X-COORD, - DIRECTION 
      call MPI_SENDRECV(fluid_point(1,1:locNY,1:locNZ),
     %     locNY*locNZ,MPI_DOUBLE_PRECISION,MPIleft,1,
     %     fluid_point(locNX+1,1:locNY,1:locNZ),locNY*locNZ,
     %     MPI_DOUBLE_PRECISION,MPIright,1,COMM_CART,stat,ierr)
      if(f.ne.7) then !Don't pass for SOR
         call MPI_SENDRECV(fluid_point(2,1:locNY,1:locNZ),
     %        locNY*locNZ,MPI_DOUBLE_PRECISION,MPIleft,1,
     %        fluid_point(locNX+2,1:locNY,1:locNZ),locNY*locNZ,
     %        MPI_DOUBLE_PRECISION,MPIright,1,COMM_CART,stat,ierr)
      endif
!---- Y-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,locNY,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     fluid_point(1:locNX,0,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.4.and.f.ne.7) then 
!G(i,-1,k) doesn't exist: don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,locNY-1,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %        fluid_point(1:locNX,-1,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %        COMM_CART,stat,ierr)
      endif
!---- Y-COORD, - DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,1,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %     fluid_point(1:locNX,locNY+1,1:locNZ),
     %     locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.7) then !Don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,2,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIlower,1,
     %        fluid_point(1:locNX,locNY+2,1:locNZ),
     %        locNX*locNZ,MPI_DOUBLE_PRECISION,MPIupper,1,
     %        COMM_CART,stat,ierr)
      endif
!---- Z-COORD, + DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,1:locNY,locNZ),
     %     locNX*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     fluid_point(1:locNX,1:locNY,0),
     %     locNX*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.5.and.f.ne.7) then !H(i,j,-1) doesn't exist, don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,1:locNY,locNZ-1),
     %        locNX*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %        fluid_point(1:locNX,1:locNY,-1),
     %        locNX*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %        COMM_CART,stat,ierr)
      endif
!---- Z-COORD, - DIRECTION
      call MPI_SENDRECV(fluid_point(1:locNX,1:locNY,1),
     %     locNX*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %     fluid_point(1:locNX,1:locNY,locNZ+1),
     %     locNX*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.7) then !Don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,1:locNY,2),
     %        locNX*locNY,MPI_DOUBLE_PRECISION,MPIbelow,1,
     %        fluid_point(1:locNX,1:locNY,locNZ+2),
     %        locNX*locNY,MPI_DOUBLE_PRECISION,MPIabove,1,
     %        COMM_CART,stat,ierr)
      endif
!----------------------
!--------y-z plane
!--------- from lowerbelow
      call MPI_SENDRECV(fluid_point(1:locNX,locNY,locNZ),
     %     locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %     fluid_point(1:locNX,0,0),
     %     locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.4.and.f.ne.7) then !-Don't pass for G or SOR 
         call MPI_SENDRECV(fluid_point(1:locNX,locNY-1,locNZ),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %        fluid_point(1:locNX,-1,0),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %        COMM_CART,stat,ierr)
      endif
      if(f.ne.5.and.f.ne.7) then !-Don't pass for H or SOR
         call MPI_SENDRECV(fluid_point(1:locNX,locNY,locNZ-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %        fluid_point(1:locNX,0,-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %        COMM_CART,stat,ierr)
      endif
      if(f.ne.4.and.f.ne.5.and.f.ne.7) then !-Don't pass for G, H, or SOR
         call MPI_SENDRECV(fluid_point(1:locNX,locNY-1,locNZ-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %        fluid_point(1:locNX,-1,-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %        COMM_CART,stat,ierr)
      endif
!--------- from upperabove
      call MPI_SENDRECV(fluid_point(1:locNX,1,1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %     fluid_point(1:locNX,locNY+1,locNZ+1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.7) then !-Don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,2,1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %        fluid_point(1:locNX,locNY+2,locNZ+1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %        COMM_CART,stat,ierr)
         call MPI_SENDRECV(fluid_point(1:locNX,1,2),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %        fluid_point(1:locNX,locNY+1,locNZ+2),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %        COMM_CART,stat,ierr)
         call MPI_SENDRECV(fluid_point(1:locNX,2,2),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerbelow,1,
     %        fluid_point(1:locNX,locNY+2,locNZ+2),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperabove,1,
     %        COMM_CART,stat,ierr)
      endif
!--------- from upperbelow
      call MPI_SENDRECV(fluid_point(1:locNX,1,locNZ),
     %     locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %     fluid_point(1:locNX,locNY+1,0),
     %     locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.5.and.f.ne.7) then !-Don't pass for H or SOR
         call MPI_SENDRECV(fluid_point(1:locNX,1,locNZ-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        fluid_point(1:locNX,locNY+1,-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        COMM_CART,stat,ierr)
         call MPI_SENDRECV(fluid_point(1:locNX,2,locNZ-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        fluid_point(1:locNX,locNY+2,-1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        COMM_CART,stat,ierr)
      endif
      if(f.ne.7) then           !-Don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,2,locNZ),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        fluid_point(1:locNX,locNY+2,0),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        COMM_CART,stat,ierr)
      endif
!--------- from lowerabove
      call MPI_SENDRECV(fluid_point(1:locNX,locNY,1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %     fluid_point(1:locNX,0,locNZ+1),
     %     locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %     COMM_CART,stat,ierr)
      if(f.ne.4.and.f.ne.7) then !-Don't pass for G or SOR
         call MPI_SENDRECV(fluid_point(1:locNX,locNY-1,1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        fluid_point(1:locNX,-1,locNZ+1),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        COMM_CART,stat,ierr)
         call MPI_SENDRECV(fluid_point(1:locNX,locNY-1,2),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        fluid_point(1:locNX,-1,locNZ+2),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        COMM_CART,stat,ierr)
      endif
      if(f.ne.7) then           !-Don't pass for SOR
         call MPI_SENDRECV(fluid_point(1:locNX,locNY,2),
     %        locNX,MPI_DOUBLE_PRECISION,MPIupperbelow,1,
     %        fluid_point(1:locNX,0,locNZ+2),
     %        locNX,MPI_DOUBLE_PRECISION,MPIlowerabove,1,
     %        COMM_CART,stat,ierr)
      endif
!-
!--info FROM plane of processors to the left (smaller proc_coord(1))
!------------ from upperleft
         call MPI_SENDRECV(fluid_point(locNX,1,1:locNZ),
     %        locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %        fluid_point(0,locNY+1,1:locNZ),
     %        locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.7) then !Don't pass for SOR
            call MPI_SENDRECV(fluid_point(locNX,2,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %           fluid_point(0,locNY+2,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.3.and.f.ne.7) then !Don't pass for V or SOR
            call MPI_SENDRECV(fluid_point(locNX-1,1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %           fluid_point(-1,locNY+1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(locNX-1,2,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %           fluid_point(-1,locNY+2,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %           COMM_CART,stat,ierr)
         endif
!------------ from lowerleft
         call MPI_SENDRECV(fluid_point(locNX,locNY,1:locNZ),
     %        locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %        fluid_point(0,0,1:locNZ),
     %        locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.4.and.f.ne.7) then !-Dont pass for G or SOR
            call MPI_SENDRECV(fluid_point(locNX,locNY-1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %           fluid_point(0,-1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.3.and.f.ne.7) then !-Dont pass for V or SOR
            call MPI_SENDRECV(fluid_point(locNX-1,locNY,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %           fluid_point(-1,0,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.3.and.f.ne.4.and.f.ne.7) then !-Dont pass for V,G, or SOR
            call MPI_SENDRECV(fluid_point(locNX-1,locNY-1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %           fluid_point(-1,-1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %           COMM_CART,stat,ierr)
         endif
!------------ from aboveleft
         call MPI_SENDRECV(fluid_point(locNX,1:locNY,1),
     %        locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %        fluid_point(0,1:locNY,locNZ+1),
     %        locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.7) then !-Dont pass for SOR
            call MPI_SENDRECV(fluid_point(locNX,1:locNY,2),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %           fluid_point(0,1:locNY,locNZ+2),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.3.and.f.ne.7) then !-Dont pass for V or SOR
            call MPI_SENDRECV(fluid_point(locNX-1,1:locNY,1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %           fluid_point(-1,1:locNY,locNZ+1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(locNX-1,1:locNY,2),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %           fluid_point(-1,1:locNY,locNZ+2),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %           COMM_CART,stat,ierr)
         endif
!------------ from belowleft
         call MPI_SENDRECV(fluid_point(locNX,1:locNY,locNZ),
     %        locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %        fluid_point(0,1:locNY,0),
     %        locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.5.and.f.ne.7) then !Dont pass for H or SOR
            call MPI_SENDRECV(fluid_point(locNX,1:locNY,locNZ-1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %           fluid_point(0,1:locNY,-1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.3.and.f.ne.7) then !Dont pass for V or SOR
            call MPI_SENDRECV(fluid_point(locNX-1,1:locNY,locNZ),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %           fluid_point(-1,1:locNY,0),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.3.and.f.ne.5.and.f.ne.7) then !Dont pass for V, H or SOR
            call MPI_SENDRECV(fluid_point(locNX-1,1:locNY,locNZ-1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %           fluid_point(-1,1:locNY,-1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %           COMM_CART,stat,ierr)
         endif
!------------ from upperaboveleft corner (single values)
         call MPI_SENDRECV(fluid_point(locNX,1,1),
     %        1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %        fluid_point(0,locNY+1,locNZ+1),
     %        1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.7) then        !Dont pass for SOR
            call MPI_SENDRECV(fluid_point(locNX,2,1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           fluid_point(0,locNY+2,locNZ+1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(locNX,1,2),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           fluid_point(0,locNY+1,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(locNX,2,2),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           fluid_point(0,locNY+2,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.3.and.f.ne.7) then        !Dont pass for v or SOR
            call MPI_SENDRECV(fluid_point(locNX-1,1,1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           fluid_point(-1,locNY+1,locNZ+1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(locNX-1,2,1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           fluid_point(-1,locNY+2,locNZ+1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(locNX-1,1,2),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           fluid_point(-1,locNY+1,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(locNX-1,2,2),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           fluid_point(-1,locNY+2,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           COMM_CART,stat,ierr)
         endif
!------------ from upperbelowleft
         call MPI_SENDRECV(fluid_point(locNX,1,locNZ),
     %        1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %        fluid_point(0,locNY+1,0),
     %        1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.7) then        !-Dont pass for SOR
            call MPI_SENDRECV(fluid_point(locNX,2,locNZ),
     %           1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %           fluid_point(0,locNY+2,0),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %           COMM_CART,stat,ierr)
            if(f.ne.5) then     !-Or H
               call MPI_SENDRECV(fluid_point(locNX,1,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %              fluid_point(0,locNY+1,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %              COMM_CART,stat,ierr)
               call MPI_SENDRECV(fluid_point(locNX,2,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %              fluid_point(0,locNY+2,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %              COMM_CART,stat,ierr)
            endif
            if(f.ne.3) then !-or V
               call MPI_SENDRECV(fluid_point(locNX-1,1,locNZ),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %              fluid_point(-1,locNY+1,0),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %              COMM_CART,stat,ierr)
               call MPI_SENDRECV(fluid_point(locNX-1,2,locNZ),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %              fluid_point(-1,locNY+2,0),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %              COMM_CART,stat,ierr)
               if(f.ne.5) then  !-Also not H
                  call MPI_SENDRECV(fluid_point(locNX-1,1,locNZ-1),
     %                 1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %                 fluid_point(-1,locNY+1,-1),
     %                 1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %                 COMM_CART,stat,ierr)
                  call MPI_SENDRECV(fluid_point(locNX-1,2,locNZ-1),
     %                 1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %                 fluid_point(-1,locNY+2,-1),
     %                 1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %                 COMM_CART,stat,ierr)
               endif
            endif
         endif
!------------ from lowerbelowleft
         call MPI_SENDRECV(fluid_point(locNX,locNY,locNZ),
     %        1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %        fluid_point(0,0,0),
     %        1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.4.and.f.ne.7) then !-Dont pass for G or SOR
            call MPI_SENDRECV(fluid_point(locNX,locNY-1,locNZ),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           fluid_point(0,-1,0),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.5.and.f.ne.7) then !-Dont pass for H or SOR
            call MPI_SENDRECV(fluid_point(locNX,locNY,locNZ-1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           fluid_point(0,0,-1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.4.and.f.ne.5.and.f.ne.7) then !-Dont pass for G,H, or SOR
            call MPI_SENDRECV(fluid_point(locNX,locNY-1,locNZ-1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           fluid_point(0,-1,-1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.3.and.f.ne.7) then !-Dont pass for V or SOR
            call MPI_SENDRECV(fluid_point(locNX-1,locNY,locNZ),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           fluid_point(-1,0,0),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           COMM_CART,stat,ierr)
            if(f.ne.4) then     !-Also not G
               call MPI_SENDRECV(fluid_point(locNX-1,locNY-1,locNZ),
     %              1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %              fluid_point(-1,-1,0),
     %              1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %              COMM_CART,stat,ierr)
            endif
            if(f.ne.5) then     !-Also not H
               call MPI_SENDRECV(fluid_point(locNX-1,locNY,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %              fluid_point(-1,0,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %              COMM_CART,stat,ierr)
            endif
            if(f.ne.4.and.f.ne.5) then !-Also not G or H
               call MPI_SENDRECV(fluid_point(locNX-1,locNY-1,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %              fluid_point(-1,-1,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %              COMM_CART,stat,ierr)
            endif
         endif
!------------ from loweraboveleft
         call MPI_SENDRECV(fluid_point(locNX,locNY,1),
     %        1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %        fluid_point(0,0,locNZ+1),1,MPI_DOUBLE_PRECISION,
     %        MPIloweraboveleft,1,COMM_CART,stat,ierr)
         if(f.ne.7.and.f.ne.4) then
            call MPI_SENDRECV(fluid_point(locNX,locNY-1,1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %           fluid_point(0,-1,locNZ+1),1,MPI_DOUBLE_PRECISION,
     %           MPIloweraboveleft,1,COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(locNX,locNY-1,2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %           fluid_point(0,-1,locNZ+2),1,MPI_DOUBLE_PRECISION,
     %           MPIloweraboveleft,1,COMM_CART,stat,ierr)
         endif
         if(f.ne.7) then            
            call MPI_SENDRECV(fluid_point(locNX,locNY,2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %           fluid_point(0,0,locNZ+2),1,MPI_DOUBLE_PRECISION,
     %           MPIloweraboveleft,1,COMM_CART,stat,ierr)
         endif
         if(f.ne.3.and.f.ne.7) then            
            call MPI_SENDRECV(fluid_point(locNX-1,locNY,1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %           fluid_point(-1,0,locNZ+1),1,MPI_DOUBLE_PRECISION,
     %           MPIloweraboveleft,1,COMM_CART,stat,ierr)
            if(f.ne.4) then
               call MPI_SENDRECV(fluid_point(locNX-1,locNY-1,1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %              fluid_point(-1,-1,locNZ+1),1,MPI_DOUBLE_PRECISION,
     %              MPIloweraboveleft,1,COMM_CART,stat,ierr)
               call MPI_SENDRECV(fluid_point(locNX-1,locNY-1,2),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %              fluid_point(-1,-1,locNZ+2),1,MPI_DOUBLE_PRECISION,
     %              MPIloweraboveleft,1,COMM_CART,stat,ierr)
            endif
            call MPI_SENDRECV(fluid_point(locNX-1,locNY,2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %           fluid_point(-1,0,locNZ+2),1,MPI_DOUBLE_PRECISION,
     %           MPIloweraboveleft,1,COMM_CART,stat,ierr)
         endif
!
!NOW DO THE SAME FOR FROM RIGHT
!--info FROM plane of processors to the right (larger proc_coord(1))
!------------ from upperright
         call MPI_SENDRECV(fluid_point(1,1,1:locNZ),
     %        locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %        fluid_point(locNX+1,locNY+1,1:locNZ),
     %        locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.7) then
            call MPI_SENDRECV(fluid_point(1,2,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %           fluid_point(locNX+1,locNY+2,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %           fluid_point(locNX+2,locNY+1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %           COMM_CART,stat,ierr)            
            call MPI_SENDRECV(fluid_point(2,2,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerleft,1,
     %           fluid_point(locNX+2,locNY+2,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperright,1,
     %           COMM_CART,stat,ierr)
         endif
!------------ from lowerright
         call MPI_SENDRECV(fluid_point(1,locNY,1:locNZ),
     %        locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %        fluid_point(locNX+1,0,1:locNZ),
     %        locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.4.and.f.ne.7) then
            call MPI_SENDRECV(fluid_point(1,locNY-1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %           fluid_point(locNX+1,-1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,locNY-1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %           fluid_point(locNX+2,-1,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.7) then
            call MPI_SENDRECV(fluid_point(2,locNY,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIupperleft,1,
     %           fluid_point(locNX+2,0,1:locNZ),
     %           locNZ,MPI_DOUBLE_PRECISION,MPIlowerright,1,
     %           COMM_CART,stat,ierr)
         endif
!------------ from aboveright
         call MPI_SENDRECV(fluid_point(1,1:locNY,1),
     %        locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %        fluid_point(locNX+1,1:locNY,locNZ+1),
     %        locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.7) then !-- Dont pass for SOR
            call MPI_SENDRECV(fluid_point(1,1:locNY,2),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %           fluid_point(locNX+1,1:locNY,locNZ+2),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,1:locNY,1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %           fluid_point(locNX+2,1:locNY,locNZ+1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,1:locNY,2),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowleft,1,
     %           fluid_point(locNX+2,1:locNY,locNZ+2),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveright,1,
     %           COMM_CART,stat,ierr)
         endif
!------------ from belowright
         call MPI_SENDRECV(fluid_point(1,1:locNY,locNZ),
     %        locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %        fluid_point(locNX+1,1:locNY,0),
     %        locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.5.and.f.ne.7) then !-- Dont pass for H or SOR
            call MPI_SENDRECV(fluid_point(1,1:locNY,locNZ-1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %           fluid_point(locNX+1,1:locNY,-1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,1:locNY,locNZ-1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %           fluid_point(locNX+2,1:locNY,-1),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.7) then !-- Dont pass for SOR
            call MPI_SENDRECV(fluid_point(2,1:locNY,locNZ),
     %           locNY,MPI_DOUBLE_PRECISION,MPIaboveleft,1,
     %           fluid_point(locNX+2,1:locNY,0),
     %           locNY,MPI_DOUBLE_PRECISION,MPIbelowright,1,
     %           COMM_CART,stat,ierr)
         endif
!- Corners.. (pass single values)
!------------ from upperaboveright
         call MPI_SENDRECV(fluid_point(1,1,1),
     %        1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %        fluid_point(locNX+1,locNY+1,locNZ+1),
     %        1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.7) then !-- Dont pass for SOR
            call MPI_SENDRECV(fluid_point(1,2,1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           fluid_point(locNX+1,locNY+2,locNZ+1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(1,1,2),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           fluid_point(locNX+1,locNY+1,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(1,2,2),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           fluid_point(locNX+1,locNY+2,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,1,1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           fluid_point(locNX+2,locNY+1,locNZ+1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,2,1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           fluid_point(locNX+2,locNY+2,locNZ+1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,1,2),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           fluid_point(locNX+2,locNY+1,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,2,2),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowleft,1,
     %           fluid_point(locNX+2,locNY+2,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveright,1,
     %           COMM_CART,stat,ierr)
         endif
!------------ from upperbelowright
         call MPI_SENDRECV(fluid_point(1,1,locNZ),
     %        1,MPI_DOUBLE_PRECISION,MPIloweraboveleft,1,
     %        fluid_point(locNX+1,locNY+1,0),
     %        1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.7) then  !-Dont pass for SOR
            call MPI_SENDRECV(fluid_point(1,2,locNZ),
     %           1,MPI_DOUBLE_PRECISION,MPIloweraboveleft,1,
     %           fluid_point(locNX+1,locNY+2,0),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %           COMM_CART,stat,ierr)
            if(f.ne.5) then     !-Also not H
               call MPI_SENDRECV(fluid_point(1,1,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveleft,1,
     %              fluid_point(locNX+1,locNY+1,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %              COMM_CART,stat,ierr)
               
               call MPI_SENDRECV(fluid_point(1,2,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveleft,1,
     %              fluid_point(locNX+1,locNY+2,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %              COMM_CART,stat,ierr)
            endif
            call MPI_SENDRECV(fluid_point(2,1,locNZ),
     %           1,MPI_DOUBLE_PRECISION,MPIloweraboveleft,1,
     %           fluid_point(locNX+2,locNY+1,0),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,2,locNZ),
     %           1,MPI_DOUBLE_PRECISION,MPIloweraboveleft,1,
     %           fluid_point(locNX+2,locNY+2,0),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %           COMM_CART,stat,ierr)
            if(f.ne.5) then     !-Also not H
               call MPI_SENDRECV(fluid_point(2,1,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveleft,1,
     %              fluid_point(locNX+2,locNY+1,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %              COMM_CART,stat,ierr)            
               call MPI_SENDRECV(fluid_point(2,2,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveleft,1,
     %              fluid_point(locNX+2,locNY+2,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowright,1,
     %              COMM_CART,stat,ierr)
            endif
         endif
!------------ from lowerbelowright
         call MPI_SENDRECV(fluid_point(1,locNY,locNZ),
     %        1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %        fluid_point(locNX+1,0,0),
     %        1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.4.and.f.ne.7) then
            call MPI_SENDRECV(fluid_point(1,locNY-1,locNZ),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           fluid_point(locNX+1,-1,0),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.5.and.f.ne.7) then
            call MPI_SENDRECV(fluid_point(1,locNY,locNZ-1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           fluid_point(locNX+1,0,-1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.4.and.f.ne.5.and.f.ne.7) then
            call MPI_SENDRECV(fluid_point(1,locNY-1,locNZ-1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           fluid_point(locNX+1,-1,-1),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.7) then
            call MPI_SENDRECV(fluid_point(2,locNY,locNZ),
     %           1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %           fluid_point(locNX+2,0,0),
     %           1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %           COMM_CART,stat,ierr)
            if(f.ne.4) then
               call MPI_SENDRECV(fluid_point(2,locNY-1,locNZ),
     %              1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %              fluid_point(locNX+2,-1,0),
     %              1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %              COMM_CART,stat,ierr)
            endif
            if(f.ne.5) then
               call MPI_SENDRECV(fluid_point(2,locNY,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %              fluid_point(locNX+2,0,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %              COMM_CART,stat,ierr)
            endif
            if(f.ne.4.and.f.ne.5) then
               call MPI_SENDRECV(fluid_point(2,locNY-1,locNZ-1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperaboveleft,1,
     %              fluid_point(locNX+2,-1,-1),
     %              1,MPI_DOUBLE_PRECISION,MPIlowerbelowright,1,
     %              COMM_CART,stat,ierr)
            endif
         endif
!------------ from loweaboveright
         call MPI_SENDRECV(fluid_point(1,locNY,1),
     %        1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %        fluid_point(locNX+1,0,locNZ+1),
     %        1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %        COMM_CART,stat,ierr)
         if(f.ne.4.and.f.ne.7) then
            call MPI_SENDRECV(fluid_point(1,locNY-1,1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %           fluid_point(locNX+1,-1,locNZ+1),
     %           1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(1,locNY-1,2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %           fluid_point(locNX+1,-1,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %           COMM_CART,stat,ierr)
         endif
         if(f.ne.7) then
            call MPI_SENDRECV(fluid_point(1,locNY,2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %           fluid_point(locNX+1,0,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %           COMM_CART,stat,ierr)
            call MPI_SENDRECV(fluid_point(2,locNY,1),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %           fluid_point(locNX+2,0,locNZ+1),
     %           1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %           COMM_CART,stat,ierr)
            if(f.ne.4) then
               call MPI_SENDRECV(fluid_point(2,locNY-1,1),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %              fluid_point(locNX+2,-1,locNZ+1),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %              COMM_CART,stat,ierr)
               call MPI_SENDRECV(fluid_point(2,locNY-1,2),
     %              1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %              fluid_point(locNX+2,-1,locNZ+2),
     %              1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %              COMM_CART,stat,ierr)
            endif
            call MPI_SENDRECV(fluid_point(2,locNY,2),
     %           1,MPI_DOUBLE_PRECISION,MPIupperbelowleft,1,
     %           fluid_point(locNX+2,0,locNZ+2),
     %           1,MPI_DOUBLE_PRECISION,MPIloweraboveright,1,
     %           COMM_CART,stat,ierr)
         endif
!-----DONE SHIFTING
!--- RELEASE THE POINTER 
      NULLIFY(fluid_point)
      RETURN
      END SUBROUTINE NEWSHIFTVAR3D


      SUBROUTINE NORTH_POLE_SHIFT(f,scalar_num)
!--- SHIFT across the north polar cap. Only utilize the edge near k=NZ and MPIabove
!    f=1-5 => T,RH,V,G,H
!    f=6   => ER
!    f=7   => not availble
!    f=8   => delta_ER
!    f=9   => pass_sc
!    f=10   => DIFF_PTR
!    f=11   => GD
!-------------------
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_grid_init
      USE dif_var_init
      USE scalar_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: f,k,stat(MPI_STATUS_SIZE),scalar_num
      integer :: low_indx,add_indx,g_low_indx,g_add_indx
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: fluid_point
      if(f.eq.1) fluid_point => T
      if(f.eq.2) fluid_point => RH
      if(f.eq.3) fluid_point => V
      if(f.eq.4) fluid_point => G
      if(f.eq.5) fluid_point => H
      if(f.eq.6) fluid_point => ER
      if(f.eq.7) then
         print *, 'f=7 (for ER) should never use north_pole_shift'
         call clean_stop
      endif
      if(f.eq.8) fluid_point => delta_ER
c      if(f.eq.9) fluid_point => scalar(scalar_num,:,:,:)
      if(f.eq.9) fluid_point => pass_sc
      if(f.eq.10) fluid_point => DIFF_PTR
      if(f.eq.11) fluid_point => GD
      if(f.ge.12) then
         print *,'not avalible:north_pole_shift'
         call clean_stop
      endif
!--set array limits for radial velocity
 !- for V(0:locNX+2,-1:locNY+2,-1:locNZ+2)
      if(f.eq.3) then
         low_indx = 0
         add_indx = 3
      else !- for everything else
         low_indx = -1
         add_indx = 4
      endif
!--set array limits for azimuthal velocity
 !--for G/GD(-1:locNX+2,0:locNY+2,-1:locNZ+2)
      if(f.eq.4.or.f.eq.11) then
         g_low_indx=0
         g_add_indx=3
      else !- for everything else
         g_low_indx=-1
         g_add_indx=4
      endif
!-everything except G and H, and GD
      if(f.ne.4.and.f.ne.5.and.f.ne.11) then
         call MPI_SENDRECV(
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ-1),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,COMM_CART,stat,ierr)
         call MPI_SENDRECV(
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ-2),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ+1),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,COMM_CART,stat,ierr)
         call MPI_SENDRECV(
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ-3),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ+2),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,COMM_CART,stat,ierr)
      endif
!-reverse sign for G, GD
c      if(f.eq.4.or.f.eq.11) then
c         call MPI_SENDRECV(
c     %        -fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ-1),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIabove_pole,1,
c     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIabove_pole,1,COMM_CART,stat,ierr)
c         call MPI_SENDRECV(
c     %        -fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ-2),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIabove_pole,1,
c     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ+1),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIabove_pole,1,COMM_CART,stat,ierr)
c         call MPI_SENDRECV(
c     %        -fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ-3),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIabove_pole,1,
c     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,locNZ+2),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIabove_pole,1,COMM_CART,stat,ierr)
c      endif
!-now do H with sign flip and defined on the xxa grid
      if(f.eq.5) then
         call MPI_SENDRECV(
     %        -fluid_point(low_indx:locNX+2,-1:locNY+2,locNZ-1),
     %        (locNX+add_indx)*(locNY+4),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,
     %        fluid_point(low_indx:locNX+2,-1:locNY+2,locNZ+1),
     %        (locNX+add_indx)*(locNY+4),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,COMM_CART,stat,ierr)
         call MPI_SENDRECV(
     %        -fluid_point(low_indx:locNX+2,-1:locNY+2,locNZ-2),
     %        (locNX+add_indx)*(locNY+4),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,
     %        fluid_point(low_indx:locNX+2,-1:locNY+2,locNZ+2),
     %        (locNX+add_indx)*(locNY+4),MPI_DOUBLE_PRECISION,
     %        MPIabove_pole,1,COMM_CART,stat,ierr)
      endif
!-----DONE SHIFTING
!--- RELEASE THE POINTER 
      NULLIFY(fluid_point)
      RETURN
      END SUBROUTINE NORTH_POLE_SHIFT

      SUBROUTINE SOUTH_POLE_SHIFT(f,scalar_num)
!--- SHIFT across the SOUTH polar cap. Only utilize the edge near k=1 and MPIbelow
!    f=1-5 => T,RH,V,G,H
!    f=6   => ER
!    f=7   => not availble
!    f=8   => delta_ER
!    f=9   => pass_sc
!    f=10   => DIFF_PTR
!    f=11   => GD
!-------------------
      USE input_init
      USE mpi_var_init
      USE fluid_var_init
      USE rad_var_init
      USE mpi_grid_init
      USE dif_var_init
      USE scalar_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: f,k,stat(MPI_STATUS_SIZE),scalar_num
      integer :: low_indx,add_indx,g_low_indx,g_add_indx
      DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: fluid_point
      if(f.eq.1) fluid_point => T
      if(f.eq.2) fluid_point => RH
      if(f.eq.3) fluid_point => V
      if(f.eq.4) fluid_point => G
      if(f.eq.5) fluid_point => H
      if(f.eq.6) fluid_point => ER
      if(f.eq.7) then
         print *, 'f=7 (for ER) should never use south_pole_shift'
         call clean_stop
      endif
      if(f.eq.8) fluid_point => delta_ER
c      if(f.eq.9) fluid_point => scalar(scalar_num,:,:,:)
      if(f.eq.9) fluid_point => pass_sc
      if(f.eq.10) fluid_point => DIFF_PTR
      if(f.eq.11) fluid_point => GD
      if(f.ge.12) then
         print *,'not avalible:south_pole_shift'
         call clean_stop
      endif
!--set array limits for radial velocity
 !- for V(0:locNX+2,-1:locNY+2,-1:locNZ+2)
      if(f.eq.3) then
         low_indx = 0
         add_indx = 3
      else !- for everything else
         low_indx = -1
         add_indx = 4
      endif
!--set array limits for azimuthal velocity
 !--for G/GD(-1:locNX+2,0:locNY+2,-1:locNZ+2)
      if(f.eq.4.or.f.eq.11) then
         g_low_indx=0
         g_add_indx=3
      else !- for everything else
         g_low_indx=-1
         g_add_indx=4
      endif

!-- from below
      if(f.ne.4.and.f.ne.5.and.f.ne.11) then !-everything except G and H,GD
         call MPI_SENDRECV(
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,1),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIbelow_pole,1,
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,0),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIbelow_pole,1,COMM_CART,stat,ierr)
         call MPI_SENDRECV(
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,2),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIbelow_pole,1,
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,-1),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIbelow_pole,1,COMM_CART,stat,ierr)
      endif
c      if(f.eq.4.or.f.eq.11) then !- G, GD gets a sign flip
c         call MPI_SENDRECV(
c     %        -fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,1),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIbelow_pole,1,
c     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,0),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIbelow_pole,1,COMM_CART,stat,ierr)
c         call MPI_SENDRECV(
c     %        -fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,2),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIbelow_pole,1,
c     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,-1),
c     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
c     %        MPIbelow_pole,1,COMM_CART,stat,ierr)
c      endif
      if(f.eq.5) then !-now do H with sign flip. There is no H(*,*,-1)
         call MPI_SENDRECV(
     %        -fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,2),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIbelow_pole,1,
     %        fluid_point(low_indx:locNX+2,g_low_indx:locNY+2,0),
     %        (locNX+add_indx)*(locNY+g_add_indx),MPI_DOUBLE_PRECISION,
     %        MPIbelow_pole,1,COMM_CART,stat,ierr)
      endif
!-----DONE SHIFTING
!--- RELEASE THE POINTER 
      NULLIFY(fluid_point)
      RETURN
      END SUBROUTINE SOUTH_POLE_SHIFT
