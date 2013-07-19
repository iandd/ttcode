!----this function is to calculate the CO fraction in equlibrium state
      double precision function XCOeq(i,j,k)
      USE global_constants
      USE fluid_var_init
      USE scalar_var_init
      IMPLICIT NONE
      integer :: i,j,k
!---a few constants
      double precision c1,c2,XH2 
!---Gibbs Free Energy
      double precision GCO             
      double precision GCH4            
      double precision GH2O
!---a few mediators            
      double precision f,a,b,c
!---chemical equlibrium constant          
      double precision K3eq
!----initialize the value of constants
      c1=5.91e-4
      c2=0.00104
      XH2=0.83
!----calculate eqilibrium state of CO fraction 
      GCO=-0.088787*T(i,j,k)-111.137747
      GCH4=0.105*T(i,j,k)-84.68434
      GH2O=0.0537546*T(i,j,k)-245.8666
      K3eq=exp(-(GCO-GCH4-GH2O)/(RGAS*1e-10*T(i,j,k)))
      f=K3eq/(PG(i,j,k)*PG(i,j,k)*XH2*XH2*XH2)*1e+12
      a=f
      b=-f*(c1+c2)-1
      c=f*c1*c2
      XCOeq=(-b-sqrt(b*b-4*a*c))/(2*a) 
      return
      end function XCOeq


!----this function is to calculate chemistry_timescale in every timestep
      double precision function chemistry_timescale(i,j,k)
      USE global_constants
      USE fluid_var_init
      USE scalar_var_init
      IMPLICIT NONE
      integer :: i,j,k
      double precision GCO,GCH3O
!----chemical eqilibrium constant
      double precision K17eq
!----a few mediators
      double precision k9_zero,k9_infinite,k9
      double precision n
      double precision XH2
      XH2=0.83 
      GCO=-0.088787*T(i,j,k)-111.137747          
      GCH3O=0.097*T(i,j,k)+8.64
      K17eq=exp(-(GCO-GCH3O)/(RGAS*1e-10*T(i,j,k)))
      k9_zero=1.4e-6*T(i,j,k)**(-1.2)*exp(-7800/T(i,j,k))
      k9_infinite=1.5e+11*T(i,j,k)*exp(-12880/T(i,j,k))
!----number density of gas
      n=PG(i,j,k)/(KBOLTZ*T(i,j,k))
      k9=k9_zero*k9_infinite/(k9_zero*n+k9_infinite)
      chemistry_timescale=K17eq/(k9*(XH2**1.5)*
     $  ((PG(i,j,k)*1e-6)**1.5)*n) 
c      if(myid.eq.0) then
c         print *,'chem',i,j,k,chemistry_timescale
c      endif
      return
      end function chemistry_timescale


!---- this subroutine is to calculate initial CO fraction of equilibrium state      
      SUBROUTINE CO_SCALAR_init
      USE global_constants
      USE fluid_var_init
      USE scalar_var_init
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,j,k
!---a few constants
      double precision :: c1,c2,XH2
!---Gibbs Free Energy  
      double precision :: GCO             
      double precision :: GCH4            
      double precision :: GH2O
!---a few mediators            
      double precision :: f,a,b,c 
!---chemical equlibrium constant        
      double precision :: K3eq            
!----initialize the value of constants
      c1=5.91e-4
      c2=0.00104
      XH2=0.83 
!----calculate CO fraction of equilibrium state 
      DO k=1,locNZ
         DO j=1,locNY
            DO i=1,locNX
               GCO=-0.088787*T(i,j,k)-111.137747
               GCH4=0.105*T(i,j,k)-84.68434
               GH2O=0.0537546*T(i,j,k)-245.8666
               K3eq=exp(-(GCO-GCH4-GH2O)/(RGAS*1e-10*T(i,j,k)))
               f=K3eq/(PG(i,j,k)*PG(i,j,k)*XH2*XH2*XH2)*1e+12
               a=f
               b=-f*(c1+c2)-1
               c=f*c1*c2
!----pass the value to SCALAR(1,i,j,k) in order to initialize SCALAR(1,*,*,*)
               SCALAR(1,i,j,k)=(-b-sqrt(b*b-4*a*c))/(2*a)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CO_SCALAR_init


    
!----this subroutine is CO relaxation scheme: update the value of SCALAR in every timestep     
      SUBROUTINE CO_relaxation
      USE fluid_var_init
      USE scalar_var_init
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,j,k
      double precision :: XCOeq,chemistry_timescale
!---update the value of SCALAR(1,*,*,*) 
      DO k=1,locNZ
         DO j=1,locNY
            DO i=1,locNX
               if(DELT.lt.chemistry_timescale(i,j,k)) then
                  SCALAR(1,i,j,k)=SCALAR(1,i,j,k)-(SCALAR(1,i,j,k)-
     $                 XCOeq(i,j,k))*DELT/chemistry_timescale(i,j,k)
                  
               else
                  SCALAR(1,i,j,k)=XCOeq(i,j,k)                 
               endif
c              if(myid.eq.0) then                                           
c                   print *, i,j,k,SCALAR(1,i,j,k)/5.91e-4   
c              endif
c               if(SCALAR(1,i,j,k)/5.91e-4>1) then
c                  print *, 'X>1 after CO relaxation at',i,j,k
c               endif
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CO_relaxation



      SUBROUTINE CALC_CARBON_MASS
      USE input_init
      USE fluid_var_init
      USE grid_var_init
      USE mpi_var_init
      USE mpi_grid_init
      USE scalar_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: i,j,k
      double precision :: loc_c_mass,global_c_mass
      loc_c_mass = 0.0
      do k=1,locNZ
         do j=1,locNY
            do i=1,locNX 
               loc_c_mass=loc_c_mass+scalar(1,i,j,k)*rh(i,j,k)*
     %              volxb(i)*volyb(j)*volzb(k)
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(loc_c_mass,global_c_mass,1,
     %     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      write(*,'(A,2(1x,I8),2(1x,e12.3))') 'carbon mass',num_iter,myid,
     $     loc_c_mass,global_c_mass

      RETURN
      END SUBROUTINE CALC_CARBON_MASS
      



      
            
      
