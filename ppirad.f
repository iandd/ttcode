

      SUBROUTINE CALC_INITIAL_MODEL
      USE saumon_var_init
      USE fluid_var_init      
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,iter,iter2
      double precision :: Score,Pcentral
      double precision :: r_interior(ngmax),m_interior(ngmax)
      double precision :: P_interior(ngmax),T_interior(ngmax)
      double precision :: m_exterior(locNX,locNY,locNZ)
      double precision :: Pbrz(locNY,locNZ)
      double precision :: Sbrz(locNY,locNZ),Rbrz(locNY,locNZ)
      double precision :: Mbrz(locNY,locNZ)
      double precision :: Mtcz(locNY,locNZ),Rtcz(locNY,locNZ)
      logical :: exterior_conv
!--make sure BCXMIN.eq.9 (ie bc using Fbottom)
!--BURROWS MODEL FOR GLOBAL_T AND GLOBAL_RH GUESS ALREADY READ IN
      
!--D_RH=?
      if(myid.eq.0) then
         print *,'INITIAL MODEL'
         do i=1,locNX
            print *,rh(i,0,0),T(i,0,0)
         enddo
      endif


!--set Fbottom (=Fcore) to some value..
      Fbottom = 1800.0

!--read in Saumon EOS for H and He
      call readtable(HTABLE,NTMP,NPRS,HTLOG,'H_TAB_I.A ')
      call readtable(HETABLE,NTMP,NPRS,HETLOG,'HE_TAB_I.A')

!--initial guess of central entropy (units of erg/g/K) and pressure
      Score = 8.45369542733284*KBOLTZ/MHYRODGEN ! 8.4785*kb/mH
      Pcentral=10.d0**(13.3134575673722) ! 10.d0**(13.31)

      do iter=1,1000
         exterior_conv = .false.
c         do iter2=1,1000
c            if(myid.eq.0) print *,'iter2=',iter2
!--call coupled_energy_eqs
            if(myid.eq.0) print *,'before coupled'
            CALL COUPLED_ENERGY_EQS
            if(myid.eq.0) print *,'after coupled'
!--call routine to calculate new density by re-establishing hydrostatic balance
            CALL APPROX_DENSITY(m_exterior)
c         enddo
!--find Pbrz, Sbrz, and Rbrz, the pressure, entropy, and 
!    radius at the base of the radiative zone
         CALL FIND_BRZ(m_exterior,Pbrz,Sbrz,Rbrz,Mbrz)
         if(myid.eq.0) print *,'BRZ',Pbrz(1,1),Sbrz(1,1),Rbrz(1,1),
     %        Mbrz(1,1)
!--call a routine to integrate out from center (using Fbottom and planetary parameters)
         CALL INTERIOR_SOLUTION(Score,Pcentral,r_interior,m_interior,
     %     P_interior,T_interior)
!--Use Pbrz to find Mtcz, Rtcz, the mass and radius at the top of the convective zone
         CALL FIND_TCZ(Pbrz,Mtcz,Rtcz,r_interior,m_interior,
     %        P_interior,T_interior)
!--compare two solutions and adjust Fbottom, Pc, and Score
!-- just do at 1,1 for now...
! 1) Sbrz = Score ?
! 2) Mbrz = Mtcz ?
! 3) Rbrz = Rtcz ?
         CALL ADJUST_STRUC(Score,Sbrz,Mbrz,Mtcz,Rbrz,Rtcz,
     %        Pcentral)

         CALL CLEAN_STOP

!--repeate until convergence
      enddo
      RETURN
      END SUBROUTINE CALC_INITIAL_MODEL


      SUBROUTINE ADJUST_STRUC(Score,Sbrz,Mbrz,Mtcz,Rbrz,Rtcz,
     %     Pcentral)
      USE fluid_var_init      
      USE mpi_var_init
      IMPLICIT NONE
      double precision :: Score,Sbrz(locNY,locNZ)
      double precision :: Mbrz(locNY,locNZ),Mtcz(locNY,locNZ)
      double precision :: Rbrz(locNY,locNZ),Rtcz(locNY,locNZ)
      double precision :: Pcentral
      print *,myid
      write(*,'(A,i8,6(1x,e12.3))') 'final',myid,Score,Sbrz(1,1),
     $     Mbrz(1,1),Mtcz(1,1),Rbrz(1,1),Rtcz(1,1)      
      END SUBROUTINE ADJUST_STRUC

      SUBROUTINE APPROX_DENSITY(mass)
      USE input_init
      USE fluid_var_init
      USE grid_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      integer :: i,j,k
      double precision :: mass(locNX,locNY,locNZ),gravity(locNX)
      double precision :: density_ln(locNX)
      call PRESSG
      do k=1,locNZ
         do j=1,locNY
!-- calculate mass 
            mass(upperbnd(j,k),j,k) = MSOL*mstar2 !--planet mass
            do i=upperbnd(j,k)-1,1,-1
               mass(i,j,k) = mass(i+1,j,k) - dxa(i+1)*4.d0*PI*
     $              (xxb(i+1)**2.d0)*rh(i+1,j,k)
            enddo
!--   calculate gravity
            do i=1,upperbnd(j,k)
               gravity(i) = GRAV*mass(i,j,k)/(xxb(i)**2.d0)
            enddo
c!--   calculate natural log of density
c            density_ln(upperbnd(j,k)) = log(10.d0**(-10.d0))
c            do i=upperbnd(j,k)-1,1,-1
c               density_ln(i) = density_ln(i+1) + (dxb(i)/T(i,j,k))*
c     %              ((mu_gas*gravity(i)/RGAS)+
c     %              (T(i+1,j,k)-T(i,j,k))/dxb(i))
c            enddo
c!--   calculate new density
            do i=1,upperbnd(j,k)
c               rh(i,j,k) = exp(density_ln(i))
               rh(i,j,k) = -(PG(i+1,j,k)-PG(i,j,k))/(dxb(i)*gravity(i))
               if(myid.eq.0.and.j.eq.1.and.k.eq.1) then
                  print *,'rh',i,rh(i,j,k),mass(i,j,k)
               endif
            enddo
            do i=upperbnd(j,k),locNX+2
               rh(i,j,k) = 10.d0**(-10.d0)
            enddo

!--another method
c            do i=1,upperbnd(j,k)
c               rh(i,j,k) = -(1.d0/gravity(i))*
c     $              (PG(i,j,k)-PG(i-1,j,k))/dxb(i)
c               if(myid.eq.0.and.j.eq.1.and.k.eq.1) then
c                  print *,'rh',i,rh(i,j,k)
c               endif
c            enddo
c            do i=upperbnd(j,k)+1,locNX+2
c               rh(i,j,k) = 10.d0**(-10.d0)
c            enddo


         enddo
      enddo
      RETURN
      END SUBROUTINE APPROX_DENSITY

      SUBROUTINE FIND_BRZ(m_exterior,Pbrz,Sbrz,Rbrz,Mbrz)
      USE fluid_var_init
      USE grid_var_init
      USE rad_var_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      include 'mpif.h'
      integer :: i,j,k
      double precision :: m_exterior(locNX,locNY,locNZ)
      double precision :: Sbrz(locNY,locNZ),Rbrz(locNY,locNZ)
      double precision :: Pbrz(locNY,locNZ),Mbrz(locNY,locNZ)
      double precision :: g_val,rho_val,S_val,del_ad,del_rad
      if(myid.eq.0) print *, 'starting find_BRZ'
      call PRESSG
      do k=1,locNZ 
         do j=1,locNY
c            do i=upperbnd(j,k),1,-1
c               if(myid.eq.0) print *,'in BRZ',PG(i,j,k),T(i,j,k),rh(i,j,k)
c               call SAUMON_EOS(PG(i,j,k),T(i,j,k),rho_val,S_val,del_ad)      
c               g_val = GRAV*m_exterior(i,j,k)/(xxb(i)**2.d0)
c               if (del_rad.gt.del_ad) then
c                  Pbrz(j,k) = PG(i,j,k)
c                  Sbrz(j,k) = S_val
c                  Rbrz(j,k) = xxb(i)
c                  Mbrz(j,k) = m_exterior(i,j,k)
c                  goto 200
c               endif
            do i=1,upperbnd(j,k)
c               call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c               write(*,'(A,i8,i8,2(1x,e12.3))') 'try',
c     %              i,myid,PG(i,j,k),T(i,j,k)
c               call MPI_BARRIER(MPI_COMM_WORLD,ierr)
               call SAUMON_EOS(PG(i,j,k),T(i,j,k),rho_val,S_val,del_ad)
               g_val = GRAV*m_exterior(i,j,k)/(xxb(i)**2.d0)
               del_rad=(xkapR(i,j,k)*FLUX_X(i,j,k)*PG(i,j,k))/
     %              (16.d0*SBCONST*(T(i,j,k)**4)*FLIMX(i,j,k)*g_val)
c               if(myid.eq.0) write(*,'(A,i8,2(1x,e12.3))') 'find',i,
c     &              del_ad,del_rad
               if (del_rad.lt.del_ad) then
                  Pbrz(j,k) = PG(i,j,k)
                  Sbrz(j,k) = S_val
                  Rbrz(j,k) = xxb(i)
                  Mbrz(j,k) = m_exterior(i,j,k)
                  goto 200
               endif
            enddo
 200        CONTINUE
         enddo
      enddo
      RETURN
      END SUBROUTINE FIND_BRZ

      SUBROUTINE FIND_TCZ(Pbrz,Mtcz,Rtcz,r_interior,m_interior,
     %     P_interior,T_interior)      
      USE fluid_var_init
      USE saumon_var_init
      IMPLICIT NONE
      integer :: i,j,k
      double precision :: Pbrz(locNY,locNZ)
      double precision :: r_interior(ngmax),m_interior(ngmax)
      double precision :: P_interior(ngmax),T_interior(ngmax)
      double precision :: Mtcz(locNY,locNZ),Rtcz(locNY,locNZ)
      do k=1,locNZ 
         do j=1,locNY
            do i=1,ngmax
               if(Pbrz(j,k).le.P_interior(i)) then
                  Mtcz(j,k) = m_interior(i)
                  Rtcz(j,k) = r_interior(i)
                  goto 202
               endif
            enddo
 202        CONTINUE
         enddo
      enddo
      RETURN
      END SUBROUTINE FIND_TCZ

      SUBROUTINE INTERIOR_SOLUTION(Score,Pc,r_interior,m_interior,
     %     P_interior,T_interior)
      USE mpi_var_init
      USE saumon_var_init
      USE global_constants
      implicit none
      integer :: ni
      double precision :: r_interior(ngmax),m_interior(ngmax)
      double precision :: P_interior(ngmax),T_interior(ngmax)
      double precision :: Score,Pc,Tc
      double precision :: r_int,m_int,P_int,T_int,rho_int
      double precision :: S_test,del_ad,dr
      double precision :: y(2),dydr(2),eps
      double precision :: dy1(2),dy2(2),dy3(2),dy4(2)
!--central temperature based on central Pressure and entropy guess
      call find_saumonT(Pc,Score,Tc)
!--find central density
      call SAUMON_EOS(Pc,Tc,rho_int,S_test,del_ad)
      print *,'entropy test: should be equal',Score,S_test,myid
!--use taylor expansion to find values just off center
      r_int=RJUP/100.d0
      m_int = 4.d0*PI/3.d0 * rho_int * r_int**3
      P_int = Pc-2.d0*PI/3.d0 * GRAV*(rho_int**2)*(r_int**2)
!--fill runga-kutta vector
      y(1)=m_int
      y(2)=P_int
!--small factor for stepsize      
      eps=2.d0**(-4)
!--integrate outward
      ni=0
      do
         ni=ni+1

         if(myid.eq.0) print*,ni,y(1),y(2)

!-- calculate stepsize
         call derivs2(r_int,y,dydr,Score)
         dr=eps * min( r_int, abs(y(1)/dydr(1)), abs(y(2)/dydr(2)) )
!-- 4th order rk step for radius
         dy1=dydr*dr
         call derivs2(r_int+dr/2,y+dy1/2,dydr,Score)
         dy2=dydr*dr
         call derivs2(r_int+dr/2,y+dy2/2,dydr,Score)
         dy3=dydr*dr
         call derivs2(r_int+dr,y+dy3,dydr,Score)
         dy4=dydr*dr
         y=y+dy1/6.0+dy2/3.0+dy3/3.0+dy4/6.0
         r_int=r_int+dr

         m_int=y(1)
         P_int=y(2)
         call find_saumonT(P_int,Score,T_int)
         call SAUMON_EOS(P_int,T_int,rho_int,S_test,del_ad)
!--save the solution for r, m, p, and T
         r_interior(ni)=r_int
         m_interior(ni)=m_int
         P_interior(ni)=P_int
         T_interior(ni)=T_int
         !--don't let the solution wander off SAUMON EOS tables
         if (P_interior(ni).lt.1e5 .or. T_interior(ni).lt.350.0) then
            print *,'OFF THE TABLE IN INTERIOR SOLUTION'
            write(*,'(A,i8,4(1x,e12.3))') 'vals:',ni,r_interior(ni),
     $           m_interior(ni),P_interior(ni),T_interior(ni)
         endif
      enddo
      RETURN
      END SUBROUTINE INTERIOR_SOLUTION


      SUBROUTINE DERIVS2(r,y,dydr,Score)
      USE global_constants
      implicit none
      double precision :: r,y(2),dydr(2),Score
      double precision :: mr,P,T,rho,Stest,del_ad      
      mr=y(1)
      P=y(2)
      call find_saumonT(P,Score,T)
      call SAUMON_EOS(P,T,rho,Stest,del_ad)
      dydr(1)=4.0*pi*r**2*rho
      dydr(2)=-GRAV*mr*rho/r**2
      RETURN
      END SUBROUTINE DERIVS2
