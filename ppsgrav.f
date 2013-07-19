!---------------------
!
!---- SELFGRAV VARIABLES
      MODULE sgrav_init
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: ellipticK,ellipticE
      DOUBLE PRECISION :: TVAL(60),RVAL(14)
      DOUBLE PRECISION :: coeffA(5),coeffB(5),coeffC(2),coeffD(2)
      DATA coeffA/
     +     1.38629436112,0.09666344259,0.03590092383,
     +     0.03742563713,0.01451196212 /
      DATA coeffB/
     +     0.5,0.12498593597,0.06880248576,0.03328355346,
     +     0.00441787012/
      DATA coeffC/
     +     0.4630151,0.1077812/
      DATA coeffD/
     +     0.2452727,0.0412496/
      END MODULE sgrav_init

      MODULE sgrav_alloc
      USE rad_var_init
      USE fluid_var_init
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE allocatesgravvar
          IMPLICIT NONE
          ALLOCATE(ellipticK(locNX))
          ALLOCATE(ellipticE(locNX))
      END SUBROUTINE allocatesgravvar
      END MODULE sgrav_alloc




      subroutine diskgrav(Nr,Mp,sigmaAVG,rdisk,dr,dphidr,vcirc)
      implicit none
      integer i,j,Nr
      double precision Mp,sigmaAVG(Nr),rdisk(Nr),dr
      double precision rSG,EllipticK(Nr),EllipticE(Nr),rk(Nr)
      double precision dphidr(Nr),vcirc(Nr),rroche,sgsoften
      sgsoften = 0.5
      do j=1,Nr
         rSG = rdisk(j)
         rroche = ((Mp/3.)**(1./3.))*rSG
         do i=1,Nr
            rk(i) = (4.d0*rSG*rdisk(i)/((rSG+rdisk(i))**2.d0))**(0.5)
         enddo
         call calcEllipticK(Nr,ellipticK,rk)
         call calcEllipticE(Nr,ellipticE,rk)
         dphidr(j) = 0.d0
         do i=1,Nr
            if(abs(rdisk(i)-rSG).ge.sgsoften*rroche) then
               dphidr(j) = dphidr(j) + (rk(i)*EllipticK(i)-
     %             (EllipticE(i)/(1-(rk(i)**2.d0)))*
     %             ((rk(i)**3.d0)/4.d0)*((rdisk(i)/rSG)-(rSG/rdisk(i))))
     %             *sigmaAVG(i)*(rdisk(i)**0.5)
            endif
         enddo
         dphidr(j) = dphidr(j)*dr/(rSG**(1.5))
         vcirc(j) = (dphidr(j)*rSG)**(0.5)
      enddo
      return
      end


      subroutine calcEllipticK(Nr,ellipticK,rk)
      implicit none
      integer :: i,Nr
      double precision :: rk(Nr)
      double precision :: m1
      do i=1,Nr
         m1 = 1.d0-rk(i)
         ellipticK(i) = coeffA(1)+(coeffA(2)*m1)+(coeffA(3)*(m1**2.d0))+
     %        (coeffA(4)*(m1**3.d0))+(coeffA(5)*(m1**4.d0))+
     %        (coeffB(1)+(coeffB(2)*m1)+(coeffB(3)*(m1**2.d0))+
     %        (coeffB(4)*(m1**3.d0))+
     %        (coeffB(5)*(m1**4.d0)))*log(1.d0/m1)
      enddo
      return
      end

      subroutine calcEllipticE(Nr,ellipticE,rk)
      implicit none
      integer :: i,Nr
      double precision :: rk(Nr)
      double precision :: m1
      do i=1,Nr
         m1 = 1.d0-rk(i)
         ellipticE(i) = 1.d0+coeffA(1)*m1+(coeffA(2)*(m1**2.d0))+
     %        (coeffB(1)*m1+(coeffB(2)*(m1**2.d0)))*log(1.d0/m1)
      enddo
      return
      end



      SUBROUTINE INITRK
      IMPLICIT NONE
      double precision :: rk(locNX,locNX) 
      do i=1,locNX
         do j=1,locNX
            rk(i,j)=(4.d0*XXA(i)*XXA(j)/((XXA(i)+XXA(j))**2.d0))**(0.5)
         enddo
      enddo
      RETURN
      END SUBROUTINE INITRK


      SUBROUTINE AVERAGESIGMA
      IMPLICIT NONE

      double precision :: locAVG(locNX)
      double precision :: gatherSIG(locNX*numprocs)
      double precision :: avgSIMGA(global_NX)

      k=1
      DO i=1,locNX
         locAVG(i) = 0.d0
         DO j=1,locNY
            locAVG(i) = locAVG(i)+RH(i,j,k)
         ENDDO
         locAVG(i) = locAVG(i)/locNY
      ENDDO



      




      do i=1,global_NX
         DO r=1,proc_dims(1)
            if(proc_coods(1).eq.r) then
               sendINFO = locAVG(i)
            else
               sendINFO = 0.d0
            endif
            
            call MPI_ALLREDUCE(sendINFO,recINFO,1,
     %           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
            avgSigma(i)
            
         ENDDO
      ENDDO
      avgSIMGA (i) = locAVG













      END SUBROUTINE AVERAGESIGMA
