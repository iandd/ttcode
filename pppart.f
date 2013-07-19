      
      SUBROUTINE INITPARTICLES
      USE particle_var_init
      USE input_init
      USE global_constants
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer indx
      double precision :: XMAX1,XMIN1
      double precision :: YMAX1,YMIN1
      double precision :: ZMAX1,ZMIN1


      if (ncosys.eq.2) then
         XMAX1=XMAX
         XMIN1=XMIN
**  If Spherical Polars: Y and Z-Coordinates in Degrees!!
         YMAX1=YMAX*2.d0*PI
         YMIN1=YMIN*2.d0*PI
         ZMAX1=(ZMAX/2.d0)*pi
         ZMIN1=(ZMIN/2.d0)*pi
      else
         write(*,*) ' Coordinate System AVALIBLE WITH PARTICLES'
         print *,'NCOSYS=',NCOSYS
         STOP
      end if

c      do indx=1,nparticles
c      part_posx(indx) = (XMAX1-XMIN1)/2.d0
c      part_posy(indx) = (YMAX1-YMIN1)/2.d0
c      part_posz(indx) = (ZMAX1-ZMIN1)/2.d0
c     enddo
c----
      part_posx(1) = 8.40*(10.**9.)
      part_posy(1) = (YMAX1-YMIN1)/2.d0
      part_posz(1) = 0.5*ZMIN1
c----
      part_posx(2) = 9.10*(10.**9.)
      part_posy(2) = YMAX1/4.d0
      part_posz(2) = 0.0
c----
      part_posx(3) = 10.10*(10.**9.)
      part_posy(3) = 3.d0*YMAX1/4.d0
      part_posz(3) = 0.5*ZMAX1

      if(myid.eq.3) then
         do indx=1,nparticles
            print *,'init',indx,part_posx(indx),part_posy(indx),
     %           part_posz(indx)
         enddo
      endif

      END SUBROUTINE INITPARTICLES

      SUBROUTINE READPOSPARTICLES
      USE particle_var_init
      USE input_init
      IMPLICIT NONE
      integer indx
      END SUBROUTINE READPOSPARTICLES



      SUBROUTINE INTERPOLATEVELOCITIES
      USE particle_var_init
      USE input_init
      USE grid_var_init
      USE fluid_var_init
      USE mpi_var_init
      IMPLICIT NONE
      include "mpif.h"
      integer indx,i,j,k
      integer x_hindex,x_lindex,y_hindex,y_lindex,z_hindex,z_lindex


c      print *,'bnds',myid,xyb(0),xyb(locNY),part_posy(1)


      do indx=1,nparticles 
!-- start by finding the processor of the particle (Using Y and Z coord)
         if(part_posy(indx).ge.xyb(0).and.part_posy(indx).lt.
     %        xyb(locNY)) then
            if(part_posz(indx).ge.xzb(0).and.part_posz(indx).lt.
     %           xzb(locNZ)) then

               print *,'XX',part_posx(indx),xxb(0),xxb(locNX)
!---- find the x-position of the particle
               do i=0,locNX
                  if(part_posx(indx).gt.xxb(i)) then
                     x_hindex = i
                     x_lindex = i-1
                     goto 100
                  endif
               enddo
 100           CONTINUE
!---- find the y-position of the particle
               print *,'YY',part_posy(indx),xyb(0),xyb(locNX)
               do j=0,locNY
                  if(part_posy(indx).gt.xyb(j)) then
                     y_hindex = j
                     y_lindex = j-1
                     goto 101
                  endif
               enddo
 101           CONTINUE
!---- find the z-position of the particle
               print *,'ZZ',part_posz(indx),xzb(0),xzb(locNX)
               do k=0,locNZ
                  if(part_posz(indx).gt.xzb(k)) then
                     z_hindex = k
                     z_lindex = k-1
                     goto 102
                  endif
               enddo
 102           CONTINUE


               print *,'INSIDE PARTICLE ROUTINE:indx=',indx,myid,
     %              part_posx(indx),
     %              part_posy(indx),part_posz(indx)
               print *,'indexs',indx,x_hindex,x_lindex,y_hindex,
     %              y_lindex,z_hindex,z_lindex

               
!            DOSOMTHING LIKE THIS 
!               part_vx(indx)=  V(x_lindex) V(x_hindex)
!               RHQX(I,J,K)=DDX1(I)*RH(I,J,K)+DDX0(I)*RH(I-1,J,K)


            endif
         endif
!         BROADCAST THE INFO TO ALL PROCESSORS
      enddo
      END SUBROUTINE INTERPOLATEVELOCITIES

      SUBROUTINE ADVANCEEULER
      USE particle_var_init
      USE input_init
      USE fluid_var_init
      IMPLICIT NONE
      integer indx
      do indx=1,nparticles 
         part_posx(indx) = part_posx(indx) + part_vx(indx)*delt
         part_posy(indx) = part_posy(indx) + part_vy(indx)*delt
         part_posz(indx) = part_posz(indx) + part_vz(indx)*delt
      enddo
      END SUBROUTINE ADVANCEEULER
