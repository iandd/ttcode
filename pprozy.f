
      SUBROUTINE RozyArtvis2
* correctly volume weighted!!!
      include 'ttdim.inc'
      include 'tt.inc'

      IMPLICIT NONE
      integer :: i,j,k,NXbou
      double precision :: rpower,fac,facm
      double precision :: rhqtemp(ny0)
      double precision :: sum,vol,dsym,Vkepfrac,vsym,Geq
      
c      rpower=0.d0
c      rpower=-1.d0
      rpower=-0.5d0
              
      NXbou=5
      do i=1,locNX        
         fac = dmax1(0.d0,dfloat(NXbou+1-i)/dfloat(NXbou))
         fac = dmax1(fac,dfloat(NXbou+i-NX)/dfloat(NXbou))**4.d0
         facm = 1.d0 - fac

        if (fac.ne.0.d0) then
          sum = 0.d0
          vol = 0.d0
          do j=1,locNY-1
            rhqtemp(j) = (RH(i,j,1)+RH(i,j+1,1))/2.d0*volya(j+1)
            sum = sum + rhqtemp(j)
            vol = vol + volya(j+1)
          end do
          dsym = sum/dfloat(locNY-1)/vol

          do j=1,NY-1
             V(i,j,1) = V(i,j,1)*facm
          end do

          Vkepfrac=1.d0 - 2.d0*PG0/2.d0          
          vsym = Vkepfrac*sqrt(GM/XXB(I))/geoxg(i)
          
          do j=1,NY-1
            Geq=G(i,j,1)+omrot
            G(i,j,1) = (vsym*fac + Geq*facm) - omrot
          end do

          dsym=RH0/dzb(1)/geoxh(i)
* add a potential slope to the density
     &         *(xxb(i)/1.5e13)**rpower
          do j=1,NY-1
            RH(i,j,1) = dsym*fac + RH(i,j,1)*facm
          enddo
        endif
      enddo
      RETURN
      END

