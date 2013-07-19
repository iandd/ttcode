 

      FUNCTION SETPLANETXGRAV(pnum)
      IMPLICIT NONE
      INTEGER :: pnum
      double precision :: setplanetxgrav
c         gravx=-1262.7538 !-- used for /cooler simulations (too high for 209458!)
c         gravx=-934.58 !-- correct value for HD209458
c         gravx=-2140.d0 !-- correct value for HD189733b
      if(pnum.eq.0) then  !--earths value
         setplanetxgrav = -980.d0
      elseif(pnum.eq.1) then !--HD209458b
         setplanetxgrav = -910.0
      elseif(pnum.eq.2) then !--HD189733b
         setplanetxgrav = -2189.5
      elseif(pnum.eq.3) then !--Tres-1
         setplanetxgrav = -1594.d0
      elseif(pnum.eq.4) then !--HD149026
         setplanetxgrav = -1674.d0
      elseif(pnum.eq.5) then !--WASP-12 (approx and large R)
         setplanetxgrav = -1044.d0
      elseif(pnum.eq.6) then !--HD17156b
         setplanetxgrav = -7607.1
      elseif(pnum.eq.7) then !--HAT-P-2
         setplanetxgrav = -23623.5
      elseif(pnum.eq.8) then !--HAT-P-7
         setplanetxgrav = -2209.44
      else
         print *,'planet_num not avail:setplanetgrav',pnum
         stop
      endif
      RETURN
      END FUNCTION SETPLANETXGRAV


      FUNCTION SETPLANETKAPEXTRA(pnum)
      IMPLICIT NONE
      INTEGER :: pnum
      double precision :: setplanetkapextra
      if(pnum.eq.1) then  !--HD209458b
c         setplanetkapextra = 0.5
c         setplanetkapextra = 0.4
c         setplanetkapextra = 0.3
c         setplanetkapextra = 0.12
         setplanetkapextra = 0.1
      elseif(pnum.eq.2) then !--HD189733b
         setplanetkapextra = 0.035
      elseif(pnum.eq.3) then !--Tres-1
         setplanetkapextra = 0.0
      elseif(pnum.eq.4) then !--HD149026
         setplanetkapextra = 0.0
      elseif(pnum.eq.5) then !--WASP12b
         setplanetkapextra = 0.0
      else
         print *,'planet_num notavail:readbur,extra',pnum
         call clean_stop
      endif
      RETURN
      END FUNCTION SETPLANETKAPEXTRA


      SUBROUTINE PRINTPLANETNAME(PLANET_NUM)
      USE MPI_VAR_INIT
      IMPLICIT NONE
      integer :: planet_num
      if(planet_num.eq.1) then  !--HD209458b
         if(myid.eq.0) print *,planet_num,'= HD209458b'
      elseif(planet_num.eq.2) then !--
         if(myid.eq.0) print *,planet_num,'= HD189733b'
      elseif(planet_num.eq.3) then !--TRES-1
         if(myid.eq.0) print *,planet_num,'= TRES-1'
      elseif(planet_num.eq.4) then !--HD149026
         if(myid.eq.0) print *,planet_num,'= HD149026'
      elseif(planet_num.eq.5) then !--WASP12b
         if(myid.eq.0) print *,planet_num,'= WASP12b'
      elseif(planet_num.eq.6) then !--HDHD17156
         if(myid.eq.0) print *,planet_num,'= HD17156'
      elseif(planet_num.eq.7) then !--HAT-P-2
         if(myid.eq.0) print *,planet_num,'= HAT-P-2'
      elseif(planet_num.eq.8) then !--HAT-P-7
         if(myid.eq.0) print *,planet_num,'= HAT-P-7'
      else
         print *,'planet_num not avail:init(2)',planet_num
         call clean_stop
      endif
      RETURN
      END SUBROUTINE PRINTPLANETNAME


      SUBROUTINE OPENOPACITYFILE_WOTIOVO
      USE input_init
      USE MPI_VAR_INIT
      IMPLICIT NONE
      if(planet_num.eq.1) then  !--HD209458b
         if(ikaptyp.ge.15.and.ikaptyp.le.19) then !-eliza...
            open(12,file="HD209_elizaopc.dat",status='old')
            if(myid.eq.0) print *,'reading HD209_elizaopc.dat'
         elseif(ikaptyp.eq.21) then
            open(12,file="BURROWS6070_woTiVO_wavelth.dat",status='old')
            if(myid.eq.0) 
     $           print *,'reading BURROWS6070_woTiVO_wavelth.dat'
         elseif(ikaptyp.eq.23.or.ikaptyp.eq.24) then
            open(12,file="BURROWS6070_woTiVO_radtransfer.dat",
     $           status='old')
            if(myid.eq.0) 
     $           print *,'BURROWS6070_woTiVO_radtransfer.dat'
         else
            open(12,file="BURROWS6117_woutTiVO.dat",status='old')
            if(myid.eq.0) print *,'reading BURROWS6117_woutTiVO.dat'
         endif
      elseif(planet_num.eq.2) then !--HD189733b
         if(ikaptyp.eq.21) then
            open(12,file="BURROWS5040_woTiVO_wavelth.dat",status='old')
            if(myid.eq.0) 
     $           print *,'reading BURROWS5040_woTiVO_wavelth.dat'
         elseif(ikaptyp.eq.23.or.ikaptyp.eq.24) then
            open(12,file="BURROWS5040_woTiVO_radtransfer.dat",
     $           status='old')
            if(myid.eq.0)
     $           print *,'BURROWS5040_woTiVO_radtransfer.dat'
         elseif(ikaptyp.eq.25) then
            open(12,file="BURROWS5040_Ray_woTiVO_radtransfer.dat",
     $           status='old')
            if(myid.eq.0) 
     $           print *,'BURROWS5040_Ray_woTiVO_radtransfer.dat'
         else
            open(12,file="BURROWS5040_woutTiVO.dat",status='old')
            if(myid.eq.0) print *,'reading BURROWS5040_woutTiVO.dat'
         endif
      elseif(planet_num.eq.3) then !--Tres-1
         open(12,file="BURROWS5226_woutTiVO.dat",status='old')
         if(myid.eq.0) print *,'reading BURROWS5226_woutTiVO.dat'
      elseif(planet_num.eq.4) then !--HD149026
         open(12,file="BURROWS6147_woutTiVO.dat",status='old')
         if(myid.eq.0) print *,'reading BURROWS6147_woutTiVO.dat'
      elseif(planet_num.eq.5) then !--WASP-12
         if(ikaptyp.eq.21) then
            open(12,file="BURROWS6300_woTiVO_wavelth.dat",status='old')
            if(myid.eq.0) 
     $           print *,'reading BURROWS6300_woTiVO_wavelth.dat'
         elseif(ikaptyp.eq.23.or.ikaptyp.eq.24) then
            open(12,file="BURROWS6300_woTiVO_radtransfer.dat",
     $           status='old')
            if(myid.eq.0) 
     $           print *,'BURROWS6300_woTiVO_radtransfer.dat'
         else
            print *,'check opacity file for WASP-12b'
            call clean_stop
         endif
      elseif(planet_num.eq.6) then !--HD17156
         open(12,file="BURROWS6079_woutTiVO.dat")
         if(myid.eq.0) print *,'reading BURROWS6079_woutTiVO.dat'
      elseif(planet_num.eq.7) then !--HAT-P-2
         open(12,file="BURROWS6290_woutTiVO.dat")
         if(myid.eq.0) print *,'reading BURROWS6290_woutTiVO.dat'
      elseif(planet_num.eq.8) then !--HAT-P-7
         open(12,file="BURROWS6350_woutTiVO.dat")
         if(myid.eq.0) print *,'reading BURROWS6350_woutTiVO.dat'
      else
         print *,'planet_num not avail',planet_num
         call clean_stop
      endif
      RETURN
      END SUBROUTINE OPENOPACITYFILE_WOTIOVO



      SUBROUTINE OPENOPACITYFILE_WITH_TIOVO
      USE MPI_VAR_INIT
      USE input_init
      IMPLICIT NONE
      if(planet_num.eq.1) then  !--HD209458b
         if(ikaptyp.eq.20) then
            open(12,file="BURROWS6070_wTiVO_wavelth.dat",status='old')
            if(myid.eq.0) 
     $           print *,'reading BURROWS6070_wTiVO_wavelth.dat'
         else
            open(12,file="BURROWS6117_wTiVO.dat",status='old')
            if(myid.eq.0) print *,'reading BURROWS6117_wTiVO.dat'
         endif
      elseif(planet_num.eq.2) then !--HD189733b
         open(12,file="BURROWS5040_wTiVO.dat",status='old')
         if(myid.eq.0) print *,'reading BURROWS4980_wTiVO.dat'
      elseif(planet_num.eq.3) then !--Tres-1
         print *,'no w/TiVO opacity file created for this planet',
     $        planet_num
         call clean_stop
      elseif(planet_num.eq.4) then !--HD149026
          print *,'no w/TiVO opacity file created for this planet',
     $        planet_num
         call clean_stop
      elseif(planet_num.eq.5) then !--WASP-12
          print *,'no w/TiVO opacity file created for this planet',
     $        planet_num
         call clean_stop
      elseif(planet_num.eq.6) then !--HD17156
          print *,'no w/TiVO opacity file created for this planet',
     $        planet_num
         call clean_stop
      elseif(planet_num.eq.7) then !--HAT-P-2
          print *,'no w/TiVO opacity file created for this planet',
     $        planet_num
         call clean_stop
      elseif(planet_num.eq.8) then !--HAT-P-7
          print *,'no w/TiVO opacity file created for this planet',
     $        planet_num
         call clean_stop
      else
         print *,'planet_num not avail',planet_num
         call clean_stop
      endif
      RETURN
      END SUBROUTINE OPENOPACITYFILE_WITH_TIOVO




