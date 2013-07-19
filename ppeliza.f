
      SUBROUTINE ElizasOPC(T_IN,PG_IN,XCO_IN,BKAP_P,BKAP_R,BKAP_A,
     $     BDXKAPP_DT,BDXKAPA_DT,i_indx,j_indx,k_indx)
!---   Routine to interpolate in Eliza Kemptons opacities tables
!---  NOTE:  Minimum pressure     = 1.00102e-4 dyn/cc
!---         Maximum pressure    = 1.0e10 dyn/cc
!---         Minimum Temperature  = 100.0
!---         Maximum Temperature  = 3000.0
!-----------------
      USE opc_init
      USE mpi_var_init
      USE fluid_var_init
      USE input_init
      IMPLICIT NONE
      include "mpif.h"
      integer :: K,IT,IT1,IW,IW1,IC,IC1
      double precision :: T_IN,PG_IN,XCO_IN,XCO1,XCH41
      double precision :: TM1,PG1,BKAP_P,BKAP_R,BKAP_A
      double precision :: BDXKAPP_DT,BDXKAPA_DT
      double precision :: WS1,WS,WS2,WS3
      double precision :: R00,R10,R11,R01,R_WS4,R_WS5
      double precision :: P00,P10,P11,P01,P_WS4,P_WS5
      double precision :: A00,A10,A11,A01,A_WS4,A_WS5
      double precision :: c1
      integer :: i_indx,j_indx,k_indx
      character proc_num*3
      character filename*12
      IF(T_IN.LE. 0.) then
         print *,' NEG. T. IN Elizas OPC=',T_IN,myid
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         STOP
      endif
      PG1 = PG_IN
      TM1 = T_IN

!--scale X_CO up to 0-1
      c1=5.91e-4
      XCO1 = XCO_IN/c1
      XCH41 = 1.d0-XCO1

!--- ikaptyp = 17 : Eliza Kempton's solar + purly CO opacities + extra absorb
!--- ikaptyp = 18 : Eliza Kempton's solar + purly CH4 opacities + extra absorb
      if(ikaptyp.eq.17) then 
         XCO1 = 1.0
         XCH41 = 0.0
      endif
      if(ikaptyp.eq.18) then 
         XCO1 = 0.0
         XCH41 = 1.0
      endif

      IF(TM1 .LT. 100.0)  TM1 = 100.001
      IF(TM1 .GT. 3000.0) TM1 = 2999.001
      IF(PG1 .LT. 1.0e-4) PG1 = 1.0001e-4
      IF(PG1 .GT. 1.0e10) PG1 = 9.999999e9
      IF(XCH41 .GT. 1.d0) XCH41=1.d0
      IF(XCH41 .LT. 0.d0) XCH41=0.d0

      DO K = 1,39
         IF(ElizaT(K) .GE. TM1) GO TO 115
      ENDDO

      IF(K.EQ.40) THEN
         PRINT *,'Overflow in opacity temperature array'
         print *,'iteration',num_iter,III
         PRINT *,'Elizas opacities: T/myid=',TM1,myid
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd(j,k)',upperbnd(j_indx,k_indx)
         print *,'upperbnd(j-1,k)',upperbnd(j_indx-1,k_indx)
         print *,'upperbnd(j+1,k)',upperbnd(j_indx+1,k_indx)
         print *,'upperbnd(j,k-1)',upperbnd(j_indx,k_indx-1)
         print *,'upperbnd(j,k+1)',upperbnd(j_indx,k_indx+1)
         print *,'v',v(i_indx,j_indx,k_indx)
         print *,'g',g(i_indx,j_indx,k_indx)
         print *,'h',h(i_indx,j_indx,k_indx)
         print *,'rh',rh(i_indx,j_indx,k_indx)
         CALL CLEAN_STOP
      ENDIF
 115  IT = K
      IT1 = K-1
      IF(IT1.eq.0) THEN
         PRINT *,'underflow in opacity temperature array'
         PRINT *,'Elizas opacities',TM1
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'iteration',num_iter,III
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         CALL CLEAN_STOP
      ENDIF
      WS1 = ElizaT(IT1) - TM1
      WS =  ElizaT(IT1) - ElizaT(IT)
      WS1 = WS1 / WS
      WS = 1.d0 - WS1
c-pressure search
      DO K = 1,57
         IF(ElizaPG(K) .GE. PG1) GO TO 125
      ENDDO
      IF(K.EQ.58) THEN
         PRINT *,'Overflow in elizas opacity pressure array',PG1,myid
         print *,'iteration',num_iter,III
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         CALL CLEAN_STOP
      ENDIF
 125  IW = K
      IW1 = K-1
      IF(IW1.EQ.0) THEN
         PRINT *,'underflow in elizas opacity pressure array'
         print *,'iteration',num_iter,III
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         CALL CLEAN_STOP
      ENDIF
      WS2 = PG1 - ElizaPG(IW1)
      WS3 = ElizaPG(IW) - ElizaPG(IW1)
      WS2 = WS2 / WS3
      WS3 = 1. - WS2

!-- FIND THE BOUNDING VALUES OF XCO
      DO K = 1,15
         IF(ElizaCH4(K) .GE. XCH41) GO TO 135
      ENDDO
 135  IC = K
      IC1 = K-1
      IF(IC1.ne.0) THEN !-only do interpolation if XCH4>0 (avoids 1/0=inf)
!-- ROSSLAND MEAN OPACITIES. EACH RXX IS INTERPOLATED USING XCH4 VALUE
         R00 = (((ElizaXKAP_ross(IT,IW,IC)-ElizaXKAP_ross(IT,IW,IC1))/
     %        (ElizaCH4(IC)-ElizaCH4(IC1)))*(XCH41-ElizaCH4(IC1))) +
     %        ElizaXKAP_ross(IT,IW,IC1)
         R10 = (((ElizaXKAP_ross(IT,IW1,IC)-ElizaXKAP_ross(IT,IW1,IC1))/
     %        (ElizaCH4(IC)-ElizaCH4(IC1)))*(XCH41-ElizaCH4(IC1))) +
     %        ElizaXKAP_ross(IT,IW1,IC1)
         R11 = (((ElizaXKAP_ross(IT1,IW1,IC)-
     %        ElizaXKAP_ross(IT1,IW1,IC1))/
     %        (ElizaCH4(IC)-ElizaCH4(IC1)))*(XCH41-ElizaCH4(IC1))) +
     %        ElizaXKAP_ross(IT1,IW1,IC1)
         R01 = (((ElizaXKAP_ross(IT1,IW,IC)-ElizaXKAP_ross(IT1,IW,IC1))/
     %        (ElizaCH4(IC)-ElizaCH4(IC1)))*(XCH41-ElizaCH4(IC1))) +
     %        ElizaXKAP_ross(IT1,IW,IC1)
      ELSE !-use the XCH4=0 value
         R00 = ElizaXKAP_ross(IT,IW,1)
         R10 = ElizaXKAP_ross(IT,IW1,1)
         R11 = ElizaXKAP_ross(IT1,IW1,1)
         R01 = ElizaXKAP_ross(IT1,IW,1)      
      ENDIF
      R_WS4 = WS2*(WS1*R00 + WS*R01)
      R_WS5 = WS3*(WS1*R10 + WS*R11)
      BKAP_R = (R_WS4+R_WS5)

!-- PLANCK MEAN OPACITIES WITH SEPERATE CO AND CH4
      P00 = ElizaXKAP_plnk_woC(IT,IW)+XCO1*ElizaXKAP_plnk_CO(IT,IW)+
     $     XCH41*ElizaXKAP_plnk_CH4(IT,IW)
      P10 = ElizaXKAP_plnk_woC(IT,IW1)+XCO1*ElizaXKAP_plnk_CO(IT,IW1)+
     $     XCH41*ElizaXKAP_plnk_CH4(IT,IW1)
      P11 = ElizaXKAP_plnk_woC(IT1,IW1)+XCO1*ElizaXKAP_plnk_CO(IT1,IW1)+
     $     XCH41*ElizaXKAP_plnk_CH4(IT1,IW1)
      P01 = ElizaXKAP_plnk_woC(IT1,IW)+XCO1*ElizaXKAP_plnk_CO(IT1,IW)+
     $     XCH41*ElizaXKAP_plnk_CH4(IT1,IW)
      P_WS4 = WS2*(WS1*P00 + WS*P01)
      P_WS5 = WS3*(WS1*P10 + WS*P11)
      BKAP_P = (P_WS4+P_WS5)

!-- ABSORPTION OPACITIES WITH SEPERATE CO AND CH4
      A00 = ElizaXKAP_abs_woC(IT,IW)+XCO1*ElizaXKAP_abs_CO(IT,IW)+
     $     XCH41*ElizaXKAP_abs_CH4(IT,IW)
      A10 = ElizaXKAP_abs_woC(IT,IW1)+XCO1*ElizaXKAP_abs_CO(IT,IW1)+
     $     XCH41*ElizaXKAP_abs_CH4(IT,IW1)
      A11 = ElizaXKAP_abs_woC(IT1,IW1)+XCO1*ElizaXKAP_abs_CO(IT1,IW1)+
     $     XCH41*ElizaXKAP_abs_CH4(IT1,IW1) 
      A01 = ElizaXKAP_abs_woC(IT1,IW)+XCO1*ElizaXKAP_abs_CO(IT1,IW)+
     $     XCH41*ElizaXKAP_abs_CH4(IT1,IW) 
      A_WS4 = WS2*(WS1*A00 + WS*A01)
      A_WS5 = WS3*(WS1*A10 + WS*A11)
      BKAP_A = (A_WS4+A_WS5)

!-- CLUDGE FOR NOW
      BDXKAPP_DT = (P00-P01)/(ElizaT(IT)-ElizaT(IT1))
      BDXKAPA_DT = (A00-A01)/(ElizaT(IT)-ElizaT(IT1))

      IF(BKAP_R.LE.0.d0.OR.BKAP_P.LE.0.d0.OR.BKAP_A.LE.0.d0) THEN
         print *,'in elizas opc...'
         print *,'neg/zero BKAP_R,BKAP_P,BKAP_A',BKAP_R,BKAP_P,BKAP_A
         PRINT *,'TM1,PG1=',TM1,PG1
         print *,'tval=',ElizaT(IT),ElizaT(IT1)
         print *,'PGval=',ElizaPG(IW),ElizaPG(IW1)
         print *,'init Ws',WS1,WS,WS2,WS3
         print *,'ZXX_R=',R00,R10,R11,R01
         print *,'ZXX_P=',P00,P10,P11,P01
         print *,'ZXX_A=',A00,A10,A11,A01
         print *,'final Ws_R',R_WS4,R_WS5
         print *,'final Ws_P',P_WS4,P_WS5
         print *,'final Ws_A',A_WS4,A_WS5
         print *,'Rosselands:',ElizaXKAP_ross(IT,IW,IC1),
     $        ElizaXKAP_ross(IT1,IW,IC1),ElizaXKAP_ross(IT,IW1,IC1),
     $        ElizaXKAP_ross(IT1,IW1,IC1)
         write(*,'(A,10(1x,e12.3))') 'R00 =',R00,
     $        (ElizaXKAP_ross(IT,IW,IC)-ElizaXKAP_ross(IT,IW,IC1)),
     $        ElizaXKAP_ross(IT,IW,IC),ElizaXKAP_ross(IT,IW,IC1),
     %        (ElizaCH4(IC)-ElizaCH4(IC1)),ElizaCH4(IC),ElizaCH4(IC1),
     $        (XCH41-ElizaCH4(IC1)),XCH41,ElizaCH4(IC1)
         print *,'XCO1,XCH41',XCO1,XCH41
         print *, IT,IT1,IW,IW1,IC,IC1
         print *,'i,j,k',i_indx,j_indx,k_indx
         print *,'upperbnd',upperbnd(j_indx,k_indx)
         write(proc_num,'(i3.3)') myid
         filename = 'SNAP_'//TRIM(ADJUSTL(proc_num))//'.'//
     $        '100'
         OPEN (12,FILE=filename,STATUS='UNKNOWN',
     $        FORM='UNFORMATTED')
         REWIND(12)
         CALL PRINTOUT(1)
         CLOSE(12)
         stop
      ENDIF
      RETURN
      END SUBROUTINE ElizasOPC



      SUBROUTINE READELIZAOPCTABLE
      USE input_init
      USE opc_init
      USE mpi_var_init
      IMPLICIT NONE
      integer :: PGindex,Tindex
      double precision :: readP,readT,readkapP_woC,readkapP_CO
      double precision :: readkapP_CH4,readkapP_s_woC,readkapP_s_CO
      double precision :: readkapP_s_CH4,readkapR_wCperc1
      double precision :: readkapR_wCperc2,readkapR_wCperc3
      double precision :: readkapR_wCperc4,readkapR_wCperc5
      double precision :: readkapR_wCperc6,readkapR_wCperc7
      double precision :: readkapR_wCperc8,readkapR_wCperc9
      double precision :: readkapR_wCperc10,readkapR_wCperc11
      double precision :: readkapR_wCperc12,readkapR_wCperc13
      double precision :: readkapR_wCperc14,readkapR_wCperc15
      double precision :: kap_ext,SETPLANETKAPEXTRA
      call OPENOPACITYFILE_WOTIOVO
      DO PGindex=1,57
         DO Tindex=1,39
            READ(12,*) readP,readT,readkapP_woC,readkapP_CO,
     $           readkapP_CH4,readkapP_s_woC,readkapP_s_CO,
     $           readkapP_s_CH4,readkapR_wCperc1,
     $           readkapR_wCperc2,readkapR_wCperc3,readkapR_wCperc4,
     $           readkapR_wCperc5,readkapR_wCperc6,readkapR_wCperc7,
     $           readkapR_wCperc8,readkapR_wCperc9,readkapR_wCperc10,
     $           readkapR_wCperc11,readkapR_wCperc12,readkapR_wCperc13,
     $           readkapR_wCperc14,readkapR_wCperc15
            ElizaXKAP_plnk_woC(Tindex,PGindex)=readkapP_woC
            ElizaXKAP_plnk_CO(Tindex,PGindex)=readkapP_CO
            ElizaXKAP_plnk_CH4(Tindex,PGindex)=readkapP_CH4
            ElizaXKAP_abs_woC(Tindex,PGindex)=readkapP_s_woC
            ElizaXKAP_abs_CO(Tindex,PGindex)=readkapP_s_CO
            ElizaXKAP_abs_CH4(Tindex,PGindex)=readkapP_s_CH4
            ElizaXKAP_ross(Tindex,PGindex,1)=readkapR_wCperc1
            ElizaXKAP_ross(Tindex,PGindex,2)=readkapR_wCperc2
            ElizaXKAP_ross(Tindex,PGindex,3)=readkapR_wCperc3
            ElizaXKAP_ross(Tindex,PGindex,4)=readkapR_wCperc4
            ElizaXKAP_ross(Tindex,PGindex,5)=readkapR_wCperc5
            ElizaXKAP_ross(Tindex,PGindex,6)=readkapR_wCperc6
            ElizaXKAP_ross(Tindex,PGindex,7)=readkapR_wCperc7
            ElizaXKAP_ross(Tindex,PGindex,8)=readkapR_wCperc8
            ElizaXKAP_ross(Tindex,PGindex,9)=readkapR_wCperc9
            ElizaXKAP_ross(Tindex,PGindex,10)=readkapR_wCperc10
            ElizaXKAP_ross(Tindex,PGindex,11)=readkapR_wCperc11
            ElizaXKAP_ross(Tindex,PGindex,12)=readkapR_wCperc12
            ElizaXKAP_ross(Tindex,PGindex,13)=readkapR_wCperc13
            ElizaXKAP_ross(Tindex,PGindex,14)=readkapR_wCperc14
            ElizaXKAP_ross(Tindex,PGindex,15)=readkapR_wCperc15

!-- add in an extra absorber for ikaptyp=16:19
            if(ikaptyp.ge.16.and.ikaptyp.le.19) then
               kap_ext = SETPLANETKAPEXTRA(planet_num)
               if(Tindex.eq.1.and.PGindex.eq.1.and.myid.eq.0) then
                  print *,'Adding extra opacity to xkapA'
                  print *,' kappa_extra =',kap_ext,'cm^2/g'
                  print *,' For P <0.03 bars'
               endif              
               if(ElizaPG(PGindex)*(10.d0**(-6.d0)).le.0.03) then !-this was used in paper2
                  ElizaXKAP_abs_woC(Tindex,PGindex)=
     $                 ElizaXKAP_abs_woC(Tindex,PGindex)+kap_ext
               endif
            endif

         ENDDO
      ENDDO
      CLOSE(12)
      RETURN
      END SUBROUTINE READELIZAOPCTABLE
