
      SUBROUTINE SAUMON_EOS(P,T,rho,S,del_ad)
      USE input_init
      USE global_constants
      USE saumon_var_init
      USE mpi_var_init
      IMPLICIT NONE
      double precision :: P,T,rho,S,del_ad
      double precision :: HDENSITY,HEDENSITY
      double precision :: S_H,S_HE,ST,SP
      double precision :: SMIX,dSMIXdT,dSMIXdP
      IF(log10(T).gt.2.10.and.log10(P).gt.4.00) then
!-- density
         call TABLEINTERPOLATE(T,P,HTLOG,HTABLE,NTMP,
     %        NPRS,4,HDENSITY)
         call TABLEINTERPOLATE(T,P,HETLOG,HETABLE,NTMP,
     %        NPRS,4,HEDENSITY)
         rho = ((XX/HDENSITY) + (YY/HEDENSITY))**(-1.d0)
!-- entropy
         call TABLEINTERPOLATE(T,P,HTLOG,HTABLE,NTMP,
     %        NPRS,5,S_H)
         call TABLEINTERPOLATE(T,P,HETLOG,HETABLE,NTMP,
     %        NPRS,5,S_HE)
         call smix_and_dsmix(T,P,SMIX,dSMIXdT,dSMIXdP)
         S = (XX*S_H) + (YY*S_HE) + SMIX
!--adiabatic gradient
         call calc_ST(T,P,S,S_H,S_HE,SMIX,dSMIXdT,ST)
         call calc_SP(T,P,S,S_H,S_HE,SMIX,dSMIXdP,SP)
         del_ad = -SP/ST
      else
         print *,'outside Saumon EOS',P,T
c         CALL CLEAN_STOP
         rho = P*mu_gas/(Rgas*T)
         del_ad = 0.4
         S = (KBOLTZ/MPROTON)*(2.5 + 
     %        log( ((KBOLTZ/(2.d0*PI*PLANCK))**1.5)*
     %        ((MPROTON*T)**1.5)* (mu_gas/rho) ) )
      endif
      return
      end

      SUBROUTINE find_SAUMONT(P,S,T)
      USE saumon_var_init
      IMPLICIT NONE
      integer :: iter
      double precision :: P,S,T,Tguess,Sguess
      double precision :: S_H,S_HE
      double precision :: SMIX,dSMIXdT,dSMIXdP
      double precision :: dSdT,ST
!-- entropy
      do iter=1,1000
         call TABLEINTERPOLATE(Tguess,P,HTLOG,HTABLE,NTMP,NPRS,5,S_H)
         call TABLEINTERPOLATE(Tguess,P,HETLOG,HETABLE,NTMP,NPRS,5,S_HE)
         call smix_and_dsmix(Tguess,P,SMIX,dSMIXdT,dSMIXdP)
         Sguess = (XX*S_H) + (YY*S_HE) + SMIX
!-- entropy derivative
         call calc_ST(Tguess,P,Sguess,S_H,S_HE,SMIX,dSMIXdT,ST)
         dSdT = (Sguess/Tguess)*ST
!-- new temperature guess
         Tguess = Tguess + 0.1*(S-Sguess)/dSdT
      enddo
      T = Tguess
      print *,'find_saumonT test, should be =:',S,Sguess
!-- THIS IS FOR A SINGLE ELEMENT... H OR He.. NOT A MIXTURE
c      PLOG = LOG10(P)
c      SLOG = LOG10(S)
c      if(PLOG .LT. 4.00 .or.PLOG .GT. 19.00) then
c         print *,'outside EOS table in find_SAUMONT',T,P
c         CALL CLEAN_STOP
c      endif
c!--find pressure using last table entry (with all pressures)
c      DO K = 1,NPRS
c         IF(TABLE(NTMP,K,1) .GT. PLOG) GO TO 125
c      ENDDO
c 125  IP = K
c      IP1 = K-1
c      WP1 = (PLOG-TABLE(NTMPS,IP1,1))/
c     %     (TABLE(NTMPS,IP,1)-TABLE(NTMPS,IP1,1))
c      WP2 = 1.d0-WP1
c!-- save entropies from table (2 lists, one for each pressure)
c      DO K = 1,NTMPS
c         SLIST_IP(k)  = TABLE(k,IP,5)
c         SLIST_IP1(k) = TABLE(k,IP1,5)
c      ENDDO
c!-- find entropy in P(IP) list
c      DO K = 1,NTMPS
c         if(SLOG.GT.SLIST_IP(k)) goto 126
c      ENDDO
c 126  IS2 = k
c      IS1 = k-1
c      WS1 = (SLOG-SLIST_IP(IS1))/(SLIST_IP(IS2)-SLIST_IP(IS1))
c      WS2 = 1.d0 - WS1
c!-- find entropy in P(IP1) list
c      DO K = 1,NTMPS
c         if(SLOG.GT.SLIST_IP1(k)) goto 127
c      ENDDO
c 127  IS4 = k
c      IS3 = k-1
c      WS3 = (SLOG-SLIST_IP1(IS3))/(SLIST_IP1(IS4)-SLIST_IP1(IS4))
c      WS4 = 1.d0 - WS3
c!-- save temperatures 
c      T00=TLOG(IS2)
c      T01=TLOG(IS1)
c      T10=TLOG(IS4)
c      T11=TLOG(IS3)
c      IF(T00.EQ.-100..OR.T01.EQ.-100..OR.T11.EQ.-100..OR.
c     %     T10.EQ.-100.) THEN
c         PRINT *,'BAD TEMPERATURES'
c         CALL CLEAN_STOP
c      ENDIF
c      WS5 = WP1*(WS1*T00 + WS2*T01)
c      WS6 = WP2*(WS3*T10 + WS4*T11)
c      T = 10.d0**(WS5+WS6)
      RETURN
      END SUBROUTINE FIND_SAUMONT


      subroutine calc_ST(T,P,S,S_H,S_HE,SMIX,dSMIXdT,ST)
      USE saumon_var_init
      IMPLICIT NONE
      double precision :: T,P,S,S_H,S_HE,SMIX,dSMIXdT,ST
      double precision :: dlogSHdlogT,dlogSHEdlogT
      call TABLEINTERPOLATE(T,P,HTLOG,HTABLE,NTMP,
     %        NPRS,9,dlogSHdlogT)
      call TABLEINTERPOLATE(T,P,HETLOG,HETABLE,NTMP,
     %        NPRS,9,dlogSHEdlogT)
      ST = (1.d0/S)*( (XX/S_H)*dlogSHdlogT + (YY/S_HE)*dlogSHEdlogT + 
     $     T*DSMIXdT )
      RETURN
      END subroutine calc_ST

      subroutine calc_SP(T,P,S,S_H,S_HE,SMIX,dSMIXdP,SP)
      USE saumon_var_init
      IMPLICIT NONE
      double precision :: T,P,S,S_H,S_HE,SMIX,dSMIXdP,SP
      double precision :: dlogSHdlogP,dlogSHEdlogP
      call TABLEINTERPOLATE(T,P,HTLOG,HTABLE,NTMP,
     %     NPRS,10,dlogSHdlogP)
      call TABLEINTERPOLATE(T,P,HETLOG,HETABLE,NTMP,
     %     NPRS,10,dlogSHEdlogP)
      SP = (1.d0/S)*( (XX/S_H)*dlogSHdlogP + (YY/S_HE)*dlogSHEdlogP + 
     $     P*DSMIXdP )
      RETURN
      END subroutine calc_SP


      subroutine smix_and_dsmix(T,P,SMIX,dSMIXdT,dSMIXdP)
      USE global_constants
      USE saumon_var_init
      implicit none
      double precision :: T,P,smix,dSMIXdT,dSMIXdP
      double precision :: XH2,XH,XHE,XHEp
      double precision :: XeH,XeHE
      double precision :: c1,c2,c3,c4
      double precision :: beta,gamma,delta
      double precision :: d1,d2,d3,d4,d5,d6
      double precision :: dXHEdT,dXHEpdT,dXH2dT,dXHdT
      double precision :: dXHEdP,dXHEpdP,dXH2dP,dXHdP
      double precision :: dGAMMAdT,dDELTAdT
      double precision :: dGAMMAdP,dDELTAdP
      double precision :: dXeHdT,dXeHEdT
      double precision :: dXeHdP,dXeHEdP
      double precision :: calc_dXHEdT,calc_dXHEpdT
      double precision :: calc_dXH2dT,calc_dXHdT
      double precision :: calc_dXHEdP,calc_dXHEpdP
      double precision :: calc_dXH2dP,calc_dXHdP
!--table values
      call TABLEINTERPOLATE(T,P,HTLOG,HTABLE,NTMP,
     %     NPRS,2,XH2)
      call TABLEINTERPOLATE(T,P,HTLOG,HTABLE,NTMP,
     %     NPRS,3,XH)
      call TABLEINTERPOLATE(T,P,HETLOG,HETABLE,NTMP,
     %     NPRS,2,XHE)
      call TABLEINTERPOLATE(T,P,HETLOG,HETABLE,NTMP,
     %     NPRS,3,XHEp)      
!--electron fractions
      XeH = 0.5*(1.d0-XH2-XH)
      XeHE = (1.d0/3.d0)*(2.d0-2.d0*XHE-XHEp)
!--constants from beta,gamma and delta
      c1 = 1.d0+XH+3.d0*XH2
      c2 = 1.d0+2.d0*XHE+XHEp
      c3 = 1.d0-XH2-XH
      c4 = 2.d0-2.d0*XHE-XHEp
      beta = MHYRODGEN*YY/(MHELIUM*XX)
      gamma = 1.5*c1/c2
      delta = (2.d0/3.d0)*(c4/c3)*beta*gamma
!--constants from smix
      d1 = 1.d0 + beta*gamma
      d2 = 1.d0 + delta
      d3 = 1.d0 + (1.d0/(beta*gamma))
      d4 = 1.d0 + (1.d0/delta)
      d5 = gamma*(beta*gamma+1.d0)
      d6 = delta*(delta+1.d0)
!--smix value
      smix = (XX/MHYRODGEN)*(2.d0/c1)*( log(d1)-XeH*log(d2)+
     +     beta*gamma*(log(D3)-XeHE*log(d4)) )
      smix = smix*KBOLTZ
!--numerical derivatives wrt T
      dXHEdT  = calc_dXHEdT(T,P)
      dXHEpdT = calc_dXHEpdT(T,P)
      dXH2dT  = calc_dXH2dT(T,P)
      dXHdT   = calc_dXHdT(T,P)
!--numerical derivatives wrt P
      dXHEdP  = calc_dXHEdP(T,P)
      dXHEpdP = calc_dXHEpdP(T,P)
      dXH2dP  = calc_dXH2dP(T,P)
      dXHdP   = calc_dXHdP(T,P)
!--gamma and delta derivatives wrt T
      dGAMMAdT = (1.5/c2)*(dXHdT+3.d0*dXH2dT)-
     +     (1.5*c1/(c2*c2))*(2.d0*dXHEdT+dXHEpdT)
      dDELTAdT = (delta/gamma)*dGAMMAdT - 
     +     beta*gamma*(2.d0/(3.d0*c3))*(2.d0*dXHEdT+dXHEpdT) +
     +     beta*gamma*(2.d0*c4/(3.d0*c3*c3))*(dXH2dT+dXHdT)
!--gamma and delta derivatives wrt P
      dGAMMAdP = (1.5/c2)*(dXHdP+3.d0*dXH2dP)-
     +     (1.5*c1/(c2*c2))*(2.d0*dXHEdP+dXHEpdP)
      dDELTAdP = (delta/gamma)*dGAMMAdP - 
     +     beta*gamma*(1.5/c3)*(2.d0*dXHEdP+dXHEpdP) +
     +     beta*gamma*(1.5*c4/(c3*c3))*(dXH2dP+dXHdP)
!--derivatives of electron fractions wrt T
      dXeHdT  = -0.5*(dXH2dT+dXHdT)
      dXeHEdT = (-1.d0/3.d0)*(2.d0*dXHEdT+dXHEpdT)
!--derivatives of electron fractions wrt P
      dXeHdP  = -0.5*(dXH2dP+dXHdP)
      dXeHEdP = (-1.d0/3.d0)*(2.d0*dXHEdP+dXHEpdP)
!--derivative of smix wrt temperature
      dSMIXdT =((-1.d0/(c1*c1))*(dXHdT+3.d0*dXH2dT)*
     +     (log(d1)-XeH*log(d2)+
     +     beta*gamma*(log(d3)-XeHE*log(d4)))) + (1.d0/c1)*( 
     +     (beta*dGAMMAdT/d1)-(dXeHdT*log(d2))-(XeH*dDELTAdT/d2)+
     +     (beta*dGAMMAdT*(log(d3)-XeHE*log(d4))) +
     +     beta*gamma*(-(dGAMMAdT/d5)-dXeHEdT*log(d4)+XeHE*dDELTAdT/d6))
      dSMIXdT = dSMIXdT*KBOLTZ*2.d0*XX/MHYRODGEN
!--derivative of smix wrt pressure
      dSMIXdP =((-1.d0/(c1*c1))*(dXHdP+3.d0*dXH2dP)*
     +     (log(d1)-XeH*log(d2)+
     +     beta*gamma*(log(d3)-XeHE*log(d4)))) + (1.d0/c1)*( 
     +     (beta*dGAMMAdP/d1)-(dXeHdP*log(d2))-(XeH*dDELTAdP/d2)+
     +     (beta*dGAMMAdP*(log(d3)-XeHE*log(d4))) +
     +     beta*gamma*(-(dGAMMAdP/d5)-dXeHEdP*log(d4)+XeHE*dDELTAdP/d6))
      dSMIXdP = dSMIXdP*KBOLTZ*2.d0*XX/MHYRODGEN
      RETURN
      END SUBROUTINE SMIX_AND_DSMIX

      function calc_dXHEdT(T,P)
      USE saumon_var_init
      implicit none
      double precision :: T,P
      double precision :: XHE1,XHE2,T1,T2
      double precision :: calc_dXHEdT
      T1 = T
      T2 = T+0.08
      call TABLEINTERPOLATE(T1,P,HETLOG,HETABLE,NTMP,
     %     NPRS,2,XHE1)
      call TABLEINTERPOLATE(T2,P,HETLOG,HETABLE,NTMP,
     %     NPRS,2,XHE2)
      calc_dXHEdT= (XHE2-XHE1)/(T2-T1)
      end

      function calc_dXHEdP(T,P)
      USE saumon_var_init
      implicit none
      double precision :: T,P
      double precision :: XHE1,XHE2,P1,P2
      double precision :: calc_dXHEdP
      P1 = P
      P2 = P+0.2
      call TABLEINTERPOLATE(T,P1,HETLOG,HETABLE,NTMP,
     %     NPRS,2,XHE1)
      call TABLEINTERPOLATE(T,P2,HETLOG,HETABLE,NTMP,
     %     NPRS,2,XHE2)
      calc_dXHEdP= (XHE2-XHE1)/(P2-P1)
      end


      function calc_dXHEpdT(T,P)
      USE saumon_var_init
      implicit none
      double precision :: T,P
      double precision :: XHEp1,XHEp2,T1,T2
      double precision :: calc_dXHEpdT
      T1 = T
      T2 = T+0.08
      call TABLEINTERPOLATE(T1,P,HETLOG,HETABLE,NTMP,
     %     NPRS,3,XHEp1)
      call TABLEINTERPOLATE(T2,P,HETLOG,HETABLE,NTMP,
     %     NPRS,3,XHEp2)
      calc_dXHEpdT= (XHEp2-XHEp1)/(T2-T1)
      end

      function calc_dXHEpdP(T,P)
      USE saumon_var_init
      implicit none
      double precision :: T,P
      double precision :: XHEp1,XHEp2,P1,P2
      double precision :: calc_dXHEpdP
      P1 = P
      P2 = P+0.2
      call TABLEINTERPOLATE(T,P1,HETLOG,HETABLE,NTMP,
     %     NPRS,3,XHEp1)
      call TABLEINTERPOLATE(T,P2,HETLOG,HETABLE,NTMP,
     %     NPRS,3,XHEp2)
      calc_dXHEpdP= (XHEp2-XHEp1)/(P2-P1)
      end

      function calc_dXH2dT(T,P)
      USE saumon_var_init
      implicit none
      double precision :: T,P
      double precision :: XH21,XH22,T1,T2
      double precision :: calc_dXH2dT
      T1 = T
      T2 = T+0.08
      call TABLEINTERPOLATE(T1,P,HTLOG,HTABLE,NTMP,
     %     NPRS,2,XH21)
      call TABLEINTERPOLATE(T2,P,HTLOG,HTABLE,NTMP,
     %     NPRS,2,XH22)
      calc_dXH2dT= (XH22-XH21)/(T2-T1)
      end

      function calc_dXH2dP(T,P)
      USE saumon_var_init
      implicit none
      double precision :: T,P
      double precision :: XH21,XH22,P1,P2
      double precision :: calc_dXH2dP
      P1 = P
      P2 = P+0.2
      call TABLEINTERPOLATE(T,P1,HTLOG,HTABLE,NTMP,
     %     NPRS,2,XH21)
      call TABLEINTERPOLATE(T,P2,HTLOG,HTABLE,NTMP,
     %     NPRS,2,XH22)
      calc_dXH2dP= (XH22-XH21)/(P2-P1)
      end

      function calc_dXHdT(T,P)
      USE saumon_var_init
      implicit none
      double precision :: T,P
      double precision :: XH1,XH2,T1,T2
      double precision :: calc_dXHdT
      T1 = T
      T2 = T+0.08
      call TABLEINTERPOLATE(T1,P,HTLOG,HTABLE,NTMP,
     %     NPRS,3,XH1)
      call TABLEINTERPOLATE(T2,P,HTLOG,HTABLE,NTMP,
     %     NPRS,3,XH2)
      calc_dXHdT= (XH2-XH1)/(T2-T1)
      end

      function calc_dXHdP(T,P)
      USE saumon_var_init
      implicit none
      double precision :: T,P
      double precision :: XH1,XH2,P1,P2
      double precision :: calc_dXHdP
      P1 = P
      P2 = P+0.2
      call TABLEINTERPOLATE(T,P1,HTLOG,HTABLE,NTMP,
     %     NPRS,3,XH1)
      call TABLEINTERPOLATE(T,P2,HTLOG,HTABLE,NTMP,
     %     NPRS,3,XH2)
      calc_dXHdP= (XH2-XH1)/(P2-P1)
      end


      SUBROUTINE TABLEINTERPOLATE(TT,PP,TVAL,TABLE,NT,NP,COL,VALUE)
      IMPLICIT NONE
c     COL tells which value to read from the table
c     4=density, 5=entropy
      INTEGER :: NT,NP,COL,IT,IT1,IW,IW1,K
      DOUBLE PRECISION :: TT,PP
      DOUBLE PRECISION :: TVAL(NT),TABLE(NT,NP,11)
      DOUBLE PRECISION :: TLOG,PLOG,VALUE,Z00,Z01,Z10,Z11
      DOUBLE PRECISION :: WS1,WS,WS2,WS3,WS4,WS5
c..   convert to actual logT, logP
      TLOG = LOG10(TT)
      PLOG = LOG10(PP)
c.....Start here with logT,logP, Reference logTs logPS
      IF(TLOG .LT. 2.10) GO TO 200 
      IF(TLOG .GT. 7.06) GO TO 200 
      IF(PLOG .LT. 4.00) GO TO 200 
      IF(PLOG .GT. 19.00) GO TO 200 
      DO K = 1,NT
         IF(TVAL(K) .GT. TLOG) GO TO 15 
      ENDDO
 15   IT = K
      IT1 = K-1
      WS1 = TVAL(IT1) - TLOG
      WS = TVAL(IT1) - TVAL(IT)
      WS1 = WS1 / WS
      WS = 1. - WS1
      DO K = 1,NP
         IF(TABLE(NT,K,1) .GT. PLOG) GO TO 25
      ENDDO
 25   IW = K
      IW1 = K-1
      WS2 = PLOG - TABLE(NT,IW1,1)
      WS3 = TABLE(NT,IW,1) - TABLE(NT,IW1,1)
      WS2 = WS2 / WS3
      WS3 = 1. - WS2
      Z00 = TABLE(IT,IW,COL)
      Z10 = TABLE(IT,IW1,COL)
      Z11 = TABLE(IT1,IW1,COL)
      Z01 = TABLE(IT1,IW,COL)
      IF(Z00.EQ.-100..OR.Z01.EQ.-100..OR.Z11.EQ.-100..OR.Z10.EQ.-100.)
     %     GO TO 200
      WS4 = WS2*(WS1*Z00 + WS*Z01)
      WS5 = WS3*(WS1*Z10 + WS*Z11)
      VALUE = 10.d0**(WS4+WS5)
      RETURN
 200  PRINT 201,TT,PP,TLOG,PLOG
 201  FORMAT(18H OUTSIDE TABLE,T =, 1PE12.3, 6H, P = , E12.3,
     X     7H, TL = ,0PF8.4,7H, PL = , F8.4)
      CALL CLEAN_STOP
      END


      SUBROUTINE READTABLE(TABLE,NTMP,NPRS,TLOG,file)
      USE mpi_var_init
      IMPLICIT NONE
      integer NTMP,NPRS,I,J,K
      double precision NPOINT(NTMP)
      DOUBLE PRECISION TLOG(NTMP),TABLE(NTMP,NPRS,11)
      character*16 file
      if(myid.eq.0) WRITE(6,*) "Reading ", file
      DO I=1,NTMP
         DO J=1,NPRS
            DO K=1,11
               TABLE(I,J,K) = -100.0
            ENDDO
         ENDDO
      ENDDO
      OPEN(UNIT=41,FILE=file,status="old")
      DO I=1,NTMP
         READ(41,*) TLOG(I),NPOINT(I)
         READ(41,*) ((TABLE(I,J,K),K=1,11),J=1,NPOINT(I))
      ENDDO
      close(41)
      return
      END
