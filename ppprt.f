      SUBROUTINE PRINTOUT(PRINTVAR)
      USE fluid_var_init
      USE grid_var_init
      USE rad_var_init
      USE mpi_var_init
      USE input_init
      USE force_var_init
      USE hypre_var_init
      USE deltaER_var_init
      IMPLICIT NONE
      integer :: i,j,k,printvar
      SELECT CASE(PRINTVAR)
      CASE (1)
!---  PRINTVAR=1 :  Recent Dataset  --> AKTTD or SNAP.xxx
         WRITE (12) global_NX,global_NY,global_NZ
         WRITE (12) locNX,locNY,locNZ,numprocs
         WRITE (12) myid,proc_coords(1),proc_coords(2),proc_coords(3)
         WRITE (12) proc_dims(1),proc_dims(2),proc_dims(3)
         WRITE (12) ZEIT,DELT,NUM_ITER
         WRITE (12) (XXA(I),I=0,locNX+2)
         WRITE (12) (XXB(I),I=-1,locNX+2)
         WRITE (12) (DXA(I),I=1,locNX+2)
         WRITE (12) (XYA(J),J=0,locNY+2)
         WRITE (12) (XYB(J),J=-1,locNY+2)
         WRITE (12) (DYA(J),J=1,locNY+2)
         WRITE (12) (XZA(K),K=0,locNZ+2)
         WRITE (12) (XZB(K),K=-1,locNZ+2)
         WRITE (12) (DZA(K),K=1,locNZ+2)
         WRITE (12) (((V(I,J,K),I=0,locNX+2),J=-1,locNY+2),K=-1,locNZ+2)
         WRITE (12) (((G(I,J,k),I=-1,locNX+2),J=0,locNY+2),K=-1,locNZ+2)
         WRITE (12) (((H(I,J,k),I=-1,locNX+2),J=-1,locNY+2),K=0,locNZ+2)
         WRITE (12) (((RH(I,J,K),I=-1,locNX+2),J=-1,locNY+2),
     %        K=-1,locNZ+2)
         WRITE (12) (((T(I,J,K),I=-1,locNX+2),J=-1,locNY+2),
     %        K=-1,locNZ+2)
         WRITE (12) (((ER(I,J,K),I=-1,locNX+2),J=-1,locNY+2),
     %        K=-1,locNZ+2)
         WRITE (12) (((PG(I,J,K),I=-1,locNX+2),J=-1,locNY+2),
     %        K=-1,locNZ+2)
         WRITE (12) (((CS(I,J,K),I=0,locNX+1),J=0,locNY+1),K=0,locNZ+1)
         WRITE (12) ((UPPERBND(J,K),J=1,locNY),K=1,locNZ)
         WRITE (12) TinitBottom,Tinittop,rhinittop,rhinitbottom,
     %        hydrotime,Maccreated
!--- NOW COMES STUFF NOT READ BACK INTO A RESTART
         WRITE (12) NCOSYS
         WRITE (12) (((xkapR(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
         WRITE (12) (((flimx(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)

c         WRITE (12) (((EdepFor(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((flux_Y(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
         WRITE (12) (((conv_flux(1,I,J,K),I=1,locNX),J=1,locNY),
     $        K=1,locNZ)
c         WRITE (12) (((test_roche(1,I,J,K),I=1,locNX),J=1,locNY),
c     $        K=1,locNZ)
c         WRITE (12) (((gravx(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((test_vadv(1,I,J,K),I=1,locNX),J=1,locNY),
c     $        K=1,locNZ)
c         WRITE (12) (((divtz(I,J,K),I=1,locNX),J=1,locNY),
c     %        K=1,locNZ)
c         WRITE (12) (((test_xcotime(I,J,K),I=1,locNX),J=1,locNY),
c     $        K=1,locNZ)
c         WRITE (12) (((ARTVIS_T(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)



c         WRITE (12) (((EdepAdv(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((flux_X(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
C         WRITE (12) (((test_roche(2,I,J,K),I=1,locNX),J=1,locNY),
C     $        K=1,locNZ)
c         WRITE (12) (((gravy(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((flux(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((flimy(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((flimz(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((centcorx(i,j,k),I=1,locNX),J=1,locNY),
c     %        K=1,locNZ)
c         WRITE (12) (((test_vadv(2,I,J,K),I=1,locNX),J=1,locNY),
c     $        K=1,locNZ)
         WRITE (12) (((conv_flux(2,I,J,K),I=1,locNX),J=1,locNY),
     $        K=1,locNZ)
c         WRITE (12) (((ARTVIS_V(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)


c         WRITE (12) (((EdepRad(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((flux_Z(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((difrx(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((gravz(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((divtx(I,J,K),I=1,locNX),J=1,locNY),
c     %        K=1,locNZ)
c         WRITE (12) (((gravz(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)
c         WRITE (12) (((test_roche(3,I,J,K),I=1,locNX),J=1,locNY),
c     $        K=1,locNZ)
c         WRITE (12) (((test_vadv(3,I,J,K),I=1,locNX),J=1,locNY),
c     $        K=1,locNZ)
c         WRITE (12) (((test_xcoeq(I,J,K),I=1,locNX),J=1,locNY),
c     $        K=1,locNZ)
         WRITE (12) (((conv_flux(3,I,J,K),I=1,locNX),J=1,locNY),
     $        K=1,locNZ)
c         WRITE (12) (((ARTVIS_G(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)

c-------
         WRITE (12) (((xkapA(I,J,K),I=1,locNX),J=1,locNY),
     $        K=1,locNZ)
c---- new stuff
         WRITE (12) (((STELLARINPUT(I,J,K),I=1,locNX),J=1,locNY),
     %        K=1,locNZ)
c         WRITE (12) (((divty(I,J,K),I=1,locNX),J=1,locNY),
c     %        K=1,locNZ)
c         WRITE (12) (((gravy(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)


c         WRITE (12) (((delta_ER(I,J,K),I=1,locNX),J=1,locNY),
c     %        K=1,locNZ)
         WRITE (12) (((disfn(I,J,K),I=1,locNX),J=1,locNY),
     %        K=1,locNZ)
c         WRITE (12) (((stellartau(I,J,K),I=1,locNX),J=1,locNY),
c     %        K=1,locNZ)
c         WRITE (12) (((saveRRX(I,J,K),I=1,locNX),J=1,locNY),
c     %        K=1,locNZ)
c         WRITE (12) (((visc_flux(1,i,j,k),I=1,locNX),J=1,locNY),
c     %        K=1,locNZ)
c         WRITE (12) (((ARTVIS_H(I,J,K),I=1,locNX),J=1,locNY),K=1,locNZ)         

         WRITE (12) (((xkapP(I,J,K),I=1,locNX),J=1,locNY),
     %        K=1,locNZ)


      CASE (2:)
         print *,'printvar=',printvar,'not yet defined: PRINTOUT' 
      CASE DEFAULT
         print *,'printvar=',printvar,'also not defined: PRINTOUT' 
      END SELECT
      RETURN
      END SUBROUTINE PRINTOUT


      SUBROUTINE PRINTSCALAR
      USE scalar_var_init
      USE fluid_var_init
      IMPLICIT NONE
      integer :: i,j,k,n
      WRITE (12) nscalar
      do n=1,nscalar
         WRITE (12) (((scalar(N,I,J,K),I=-1,locNX+2),J=-1,locNY+2),
     %        K=-1,locNZ+2)
      enddo      
      RETURN
      END SUBROUTINE PRINTSCALAR





      SUBROUTINE SCREEN_TIME_PRINTOUT
      USE input_init
      USE fluid_var_init
      USE global_constants
      IMPLICIT NONE
      if(mstar2.eq.0) then
         write(*,'(A,I8,3(1x,e12.6))') 
     &        'N,dt,t,hydrot(d)',num_iter,delt,zeit,
     $        hydrotime/DAY
      else
c         write(*,'(A,I8,4(1x,e12.6))') 
c     &        'N,dt,t,hydrot(d),Mac',num_iter,delt,zeit,
c     %        hydrotime/DAY,Maccreated
         write(*,'(A,I8,4(1x,e12.6))') 
     &        'N,dt/fdelt,t,hydrot(d),Mac',num_iter,delt/fdelt,zeit,
     %        hydrotime/DAY,Maccreated
      endif
      RETURN
      END SUBROUTINE SCREEN_TIME_PRINTOUT



      SUBROUTINE MATRIX_ITERATIONS_PRINTOUT
      USE input_init
      USE rad_var_init
      USE sor_var_init
      IMPLICIT NONE      
      write(*,'(A,2(I8,1x),2(e12.6,1x))') 
     $     ' coupledenergy/SOR iterations,conv/SOR param,: ',
     %     energy_iter,SOR_ITER,conv_parm,SORPARAM

      RETURN
      END SUBROUTINE MATRIX_ITERATIONS_PRINTOUT
