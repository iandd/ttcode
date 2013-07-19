1000000,10001,0.98,1.e12		ISTEPMAX,IREA,FDELT,DTMAX
1000,0				iprstep,FILESTART
0,0,0                   	ietot,isoth,isoden
1,1,1,1,3,0,0,0             	IADV,IFOR,IGRAV,IVIS,IRAD,INUC,ISCAL,IDIF
1,1.e-11			IRHO,RHMIN
0,1.e+4,1.e1	          	ITEMP,TEMPMAX,TEMPMIN
0,6.e5,0.0,100000		IVEL,VELMAX,VELDAMP,DAMP_STEPS
1,1.e9,1,1.0e-4,4,1,0.0		IVISTYP,XNUE,IDISS,alpha,IETA,IARTVIS,PR_NUM
1,1.4,2.3,0,0.04 		IGASTYP,GAMMA,mu_gas,ICSTYP,aspect
21,0,1.0,0.1		        IKAPTYP,ISIGTYP,ZKAPMULT,gamma_OPC
1,3.52,0.0,0.0,1.0e+0		IROTATE,PROT,ramp_start,ramp_full,windAMP
9,2,2,9,15844.230		IBDRAD,IFLIMTYP,F_INCIDENT,BCXMIN,FBOT
1.71,500,1.0e-7			SORPARAM,ITMAX,EPSMAT
***********************  More Switches ******************
0,0,100			        iperturb,perturb_typ,perturb_iter
3.52,0.0,0.0			Porb,ecc,ecc_start(days)
0,3000,2000.0			iFbotramp,Framp_steps,Fbot_final
-1,-1.0,-1.0,-1.0		tmp,tmp,tmp,tmp
***********************  Fixed Modelparameter ***********
2,1                            	NCOSYS,inc_poles
2,3,1				MODTYP,MODVER,PLANET_NUM
0.,1.155,6070.0,1.119,2048.0,0.	STARTYP,RSTAR,TSTAR,MSTAR1,TIRR,RH0
0.e0,0.000654,0.0,1.0,1.0	RSTAR2,mstar2,rp,phip,psoft
3      				numdim
200,8.50e9,1.03e10,0.0		global_NX,XMIN,XMAX,DXMIN
160,0.00e0,1.0,0.0		global_NY,YMIN,YMAX,DYMIN
64,-0.95,0.95,0.0		global_NZ,ZMIN,ZMAX,DZMIN