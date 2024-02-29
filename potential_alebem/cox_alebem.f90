!===============================================================================
!	SOLVE 2D WATER WAVE PROBLEM WITH ONE SIDE WAVEMAKER AND ONESIDE FREE BY BEM
!	WEN-HAUI TSAO 2021 LSU
!===============================================================================
      PROGRAM WATER_WAVE
      IMPLICIT NONE
      INTEGER SMTYP,NNRL,MDEG
	  INTEGER I,J,NT,NGA,NARC,NTIM,OUTSTEP,NFIELD,NPL,WTYP,NNODE,NELM,ICON,NITER,NWG,E1LOC,OUTYP,ALETYP      
	  INTEGER,ALLOCATABLE::NELEM(:),ME(:),NS(:),BUNTYP(:),LN(:,:),IPIV(:)
	  REAL*8 DEP,CLEN,D_OUT,GRAV,MU,ARDZONE,WIDTH,THO,DELTTIME,WGX(10)
	  REAL*8 TIME,DIS,VEL,ACC,P_ATM,FOR,DDIS,DTEMP,WC,E1,E2,ETOL
      REAL*8 WAVE_PERIOD,WAVE_HEIGHT,PSI,AMP,OMEGA,WAVE_K,WAVE_C
	  REAL*8,ALLOCATABLE::COOR(:,:),SIDE_L(:)
	  REAL*8,ALLOCATABLE::NODE(:,:),NORM(:,:),JCB(:),LENG(:),PHI(:),PPHI(:),PHIT(:),PPHIT(:)
	  REAL*8,ALLOCATABLE::KER1(:,:),KER2(:,:),DP(:,:),DPDS(:),DPDSS(:),PR(:),DPDT(:),ACCMO(:)
	  REAL*8,ALLOCATABLE::PHIT_TEMP(:)
	  REAL*8,ALLOCATABLE::WT(:),RT(:),SHA1(:),SHA2(:),SH(:,:),AWT(:),ART(:)
      REAL*8,ALLOCATABLE::G1(:,:),H1(:,:),EYE(:)
	  REAL*8,ALLOCATABLE::DI(:),VE(:),AC(:)
      REAL*8,ALLOCATABLE::DPM(:,:),DPDTM(:),CV(:,:),D2PM(:,:),D2PDTM(:)

      OPEN(UNIT=1,FILE='1.ipt',STATUS='OLD')
      OPEN(UNIT=2,FILE='2.ipt',STATUS='OLD')
      OPEN(UNIT=3,FILE='3.ipt',STATUS='OLD')
      OPEN(UNIT=5,FILE='IO.DAT')
      OPEN(UNIT=6,FILE='S.DAT')
      OPEN(UNIT=7,FILE='WG.DAT')
!      OPEN(UNIT=8,FILE='P.DAT')
!      OPEN(UNIT=9,FILE='F.DAT')
!      OPEN(UNIT=10,FILE='E.DAT')
!      OPEN(UNIT=11,FILE='DOMAIN.DAT')
	  OPEN(UNIT=21,FILE='ERR.DAT')
	  OPEN(UNIT=22,FILE='ABORT.TXT')
	  OPEN(UNIT=23,FILE='CFL.DAT')
      OPEN(UNIT=99,FILE='TEST.TXT')

!---SMOOTHING PARAMETERS
    CALL INPUT_3(SMTYP,NNRL,MDEG,NARC,ALETYP)
    ALLOCATE(AWT(NARC),ART(NARC))
    AWT=0.D0
    ART=0.D0
    CALL GAUSS(AWT,ART,NARC)
    
!---TOPOGRAGHY AND WAVE TYPE
	CALL INPUT_2(NPL,WTYP,OUTYP,NWG,WGX,WAVE_PERIOD,WAVE_HEIGHT,PSI)
	ALLOCATE(NELEM(NPL),ME(NPL),NS(NPL),BUNTYP(NPL),COOR(NPL,2),SIDE_L(NPL))
	NELEM=0
	NS=0
	BUNTYP=0
	COOR=0.D0
	SIDE_L=0.D0

!---INPUT ALL KINDS OF PARAMETERS
	CALL INPUT_1(NPL,COOR,NFIELD,NNODE,NELM,NELEM,ME,NS,BUNTYP,DEP,CLEN,NGA,GRAV,MU,ARDZONE,WIDTH,THO,&
			  &NTIM,DELTTIME,OUTSTEP,ICON,NITER,ETOL)
	ALLOCATE(LN(NELM,2),NODE(NNODE,2),NORM(NELM,2),JCB(NELM),LENG(NELM))
	ALLOCATE(PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),PHIT_TEMP(NNODE))
	ALLOCATE(KER1(NNODE,NNODE),KER2(NNODE,NNODE),DP(NNODE,2))
	ALLOCATE(DPDS(NNODE),DPDSS(NS(1)),DPDT(NNODE),PR(NNODE),ACCMO(NNODE))
	ALLOCATE(WT(NGA),RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA))
	ALLOCATE(DI(NTIM),VE(NTIM),AC(NTIM))
    ALLOCATE(G1(NNODE,NNODE),H1(NNODE,NNODE),EYE(NNODE),IPIV(NNODE))
    ALLOCATE(DPM(NS(1),2),DPDTM(NS(1)),CV(NS(1),2),D2PM(NS(1),2),D2PDTM(NS(1)))
	LN=0
	NODE=0.D0
	NORM=0.D0
	JCB=0.D0
	LENG=0.D0
	PHI=0.D0
	PPHI=0.D0
	PHIT=0.D0
	PPHIT=0.D0
	KER1=0.D0
	KER2=0.D0
	DP=0.D0
	DPDS=0.D0
    DPDSS=0.D0
	DPDT=0.D0
	PR=0.D0
	ACCMO=0.D0
	PHIT_TEMP=0.D0
	WT=0.D0
	RT=0.D0
	SHA1=0.D0
	SHA2=0.D0
	SH=0.D0
	DI=0.D0
	VE=0.D0
	AC=0.D0
    G1=0.D0
    H1=0.D0
    EYE=1.D0
    IPIV=0
    DPM=0.D0
    DPDTM=0.D0
    CV=0.D0
    D2PM=0.D0
    D2PDTM=0.D0
    
!---PREPARE IF OUTLET IS A WALL (IF A WALL, NO ITERATION NEEDED)
    IF (OUTYP==0)THEN
        NITER=1
    END IF
        
!---SHAPE FUNCTION AND MESH
    CALL GAUSS(WT,RT,NGA)
    CALL SHAP(SHA1,SHA2,SH,NGA,RT)
	CALL LENGTH(NPL,COOR,SIDE_L)
	CALL MESH(NPL,NNODE,NELM,NELEM,LN,COOR,SIDE_L,NODE)
    WRITE(*,*) 'PASS MESH'

!---GENERATE WAVES
	SELECT CASE (WTYP)
    CASE(1)
        CALL PERIODIC(DEP,WAVE_PERIOD,WAVE_HEIGHT,AMP)
        OMEGA=2.D0*DACOS(-1.D0)/WAVE_PERIOD
		DO I=1,NTIM
		 TIME=(I-1)*DELTTIME
		 DI(I)=AMP*DCOS(OMEGA*TIME+PSI)-AMP !AMP*SIN(OMEGA*TIME+PSI)
		 VE(I)=-AMP*OMEGA*DSIN(OMEGA*TIME+PSI) !AMP*OMEGA*COS(OMEGA*TIME+PSI)
		 AC(I)=-AMP*OMEGA**2*DCOS(OMEGA*TIME+PSI) !-AMP*OMEGA**2*SIN(OMEGA*TIME+PSI)
		 END DO
	CASE(2)
		NT=OMEGA/DELTTIME
		DO I=1,NT+1
		TIME=(I-1)*DELTTIME
		DI(I)=0.5D0*AMP-0.5D0*AMP*DCOS(DACOS(-1.D0)/OMEGA*TIME) !STR/DUR*TIME
		VE(I)=DACOS(-1.D0)/OMEGA*0.5D0*AMP*DSIN(DACOS(-1.D0)/OMEGA*TIME) !STR/DUR
		AC(I)=(DACOS(-1.D0)/OMEGA)**2*0.5D0*AMP*DCOS(DACOS(-1.D0)/OMEGA*TIME) !0.D0
		END DO
		DI(NT+2:NTIM)=DI(NT+1)
    CASE(3)
        NT=10.D0/DELTTIME
        CALL SOLITARY(DEP,WAVE_HEIGHT,GRAV,1,1,NT,DELTTIME,DI,VE,AC)
        DI(NT+1:NTIM)=DI(NT)
        VE(NT+1:NTIM)=VE(NT)
        AC(NT+1:NTIM)=AC(NT)

    END SELECT  
    
DO NT=1,NTIM
	WRITE(*,*) NT,'TH'
     TIME=(NT-1)*DELTTIME
!---WAVEMAKER DIS, VEL, ACC, PHASE LAG (IN RADIUS)
	 DIS=DI(NT)
     VEL=VE(NT)
     ACC=AC(NT)
	 WRITE(5,"(7(E15.8,1X))") TIME,DIS,VEL,ACC

!---CALCULATE WAVE SPEED FOR RADIATION CONDITION
	D_OUT=NODE(NS(1),2)-COOR(NPL-2,2)
	CALL WAVE_SPD(GRAV,OMEGA,D_OUT,WC)

!---REMESH LOCATION AND POTENTIAL OF A FREE-SURFACE NODE AT CURRENT (NT-1)TH TIME
!---THE FOLLOWING CALCULATION IS FOR THE NEXT TIME-STEP)
      CALL REMESH(NPL,NNODE,NELEM,NODE,NS,TIME,DEP,AMP,NWG,WGX)
      WRITE(*,*) 'PASS REMESH'
!---OUTPUT BOUNDARY
      IF(MOD(NT,OUTSTEP)==1)THEN
          WRITE(6,"(5000(E15.8,1X))") NODE(:,1) 
          WRITE(6,"(5000(E15.8,1X))") NODE(:,2)
      END IF
      
!---BUILD KERNEL OF ZONE1, ZONE2, ZONE3 AND ASSEMBLE THEM INTO BIG KERNEL
	   CALL KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT,EYE)
       WRITE(*,*) 'PASS KERNEL'

!---CALCULATE PRESSURE ON THE BOUNDARY
	  CALL PRESSURE(ICON,TIME,THO,GRAV,DEP,NPL,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)

!---CALCULATE PRESSURE AND VELOCITY IN THE DOMAIN
!	  CALL DOMAIN(NPL,NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,&
!				 &SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)

!**************IMPLICIT SOLVER FOR PHIT AND PPHIT ON THE ABSORPTION SIDE**************
DO J=1,NITER

	  PHIT_TEMP=PHIT
!---APPLY BC FOR SOLVING FREE-SURFACE PPHI
       CALL BOUND(NPL,NNODE,NELEM,NS,BUNTYP,OUTYP,PHI,PPHI,PHIT,VEL,WC)
!	   CALL SOLVE_LAPACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!	   CALL SOLVE_BACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
        CALL SOLVE_LAPACK2_1(NPL,PHI,PPHI,KER1,KER2,G1,H1,NNODE,NELEM,BUNTYP,IPIV)

!---FIRST-ORDER TAYLOR SERIES EXPANSION (ALSO GET TANGENTIAL VELOCITY ON ZONE BOUNDARY)
!	  CALL TAYES1(NPL,NNODE,NELM,NELEM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,DPDS,DP,PHIT,DPDT,&
!				 &DEP,CLEN,GRAV,MU,ARDZONE,VEL,WC)
	  CALL TAYES1_CS(NPL,NNODE,NELM,NELEM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,DPDS,DPDSS,DP,PHIT,DPDT,&
				 &DEP,CLEN,GRAV,MU,ARDZONE,VEL,WC,AWT,ART,NARC)

!---APPLY BC FOR SOLVING FREE-SURFACE PPHIT
      CALL ACCBC(NPL,NNODE,NELM,NELEM,ME,LN,NORM,JCB,PPHI,DP,ACCMO)
      CALL BOUNDT(NPL,NNODE,NELM,NS,NELEM,BUNTYP,OUTYP,LN,PHIT,PPHIT,DPDS,JCB,WC,ACC,ACCMO)
!	  CALL SOLVE_LAPACK(NPL,PHIT,PPHIT,KER1,KER2,NNODE,NELEM,BUNTYP)
!	  CALL SOLVE_BACK(NPL,PHIT,PPHIT,KER1,KER2,NNODE,NELEM,BUNTYP)
      CALL SOLVE_LAPACK2_2(NPL,PHIT,PPHIT,G1,H1,NNODE,NELEM,BUNTYP,IPIV)

IF (OUTYP==1)THEN
CALL CONVERGE(NELEM(2)+1,PHIT_TEMP(NS(1)+1:NS(2)),PHIT(NS(1)+1:NS(2)),E1LOC,E1,E2)

	IF(E1<=ETOL.AND.E2<ETOL)THEN
		WRITE(*,*) 'PASS CONVERGED',J
		WRITE(21,*) TIME,J,E1,E2
		GOTO 205
	ELSE IF(J>=NITER)THEN
		WRITE(22,*) 'CG FAIL',TIME,E1LOC,E1,E2
		WRITE(*,*) 'CONVERGE FAIL'
		STOP
	END IF
END IF

END DO
!**************************************************************************************

205 CONTINUE

!---CHECK THE CFL NUMBER	
!	  CALL COURANT(TIME,DELTTIME,NNODE,NELM,LN,NODE,DP,JCB)
    
!**************FREE SURFACE UPDATING IN EL(ALETYP=0) OR ALE(ALETYP=1) FRAME**************    
    IF (ALETYP==0)THEN
!---SECOND-ORDER TAYLOR SERIES EXPANSION IN EL FRAME
        CALL TAYES2(PHI,PPHI,PHIT,PPHIT,DPDS,DPDT,DP,NPL,NNODE,NELEM,NODE,NELM,NORM,JCB,LN,DELTTIME,GRAV,ACC)
        
    ELSE IF (ALETYP==1)THEN
!---CALCULATE MESH VELOCITY
        CALL ALE_VEL(NPL,NNODE,NELM,NELEM,NS,LN,NODE,NORM,DP,DPDT,DPM,DPDTM,CV)
    
!---CALCULATE MESH ACCELERATION
        CALL ALE_ACC(DELTTIME,GRAV,MU,CLEN,ARDZONE,NPL,NNODE,NELEM,NS,NELM,LN,NODE,NORM,JCB,&
                    &PHI,PPHI,PHIT,PPHIT,DP,DPDS,DPDSS,DPDT,CV,D2PM,D2PDTM)
        WRITE(*,*) 'PASS ALE'        
        
!---TSE SUM IN ALE FRAME
        CALL ALE_TSE_SUM(DELTTIME,NPL,NNODE,NS,NODE,PHI,DPM,DPDTM,D2PM,D2PDTM)
    END IF
    
WRITE(*,*) 'PASS TAYES'

    
    
    
    
!---SMOOTHING FOR FREE-SURFACE NODE AND PHI
      IF (SMTYP==1)THEN
          CALL FS_SMOOTH(NNRL,MDEG,NS(1),NODE(1:NS(1),2))
          CALL FS_SMOOTH(NNRL,MDEG,NS(1),PHI(1:NS(1)))
          WRITE(*,*) 'PASS SMOOTH'
      END IF
END DO

STOP
    END


!**********************************************************************
      SUBROUTINE HEADLINE(ID,IREAD)
!**********************************************************************
      CHARACTER*2 ID
      ID =  '*'
      DO WHILE (ID .EQ. '*')
      READ(IREAD,'(A1)') ID
      END DO
      RETURN
      END
!**********************************************************************
SUBROUTINE INPUT_1(NPL,COOR,NFIELD,NNODE,NELM,NELEM,ME,NS,BUNTYP,DEP,CLEN,NGA,GRAV,MU,ARDZONE,WIDTH,THO,&
				&NTIM,DELTTIME,OUTSTEP,ICON,NITER,ETOL)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER I,J,NPL,NFIELD,NNODE,NELM,NGA,NTIM,ICON,NITER,OUTSTEP
	  INTEGER NELEM(NPL),ME(NPL),NS(NPL),BUNTYP(NPL)
	  REAL*8 COOR(NPL,2),DEP,CLEN,ENDTIME,DELTTIME,GRAV,MU,ARDZONE,WIDTH,THO,ETOL
      CHARACTER*2 ID
         ID = '*'
!---NODAL COORDINATES
         CALL HEADLINE(ID,1)
         READ(1,*)  ((COOR(I,J),J=1,2),I=1,NPL)
!---ELEMENT MESH NUMBER
         CALL HEADLINE(ID,1)
        READ(1,*) (NELEM(I),I=1,NPL)
        NNODE = 0
        NELM  = 0
		ME(1) = NELEM(1)
		NS	  = 0
        DO I=1,NPL
         NELM = NELM+NELEM(I)
         NNODE = NNODE+(NELEM(I)+1)
        END DO

		DO I=2,NPL
		 ME(I)=ME(I-1)+NELEM(I)
		END DO

      NS(1)=NELEM(1)+1
	  DO I=2,NPL
      NS(I)=NS(I-1)+NELEM(I)+1
	  END DO

	  NFIELD=(NELEM(1)-1)*(NELEM(NPL)-1)

!---BOUNDARY TYPE
         CALL HEADLINE(ID,1)
        READ(1,*) (BUNTYP(I),I=1,NPL)
!---READ THE GAUSSIAN INTEGRATION POINT NO.
         CALL HEADLINE(ID,1)
         READ(1,*)  NGA
!---READ GRAV ACC, MAX OF WAVE DAMPING COEF, ZONE OF WAVE DAMPING, TANK WIDTH, FLUID DENSITY
         CALL HEADLINE(ID,1)
         READ(1,*)  GRAV,MU,ARDZONE,WIDTH,THO
!---READ TIME,TIME INTERVAL,OUTPUT STEP,INDEX OF BERNOULLI CONSTANT, NUMBER IF ITERATION,TOLERANCE ERROR
         CALL HEADLINE(ID,1)
         READ(1,*)  ENDTIME,DELTTIME,OUTSTEP,ICON,NITER,ETOL
		NTIM=ENDTIME/DELTTIME+1
		DEP=MAXVAL(COOR(:,2))
        CLEN=MAXVAL(COOR(:,1))

      RETURN
      END
!**********************************************************************
      SUBROUTINE INPUT_2(NPL,WTYP,OUTYP,NWG,WGX,WAVE_PERIOD,WAVE_HEIGHT,PSI)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER NPL,WTYP,OUTYP,NWG
	  REAL*8 WAVE_PERIOD,WAVE_HEIGHT,PSI,WGX(10)
      CHARACTER*2 ID
         ID = '*'
!---NUMBER OF PLANES
         CALL HEADLINE(ID,2)
         READ(2,*) NPL
!---WAVE GENERATION: 1=PERIODIC; 2=SOLITARY
         CALL HEADLINE(ID,2)
         READ(2,*) WTYP
!---READ WAVE_PERIOD,WAVE_HEIGHT
         CALL HEADLINE(ID,2)
         READ(2,*)  WAVE_PERIOD,WAVE_HEIGHT,PSI
!---WAVE GENERATION: 0=WALL; 1=RADIATION
         CALL HEADLINE(ID,2)
         READ(2,*) OUTYP
!---WAVE GAUGE NUMBER
         CALL HEADLINE(ID,2)
         READ(2,*) NWG
		 IF (NWG>10)THEN
		 WRITE(*,*) "NEED MORE ALLOCATION FOR WAVE GAUGE"
		 WRITE(22,*) "NEED MORE ALLOCATION FOR WAVE GAUGE"
		 STOP
		 END IF
!---WAVE GAUGE LOCATION
         CALL HEADLINE(ID,2)
         READ(2,*) WGX(1:NWG)

      RETURN
    END
!**********************************************************************
      SUBROUTINE INPUT_3(SMTYP,NNRL,MDEG,NARC,ALETYP)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER SMTYP,NNRL,MDEG,NARC,ALETYP
      CHARACTER*2 ID
         ID = '*'
!---Do you need free-surface smoothing: 1 = yes; 0 = no
         CALL HEADLINE(ID,3)
         READ(3,*) SMTYP
!---number of neighboring node and degree of polynomial for SG filter
!---note that nl = nr and nl + nr < m
         CALL HEADLINE(ID,3)
         READ(3,*) NNRL,MDEG

!---number of gaussian quadrature for calculating arc length
         CALL HEADLINE(ID,3)
         READ(3,*) NARC
         
!---do you use ALE approach? 0 = NO; 1 = YES
         CALL HEADLINE(ID,3)
         READ(3,*) ALETYP
         
      RETURN
    END
!**********************************************************************
      SUBROUTINE LENGTH(NPL,COOR1,SIDE_L1)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
	  INTEGER NPL
      REAL*8  COOR1(NPL,2),SIDE_L1(NPL)
      DO I=1,NPL-1
        SIDE_L1(I)=DSQRT((COOR1(I+1,1)-COOR1(I,1))**2+(COOR1(I+1,2)-COOR1(I,2))**2)
      END DO
        SIDE_L1(NPL)=DSQRT((COOR1(NPL,1)-COOR1(1,1))**2+(COOR1(NPL,2)-COOR1(1,2))**2)

      RETURN
      END
!**********************************************************************
      SUBROUTINE MESH(NPL,NNODE,NELM,NELEM,LN,COOR,LENG,NODE)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NPL,NEL(NPL),NNODE,NELM,NELEM(NPL),LN(NELM,2)
      REAL*8 SX,SY,NORM,DELT,LENG(NPL),COOR(NPL,2),NODE(NNODE,2)

K=0
DO I=1,NPL-1
	J=NPL-I
    NEL(I) = NELEM(I)+1
    DELT=LENG(J)/NELEM(I)
	SX=COOR(J,1)-COOR(J+1,1)
	SY=COOR(J,2)-COOR(J+1,2)
	NORM=DSQRT(SX**2+SY**2)
	SX=SX/NORM
	SY=SY/NORM
    DO L=1,NELEM(I)+1
       NODE(L+K,1)=COOR(J+1,1)+(L-1)*DELT*SX
       NODE(L+K,2)=COOR(J+1,2)+(L-1)*DELT*SY
    END DO
    K=K+NEL(I)
END DO

    NEL(NPL) = NELEM(NPL)+1
    DELT=LENG(NPL)/NELEM(NPL)
	SX=COOR(NPL,1)-COOR(1,1)
	SY=COOR(NPL,2)-COOR(1,2)
	NORM=DSQRT(SX**2+SY**2)
	SX=SX/NORM
	SY=SY/NORM
    DO I=1,NELEM(NPL)+1
      NODE(I+K,1)=COOR(1,1)+(I-1)*DELT*SX
      NODE(I+K,2)=COOR(1,2)+(I-1)*DELT*SY
    END DO

!----TO GIVE THE LOCAL ELEMENT NODE NUMBER
      L=1
	  N=1
      DO I=1,NPL
       DO J=1,NELEM(I)
        LN(N,1)=L
        LN(N,2)=L+1
        L=L+1
		N=N+1
       END DO
       L=L+1
      END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE SHAP(SHA1,SHA2,SH,NGA,RT)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NGA
      REAL*8 RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
      DO M=1,NGA
        SHA1(M)=0.5D0*(1-RT(M))
        SHA2(M)=0.5D0*(1+RT(M))

        SH(1,M)=SHA1(M)
        SH(2,M)=SHA2(M)
      END DO
	RETURN
	END 
!**********************************************************************
      SUBROUTINE REMESH(NPL,NNODE,NELEM,NODE,NS,TIME,DEP,AMP,NWG,WGX)
!**********************************************************************
      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NPL,NNODE,NELEM(NPL),NS(NPL),NWG,IL,IR
      REAL*8 TIME,DEP,AMP,LENG,NODE(NNODE,2),SX,SY,NORM,WGX(10),WGY(10)
      
!------ENSURE DUPLICATE POINT ON END NODE OF FREE SURFACE-----
    NODE(NNODE,:)=NODE(1,:)
	NODE(NS(1)+1,:)=NODE(NS(1),:)

!------BOTTOM END NODE GOES WITH FREE SURFACE -----
	NODE(NS(NPL-1),1)=NODE(1,1)
	NODE(NS(2),1)=NODE(NS(1),1)

!----- REMESH FOR ALL PLANES EXCEPT THE FREE SURFACE
DO I=2,NPL
	SX=NODE(NS(I),1)-NODE(NS(I-1),1)
	SY=NODE(NS(I),2)-NODE(NS(I-1),2)
	NORM=DSQRT(SX**2+SY**2)
	SX=SX/NORM
	SY=SY/NORM
	DELT=NORM/NELEM(I)
	  K=1
      DO L=NS(I-1)+1,NS(I)
        NODE(L,1)=NODE(NS(I-1),1)+DELT*(K-1)*SX
        NODE(L,2)=NODE(NS(I-1),2)+DELT*(K-1)*SY
		K=K+1
      END DO
END DO

!==========OUTPUT WAVE ELEVATION AT THE WAVE GAUGES (IN CM)==========
DO I=1,NWG
CALL BWLOC(-WGX(I),NS(1),-NODE(1:NS(1),1),0,IR,IL)
TEMP=(WGX(I)-NODE(IL,1))/(NODE(IR,1)-NODE(IL,1))
WGY(I)=NODE(IL,2)+TEMP*(NODE(IR,2)-NODE(IL,2))-DEP
END DO
WRITE(7,"(11(E15.8,1X))") TIME,WGY(1:NWG)*100.D0

      RETURN
    END
!********************************************************************
      SUBROUTINE BOUND(NPL,NNODE,NELEM,NS,BUNTYP,OUTYP,PHI,PPHI,PHIT,VEL,WC)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8 (A-H,O-Z)
       INTEGER NPL,NNODE,NELEM(NPL),NS(NPL),BUNTYP(NPL),OUTYP
	   REAL*8 R,VEL,WC
       REAL*8 PHI(NNODE),PPHI(NNODE),PHIT(NNODE)

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
			   PHI(J)=PHI(J)
			   PPHI(J)=0.D0
            ELSE
			  IF (I==2)THEN
                  IF(OUTYP==0)THEN
                      PPHI(J)=0.D0
                      PHI(J)=0.D0
                  ELSE
                      PPHI(J)=-PHIT(J)/WC
                      PHI(J)=0.D0
                  END IF
			  ELSE IF (I==NPL)THEN
			   PPHI(J)=-VEL
			   PHI(J)=0.D0
			  ELSE
			   PPHI(J)=0.D0
			   PHI(J)=0.D0
			  END IF
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT,EYE)
!********************************************************************
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,NNODE,NELM,NGA,LN(NELM,2)
      REAL*8 RD,SIGMA(NNODE),EYE(NNODE)
      REAL*8  KER1(NNODE,NNODE),KER2(NNODE,NNODE)
      REAL*8  NX,NY,H(2),G(2),XFUNC(10),YFUNC(10),PXI1(2)
      REAL*8  LENG(NELM),NORM(NELM,2),JCB(NELM),NODE(NNODE,2)
      REAL*8  WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)

        KER1=0.D0
        KER2=0.D0
!**** CALCULATE THE JACOBIAN
      DO J=1,NELM
      LENG(J)=DSQRT((NODE(LN(J,1),1)-NODE(LN(J,2),1))**2+(NODE(LN(J,1),2)-NODE(LN(J,2),2))**2)
		DO L=1,2
          PXI1(L)=(-0.5D0)*NODE(LN(J,1),L)+ 0.5D0*NODE(LN(J,2),L)
        END DO
        NX=PXI1(2)
        NY=PXI1(1)
        JCB(J)=DSQRT(NX**2+NY**2)
        NORM(J,1)=-NX/JCB(J)
        NORM(J,2)=NY/JCB(J)
      END DO
      
!$omp parallel do private(I,J,K,M,XFUNC,YFUNC,RD,G,H)
!***THE SURFACE KERNELS
      DO I = 1,NNODE
       DO J=1,NELM
       DO M=1,NGA
          XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
          YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
       END DO
        DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
         RD=DSQRT((NODE(I,1)-NODE(LN(J,K),1))**2+(NODE(I,2)-NODE(LN(J,K),2))**2)
        IF (RD .LE. 0.0000001D0) THEN
        H(K)=0.D0
	    G(K)=LENG(J)/2*(1.5D0-DLOG(LENG(J)))
		ELSE !---NON DSINGULER TERM
         G(K)=0.D0
         DO M=1,NGA
            H(K)=H(K)+(-1.D0)/((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2)*&
                    &((XFUNC(M)-NODE(I,1))*NORM(J,1)+(YFUNC(M)-NODE(I,2))*NORM(J,2))&
                    &*JCB(J)*SH(K,M)*WT(M)
            G(K)=G(K)+DLOG(1.D0/((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2)**0.5D0)*JCB(J)*SH(K,M)*WT(M)
         END DO
      END IF
         KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
      END DO
      END DO
      END DO
!$omp end parallel do
      
!***DSINGULAR OF KER1
    DO I=1,NNODE
    KER1(I,I)=0.D0
    END DO      
      CALL DGEMM('N','N',NNODE,1,NNODE,1.D0,KER1,NNODE,EYE,NNODE,0.D0,SIGMA,NNODE)
    DO I=1,NNODE
    KER1(I,I)=-SIGMA(I)
    END DO  
!       DO I=1,NNODE
!          SIGMA=0.D0
!          DO J=1,NNODE
!             SIGMA=KER1(I,J)+SIGMA
!          END DO
!          KER1(I,I)=-SIGMA
!       END DO

      RETURN
      END
!**********************************************************************
       SUBROUTINE SOLVE_BACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!**********************************************************************
!      TO SOLVE KER1*PHI=KER2*PPHI
!      PPHI=PARTIAL PHI OVER PARTIAL N
!      USING GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION
!======================================================================
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8    (A-H,O-Z)
       INTEGER NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
       REAL*8    PHI(NNODE),PPHI(NNODE)
       REAL*8    KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8    H1(NNODE,NNODE),Q1(NNODE),TEMP(NNODE)
       REAL*8  SUM,A,SIG,G1(NNODE,NNODE),P1(NNODE)

!********** TO MOVE THE KER1 AND KER2**********************************
!-----PHI PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO
!-------------------
         DO I=1,NNODE
           N=0
           DO L=1,NPL
             DO J=K+N,(NELEM(L)+1)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+(NELEM(L)+1)
          END DO
        END DO

       DO I=1,NNODE
          TEMP(I)=0.D0
          DO J=1,NNODE
          TEMP(I)=TEMP(I)+H1(I,J)*Q1(J)
          END DO
       END DO
!*************GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION*********
      DO I=1,NNODE
        G1(I,NNODE+1)=TEMP(I)
      END DO
      DO I=1,NNODE
         SUM=0.D0
         DO K=I,NNODE
            IF (G1(K,I) .NE. 0) THEN
               IF (K .NE. I) THEN
               IF (G1(I,I) .EQ. 0.D0) THEN
                 WRITE(*,*) 'SOME OF THE DIAG-TERMS ARE ZERO'
                 WRITE(22,*) 'SOME OF THE DIAG-TERMS ARE ZERO'
                 STOP
               END IF
               A=G1(K,I)/G1(I,I)
               DO J=I,NNODE+1
                  G1(K,J)=G1(K,J)-A*G1(I,J)
               END DO
               END IF
            END IF
            SUM=SUM+G1(K,I)
          END DO
          IF (SUM .EQ. 0.D0) THEN
          WRITE(*,*) 'NO UNIQUE SOLUTION EXISTS  STOP AT GAUSSELI'
          WRITE(22,*) 'NO UNIQUE SOLUTION EXISTS  STOP AT GAUSSELI'
          STOP
          END IF
      END DO

      P1(NNODE)=G1(NNODE,NNODE+1)/G1(NNODE,NNODE)
      DO I=NNODE-1,1,-1
         SIG=0.D0
         DO J=I+1,NNODE
            SIG=G1(I,J)*P1(J)+SIG
          END DO
         P1(I)=(G1(I,NNODE+1)-SIG)/G1(I,I)
      END DO
!=================================
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO

      RETURN
      END
!**********************************************************************
       SUBROUTINE SOLVE_LAPACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!**********************************************************************
!      TO SOLVE KER1*PHI=KER2*PPHI
!      PPHI=PARTIAL PHI OVER PARTIAL N
!======================================================================
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8    (A-H,O-Z)
       INTEGER NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8    PHI(NNODE),PPHI(NNODE)
       REAL*8    KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8    H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)

!********** TO MOVE THE KER1 AND KER2**********************************
!-----PHI PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO
!-------------------
         DO I=1,NNODE
           N=0
           DO L=1,NPL
             DO J=K+N,(NELEM(L)+1)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+(NELEM(L)+1)
          END DO
        END DO

       DO I=1,NNODE
          P1(I)=0.D0
          DO J=1,NNODE
          P1(I)=P1(I)+H1(I,J)*Q1(J)
          END DO
       END DO

!*************SOLVE BY CALLING LAPACK*********
CALL DGESV(NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)

!=================================
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO

      RETURN
    END
!**********************************************************************
       SUBROUTINE SOLVE_LAPACK2_1(NPL,PHI,PPHI,KER1,KER2,G1,H1,NNODE,NELEM,BUNTYP,IPIV)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE),KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)
	   CHARACTER*1 TRANS
       TRANS = 'N'
       
!-----MOVE PHI AND PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO
       
!-----MOVE KER1 AND KER2---- 
         DO I=1,NNODE
           N=0
           DO L=1,NPL
             DO J=K+N,(NELEM(L)+1)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+(NELEM(L)+1)
          END DO
        END DO

       P1=MATMUL(H1,Q1)
       
!*************SOLVE BY CALLING LAPACK*********
CALL DGETRF(NNODE,NNODE,G1,NNODE,IPIV,INFO)
CALL DGETRS(TRANS,NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)
       
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO
                     
      RETURN
    END
!**********************************************************************
       SUBROUTINE SOLVE_LAPACK2_2(NPL,PHI,PPHI,G1,H1,NNODE,NELEM,BUNTYP,IPIV)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)
	   CHARACTER*1 TRANS
       TRANS = 'N'
       
!-----MOVE PHI AND PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO

       P1=MATMUL(H1,Q1)

!*************SOLVE BY CALLING LAPACK*********
CALL DGETRS(TRANS,NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO

      RETURN
      END
!********************************************************************
SUBROUTINE TAYES1(NPL,NNODE,NELM,NELEM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,&
	  &DPDS,DP,PHIT,DPDT,DEP,CLEN,GRAV,MU,ARDZONE,VEL,WC)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NPL,NNODE,NELM,NELEM(NPL),ME(NPL),NS(NPL),LN(NELM,2)
    REAL*8 WC,DEP,CLEN,GRAV,MU,ARDZONE,VEL
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),DP(NNODE,2)
	REAL*8 DPDS(NNODE),PHIT(NNODE),DPDT(NNODE)

	DPDS=0.D0
!*********************ON FREE SURFACE*********************
    DO I=1,ME(1)
    DO J=1,2
		IF(LN(I,J).EQ.1) THEN
            DP(LN(I,J),1)=-PPHI(NNODE)
            DPDS(LN(I,J))=(DP(LN(I,J),1)-PPHI(LN(I,J))*NORM(I,1))/NORM(I,2)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(1)) THEN
		    DP(LN(I,J),1)=PPHI(NS(1)+1)
            DPDS(LN(I,J))=(DP(LN(I,J),1)-PPHI(LN(I,J))*NORM(I,1))/NORM(I,2)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
        ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO
        
    DO I=1,NS(1)
        IF (NODE(I,1)>=ARDZONE)THEN
          DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)-MU*(NODE(I,1)-ARDZONE)/(CLEN-ARDZONE)*PHI(I)
        ELSE
          DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)
        END IF
        PHIT(I)=DPDT(I)-(DP(I,1)**2+DP(I,2)**2)
    END DO

!*********************ON RIGHT WALL*********************
    DO I=ME(1)+1,ME(2)
    DO J=1,2
		IF(LN(I,J).EQ.NS(1)+1) THEN
			DPDS(LN(I,J))=-DP(NS(1),2)
            DP(LN(I,J),1)=PPHI(LN(I,J))
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(2)) THEN
			DPDS(LN(I,J))=0.D0
            DP(LN(I,J),1)=PPHI(LN(I,J))
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=PPHI(LN(I,J))
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

!*********************LEFT WALL*********************
    DO I=ME(NPL-1)+1,ME(NPL)
    DO J=1,2
		IF(LN(I,J).EQ.NS(NPL-1)+1) THEN
			DPDS(LN(I,J))=0.D0
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NNODE) THEN
			DPDS(LN(I,J))=DP(1,2)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

!*********************BOTTOM*********************
    DO I=ME(2)+1,ME(NPL-1)
    DO J=1,2
		IF(LN(I,J).EQ.NS(2)+1) THEN
          DPDS(LN(I,J))=-PPHI(NS(2))
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
          DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(NPL-1)) THEN
          DPDS(LN(I,J))=PPHI(NS(NPL-1)+1)
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
          DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
		  DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
		  DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

      RETURN
    END
!********************************************************************
SUBROUTINE TAYES1_CS(NPL,NNODE,NELM,NELEM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,&
	  &DPDS,DPDSS,DP,PHIT,DPDT,DEP,CLEN,GRAV,MU,ARDZONE,VEL,WC,AWT,ART,NARC)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NARC,NPL,NNODE,NELM,NELEM(NPL),ME(NPL),NS(NPL),LN(NELM,2)
    REAL*8 WC,DEP,CLEN,GRAV,MU,ARDZONE,VEL
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),DP(NNODE,2)
	REAL*8 DPDS(NNODE),PHIT(NNODE),DPDT(NNODE),AWT(NARC),ART(NARC),DPDSS(NS(1))

	DPDS=0.D0
!*********************ON FREE SURFACE*********************
    CALL FS_CSDIFF(NS(1),NODE(1:NS(1),1),NODE(1:NS(1),2),PHI(1:NS(1)),AWT,ART,NARC,DPDS(1:NS(1)),DPDSS)
    
    DO I=1,ME(1)
    DO J=1,2
		IF(LN(I,J).EQ.1) THEN
            DP(LN(I,J),1)=-PPHI(NNODE)
            DPDS(LN(I,J))=(DP(LN(I,J),1)-PPHI(LN(I,J))*NORM(I,1))/NORM(I,2)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(1)) THEN
		    DP(LN(I,J),1)=PPHI(NS(1)+1)
            DPDS(LN(I,J))=(DP(LN(I,J),1)-PPHI(LN(I,J))*NORM(I,1))/NORM(I,2)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
        ELSE
			!DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I) already got from cubic spline interpolation
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO
        
    DO I=1,NS(1)
        IF (NODE(I,1)>=ARDZONE)THEN
          DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)-MU*(NODE(I,1)-ARDZONE)/(CLEN-ARDZONE)*PHI(I)
        ELSE
          DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)
        END IF
        PHIT(I)=DPDT(I)-(DP(I,1)**2+DP(I,2)**2)
    END DO

!*********************ON RIGHT WALL*********************
    DO I=ME(1)+1,ME(2)
    DO J=1,2
		IF(LN(I,J).EQ.NS(1)+1) THEN
			DPDS(LN(I,J))=-DP(NS(1),2)
            DP(LN(I,J),1)=PPHI(LN(I,J))
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(2)) THEN
			DPDS(LN(I,J))=0.D0
            DP(LN(I,J),1)=PPHI(LN(I,J))
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=PPHI(LN(I,J))
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

!*********************LEFT WALL*********************
    DO I=ME(NPL-1)+1,ME(NPL)
    DO J=1,2
		IF(LN(I,J).EQ.NS(NPL-1)+1) THEN
			DPDS(LN(I,J))=0.D0
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NNODE) THEN
			DPDS(LN(I,J))=DP(1,2)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

!*********************BOTTOM*********************
    DO I=ME(2)+1,ME(NPL-1)
    DO J=1,2
		IF(LN(I,J).EQ.NS(2)+1) THEN
          DPDS(LN(I,J))=-PPHI(NS(2))
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
          DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(NPL-1)) THEN
          DPDS(LN(I,J))=PPHI(NS(NPL-1)+1)
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
          DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
		  DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
		  DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

      RETURN
    END
!**********************************************************************
      SUBROUTINE FS_CSDIFF(N,X,Y,PHI,AWT,ART,NARC,DPDS,DPDSS)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,N,NE,NARC
      REAL*8 X(N),Y(N),PHI(N),DPDS(N),DPDSS(N),AWT(NARC),ART(NARC),S(N),DARC,U
      REAL*8 B(N),C(N),D(N),Y0,Y1,Y2

!---CALCULATE THE ARC LENGTH AND TANGENTIAL COORDINATE
      CALL SPLINE(N,X,Y,B,C,D)
      S(1)=0.D0
      DO I=2,N
          DARC=0.D0
          DO J=1,NARC
              U=0.5D0*(1-ART(J))*X(I-1)+0.5D0*(1+ART(J))*X(I)
              CALL SEVAL(N,U,X,Y,B,C,D,Y0,Y1,Y2)
              DARC=DARC+DSQRT(1.D0+Y1**2)*0.5D0*(X(I)-X(I-1))*AWT(J)
          END DO
          S(I)=S(I-1)+DARC      
      END DO

!---calculate dPHI/dS
      CALL SPLINE(N,S,PHI,B,C,D)
      DO I=1,N
          CALL SEVAL(N,S(I),S,PHI,B,C,D,Y0,DPDS(I),DPDSS(I))
      END DO

      RETURN
    END
!********************************************************************
      SUBROUTINE ACCBC(NPL,NNODE,NELM,NELEM,ME,LN,NORM,JCB,PPHI,DP,ACCMO)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NPL,NNODE,NELM,NELEM(NPL),ME(NPL),LN(NELM,2)
    REAL*8 DPNDX,DPNDY,DPNDS
	REAL*8 NORM(NELM,2),JCB(NELM),PPHI(NNODE),DP(NNODE,2),ACCMO(NNODE)

DO K=1,NPL-1
  DO I=ME(K)+1,ME(K+1)
    DO J=1,2
		DPNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
		DPNDX=DPNDS*NORM(I,2)	! PHINN=PHISS=0 FOR LINEAR ELEMENT
		DPNDY=-DPNDS*NORM(I,1)
		ACCMO(LN(I,J))=DP(LN(I,J),1)*DPNDX+DP(LN(I,J),2)*DPNDY
    END DO
  END DO
END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE BOUNDT(NPL,NNODE,NELM,NS,NELEM,BUNTYP,OUTYP,LN,PHIT,PPHIT,DPDS,JCB,WC,ACC,ACCMO)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8    (A-H,O-Z)
       INTEGER NNODE,NELM,NS(NPL),NELEM(NPL),BUNTYP(NPL),LN(NELM,2),OUTYP
	   REAL*8 R,ACC,WC
       REAL*8 JCB(NELM),PHIT(NNODE),PPHIT(NNODE),DPDS(NNODE),ACCMO(NNODE),DPDSS(NNODE)

	PPHIT=0.D0

    DO I=NELEM(1)+1,NELEM(1)+NELEM(2)
    DO J=1,2
		IF(LN(I,J).EQ.NS(2)) THEN
		  DPDSS(LN(I,J))=0.D0 ! ASSUME A FLAT GROUND SO VN=0
		ELSE
		  DPDSS(LN(I,J))=0.5D0*(-DPDS(LN(I,1))+DPDS(LN(I,2)))/JCB(I)
		END IF
    END DO
    END DO

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
			 PHIT(J)=PHIT(J)
			 PPHIT(J)=0.D0
            ELSE
			  IF (I==2)THEN
                  IF (OUTYP==0)THEN
                      PPHIT(J)=0.D0
                      PHIT(J)=0.D0 
                  ELSE
                      PPHIT(J)=WC*DPDSS(J)
                      PHIT(J)=0.D0
                  END IF
			  ELSE IF (I==NPL)THEN
			    PPHIT(J)=-ACC-ACCMO(J)
			    PHIT(J)=0.D0
			  ELSE
			    PPHIT(J)=0.D0-ACCMO(J)
			    PHIT(J)=0.D0
			  END IF
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO

      RETURN
      END
!**********************************************************************
SUBROUTINE TAYES2(PHI,PPHI,PHIT,PPHIT,DPDS,DPDT,DP,NPL,NNODE,NELEM,NODE,NELM,NORM,&
				 &JCB,LN,DELTTIME,GRAV,ACC)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NPL,NNODE,NELM,NELEM(NPL),LN(NELM,2)
      REAL*8 DELTTIME,GRAV,ACC,D2PDT,DPTDS,DPDNDS
      REAL*8 NORM(NELM,2),JCB(NELM)
      REAL*8 NODE(NNODE,2),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE)
	  REAL*8 DP(NNODE,2),DPDS(NNODE),DPDT(NNODE)
      REAL*8 NEW_NODE(NELEM(1)+1,2),NEW_PHI(NELEM(1)+1),D2P(2)

	DO I=1,NELEM(1)
	  DO J=1,2
		IF(LN(I,J) .EQ. 1) THEN
				D2P(1)=-PPHIT(NNODE)
				DPDNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
				DPTDS=(D2P(1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
				D2P(2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		ELSE IF (LN(I,J) .EQ. NELEM(1)+1) THEN
				D2P(1)=PPHIT(NELEM(1)+2)
				DPDNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
				DPTDS=(D2P(1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
				D2P(2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		ELSE
			DPTDS=(-0.5D0*PHIT(LN(I,1))+0.5D0*PHIT(LN(I,2)))/JCB(I)
			DPDNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
			D2P(1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
			D2P(2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		END IF

		D2PDT=DP(LN(I,J),1)*D2P(1)+DP(LN(I,J),2)*D2P(2)-GRAV*DP(LN(I,J),2)
		NEW_PHI(LN(I,J))=PHI(LN(I,J))+DELTTIME*DPDT(LN(I,J))+D2PDT*DELTTIME**2/2
		NEW_NODE(LN(I,J),1)=NODE(LN(I,J),1)+DP(LN(I,J),1)*DELTTIME+0.5D0*D2P(1)*DELTTIME**2
		NEW_NODE(LN(I,J),2)=NODE(LN(I,J),2)+DP(LN(I,J),2)*DELTTIME+0.5D0*D2P(2)*DELTTIME**2
	  END DO
	END DO

	DO I=1,NELEM(1)+1
	NODE(I,:)=NEW_NODE(I,:)
	PHI(I)=NEW_PHI(I)
	END DO

      RETURN
    END
!**********************************************************************
      SUBROUTINE ALE_VEL(NPL,NNODE,NELM,NELEM,NS,LN,NODE,NORM,DP,DPDT,DPM,DPDTM,CV)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,NPL,NNODE,NELM,NELEM(NPL),NS(NPL),LN(NELM,2)
      REAL*8 NODE(NNODE,2),NORM(NELM,2),DP(NNODE,2),DPDT(NNODE),DPM(NS(1),2),DPDTM(NS(1)),CV(NS(1),2)
      
!---CALCULATE MESH VELOCITY      
      DO I=1,NELEM(1)
          DO J=1,2
              DPM(LN(I,J),1)=0.D0
              DPM(LN(I,J),2)=DP(LN(I,J),2)+(DP(LN(I,J),1)-DPM(LN(I,J),1))*NORM(I,1)/NORM(I,2)
          END DO
      END DO
      
!---CALCULATE CONVECTIVE VELOCITY AND DPDT IN ALE FRAME
      DO I=1,NS(1)
          CV(I,1)=DP(I,1)-DPM(I,1)
          CV(I,2)=DP(I,2)-DPM(I,2)
          DPDTM(I)=DPDT(I)-CV(I,1)*DP(I,1)-CV(I,2)*DP(I,2)
      END DO
      
      RETURN
    END
!**********************************************************************
SUBROUTINE ALE_ACC(DELTTIME,GRAV,MU,CLEN,ARDZONE,NPL,NNODE,NELEM,NS,NELM,LN,NODE,NORM,JCB,&
                  &PHI,PPHI,PHIT,PPHIT,DP,DPDS,DPDSS,DPDT,CV,D2PM,D2PDTM)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,NPL,NNODE,NELM,NELEM(NPL),NS(NPL),LN(NELM,2)
      REAL*8 DELTTIME,GRAV,MU,CLEN,ARDZONE,DPTDS,DPDNDS,TEMP
      REAL*8 NORM(NELM,2),JCB(NELM)
      REAL*8 NODE(NNODE,2),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE)
	  REAL*8 DP(NNODE,2),DPDS(NNODE),DPDT(NNODE),D2P(NNODE,2),D2PDT(NNODE)
      REAL*8 CV(NS(1),2),D2PM(NS(1),2),D2PDTM(NS(1)),DPDSS(NS(1))
      REAL*8 UXX(NS(1)),UXY(NS(1)),UYY(NS(1)),DPDTDX(NS(1),2)

!---CALCULATE THE SECOND ORDER TERM OF TSE IN EL FRAME (DPDSS=DPDNN=0 FOR LINEAR ELEMENT)
	DO I=1,NELEM(1)
	  DO J=1,2
		IF(LN(I,J) .EQ. 1) THEN
            D2P(LN(I,J),1)=-PPHIT(NNODE)
            DPDNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
            DPTDS=(D2P(LN(I,J),1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
            D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
        ELSE IF (LN(I,J) .EQ. NELEM(1)+1) THEN
            D2P(LN(I,J),1)=PPHIT(NELEM(1)+2)
            DPDNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
            DPTDS=(D2P(LN(I,J),1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
            D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		ELSE
			DPTDS=(-0.5D0*PHIT(LN(I,1))+0.5D0*PHIT(LN(I,2)))/JCB(I)
			DPDNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
			D2P(LN(I,J),1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
			D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
        END IF
        UXX(LN(I,J))=2.D0*DPDNDS*NORM(I,1)*NORM(I,2)+DPDSS(LN(I,J))*(NORM(I,2)**2-NORM(I,1)**2) !2.D0*DPDNDS*NORM(I,1)*NORM(I,2)
        UXY(LN(I,J))=DPDNDS*(NORM(I,2)**2-NORM(I,1)*2)-2.D0*DPDSS(LN(I,J))*NORM(I,2)*NORM(I,1) !DPDNDS*(NORM(I,2)*NORM(I,2)-NORM(I,1)*NORM(I,1))
        UYY(LN(I,J))=-2.D0*DPDNDS*NORM(I,1)*NORM(I,2)-DPDSS(LN(I,J))*(NORM(I,2)**2-NORM(I,1)**2) !-2.D0*DPDNDS*NORM(I,1)*NORM(I,2)
		D2PDT(LN(I,J))=DP(LN(I,J),1)*D2P(LN(I,J),1)+DP(LN(I,J),2)*D2P(LN(I,J),2)-GRAV*DP(LN(I,J),2)
        DPDTDX(LN(I,J),1)=DP(LN(I,J),1)*UXX(LN(I,J))+DP(LN(I,J),2)*UXY(LN(I,J))+GRAV*NORM(I,1)/NORM(I,2)
        DPDTDX(LN(I,J),2)=DP(LN(I,J),1)*UXY(LN(I,J))+DP(LN(I,J),2)*UYY(LN(I,J))
	  END DO
    END DO
    
    DO I=1,NS(1)
        IF (NODE(I,1)>=ARDZONE)THEN
            TEMP=MU*(NODE(I,1)-ARDZONE)/(CLEN-ARDZONE)
            D2PDT(I)=D2PDT(I)-TEMP*DPDT(I)
            DPDTDX(I,1)=DPDTDX(I,1)-TEMP*DP(I,1)
            DPDTDX(I,2)=DPDTDX(I,2)-TEMP*DP(I,2)
        END IF
    END DO
    
!---calculate the mesh acceleration
      DO I=1,NELEM(1)+1
          D2PM(I,1)=D2P(I,1)-CV(I,1)*UXX(I)-CV(I,2)*UXY(I)
          D2PM(I,2)=D2P(I,2)-CV(I,1)*UXY(I)-CV(I,2)*UYY(I)
          D2PDTM(I)=D2PDT(I)-CV(I,1)*DPDTDX(I,1)-CV(I,2)*DPDTDX(I,2)
      END DO
    
      RETURN
    END
!**********************************************************************
SUBROUTINE ALE_TSE_SUM(DELTTIME,NPL,NNODE,NS,NODE,PHI,DPM,DPDTM,D2PM,D2PDTM)
!**********************************************************************
      IMPLICIT NONE 
      INTEGER I,NPL,NNODE,NS(NPL)
      REAL*8 DELTTIME
      REAL*8 NODE(NNODE,2),PHI(NNODE),NEW_NODE(NS(1),2),NEW_PHI(NS(1))
	  REAL*8 DPM(NS(1),2),DPDTM(NS(1)),D2PM(NS(1),2),D2PDTM(NS(1))

	DO I=1,NS(1)
		NEW_PHI(I)=PHI(I)+DELTTIME*DPDTM(I)+0.5D0*DELTTIME**2*D2PDTM(I)
		NEW_NODE(I,1)=NODE(I,1)+DELTTIME*DPM(I,1)+0.5D0*DELTTIME**2*D2PM(I,1)
		NEW_NODE(I,2)=NODE(I,2)+DELTTIME*DPM(I,2)+0.5D0*DELTTIME**2*D2PM(I,2)
	END DO

	DO I=1,NS(1)
	    NODE(I,:)=NEW_NODE(I,:)
	    PHI(I)=NEW_PHI(I)
	END DO

      RETURN
    END
!**********************************************************************
      SUBROUTINE FS_SMOOTH(NNRL,MDEG,N,Y)
!**********************************************************************
      IMPLICIT NONE
      INTEGER N,NNRL,MDEG,FLAG
      REAL*8 Y(N)
      
      CALL savgol_filter(NNRL,NNRL,0,MDEG,N,Y,flag)
      
      RETURN
    END
!********************************************************************
SUBROUTINE PRESSURE(ICON,TIME,THO,GRAV,DEP,NPL,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,ICON,NPL,NNODE,NS(NPL)
      REAL*8 TIME,DEP,THO,GRAV,P1,P2,P3,P_ATM
      REAL*8 NODE(NNODE,2),PHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL*8 CP1(NS(1)),CP2(NS(1)),CP3(NS(1)),CP(NS(1))

!----ATOM PRESSURE (BERNOULLI CONSTANT) ON THE FREE SURFACE
	DO I=ICON,ICON ! NS(1) !
	CP1(I)=THO*PHIT(I)
	CP2(I)=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	CP3(I)=THO*GRAV*(NODE(I,2)-DEP)
	CP(I)=CP1(I)+CP2(I)+CP3(I)
	ENDDO
	P_ATM=CP(ICON)

!----PRESSURE ON ZONE 1 BOUNDARY
	DO I=NS(1)+2,NNODE-1
	P1=THO*PHIT(I)
	P2=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	P3=THO*GRAV*(NODE(I,2)-DEP)
	PR(I)=P_ATM-(P1+P2+P3)
	END DO

      RETURN
      END
!********************************************************************
SUBROUTINE DOMAIN(NPL,NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,&
				 &SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,L,M,IL,IR,NPL,NGA,NFIELD,NNODE,NELM,NELEM(NPL),NS(NPL),LN(NELM,2)
	  REAL*8 THO,GRAV,DEP,P_ATM
	  REAL*8 HB,DX,DY,PI2,TEMP,P1,P2,P3
	  REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM)
	  REAL*8 PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL*8 DNODE(NFIELD,2),DVX(NFIELD),DVY(NFIELD),DPHIT(NFIELD),DPR(NFIELD)
	  REAL*8 KER1(NFIELD,NNODE),KER2(NFIELD,NNODE)
      REAL*8 H(2),G(2),XFUNC(10),YFUNC(10),PXI1(2)
	  REAL*8 WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
	
	PI2=2.D0*DACOS(-1.D0)

!----CREATE DOMAIN POINT
	L=1
	DO I=2,NS(1)-1
	  CALL BWLOC(NODE(I,1),NS(NPL-1)-NS(2),NODE(NS(2)+1:NS(NPL-1),1),NS(2),IL,IR)
	  HB=NODE(IL,2)+(NODE(I,1)-NODE(IL,1))/(NODE(IR,1)-NODE(IL,1))*(NODE(IR,2)-NODE(IL,2))+0.01D0 ! keep it a little far away from the boundary
	  DY=-NODE(I,2)/NELEM(NPL)
		DO J=2,NELEM(NPL)
		  TEMP=NODE(I,2)+DY*(J-1)
			IF(TEMP>HB)THEN
			DNODE(L,1)=NODE(I,1)
			DNODE(L,2)=TEMP !NODE(I,2)+DY*(J-1)
			L=L+1
			END IF
		END DO
	END DO

!--- SET A DUMMY NODE TO USELESS DNODE
	DO I=L,NFIELD
	DNODE(I,:)=DNODE(1,:)
	END DO

	KER1=0.D0
	KER2=0.D0
	DVX=0.D0
!----CALCULATE X VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((YFUNC(M)-DNODE(I,2))**2-(XFUNC(M)-DNODE(I,1))**2)*NORM(J,1)-&
					&2.D0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,2))/&
					&((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(XFUNC(M)-DNODE(I,1))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVX=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2

	KER1=0.D0
	KER2=0.D0
	DVY=0.D0
!----CALCULATE Y VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((XFUNC(M)-DNODE(I,1))**2-(YFUNC(M)-DNODE(I,2))**2)*NORM(J,2)-&
					&2.D0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,1))/&
					&((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(YFUNC(M)-DNODE(I,2))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVY=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2
	
	KER1=0.D0
	KER2=0.D0
	DPHIT=0.D0
!----CALCULATE PARTIAL POTENTIAL OVER TIME BY BIE
      DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(-1.D0)/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*&
					&((XFUNC(M)-DNODE(I,1))*NORM(J,1)+(YFUNC(M)-DNODE(I,2))*NORM(J,2))*TEMP
            G(K)=G(K)+DLOG(1.D0/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**0.5D0)*TEMP	 
		 END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DPHIT=(MATMUL(KER2,PPHIT)-MATMUL(KER1,PHIT))/PI2

!----CALCULATE PRESSURE DISTRIBUTION IN DOMAIN
	DO I=1,NFIELD
	P1=THO*DPHIT(I)
	P2=THO*0.5D0*(DVX(I)**2+DVY(I)**2)
	P3=THO*GRAV*(DNODE(I,2)-DEP)
	DPR(I)=P_ATM-(P1+P2+P3)
	END DO

	WRITE(11,'(3000(1X,F15.7))') NODE(:,1),DNODE(:,1)
	WRITE(11,'(3000(1X,F15.7))') NODE(:,2),DNODE(:,2)
	WRITE(11,'(3000(1X,F15.7))') DP(:,1),DVX
	WRITE(11,'(3000(1X,F15.7))') DP(:,2),DVY
	WRITE(11,'(3000(1X,F15.7))') PR,DPR

      RETURN
      END
!********************************************************************
      SUBROUTINE GAUSS(WT,RT,NGA)
!********************************************************************
      INTEGER NGA
      REAL*8   WT(NGA),RT(NGA)

      SELECT CASE(NGA)
       CASE(3)
        WT(1)=0.55555555
        WT(2)=0.88888889
        WT(3)=0.55555555
        RT(1)=0.77459667
        RT(2)=0.D0
        RT(3)=-0.77459667
       CASE(4)
        WT(1)=0.65214515
        WT(2)=0.34785484
        WT(3)=0.34785484
        WT(4)=0.65214515
        RT(1)=0.33998104
        RT(2)=0.86113631
        RT(3)=-0.86113631
        RT(4)=-0.33998104
       CASE(5)
        WT(1)=0.23692689
        WT(2)=0.47862867
        WT(3)=0.56888889
        WT(4)=0.47862867
        WT(5)=0.23692689
        RT(1)=0.90617985
        RT(2)=0.53846931
        RT(3)=0.D0
        RT(4)=-0.53846931
        RT(5)=-0.90617985
	 CASE(6)
	  WT(1)=0.17132449
	  WT(2)=0.36076157
	  WT(3)=0.46791393
	  WT(4)=0.46791393
	  WT(5)=0.36076157
	  WT(6)=0.17132449
	  RT(1)=0.93246951
	  RT(2)=0.66120938
	  RT(3)=0.23861918
	  RT(4)=-0.23861918
	  RT(5)=-0.66120938
	  RT(6)=-0.9346951
       CASE(8)
        WT(1)=0.1012285362903763D0
        WT(2)=0.2223810344533745D0
        WT(3)=0.3137066458778873D0
        WT(4)=0.3626837833783620D0
        WT(8)=0.1012285362903763D0
        WT(7)=0.2223810344533745D0
        WT(6)=0.3137066458778873D0
        WT(5)=0.3626837833783620D0
        RT(1)=0.9602898564975363D0
        RT(2)=0.7966664774136267D0
        RT(3)=0.5255324099163290D0
        RT(4)=0.1834346424956498D0
        RT(8)=-0.9602898564975363D0
        RT(7)=-0.7966664774136267D0
        RT(6)=-0.5255324099163290D0
        RT(5)=-0.1834346424956498D0
       CASE(10)
        WT(1)=0.D06667134
        WT(2)=0.14945134
        WT(3)=0.21908636
        WT(4)=0.26926671
        WT(5)=0.29552422
        WT(10)=0.D06667134
        WT(9)=0.14945134
        WT(8)=0.21908636
        WT(7)=0.26926671
        WT(6)=0.29552422
        RT(1)=0.97390652
        RT(2)=0.86506336
        RT(3)=0.67940956
        RT(4)=0.43339539
        RT(5)=0.14887433
        RT(10)=-0.97390652
        RT(9)=-0.86506336
        RT(8)=-0.67940956
        RT(7)=-0.43339539
        RT(6)=-0.14887433
      END SELECT

      RETURN
      END
!********************************************************************
	SUBROUTINE CONVERGE(N,P1,P2,E1LOC,E1,E2)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N,K(1),E1LOC
	REAL*8 P1(N),P2(N),E1,E2

	E1=MAXVAL(DABS(P1-P2))
	E1LOC=MAXLOC(DABS(P1-P2),1)

	E2=0.D0
	DO I=1,N
	E2=E2+(P1(I)-P2(I))**2
	END DO
	E2=DSQRT(E2/N)

	RETURN
	END
!********************************************************************
	SUBROUTINE BWLOC(PX,N,X,IST,IL,IR)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N,IST,IL,IR
	REAL*8 PX,X(N)

	DO I=1,N-1
		IF(X(I)>PX.AND.X(I+1)<=PX)THEN
		IR=IST+I
		IL=IST+I+1
		GOTO 777
		END IF
	END DO
777 CONTINUE

	RETURN
	END
!********************************************************************
	SUBROUTINE WAVE_SPD(GRAV,OMEGA,D,C)
!********************************************************************
	IMPLICIT NONE
	INTEGER I
	REAL*8 K,K2,GRAV,OMEGA,D,C,PI,F0,F1
	PI=DACOS(-1.D0)

	K=1.D0
	DO I=1,100
	F0=K*DTANH(K*D)-OMEGA**2/GRAV
	F1=DTANH(K*D)+K-K*(DTANH(K*D)**2)
	K2=K-(F0)/(F1)
		IF((K2-K)/K<=0.000001D0) THEN
		GOTO 717
		END IF
	K=K2
	END DO
	717 CONTINUE

	C=DSQRT(GRAV*DTANH(K*D)/K)

	RETURN
	END
!********************************************************************
SUBROUTINE COURANT(TIME,DELTTIME,NNODE,NELM,LN,NODE,DP,JCB)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,CFLOC,NNODE,NELM,LN(NELM,2)
    REAL*8 TIME,DELTTIME,U,V,VE
	REAL*8 NODE(NNODE,2),DP(NNODE,2),JCB(NELM),CN(NELM),CFL

  DO I=1,NELM
    U=DSQRT(DP(LN(I,1),1)**2+DP(LN(I,1),2)**2)
    V=DSQRT(DP(LN(I,2),1)**2+DP(LN(I,2),2)**2)
    VE=MAX(U,V)
    CN(I)=0.5D0*VE*DELTTIME/JCB(I)
  END DO
  CFL=MAXVAL(CN)
  CFLOC=MAXLOC(CN,1)
  
	WRITE(23,*) TIME,CFL

    IF (CFL>=250.D0)THEN
    WRITE(22,*) TIME,"CFL=",CFL,"@ ELEMENT",CFLOC
    STOP
    END IF    

	RETURN
	END
!**********************************************************************
SUBROUTINE LOBATTO(N,X,W,A,B)
!**********************************************************************
IMPLICIT NONE
INTEGER I,ITER,N
REAL A,B,C,D,X(N),W(N)
REAL X0,X1,E,P1,P2,P3
REAL,PARAMETER::PI=3.1415926
REAL,PARAMETER::EMAX=1.0E-5
REAL,EXTERNAL::LEGENDRE

X(1)=-1.0
W(1)=2.0/N/(N-1)
X(N)=1.0
W(N)=W(1)

DO I=2,N-1
ITER=1
E=100.0
X0=(1.0-3.0*(N-2)/8.0/(N-1)**3)*COS((4.0*I-3)*PI/(4.0*N-3))
    DO WHILE (E>=EMAX.AND.ITER<=1000)
    P1=(LEGENDRE(N-2,X0)-X0*LEGENDRE(N-1,X0))*(N-1)/(1.0-X0**2)
    P2=(2.0*X0*P1-N*(N-1)*LEGENDRE(N-1,X0))/(1.0-X0**2)
    P3=(2.0*X0*P2-(N*(N-1)-2)*P1)/(1.0-X0**2)
    X1=X0-2.0*P1*P2/(2.0*P2**2-P1*P3)
    E=ABS(X1-X0)
    ITER=ITER+1
    X0=X1
    END DO
X(N-I+1)=X1
END DO

DO I=2,N-1
W(I)=2.0/N/(N-1)/LEGENDRE(N-1,X(I))**2
END DO

C=(B-A)/2.0
D=(B+A)/2.0
DO I=1,N
X(I)=C*X(I)+D
W(I)=C*W(I)
END DO

RETURN
END
!**********************************************************************
FUNCTION LEGENDRE(N,X)
!**********************************************************************
REAL*8 FI,PI,PIM1,PIM2
REAL LEGENDRE
INTEGER I
IF (N.EQ.0) THEN
    LEGENDRE=1
    ELSEIF (N.EQ.1) THEN
    LEGENDRE=X
    ELSE
    PIM1=1
    PI=X
        DO I=2,N
        FI=I
        PIM2=PIM1
        PIM1=PI
        PI=((I+I-1)*X*PIM1-(I-1)*PIM2)/FI
        END DO
    LEGENDRE=PI
ENDIF
END
      
subroutine savgol_filter(nl,nr,ld,m,n1,y,flag)
!-----------------------------------------------------------------------------------
! This routine is used to perform the Savitzky-Golay algorithm.
!-----------------------------------------------------------------------------------
!    nl:: input, integer, the number of leftward data points used.
!    nr:: input, integer, the number of rightward data points used.
!    ld:: input, integer, the order of the derivative desired.
!     m:: input, integer, the order of the smoothing polynomial.
!    n1:: input, integer, the number of data points.
! y(n1):: input/output, real values, the data to be smoothed.
!  flag:: output, integer, error message, 0=success, 1=failure.
!-----------------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.26.
!-----------------------------------------------------------------------------------
! Dependence:: subroutine savgol.
! -----------------------------------------------------------------------------------

    implicit none
    integer(kind=4), intent(in):: nl, nr, ld, m, n1
    real   (kind=8), intent(inout):: y(n1)
    integer(kind=4), intent(out):: flag
    ! Local variables.
    integer(kind=4):: i, j, xl(nl+nr+1)
    real   (kind=8):: y0(n1), coef(nl+nr+1)

    xl(1) = 0

    y0 = y

    do i=1, nl

        xl(i+1) = -i

    end do

    do i=1, nr

        xl(1+nl+i) = nr-i+1

    end do

    call savgol(nl,nr,ld,m,coef,flag)

    if (flag/=0) return

    do i=1, n1-nr

        y(i) = 0.0

        do j=1, nl+nr+1

            if (i+xl(j) .gt. 0) then

                y(i) = y(i) + coef(j)*y0(i+xl(j))

            end if

        end do

    end do

    if (ld==0) then

        y(1:nl) = y0(1:nl)

        y(n1-nr+1:n1) = y0(n1-nr+1:n1)

    else 

        y(1:nl) = y(nl+1)

        y(n1-nr+1:n1) = y(n1-nr)
 
    end if

    return

end subroutine savgol_filter

subroutine savgol(nl,nr,ld,m,coef,flag)
!-----------------------------------------------------------------------------------
! This routine is used to calculate a set of Savitzky-Golay filter coefficients.
!-----------------------------------------------------------------------------------
!            nl:: input, integer, the number of leftward data points used.
!            nr:: input, integer, the number of rightward data points used.
!            ld:: input, integer, the order of the derivative desired.
!             m:: input, integer, the order of the smoothing polynomial.
! coef(nl+nr+1):: output, real values, calculated coefficents in wrap-around order.
!          flag:: output, integer, error message, 0=success, 1=failure.
!-----------------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.20.
!-----------------------------------------------------------------------------------
! Dependence:: subroutine ludcmp;
!              subroutine lubksb.
!-----------------------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.646 IN Press et al.
! -----------------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nl, nr, ld, m
    real   (kind=8), intent(inout):: coef(nl+nr+1)
    integer(kind=4), intent(out):: flag
    ! Local variables.
    integer(kind=4):: imj, ipj, k, kk, mm, indx(m+1)

    real   (kind=8):: d, fac, summ, a(m+1,m+1), b(m+1)

    flag = 0

    if (nl < 0 .or. nr < 0 .or. ld > m .or. nl+nr < m) then

        flag = 1

        return

    end if

    do ipj=0, 2*m

        summ = 0.0

        if (ipj .eq. 0) summ = 1.0

        do k=1, nr

            summ = summ + (float(k))**ipj

        end do

        do k=1, nl

            summ = summ + (float(-k))**ipj

        end do

        mm = min(ipj, 2*m-ipj)

        do imj=-mm, mm, 2

            a(1+(ipj+imj)/2,1+(ipj-imj)/2) = summ

        end do

    end do

    call ludcmp(a,m+1,indx,d,flag)

    if (flag .ne. 0) return

    b = 0.0

    b(ld+1) = 1.0

    call lubksb(a,m+1,indx,b)

    coef = 0.0

    do k=-nl, nr

        summ = b(1)

        fac = 1.0

        do mm=1, m

            fac = fac * k

            summ = summ + b(mm+1) * fac

        end do

        kk = mod(nl+nr+1-k, nl+nr+1) + 1

        coef(kk) = summ

    end do

    return

end subroutine savgol

subroutine lubksb(a,n,indx,b)
!-------------------------------------------------------------------------
!  a(n,n):: input, real values, the LU decomposition of a matrix.
!       n:: input, integer, the dimenstion of the matrix.
! indx(n):: input,  integer values, vector that records the row 
!           permutation effected by the partial pivoting.
!    b(n):: output, real values, the solution vector X for 
!                   linear equations A*X=B.
!-------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.18.
!-------------------------------------------------------------------------
! Dependence:: No.--------------------------------------------------------
!-------------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.39 IN Press et al.
! -------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n, indx(n)
    real   (kind=8), intent(in):: a(n,n)
    real   (kind=8), intent(inout):: b(n)
   ! Local variables.
    integer(kind=4):: i, ii, j, ll
    real   (kind=8):: summ

    ii = 0

    do i=1, n

        ll = indx(i)

        summ = b(ll)

        b(ll) = b(i)

        if (ii .ne. 0) then

            do j=ii, i-1

                summ = summ - a(i,j) * b(j)

            end do

        else if (summ .ne. 0.0) then
            
            ii = i

        end if

        b(i) = summ

    end do

    do i=n, 1, -1

        summ = b(i)

        do j=i+1, n

            summ = summ - a(i,j) * b(j)

        end do

        b(i) = summ / a(i,i)

    end do

    return

end subroutine lubksb

subroutine ludcmp(a,n,indx,d,flag)
!-------------------------------------------------------------------------
!This routine is used in combination with lubksb to solve 
!linear equations or invert a matrix.
!-------------------------------------------------------------------------
!  a(n,n):: input, real values, a matrix to be decomposed.
!       n:: input, integer, the dimension of the matrix.
! indx(n):: output, integer values, vector that records the row 
!           permutation effected by the partial pivoting.
!       d:: output, integer, output as 1 or -1 depending on whether 
!           the number of row interchanges was even or odd.
!    flag:: output, integer, error message, 0=success, 1=singular matrix.
!-------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.20.
!-------------------------------------------------------------------------
! Dependence:: No.--------------------------------------------------------
!-------------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.38 IN Press et al.
! ------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n
    integer(kind=4), intent(out):: indx(n), flag
    real   (kind=8), intent(inout):: a(n,n)
    real   (kind=8), intent(out):: d
    ! Local variables.
    integer(kind=4):: i, j, k, imax
    real   (kind=8):: aamax, dum, summ, vv(n)
   
    indx = 0

    flag = 0

    d = 1.0

    do i=1, n

        aamax = 0.0

        do j=1, n

            if (abs(a(i,j)) .gt. aamax)  aamax = abs(a(i,j))

        end do

        if (aamax .eq. 0.0) then 
 
            flag = 1

            return

        end if

        vv(i) = 1.0/aamax

    end do

    do j=1, n

        do i=1, j-1

            summ = a(i,j)

            do k=1, i-1

                summ = summ - a(i,k) * a(k,j)

            end do

            a(i,j) = summ

        end do

        aamax = 0.0

        do i=j, n

            summ = a(i,j)

            do k=1, j-1

                summ = summ - a(i,k) * a(k,j)

            end do

            a(i,j) = summ

            dum = vv(i) * abs(summ)

            if (dum .ge. aamax) then
       
                imax = i
          
                aamax = dum

            end if

        end do

        if (j .ne. imax) then

            do k=1, n

                dum = a(imax,k)

                a(imax,k) = a(j,k)

                a(j,k) = dum
  
            end do

            d = -d

            vv(imax) = vv(j)

        end if

        indx(j) = imax

        if (a(j,j) .eq. 0.0) a(j,j) = tiny(0.0D+00)

        if (j .ne. n) then

            dum = 1.0 / a(j,j)

            do i=j+1, n

                a(i,j) = a(i,j) * dum

            end do

        end if

    end do

    return

    end subroutine ludcmp    
      subroutine SEVAL (N,U,X,Y,B,C,D,V_0,V_1,V_2)
!------------------------------------------------------------------------
!     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
!     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
!     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
!     BY THE SPLINE SUBROUTINE.
!
!     INPUTS:
!     N       NUMBER OF POINTS OF CURVE Y = F(X)
!     U       ABSCISSA OF POINT TO BE INTERPOLATED
!     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
!             OF CURVE F(X)
!     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
!             CUBIC SPLINE
!
!     OUTPUTS:
!     SEVAL   INTERPOLATED VALUE
!             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
!
!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!------------------------------------------------------------------------
      REAL *8 B(N),C(N),D(N),X(N),Y(N),U,DX
      real *8 V_0,V_1,V_2
      DATA I/1/

!     BINARY SEARCH

      IF (I.GE.N) I = 1
      IF (U.LT.X(I)) GO TO 101
      IF (U.LE.X(I+1)) GO TO 301
101 I = 1
      J = N+1
201 K = (I+J)/2
      IF (U.LT.X(K)) J = K
      IF (U.GE.X(K)) I = K
      IF (J.GT.I+1) GO TO 201

!     SPLINE EVALUATION

   301 DX = U-X(I)
!      SEVAL = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      V_0 = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      V_1 = B(I)+2.D0*C(I)*DX+3.D0*D(I)*DX**2
      V_2 = 2.D0*C(I)+6.D0*D(I)*DX
      RETURN
      END

      SUBROUTINE SPLINE (N,X,Y,B,C,D)
!---------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
!     SPLINE TO BEST APPROXIMATE A DISCRETE FUNCTION GIVEN BY N POINTS
!
!     INPUTS:
!     N       NUMBER OF GIVEN POINTS
!     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
!             OF FUNCTION F(X)
!
!     OUTPUTS:
!     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
!             OF THE CUBIC SPLINE
!
!     REFERENCE:
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!---------------------------------------------------------------------
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION B(N),C(N),D(N),X(N),Y(N)
      NM1 = N-1
      IF (N.LT.2) RETURN
      IF (N.LT.3) GO TO 501

!     BUILD THE TRIDIAGONAL SYSTEM
!     B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)

      D(1) = X(2)-X(1)
      C(2) = (Y(2)-Y(1))/D(1)
      DO I = 2,NM1
      D(I) = X(I+1)-X(I)
      B(I) = 2.D0*(D(I-1)+D(I))
      C(I+1) = (Y(I+1)-Y(I))/D(I)
      C(I) = C(I+1)-C(I)
      END DO

!     CONDITIONS AT LIMITS
!     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES

      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.D0
      C(N) = 0.D0
      IF (N.EQ.3) GO TO 151
      C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))

!     FORWARD ELIMINATION

151 DO I = 2,N
      T = D(I-1)/B(I-1)
      B(I) = B(I)-T*D(I-1)
      C(I) = C(I)-T*C(I-1)
      END DO

!     BACK SUBSTITUTION

      C(N) = C(N)/B(N)
      DO L = 1,NM1
      I = N-L
      C(I) = (C(I)-D(I)*C(I+1))/B(I)
      END DO

!     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL

      B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D0*C(N))
      DO I = 1,NM1
      B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D0*C(I))
      D(I) = (C(I+1)-C(I))/D(I)
      C(I) = 3.D0*C(I)
      END DO
      C(N) = 3.D0*C(N)
      D(N) = D(NM1)
      RETURN

!     CAS N = 2

501 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.D0
      D(1) = 0.D0
      B(2) = B(1)
      C(2) = 0.D0
      D(2) = 0.D0
      RETURN
    END
!********************************************************************
 subroutine periodic(h,wave_period,wave_height,wave_amplitude)
!********************************************************************
! the subroutine is used to give the stroke of the wavemaker
 implicit none
 integer i
 real*8 pi,k0,k1,hs,s,e
 real*8 h,wave_period,wave_height,wave_amplitude
 pi=dacos(-1.d0)
    k0 = 1.d0
    e = 1.d0
    do while (e > 1.e-7)
        k1 = k0-(k0*dtanh(k0*h)-(2.d0*pi/wave_period)**2/9.81)/(dtanh(k0*h)+k0-k0*(dtanh(k0*h))**2)
        e = dabs((k1-k0)/k0)
        hs = 2.d0*(dcosh(2.d0*k1*h)-1.d0)/(dsinh(2.d0*k1*h)+2.d0*k1*h)
        s = wave_height/hs
        wave_amplitude = 0.5d0*s
        k0 = k1
    end do
    
    return
    end
!********************************************************************
    subroutine solitary(h0,wave_H,g,Npp,total_Nc,total_step,dt,dis,vb_x,acc)
!********************************************************************
    implicit none
    integer i,Npp,step,total_Nc,total_step,No(Npp)
    real*8 h0,wave_H,wave_K,wave_C,g,dt,x0
    real*8 vb_x(total_Nc,0:total_step-1),vb_z(total_Nc,0:total_step-1),x(total_Nc,0:total_step-1)
    real*8 dis(total_Nc,0:total_step-1),acc(total_Nc,0:total_step-1)
    
    call fenton_parameters(h0,wave_H,wave_K,wave_C,g)
    
    x0=-5.d0*wave_c
    x=0.d0
    do i=1,Npp
        No(i)=i
    end do

    do step = 0,total_step-1
        call wmbc(total_Nc,total_step,vb_x,vb_z,step,dt,wave_H,wave_K,wave_C,h0,x0,x,g,No,Npp)
    end do
    
    call get_wmk_dis_acc(total_Nc,total_step,vb_x,dt,dis,acc)

    return
    end
!********************************************************************
    subroutine get_wmk_dis_acc(total_Nc,total_step,vb_x,dt,dis,acc)
!********************************************************************
    implicit none
    integer i,j,total_Nc,total_step
    real*8 vb_x(total_Nc,0:total_step-1),dt,dis(total_Nc,0:total_step-1),acc(total_Nc,0:total_step-1)
    
    dis=0.d0
    acc=0.d0
    do i=1,total_Nc
        do j=1,total_step-1
            dis(i,j)=dis(i,j-1)+0.5d0*dt*(vb_x(i,j-1)+vb_x(i,j))
            if (j==total_step-1)then
                acc(i,j)=acc(i,j-1)
            else if (j==0)then
                acc(i,j)=(vb_x(i,j+1)-vb_x(i,j))/dt
            else
                acc(i,j)=0.5d0*(vb_x(i,j+1)-vb_x(i,j-1))/dt
            end if
        end do
    end do

    return
    end
!********************************************************************
    subroutine wmbc(total_Nc,total_step,vb_x,vb_z,step,dt,wave_H,wave_K,wave_C,h0,x0,x,g,No,Npp)
!********************************************************************
    ! the subroutine is used to give the velocity of the wavemaker
    ! x(total_Nc,total_step) is the coordinate of the observation point at time t
    ! x0 = initial position of wave crest
    implicit none
   integer j,total_Nc,total_step,step
   integer::zz,Npp,No(Npp)
   real*8 vb_x(total_Nc,0:total_step-1),vb_z(total_Nc,0:total_step-1),dt,wave_H,wave_K,wave_C,h0,x0
   real*8 x(total_Nc,0:total_step-1),zeta,kx,t,alpha,kk,cc,g,XX,S,xxx
   
   t=step*dt
   alpha=wave_H/h0
   kk=wave_K*h0
   cc=wave_C/(g*h0)**0.5
   
   do j=1,total_Nc
      do zz=1,Npp
          
         if(j.eq.No(zz)) then
             xxx=x(j,step) ! xxx: xi = position of the wave paddle at time t (step)
             XX=(xxx-x0-wave_C*t)/h0 ! XX: X = xi - Ct -x0
             S=1.d0/cosh(kk*XX) ! sech(K*X)
              zeta=S**2*alpha
              zeta=zeta+(-.75*S**2+.75*S**4)*alpha**2
              zeta=zeta+(.625*S**2-1.8875*S**4+1.2625*S**6)*alpha**3
           zeta=zeta+(-1.36817*S**2+3.88033*S**4-4.68304*S**6+2.17088*S**8)*alpha**4
        zeta=zeta+(1.86057*S**2-7.45136*S**4+12.7637*S**6-11.4199*S**8+4.24687*S**10)*alpha**5
    zeta=zeta+(-2.57413*S**2+13.2856*S**4-31.1191*S**6+40.1068*S**8-28.4272*S**10+8.728*S**12)*alpha**6
zeta=zeta+(3.4572*S**2-22.782*S**4+68.258*S**6-116.974*S**8+120.49*S**10-71.057*S**12+18.608*S**14)*alpha**7
zeta=zeta+(-4.6849*S**2+37.67*S**4-139.28*S**6+301.442*S**8-411.416*S**10+355.069*S**12-180.212*S**14+41.412*S**16)*alpha**8
zeta=zeta+(6.191*S**2-60.57*S**4+269.84*S**6-712.125*S**8+1217.98*S**10-1384.37*S**12+1023.07*S**14-450.29*S**16+90.279*S**18)*alpha**9
              zeta=zeta*h0
              vb_x(j,step)=wave_C*zeta/(h0+zeta) ! the velocity of the wavemaker
              goto 111
              
         else
            vb_x(j,step)=0.d0
         end if
         
      end do
     111  vb_z(j,step)=0.d0
   end do
   
    end
!********************************************************************
 subroutine fenton_parameters(h0,wave_H,wave_K,wave_C,g)
!********************************************************************
! the subroutine is used to give the parameters of the solitary wave by Fenton theory
! wave_H = wave height
! h0 = still water depth
! g = gravity
      real*8 wave_H,wave_K,wave_C,h0,g,alpha
      alpha=wave_H/h0
      wave_K=sqrt(3./4.d0*alpha)
      wave_K=wave_K*(1-.625*alpha+.554688*alpha**2-.561535*alpha**3+.567095*alpha**4-.602969*alpha**5 &
                  +.624914*alpha**6-.670850*alpha**7+.700371*alpha**8)
      wave_C=1+alpha
      wave_C=wave_C-.05*alpha**2-.0428571*alpha**3-.0342857*alpha**4-.0315195*alpha**5-.0292784*alpha**6 &
          -.0268451*alpha**7-.0302634*alpha**8-.0219347*alpha**9
      wave_K=wave_K/h0
      wave_C=wave_C**0.5*(g*h0)**0.5
    end