!-------------------------------------------------------------------------------
!------- DROPLET DISSOLUTION
!-------------------------------------------------------------------------------
!------- LATTICE CONSTANT FOR D2Q9 MODEL
!-------------------------------------------------------------------------------
MODULE D2Q9CONST
!------- D2Q9 WEIGHTS
	REAL(KIND=8),PARAMETER::W(1:9)=(/4.0D0/9.0D0,1.0D0/9.0D0, &
	& 1.0D0/9.0D0,1.0D0/9.0D0,1.0D0/9.0D0,1.0D0/36.0D0,1.0D0/36.0D0, &
	& 1.0D0/36.0D0,1.0D0/36.0D0/)
!------- D2Q9 DIRECTIONS
	INTEGER::E1(1:9)=(/0,1,0,-1,0,1,-1,-1,1/)
	INTEGER::E2(1:9)=(/0,0,1,0,-1,1,1,-1,-1/)
END MODULE D2Q9CONST
!-------------------------------------------------------------------------------
!------- SIMULATION PARAMETERS DETERMINATION
!-------------------------------------------------------------------------------
MODULE SIMPARAM
	INTEGER,PARAMETER::NI=200
	INTEGER,PARAMETER::NJ=80
	INTEGER,PARAMETER::NL=9
	INTEGER,PARAMETER::NDIM=2
	INTEGER,PARAMETER::TMAX=15000
	INTEGER,PARAMETER::NFLD=50
	INTEGER,PARAMETER::INIT=1
	REAL(KIND=8),PARAMETER::CS=1.0D0/SQRT(3.0D0)
	REAL(KIND=8),PARAMETER::A=-0.010D0
	REAL(KIND=8),PARAMETER::B= 0.010D0
	REAL(KIND=8),PARAMETER::KAPPA=-1.0D0*A
	REAL(KIND=8),PARAMETER::GAMA=0.50D0
	REAL(KIND=8),PARAMETER::TAUG=1.0D0/(3.0D0-SQRT(3.0D0))
	REAL(KIND=8),PARAMETER::RADIUS=16.0D0
	REAL(KIND=8),PARAMETER::AMURATIO=1.0D0
	REAL(KIND=8),PARAMETER::AMUO=1.0D0/6.0D0 !0.50D0
	REAL(KIND=8),PARAMETER::AMUW=AMUO*AMURATIO
END MODULE SIMPARAM	 
!-------------------------------------------------------------------------------
!------- FLOW SIMULATION USING LBE MODEL (MAIN LOOP)
!-------------------------------------------------------------------------------
	PROGRAM MAIN
	USE SIMPARAM,ONLY:NI,NJ,NL,NFLD,NDIM,TMAX,INIT
	IMPLICIT NONE
!-------
	REAL(KIND=8)::TIME1,TIME2,TIMETOTAL,M,D,AMUO,AMUW
	REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE::F,FEQ,U,G,GEQ,FRC
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE::RHO,PHI,AMU,TAUF
	INTEGER::TSTEP,I,J,L,IFLD
!-------
	ALLOCATE(F(NI,NJ,NL))
	ALLOCATE(FEQ(NI,NJ,NL))
	ALLOCATE(G(NI,NJ,NL))
	ALLOCATE(GEQ(NI,NJ,NL))
	ALLOCATE(U(NI,NJ,NDIM))
	ALLOCATE(RHO(NI,NJ))
	ALLOCATE(PHI(NI,NJ))
	ALLOCATE(AMU(NI,NJ))
	ALLOCATE(TAUF(NI,NJ))
	ALLOCATE(FRC(NI,NJ,NDIM))
!-------
100	FORMAT(E13.5)
!-------
	CALL COMPPARAM(M)
	CALL INITMACRO(RHO,U,PHI)
	CALL COMPEQ(FEQ,GEQ,RHO,U,PHI,AMU,TAUF)
!-------
	F=FEQ
	G=GEQ
!-------
	IFLD=0.0
	TIMETOTAL=0.0
	DO TSTEP=1,TMAX
	  CALL CPU_TIME(TIME1)
!-------
	  CALL COMPFRC(FRC,PHI,TSTEP)
	  CALL COMPEQ(FEQ,GEQ,RHO,U,PHI,AMU,TAUF)
	  CALL COLLISION(F,FEQ,G,GEQ,TAUF,FRC,U,RHO)
	  CALL STREAM(F)
	  CALL STREAM(G)
	  CALL COMPMACRO(F,RHO,U,G,PHI,FRC)
!-------
	  IF (MOD(TSTEP,NFLD).EQ.0) THEN
	    IFLD=IFLD+1
	    CALL UOUT(IFLD,U)
	    CALL RHOOUT(IFLD,RHO)
	    CALL PHIOUT(IFLD,PHI)
	    CALL AMUOUT(IFLD,AMU)
	    CALL TAUFOUT(IFLD,TAUF)
	    CALL PRESOUT(IFLD,U,RHO,F,PHI)
	    CALL FRCOUT(IFLD,FRC)
	  ENDIF
!-------
	  CALL CPU_TIME(TIME2)
	  TIMETOTAL=TIMETOTAL+(TIME2-TIME1)
	ENDDO
!-------
	WRITE(*,*) 'CPU TIME=',TIMETOTAL
!-------
	DEALLOCATE(F)
	DEALLOCATE(FEQ)
	DEALLOCATE(G)
	DEALLOCATE(GEQ)
	DEALLOCATE(U)
	DEALLOCATE(RHO)
	DEALLOCATE(PHI)
	DEALLOCATE(AMU)
	DEALLOCATE(TAUF)
	DEALLOCATE(FRC)
!-------
	END PROGRAM MAIN
!-------
!-------------------------------------------------------------------------------
!------- COMPUTATION OF DIMENSIONLESS PARAMETERS
!-------------------------------------------------------------------------------
	SUBROUTINE COMPPARAM(M)
	USE SIMPARAM,ONLY:CS,NI,TAUG,A,B,KAPPA,GAMA,RADIUS,TMAX
	IMPLICIT NONE
!-------
	REAL(KIND=8),INTENT(INOUT)::M
	REAL(KIND=8)::ZETA,ZETA2,DX,SIGMA,CH,CHSTAR,WE,PHI0,E,PI,FO,DI
!-------
	PI=3.141592654
!-------
	M=GAMA*(TAUG-0.50D0)
	ZETA=SQRT(-2.0D0*KAPPA/A)
	ZETA2=2.0*ZETA
	DX=1.0D0
	PHI0=SQRT(-A/B)
	SIGMA=4.0/3.0*KAPPA*PHI0**2.0D0/ZETA
!-------
	CH=ZETA/RADIUS
	CHSTAR=ZETA/DX
	WE=SIGMA/(ZETA*1.0D0*CS**2.0D0)
	E=PI*RADIUS**2.0D0/NI**2.0D0
	FO=-M*A*TMAX/RADIUS**2.0D0
	DI=CH/(E*(E-1.0)**2.0)
!-------
	WRITE(*,*) 'A=',A
	WRITE(*,*) 'B=',B
	WRITE(*,*) 'KAPPA=',KAPPA
	WRITE(*,*) 'GAMA=',GAMA
	WRITE(*,*) 'M=',M
	WRITE(*,*) 'ZETA=',ZETA
	WRITE(*,*) 'ZETA2=',ZETA2
	WRITE(*,*) 'PHI0=',PHI0
	WRITE(*,*) 'SIGMA=',SIGMA
	WRITE(*,*) 'TAUG=',TAUG
!-------
	WRITE(*,*) '-------------------------------------'
	WRITE(*,*) 'CH=',CH
	WRITE(*,*) 'CHSTAR=',CHSTAR
	WRITE(*,*) 'WE=',WE
	WRITE(*,*) 'E=',E
	WRITE(*,*) 'FO=',FO
	WRITE(*,*) 'DI=',DI
!-------
	END SUBROUTINE COMPPARAM
!-------------------------------------------------------------------------------
!------- INITIAL VALUES OF DENSITY AND VELOCITY        
!-------------------------------------------------------------------------------
	SUBROUTINE INITMACRO(RHO,U,PHI)
	USE SIMPARAM,ONLY:NI,NJ,NDIM,INIT,RADIUS,A,B
	IMPLICIT NONE
!-------
	REAL(KIND=8),INTENT(INOUT)::RHO(NI,NJ),U(NI,NJ,NDIM),PHI(NI,NJ)
	INTEGER::I,J,I0,J0,I1,J1,IP,JP
	REAL(KIND=8)::PHI0
	CHARACTER(LEN=22) DATFILE1,DATFILE2
	CHARACTER(LEN=21) DATFILE3,DATFILE4
!-------
100	FORMAT(E13.5)
!-------
	DATFILE1(1:22)='VELFIELDS/datx0100.ssv'
	DATFILE2(1:22)='VELFIELDS/daty0100.ssv'
	DATFILE3(1:21)='RHOFIELDS/dat0100.ssv'
	DATFILE4(1:21)='PHIFIELDS/dat0100.ssv'
!-------
	IF (INIT==2) THEN
	OPEN (UNIT=1310,FILE=DATFILE1,FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN (UNIT=1311,FILE=DATFILE2,FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN (UNIT=1312,FILE=DATFILE3,FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN (UNIT=1313,FILE=DATFILE4,FORM='FORMATTED',STATUS='UNKNOWN')
!-------
	DO J=1,NJ
	DO I=1,NI
	  READ(1310,100) U(I,J,1)
	  READ(1311,100) U(I,J,2)
	  READ(1312,100) RHO(I,J)
	  READ(1313,100) PHI(I,J)
	ENDDO
	ENDDO
!-------
	DO J=1,NJ
	DO I=1,INT(NI/2)
	  IF ((PHI(I,J)<1.050D0).AND.(PHI(I,J)>-0.920D0)) THEN
	    IP=I
	    JP=J
	    U(IP,JP,1)=U(IP,JP,1)+0.050D0
	  ENDIF
	ENDDO
	ENDDO
!-------
	CLOSE(UNIT=1310)
	CLOSE(UNIT=1311)
	CLOSE(UNIT=1312)
	CLOSE(UNIT=1313)
!-------
	ELSE
!-------
	U(:,:,1)=0.0D0
	U(:,:,2)=0.0D0
	RHO(:,:)=1.0D0
	PHI(:,:)=-1.0D0
!-------
	I0=INT(NI/2.0)
	J0=INT(NJ/2.0)
!-------
	DO J=1,NJ
	DO I=1,NI
	  I1=I-I0-1.5*RADIUS
	  J1=J-J0
	    IF (SQRT(I1**2.0+J1**2.0)<=RADIUS) THEN
	      PHI(I,J)=1.0D0
	    ENDIF
	  I1=I-I0+1.5*RADIUS !+INT(RADIUS/8.0D0)
	  J1=J-J0
	    IF (SQRT(I1**2.0+J1**2.0)<=RADIUS) THEN
	      PHI(I,J)=1.0D0
	    ENDIF
	ENDDO
	ENDDO
!-------
	ENDIF
!-------	  
	END SUBROUTINE INITMACRO
!-------------------------------------------------------------------------------
!-------  COMPUTATION OF EQUILIBRIUM DISTRIBUTION FUNCTION
!-------------------------------------------------------------------------------
	SUBROUTINE COMPEQ(FEQ,GEQ,RHO,U,PHI,AMU,TAUF)
	USE SIMPARAM,ONLY:NI,NJ,NL,NDIM,KAPPA,GAMA,A,B,CS,AMUO,AMUW
	USE D2Q9CONST
	IMPLICIT NONE
!-------
	REAL(KIND=8),INTENT(IN)::RHO(NI,NJ),U(NI,NJ,NDIM),PHI(NI,NJ)
	REAL(KIND=8),INTENT(INOUT)::FEQ(NI,NJ,NL),GEQ(NI,NJ,NL),AMU(NI,NJ),TAUF(NI,NJ)
	REAL(KIND=8)::UXY,UXEQ,UYEQ,USQR,EC1,EC2,Q,R,MM,N,ALPHA,SUMM,PHI0
        INTEGER::I,J,L,IM,IP,JM,JP
	REAL(KIND=8),DIMENSION(:,:)::D2PHI(NI,NJ),DPHIDX(NI,NJ),DPHIDY(NI,NJ), &
& P0(NI,NJ),PXX(NI,NJ),PYY(NI,NJ),PXY(NI,NJ),F2(NI,NJ),F3(NI,NJ),F6(NI,NJ),F7(NI,NJ), &
& MU(NI,NJ),SUMF(NI,NJ),SUMG(NI,NJ)
!-------
	ALPHA=1.0D0/4.0D0
!-------
	DO J=1,NJ
	DO I=1,NI
	DO L=2,NL
	  EC1=E1(L)
	  EC2=E2(L)
	  UXEQ=U(I,J,1)
	  UYEQ=U(I,J,2)
	  UXY=UXEQ*EC1+UYEQ*EC2
	  USQR=(UXEQ**2.0+UYEQ**2.0)
	  FEQ(I,J,L)=W(L)*RHO(I,J)*(3.0*UXY+4.5*UXY*UXY-1.5*USQR)
	  GEQ(I,J,L)=W(L)*PHI(I,J)*(3.0*UXY+4.5*UXY*UXY-1.5*USQR)
	ENDDO
	ENDDO
	ENDDO
!-------------------------------------------------------------------------------
	Q=4.0D0
	R=1.0D0
	MM=1.0D0
	N=4.0D0
!-------
	DO J=1,NJ
	  JP=J+1
	  IF (JP.GT.NJ) JP=JP-NJ
	  JM=J-1
	  IF (JM.LT.1) JM=JM+NJ
	    DO I=1,NI
	      IP=I+1
	      IF (IP.GT.NI) IP=IP-NI
	      IM=I-1
	      IF (IM.LT.1) IM=IM+NI
!-------
	  D2PHI(I,J)=(-4.0D0*(Q+R)*PHI(I,J)+Q*(PHI(I,JP)+PHI(I,JM)+PHI(IP,J)+PHI(IM,J))+ &
& R*(PHI(IP,JP)+PHI(IM,JP)+PHI(IM,JM)+PHI(IP,JM)))/6.0D0
!	  WRITE(*,*) D2PHI(I,J)
!-------
	  DPHIDX(I,J)=(MM*(PHI(IP,JP)-PHI(IM,JP)-PHI(IM,JM)+PHI(IP,JM))+N*(PHI(IP,J)-PHI(IM,J)))/12.0D0
	  DPHIDY(I,J)=(MM*(PHI(IP,JP)+PHI(IM,JP)-PHI(IM,JM)-PHI(IP,JM))+N*(PHI(I,JP)-PHI(I,JM)))/12.0D0
!	  WRITE(*,*) DPHIDY(I,J)
	ENDDO
	ENDDO
!-------------------------------------------------------------------------------
	DO J=1,NJ
	DO I=1,NI
	  P0(I,J)=A/2.0D0*PHI(I,J)*PHI(I,J)+3.0D0/4.0D0*B*(PHI(I,J)**4.0D0)- &
& KAPPA*PHI(I,J)*D2PHI(I,J)-KAPPA/2.0D0*(DPHIDX(I,J)*DPHIDX(I,J)+DPHIDY(I,J)*DPHIDY(I,J))
!	  WRITE(*,*) P0(I,J)
!-------
	  PXX(I,J)=RHO(I,J)*CS**2.0+P0(I,J)+KAPPA*DPHIDX(I,J)*DPHIDX(I,J)
	  PYY(I,J)=RHO(I,J)*CS**2.0+P0(I,J)+KAPPA*DPHIDY(I,J)*DPHIDY(I,J)
	  PXY(I,J)=KAPPA*DPHIDX(I,J)*DPHIDY(I,J)
!	  WRITE(*,*) PXY(I,J)
!-------
	  F2(I,J)=(1.0D0-ALPHA)/2.0D0*PXX(I,J)-ALPHA/2.0D0*PYY(I,J)
	  F3(I,J)=(1.0D0-ALPHA)/2.0D0*PYY(I,J)-ALPHA/2.0D0*PXX(I,J)
	  F6(I,J)=ALPHA/4.0D0*(PXX(I,J)+PYY(I,J))+PXY(I,J)/4.0D0
	  F7(I,J)=ALPHA/4.0D0*(PXX(I,J)+PYY(I,J))-PXY(I,J)/4.0D0
!	  WRITE(*,*) F7(I,J)
!-------
	  MU(I,J)=A*PHI(I,J)+B*PHI(I,J)**3.0D0-KAPPA*D2PHI(I,J)
!	  IF (MU(I,J).NE. 0.0) THEN
!	  WRITE(*,*) MU(I,J)
!	  ENDIF
	ENDDO
	ENDDO
!-------
!	DO J=1,NJ
!	  WRITE(*,*) J,(MU(I,J),I=1,NI)
!	ENDDO
!-------
!-------------------------------------------------------------------------------
	DO J=1,NJ
	DO I=1,NI
	  FEQ(I,J,2)=FEQ(I,J,2)+F2(I,J)
	  FEQ(I,J,4)=FEQ(I,J,4)+F2(I,J)
	  FEQ(I,J,3)=FEQ(I,J,3)+F3(I,J)
	  FEQ(I,J,5)=FEQ(I,J,5)+F3(I,J)
	  FEQ(I,J,6)=FEQ(I,J,6)+F6(I,J)
	  FEQ(I,J,8)=FEQ(I,J,8)+F6(I,J)
	  FEQ(I,J,7)=FEQ(I,J,7)+F7(I,J)
	  FEQ(I,J,9)=FEQ(I,J,9)+F7(I,J)
	DO L=2,NL
	  GEQ(I,J,L)=GEQ(I,J,L)+W(L)*GAMA*MU(I,J)/(CS**2.0D0)
	ENDDO
	ENDDO
	ENDDO
!-------
	SUMF=0.0
	SUMG=0.0
!-------
	DO J=1,NJ
	DO I=1,NI
	DO L=2,NL
	  SUMF(I,J)=SUMF(I,J)+FEQ(I,J,L)
	  SUMG(I,J)=SUMG(I,J)+GEQ(I,J,L)
	ENDDO
	  FEQ(I,J,1)=RHO(I,J)-SUMF(I,J)
	  GEQ(I,J,1)=PHI(I,J)-SUMG(I,J)
	ENDDO
	ENDDO
!-------
!	SUMM=0.0
!	DO L=1,NL
!	  SUMM=SUMM+FEQ(10,10,L)
!	ENDDO
!	WRITE(*,*) 'SUMM=',SUMM
!-------
	PHI0=SQRT(-A/B)
	DO J=1,NJ
	DO I=1,NI
	  AMU(I,J)=PHI0/((PHI0-PHI(I,J)/(2.0D0*AMUW))+(PHI0+PHI(I,J))/(2.0D0*AMUO))
	  TAUF(I,J)=AMU(I,J)/(CS**2.0D0)+0.50D0
	ENDDO
	ENDDO
!-------
	END SUBROUTINE COMPEQ
!-------------------------------------------------------------------------------
!------- COMPUTATION OF MACROSCOPIC DENSITY (RHO) AND VELOCITY U=(UX,UY)
!-------------------------------------------------------------------------------
	SUBROUTINE COMPMACRO(F,RHO,U,G,PHI,FRC)
	USE SIMPARAM,ONLY:NI,NJ,NL,NDIM
	USE D2Q9CONST
	IMPLICIT NONE
!-------
	REAL(KIND=8),INTENT(IN)::F(NI,NJ,NL),G(NI,NJ,NL),FRC(NI,NJ,NDIM)
	REAL(KIND=8),INTENT(INOUT)::RHO(NI,NJ),U(NI,NJ,NDIM),PHI(NI,NJ)
        INTEGER::I,J,L
!-------
        DO J=1,NJ
        DO I=1,NI
	  RHO(I,J)=F(I,J,1)+F(I,J,2)+F(I,J,3)+F(I,J,4)+F(I,J,5)+ &
 & F(I,J,6)+F(I,J,7)+F(I,J,8)+F(I,J,9)
	  U(I,J,1)=(F(I,J,2)-F(I,J,4)+F(I,J,6)-F(I,J,7)- &
 & F(I,J,8)+F(I,J,9)+0.5*FRC(I,J,1))/RHO(I,J)
	  U(I,J,2)=(F(I,J,3)-F(I,J,5)+F(I,J,6)+F(I,J,7)- &
 & F(I,J,8)-F(I,J,9)+0.5*FRC(I,J,2))/RHO(I,J)
	  PHI(I,J)=G(I,J,1)+G(I,J,2)+G(I,J,3)+G(I,J,4)+G(I,J,5)+ &
 & G(I,J,6)+G(I,J,7)+G(I,J,8)+G(I,J,9)
	ENDDO
	ENDDO
!-------
	END SUBROUTINE COMPMACRO
!-------------------------------------------------------------------------------
!-------  COLLISION STEP FOR F
!-------------------------------------------------------------------------------
	SUBROUTINE COLLISION(F,FEQ,G,GEQ,TAUF,FRC,U,RHO)
	USE SIMPARAM,ONLY:NI,NJ,NL,TAUG,NDIM,CS
	USE D2Q9CONST
	IMPLICIT NONE
!-------
	REAL(KIND=8),INTENT(IN)::FEQ(NI,NJ,NL),GEQ(NI,NJ,NL),TAUF(NI,NJ),FRC(NI,NJ,NDIM), &
 & U(NI,NJ,NDIM),RHO(NI,NJ)
	REAL(KIND=8),INTENT(INOUT)::F(NI,NJ,NL),G(NI,NJ,NL)
	REAL(KIND=8),DIMENSION(:,:,:)::FORCE(NI,NJ,NL)
        INTEGER::I,J,L
!-------
	DO J=1,NJ
	DO I=1,NI
	DO L=1,NL
	  FORCE(I,J,L)=(1.0D0-1.0D0/(2.0D0*TAUF(I,J)))*W(L)*(((E1(L)-U(I,J,1))*FRC(I,J,1)+ &
 & (E2(L)-U(I,J,2))*FRC(I,J,2))/(CS**2.0D0)+(E1(L)*U(I,J,1)+ &
 &  E2(L)*U(I,J,2))*(E1(L)*FRC(I,J,1)+E2(L)*FRC(I,J,2))/(CS**4.0D0))
	  F(I,J,L)=F(I,J,L)-(F(I,J,L)-FEQ(I,J,L))/TAUF(I,J)+FORCE(I,J,L)
	  G(I,J,L)=G(I,J,L)-(G(I,J,L)-GEQ(I,J,L))/TAUG
	ENDDO
	ENDDO
	ENDDO
!-------
	END SUBROUTINE COLLISION
!-------------------------------------------------------------------------------
!       STREAMING STEP
!-------------------------------------------------------------------------------
	SUBROUTINE STREAM(ARRAY)
	USE SIMPARAM,ONLY:NI,NJ,NL
	IMPLICIT NONE
!-------
	REAL(KIND=8),INTENT(INOUT)::ARRAY(NI,NJ,NL)
	REAL(KIND=8),DIMENSION(:,:,:)::ARRAYSTR(NI,NJ,NL)
        INTEGER::I,J,JM,JP,IM,IP
!-------
	DO J=1,NJ
          IF (J.GT.1) THEN 
           JM=J-1
          ELSE
           JM=NJ
          ENDIF
!-------      
          IF(J.LT.NJ) THEN
           JP=J+1
          ELSE
           JP=1
          ENDIF
!-------
        DO I=1,NI
          IF(I.GT.1) THEN
            IM=I-1
          ELSE 
            IM=NI 
          ENDIF
!-------
          IF(I.LT.NI) THEN 
            IP=I+1
          ELSE 
            IP=1
          ENDIF
!-------
          ARRAYSTR(I, J, 1)=ARRAY(I,J,1)
          ARRAYSTR(IP,J, 2)=ARRAY(I,J,2)
          ARRAYSTR(I, JP,3)=ARRAY(I,J,3)
          ARRAYSTR(IM,J, 4)=ARRAY(I,J,4)
          ARRAYSTR(I, JM,5)=ARRAY(I,J,5)
          ARRAYSTR(IP,JP,6)=ARRAY(I,J,6)
          ARRAYSTR(IM,JP,7)=ARRAY(I,J,7)
          ARRAYSTR(IM,JM,8)=ARRAY(I,J,8)
          ARRAYSTR(IP,JM,9)=ARRAY(I,J,9)
!-------
	ENDDO
	ENDDO
!-------
	ARRAY=ARRAYSTR
!-------
	END SUBROUTINE STREAM
!-------------------------------------------------------------------------------
!------- VELOCITY FIELD SAVING 
!-------------------------------------------------------------------------------
	SUBROUTINE UOUT(IFLD,U)
	USE SIMPARAM,ONLY:NI,NJ,NDIM
	IMPLICIT NONE
!-------
100	FORMAT(E13.5)
	REAL(KIND=8),INTENT(IN)::U(NI,NJ,NDIM)
	INTEGER,INTENT(INOUT)::IFLD
!-------
	INTEGER::I,J,D0,D1,D2,D3
	CHARACTER(LEN=22) DATFILE1,DATFILE2
	CHARACTER(LEN=21) DATFILE
!-------
	DATFILE1(1:14)='VELFIELDS/datx'
	DATFILE2(1:14)='VELFIELDS/daty'
	DATFILE(1:13)='VELFIELDS/dat'
!-------
	IF (IFLD .LT. 10) THEN
	  D0=0
	  D1=0
	  D2=0
	  D3=IFLD
	ENDIF
!-------
	IF ((IFLD .GE. 10) .AND. (IFLD .LT. 100)) THEN
	  D0=0
	  D1=0
	  D2=INT(IFLD/10)
	  D3=IFLD-D2*10
	ENDIF
!-------
	IF ((IFLD .GE. 100) .AND. (IFLD .LT. 1000)) THEN
	  D0=0
	  D1=INT(IFLD/100)
	  D2=IFLD-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D1*100-D2*10
	ENDIF
!-------
	IF (IFLD .GE. 1000) THEN
	  D0=INT(IFLD/1000)
	  D1=IFLD-D0*1000
	  D1=INT(D1/100)
	  D2=IFLD-D0*1000-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D0*1000-D1*100-D2*10
	ENDIF
!-------
	DATFILE1(15:15)=CHAR(D0+48)
	DATFILE1(16:16)=CHAR(D1+48)
	DATFILE1(17:17)=CHAR(D2+48)
	DATFILE1(18:18)=CHAR(D3+48)
	DATFILE1(19:22)='.ssv'
!-------
	DATFILE2(15:15)=CHAR(D0+48)
	DATFILE2(16:16)=CHAR(D1+48)
	DATFILE2(17:17)=CHAR(D2+48)
	DATFILE2(18:18)=CHAR(D3+48)
	DATFILE2(19:22)='.ssv'
!-------
	DATFILE(14:14)=CHAR(D0+48)
	DATFILE(15:15)=CHAR(D1+48)
	DATFILE(16:16)=CHAR(D2+48)
	DATFILE(17:17)=CHAR(D3+48)
	DATFILE(18:21)='.ssv'
!-------
	OPEN (UNIT=310,FILE=DATFILE1,FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN (UNIT=311,FILE=DATFILE2,FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN (UNIT=312,FILE=DATFILE,FORM='FORMATTED',STATUS='UNKNOWN')
!-------
	DO J=1,NJ
	DO I=1,NI
	  WRITE(310,100) U(I,J,1)
	  WRITE(311,100) U(I,J,2)
	  WRITE(312,100) SQRT(U(I,J,1)**2.0D0+U(I,J,2)**2.0D0)
	ENDDO
	ENDDO
!-------
	CLOSE(UNIT=310)
	CLOSE(UNIT=311)
	CLOSE(UNIT=312)
!-------
	END SUBROUTINE UOUT
!-------------------------------------------------------------------------------
!-------  DENSITY FIELD SAVING 
!-------------------------------------------------------------------------------
	SUBROUTINE RHOOUT(IFLD,RHO)
	USE SIMPARAM,ONLY:NI,NJ
	IMPLICIT NONE
!-------
100	FORMAT(E13.5)
	REAL(KIND=8),INTENT(IN)::RHO(NI,NJ)
	INTEGER,INTENT(INOUT)::IFLD
!-------
	INTEGER::I,J,D0,D1,D2,D3
	CHARACTER DATFILE*(21)
!-------
	DATFILE(1:13)='RHOFIELDS/dat'
!-------
	IF (IFLD .LT. 10) THEN
	  D0=0
	  D1=0
	  D2=0
	  D3=IFLD
	ENDIF
!-------
	IF ((IFLD .GE. 10) .AND. (IFLD .LT. 100)) THEN
	  D0=0
	  D1=0
	  D2=INT(IFLD/10)
	  D3=IFLD-D2*10
	ENDIF
!-------
	IF ((IFLD .GE. 100) .AND. (IFLD .LT. 1000)) THEN
	  D0=0
	  D1=INT(IFLD/100)
	  D2=IFLD-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D1*100-D2*10
	ENDIF
!-------
	IF (IFLD .GE. 1000) THEN
	  D0=INT(IFLD/1000)
	  D1=IFLD-D0*1000
	  D1=INT(D1/100)
	  D2=IFLD-D0*1000-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D0*1000-D1*100-D2*10
	ENDIF
!-------
	DATFILE(14:14)=CHAR(D0+48)
	DATFILE(15:15)=CHAR(D1+48)
	DATFILE(16:16)=CHAR(D2+48)
	DATFILE(17:17)=CHAR(D3+48)
	DATFILE(18:21)='.ssv'
!-------
	OPEN (UNIT=410,FILE=TRIM(DATFILE),FORM='FORMATTED',STATUS='UNKNOWN')
!-------
	DO J=1,NJ
	DO I=1,NI
	  WRITE(410,100) RHO(I,J)
	ENDDO
	ENDDO
!-------
	CLOSE(UNIT=410)
!-------
	END SUBROUTINE RHOOUT
!--------------------------------------------------------------------------------
!-------  PHI FIELD SAVING 
!--------------------------------------------------------------------------------
	SUBROUTINE PHIOUT(IFLD,PHI)
	USE SIMPARAM,ONLY:NI,NJ
	IMPLICIT NONE
!-------
100	FORMAT(E13.5)
	REAL(KIND=8),INTENT(IN)::PHI(NI,NJ)
	INTEGER,INTENT(INOUT)::IFLD
!-------
	INTEGER::I,J,D0,D1,D2,D3
	CHARACTER DATFILE*(21)
!-------
	DATFILE(1:13)='PHIFIELDS/dat'
!-------
	IF (IFLD .LT. 10) THEN
	  D0=0
	  D1=0
	  D2=0
	  D3=IFLD
	ENDIF
!-------
	IF ((IFLD .GE. 10) .AND. (IFLD .LT. 100)) THEN
	  D0=0
	  D1=0
	  D2=INT(IFLD/10)
	  D3=IFLD-D2*10
	ENDIF
!-------
	IF ((IFLD .GE. 100) .AND. (IFLD .LT. 1000)) THEN
	  D0=0
	  D1=INT(IFLD/100)
	  D2=IFLD-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D1*100-D2*10
	ENDIF
!-------
	IF (IFLD .GE. 1000) THEN
	  D0=INT(IFLD/1000)
	  D1=IFLD-D0*1000
	  D1=INT(D1/100)
	  D2=IFLD-D0*1000-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D0*1000-D1*100-D2*10
	ENDIF
!-------
	DATFILE(14:14)=CHAR(D0+48)
	DATFILE(15:15)=CHAR(D1+48)
	DATFILE(16:16)=CHAR(D2+48)
	DATFILE(17:17)=CHAR(D3+48)
	DATFILE(18:21)='.ssv'
!-------
	OPEN (UNIT=510,FILE=TRIM(DATFILE),FORM='FORMATTED',STATUS='UNKNOWN')
!-------
	DO J=1,NJ
	DO I=1,NI
	  WRITE(510,100) PHI(I,J)
	ENDDO
	ENDDO
!-------
	CLOSE(UNIT=510)
!-------
	END SUBROUTINE PHIOUT
!----------------------------------------------------------------------------------
!------- CHECK
!----------------------------------------------------------------------------------
	SUBROUTINE CHECKAC(F,D)
	USE SIMPARAM,ONLY:NI,NJ,NL
	IMPLICIT NONE
!-------
	REAL(KIND=8),INTENT(IN)::F(NI,NJ,NL)
	INTEGER::I,J,L
	REAL(KIND=8),INTENT(OUT)::D
!-------
	D=0.0
	J=1
	I=1
	DO L=1,NL
	  D=D+F(I,J,L)
	ENDDO
!-------	  
	END SUBROUTINE CHECKAC
!----------------------------------------------------------------------------------
!------- CHECK BEFORE COLLISION
!----------------------------------------------------------------------------------
	SUBROUTINE CHECKBC(F,D)
	USE SIMPARAM,ONLY:NI,NJ,NL,CS,NDIM
	USE D2Q9CONST
	IMPLICIT NONE
!-------
	REAL(KIND=8),INTENT(IN)::F(NI,NJ,NL)
	INTEGER::I,J,L
	REAL(KIND=8),INTENT(OUT)::D
!-------
	D=0.0
	J=1
	I=1
	DO L=1,NL
	  D=D+F(I,J,L)
	ENDDO
!-------	  
	END SUBROUTINE CHECKBC
!--------------------------------------------------------------------------------
!-------  PRES FIELD SAVING 
!--------------------------------------------------------------------------------
	SUBROUTINE PRESOUT(IFLD,U,RHO,F,PHI)
	USE D2Q9CONST
	USE SIMPARAM,ONLY:NI,NJ,NL,NDIM,CS,A,B,KAPPA
	IMPLICIT NONE
!-------
100	FORMAT(E13.5)
	REAL(KIND=8),INTENT(IN)::F(NI,NJ,NL),RHO(NI,NJ),U(NI,NJ,NDIM),PHI(NI,NJ)
	INTEGER,INTENT(INOUT)::IFLD
!-------
	INTEGER::I,J,IP,JP,IM,JM,D0,D1,D2,D3
	CHARACTER DATFILE*(21)
	REAL(KIND=8)::MM,Q,R,N
	REAL(KIND=8),DIMENSION(:,:)::DPXX(NI,NJ),DPYY(NI,NJ),P0(NI,NJ),D2PHI(NI,NJ), &
 & DPHIDX(NI,NJ),DPHIDY(NI,NJ)
!-------
	DATFILE(1:13)='PRSFIELDS/dat'
!-------
	IF (IFLD .LT. 10) THEN
	  D0=0
	  D1=0
	  D2=0
	  D3=IFLD
	ENDIF
!-------
	IF ((IFLD .GE. 10) .AND. (IFLD .LT. 100)) THEN
	  D0=0
	  D1=0
	  D2=INT(IFLD/10)
	  D3=IFLD-D2*10
	ENDIF
!-------
	IF ((IFLD .GE. 100) .AND. (IFLD .LT. 1000)) THEN
	  D0=0
	  D1=INT(IFLD/100)
	  D2=IFLD-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D1*100-D2*10
	ENDIF
!-------
	IF (IFLD .GE. 1000) THEN
	  D0=INT(IFLD/1000)
	  D1=IFLD-D0*1000
	  D1=INT(D1/100)
	  D2=IFLD-D0*1000-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D0*1000-D1*100-D2*10
	ENDIF
!-------
	DATFILE(14:14)=CHAR(D0+48)
	DATFILE(15:15)=CHAR(D1+48)
	DATFILE(16:16)=CHAR(D2+48)
	DATFILE(17:17)=CHAR(D3+48)
	DATFILE(18:21)='.ssv'
!-------
	OPEN (UNIT=610,FILE=TRIM(DATFILE),FORM='FORMATTED',STATUS='UNKNOWN')
!-------
	DO J=1,NJ
	DO I=1,NI
	  P0(I,J)=RHO(I,J)*(CS**2.0D0)+A*PHI(I,J)**2.0D0+3.0D0/4.0D0*PHI(I,J)**4.0D0
	  WRITE(610,100) P0(I,J)
	ENDDO
	ENDDO
!-------
!	DO J=1,NJ
!	DO I=1,NI
!	  DPXX(I,J)=E1(1)*E2(1)*F(I,J,1)+E1(2)*E2(2)*F(I,J,2)+E1(3)*E2(3)*F(I,J,3)+ &
! & E1(4)*E2(4)*F(I,J,4)+E1(5)*E2(5)*F(I,J,5)+E1(6)*E2(6)*F(I,J,6)+ &
! & E1(7)*E2(7)*F(I,J,7)+E1(8)*E2(8)*F(I,J,8)+E1(9)*E2(9)*F(I,J,9)-RHO(I,J)*U(I,J,1)*U(I,J,1)
!	  DPYY(I,J)=E1(1)*E2(1)*F(I,J,1)+E1(2)*E2(2)*F(I,J,2)+E1(3)*E2(3)*F(I,J,3)+ &
! & E1(4)*E2(4)*F(I,J,4)+E1(5)*E2(5)*F(I,J,5)+E1(6)*E2(6)*F(I,J,6)+ &
! & E1(7)*E2(7)*F(I,J,7)+E1(8)*E2(8)*F(I,J,8)+E1(9)*E2(9)*F(I,J,9)-RHO(I,J)*U(I,J,2)*U(I,J,2)
!	  P0(I,J)=(0.50D0*(DPXX(I,J)+DPYY(I,J))-RHO(I,J)*CS**2.0D0)/(RHO(I,J)*CS**2.0D0)
!	  WRITE(610,100) P0(I,J)
!	ENDDO
!	ENDDO
!-------
!	Q=1.0D0
!	R=0.0D0
!	MM=0.0D0
!	N=0.50D0
!-------
!	DO J=1,NJ
!	  JP=J+1
!	  IF (JP.GT.NJ) JP=JP-NJ
!	  JM=J-1
!	  IF (JM.LT.1) JM=JM+NJ
!	    DO I=1,NI
!	      IP=I+1
!	      IF (IP.GT.NI) IP=IP-NI
!	      IM=I-1
!	      IF (IM.LT.1) IM=IM+NI
!-------
!	  D2PHI(I,J)=-4.0D0*(Q+R)*PHI(I,J)+Q*(PHI(I,JP)+PHI(I,JM)+PHI(IP,J)+PHI(IM,J))+ &
!& R*(PHI(IP,JP)+PHI(IM,JP)+PHI(IM,JM)+PHI(IP,JM))
!	  WRITE(*,*) D2PHI(I,J)
!-------
!	  DPHIDX(I,J)=(MM* (PHI(IP,JP)-PHI(IM,JP)-PHI(IM,JM)+PHI(IP,JM))+N*(PHI(IP,J)-PHI(IM,J)))
!	  DPHIDY(I,J)=(MM*(-PHI(IP,JP)-PHI(IM,JP)+PHI(IM,JM)+PHI(IP,JM))+N*(PHI(I,JP)-PHI(I,JM)))
!	  WRITE(*,*) DPHIDY(I,J)
!	ENDDO
!	ENDDO
!-------------------------------------------------------------------------------
!	DO J=1,NJ
!	DO I=1,NI
!	  P0(I,J)=(A/2.0D0*PHI(I,J)*PHI(I,J)+3.0D0/4.0D0*B*(PHI(I,J)**4.0D0)- &
! & KAPPA*PHI(I,J)*D2PHI(I,J)-KAPPA/2.0D0*(DPHIDX(I,J)*DPHIDX(I,J)+DPHIDY(I,J)*DPHIDY(I,J))-RHO(I,J)*CS**2.0D0)/ &
! & (RHO(I,J)*CS**2.0D0)
!	ENDDO
!	ENDDO
!-------
	CLOSE(UNIT=610)
!-------
	END SUBROUTINE PRESOUT
!----------------------------------------------------------------------------------
!-------  PHI FIELD SAVING 
!--------------------------------------------------------------------------------
	SUBROUTINE AMUOUT(IFLD,AMU)
	USE SIMPARAM,ONLY:NI,NJ
	IMPLICIT NONE
!-------
100	FORMAT(E13.5)
	REAL(KIND=8),INTENT(IN)::AMU(NI,NJ)
	INTEGER,INTENT(INOUT)::IFLD
!-------
	INTEGER::I,J,D0,D1,D2,D3
	CHARACTER DATFILE*(21)
!-------
	DATFILE(1:13)='AMUFIELDS/dat'
!-------
	IF (IFLD .LT. 10) THEN
	  D0=0
	  D1=0
	  D2=0
	  D3=IFLD
	ENDIF
!-------
	IF ((IFLD .GE. 10) .AND. (IFLD .LT. 100)) THEN
	  D0=0
	  D1=0
	  D2=INT(IFLD/10)
	  D3=IFLD-D2*10
	ENDIF
!-------
	IF ((IFLD .GE. 100) .AND. (IFLD .LT. 1000)) THEN
	  D0=0
	  D1=INT(IFLD/100)
	  D2=IFLD-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D1*100-D2*10
	ENDIF
!-------
	IF (IFLD .GE. 1000) THEN
	  D0=INT(IFLD/1000)
	  D1=IFLD-D0*1000
	  D1=INT(D1/100)
	  D2=IFLD-D0*1000-D1*100

	  D2=INT(D2/10)
	  D3=IFLD-D0*1000-D1*100-D2*10
	ENDIF
!-------
	DATFILE(14:14)=CHAR(D0+48)
	DATFILE(15:15)=CHAR(D1+48)
	DATFILE(16:16)=CHAR(D2+48)
	DATFILE(17:17)=CHAR(D3+48)
	DATFILE(18:21)='.ssv'
!-------
	OPEN (UNIT=510,FILE=TRIM(DATFILE),FORM='FORMATTED',STATUS='UNKNOWN')
!-------
	DO J=1,NJ
	DO I=1,NI
	  WRITE(510,100) AMU(I,J)
	ENDDO
	ENDDO
!-------
	CLOSE(UNIT=510)
!-------
	END SUBROUTINE AMUOUT
!----------------------------------------------------------------------------------
!-------  PHI FIELD SAVING 
!--------------------------------------------------------------------------------
	SUBROUTINE TAUFOUT(IFLD,TAUF)
	USE SIMPARAM,ONLY:NI,NJ
	IMPLICIT NONE
!-------
100	FORMAT(E13.5)
	REAL(KIND=8),INTENT(IN)::TAUF(NI,NJ)
	INTEGER,INTENT(INOUT)::IFLD
!-------
	INTEGER::I,J,D0,D1,D2,D3
	CHARACTER DATFILE*(21)
!-------
	DATFILE(1:13)='TAUFIELDS/dat'
!-------
	IF (IFLD .LT. 10) THEN
	  D0=0
	  D1=0
	  D2=0
	  D3=IFLD
	ENDIF
!-------
	IF ((IFLD .GE. 10) .AND. (IFLD .LT. 100)) THEN
	  D0=0
	  D1=0
	  D2=INT(IFLD/10)
	  D3=IFLD-D2*10
	ENDIF
!-------
	IF ((IFLD .GE. 100) .AND. (IFLD .LT. 1000)) THEN
	  D0=0
	  D1=INT(IFLD/100)
	  D2=IFLD-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D1*100-D2*10
	ENDIF
!-------
	IF (IFLD .GE. 1000) THEN
	  D0=INT(IFLD/1000)
	  D1=IFLD-D0*1000
	  D1=INT(D1/100)
	  D2=IFLD-D0*1000-D1*100

	  D2=INT(D2/10)
	  D3=IFLD-D0*1000-D1*100-D2*10
	ENDIF
!-------
	DATFILE(14:14)=CHAR(D0+48)
	DATFILE(15:15)=CHAR(D1+48)
	DATFILE(16:16)=CHAR(D2+48)
	DATFILE(17:17)=CHAR(D3+48)
	DATFILE(18:21)='.ssv'
!-------
	OPEN (UNIT=510,FILE=TRIM(DATFILE),FORM='FORMATTED',STATUS='UNKNOWN')
!-------
	DO J=1,NJ
	DO I=1,NI
	  WRITE(510,100) TAUF(I,J)
	ENDDO
	ENDDO
!-------
	CLOSE(UNIT=510)
!-------
	END SUBROUTINE TAUFOUT
!----------------------------------------------------------------------------------
!------- COMPUTATION OF THE FORCE TERM
!----------------------------------------------------------------------------------
	SUBROUTINE COMPFRC(FRC,PHI,TSTEP)
	USE SIMPARAM,ONLY:NI,NJ,NDIM
	IMPLICIT NONE
!-------
	INTEGER,INTENT(IN)::TSTEP
	REAL(KIND=8),INTENT(IN)::PHI(NI,NJ)
	REAL(KIND=8),INTENT(OUT)::FRC(NI,NJ,NDIM)
	INTEGER::I,J,JM,JP,IM,IP
!-------
!	WRITE(*,*) TSTEP
	IF (TSTEP<2000) THEN
	  FRC(:,:,1)=0.0D0
	  FRC(:,:,2)=0.0D0
	ENDIF
!-------
	DO J=1,NJ
	DO I=1,INT(NI/2)
	  IF ((TSTEP>=2000).AND.(PHI(I,J)<1.050D0).AND.(PHI(I,J)>-0.920D0)) THEN
	    IP=I
	    JP=J
	    FRC(IP,JP,1)=0.000030D0
	    FRC(IP,JP,2)=0.0D0
	  ENDIF
	ENDDO
	ENDDO
!-------
	DO J=1,NJ
	DO I=INT(NI/2),NI
	  IF ((TSTEP>=2000).AND.(PHI(I,J)<1.050D0).AND.(PHI(I,J)>-0.920D0)) THEN
	    IP=I
	    JP=J
	    FRC(IP,JP,1)=-0.000030D0
	    FRC(IP,JP,2)=0.0D0
	  ENDIF
	ENDDO
	ENDDO
!-------
	IF (TSTEP>8500) THEN
	  FRC(:,:,1)=0.0D0
	  FRC(:,:,2)=0.0D0
	ENDIF
!	WRITE(*,*) FRC(70,70,1)
!-------
!	DO J=1,NJ
!	  WRITE(*,*) J,(FRC(I,J,2),I=1,NI)
!	ENDDO
!-------
	END SUBROUTINE COMPFRC
!-------------------------------------------------------------------------------
!------- FORCE FIELD SAVING 
!-------------------------------------------------------------------------------
	SUBROUTINE FRCOUT(IFLD,FRC)
	USE SIMPARAM,ONLY:NI,NJ,NDIM
	IMPLICIT NONE
!-------
100	FORMAT(E13.5)
	REAL(KIND=8),INTENT(IN)::FRC(NI,NJ,NDIM)
	INTEGER,INTENT(INOUT)::IFLD
!-------
	INTEGER::I,J,D0,D1,D2,D3
	CHARACTER(LEN=22) DATFILE1,DATFILE2
!-------
	DATFILE1(1:14)='FRCFIELDS/datx'
	DATFILE2(1:14)='FRCFIELDS/daty'
!-------
	IF (IFLD .LT. 10) THEN
	  D0=0
	  D1=0
	  D2=0
	  D3=IFLD
	ENDIF
!-------
	IF ((IFLD .GE. 10) .AND. (IFLD .LT. 100)) THEN
	  D0=0
	  D1=0
	  D2=INT(IFLD/10)
	  D3=IFLD-D2*10
	ENDIF
!-------
	IF ((IFLD .GE. 100) .AND. (IFLD .LT. 1000)) THEN
	  D0=0
	  D1=INT(IFLD/100)
	  D2=IFLD-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D1*100-D2*10
	ENDIF
!-------
	IF (IFLD .GE. 1000) THEN
	  D0=INT(IFLD/1000)
	  D1=IFLD-D0*1000
	  D1=INT(D1/100)
	  D2=IFLD-D0*1000-D1*100
	  D2=INT(D2/10)
	  D3=IFLD-D0*1000-D1*100-D2*10
	ENDIF
!-------
	DATFILE1(15:15)=CHAR(D0+48)
	DATFILE1(16:16)=CHAR(D1+48)
	DATFILE1(17:17)=CHAR(D2+48)
	DATFILE1(18:18)=CHAR(D3+48)
	DATFILE1(19:22)='.ssv'
!-------
	DATFILE2(15:15)=CHAR(D0+48)
	DATFILE2(16:16)=CHAR(D1+48)
	DATFILE2(17:17)=CHAR(D2+48)
	DATFILE2(18:18)=CHAR(D3+48)
	DATFILE2(19:22)='.ssv'
!-------
	OPEN (UNIT=210,FILE=DATFILE1,FORM='FORMATTED',STATUS='UNKNOWN')
	OPEN (UNIT=211,FILE=DATFILE2,FORM='FORMATTED',STATUS='UNKNOWN')
!-------
	DO J=1,NJ
	DO I=1,NI
	  WRITE(210,100) FRC(I,J,1)
	  WRITE(211,100) FRC(I,J,2)
	ENDDO
	ENDDO
!-------
	CLOSE(UNIT=210)
	CLOSE(UNIT=211)
!-------
	END SUBROUTINE FRCOUT
!----------------------------------------------------------------------------------
