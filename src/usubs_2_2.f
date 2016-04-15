C     SUBROUTINE        B  C  T  I  M  E           SUTRA VERSION 2.2     BCTIME.........100
C                                                                        BCTIME.........200
C *** PURPOSE :                                                          BCTIME.........300
C ***  USER-PROGRAMMED SUBROUTINE WHICH ALLOWS THE USER TO SPECIFY:      BCTIME.........400
C ***   (1) TIME-DEPENDENT SPECIFIED PRESSURES AND TIME-DEPENDENT        BCTIME.........500
C ***       CONCENTRATIONS OR TEMPERATURES OF INFLOWS AT THESE POINTS    BCTIME.........600
C ***   (2) TIME-DEPENDENT SPECIFIED CONCENTRATIONS OR TEMPERATURES      BCTIME.........700
C ***   (3) TIME-DEPENDENT FLUID SOURCES AND CONCENTRATIONS              BCTIME.........800
C ***       OR TEMPERATURES OF INFLOWS AT THESE POINTS                   BCTIME.........900
C ***   (4) TIME-DEPENDENT ENERGY OR SOLUTE MASS SOURCES                 BCTIME........1000
C                                                                        BCTIME........1100
      SUBROUTINE BCTIME(IPBC,PBC,IUBC,UBC,QIN,UIN,QUIN,IQSOP,IQSOU,      BCTIME........1200
     1   IPBCT,IUBCT,IQSOPT,IQSOUT,X,Y,Z,IBCPBC,IBCUBC,IBCSOP,IBCSOU,
     2   PM1,UM1,CJGNUP,CJGNUU,RCIT,SW,POR,NREG,YY,SAREA,SM)     ! Chengji 2015-03-31
      USE M_PARAMS
      USE M_ET
      USE M_TIDE
      USE M_CTRL
      USE M_AEROR
      USE M_SURFR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                BCTIME........1400
C      IMPLICIT NONE FOR BCTIME TO PREVENT TYPO
C      IMPLICIT NONE
      DIMENSION IPBC(NBCN),PBC(NBCN),IUBC(NBCN),UBC(NBCN),               BCTIME........1500
     1   QIN(NN),UIN(NN),QUIN(NN),IQSOP(NSOP),IQSOU(NSOU),               BCTIME........1600
     2   X(NN),Y(NN),Z(NN)
      DIMENSION CJGNUP(NBCN),CJGNUU(NBCN)      ! Chengji 2015-03-31
      DIMENSION PM1(NN)  ! Chengji 2015-03-31 
      DIMENSION NREG(NN),YY(NSOP),SAREA(NSOP)
      DIMENSION SW(NN),POR(NN),RCIT(NN),UM1(NN),SM(NN)
      INTEGER(1) IBCPBC(NBCN),IBCUBC(NBCN),IBCSOP(NSOP),IBCSOU(NSOU)     BCTIME........1800
      COMMON /DIMS/ NN,NE,NIN,NBI,NCBI,NB,NBHALF,NPBC,NUBC,              BCTIME........1900
     1   NSOP,NSOU,NBCN,NCIDB                                            BCTIME........2000
      COMMON /FUNITS/ K00,K0,K1,K2,K3,K4,K5,K6,K7,K8,K9,                 BCTIME........2100
     1   K10,K11,K12,K13                                                 BCTIME........2200
      COMMON /GRAVEC/ GRAVX,GRAVY,GRAVZ                                  BCTIME........2300
      COMMON /TIMES/ DELT,TSEC,TMIN,THOUR,TDAY,TWEEK,TMONTH,TYEAR,       BCTIME........2400
     1   TMAX,DELTP,DELTU,DLTPM1,DLTUM1,IT,ITBCS,ITRST,ITMAX,TSTART      BCTIME........2500
      COMMON /CONTRL/ GNUP,GNUU,UP,DTMULT,DTMAX,ME,ISSFLO,ISSTRA,ITCYC, 
     1   NPCYC,NUCYC,NPRINT,NBCFPR,NBCSPR,NBCPPR,NBCUPR,IREAD,          
     2   ISTORE,NOUMAT,IUNSAT,KTYPE                                     
      DIMENSION KTYPE(2)                                                
      DOUBLE PRECISION TIDE
      REAL SEEPX,SEEPY ! The x- and y-coordinate of seepage node
      
      SEEPX=0
      SEEPY=0
C      DIMENSION TIDE
C                                                                        BCTIME........2600
C.....DEFINITION OF REQUIRED VARIABLES                                   BCTIME........2700
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........2800
C     NN = EXACT NUMBER OF NODES IN MESH                                 BCTIME........2900
C     NPBC = EXACT NUMBER OF SPECIFIED PRESSURE NODES                    BCTIME........3000
C     NUBC = EXACT NUMBER OF SPECIFIED CONCENTRATION                     BCTIME........3100
C            OR TEMPERATURE NODES                                        BCTIME........3200
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........3300
C     IT = NUMBER OF CURRENT TIME STEP                                   BCTIME........3400
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........3500
C     TSEC = TIME AT END OF CURRENT TIME STEP IN SECONDS                 BCTIME........3600
C     TMIN = TIME AT END OF CURRENT TIME STEP IN MINUTES                 BCTIME........3700
C     THOUR = TIME AT END OF CURRENT TIME STEP IN HOURS                  BCTIME........3800
C     TDAY = TIME AT END OF CURRENT TIME STEP IN DAYS                    BCTIME........3900
C     TWEEK = TIME AT END OF CURRENT TIME STEP IN WEEKS                  BCTIME........4000
C     TMONTH = TIME AT END OF CURRENT TIME STEP IN MONTHS                BCTIME........4100
C     TYEAR = TIME AT END OF CURRENT TIME STEP IN YEARS                  BCTIME........4200
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........4300
C     PBC(IP) = SPECIFIED PRESSURE VALUE AT IP(TH) SPECIFIED             BCTIME........4400
C               PRESSURE NODE                                            BCTIME........4500
C     UBC(IP) = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE OF ANY      BCTIME........4600
C               INFLOW OCCURRING AT IP(TH) SPECIFIED PRESSURE NODE       BCTIME........4700
C     IPBC(IP) = ACTUAL NODE NUMBER OF IP(TH) SPECIFIED PRESSURE NODE    BCTIME........4800
C                {WHEN NODE NUMBER I=IPBC(IP) IS NEGATIVE (I<0),         BCTIME........4900
C                VALUES MUST BE SPECIFIED FOR PBC AND UBC.}              BCTIME........5000
C     IBCPBC(IP) = INDICATOR OF WHERE THIS PRESSURE SPECIFICATION        BCTIME........5100
C                  WAS MADE. MUST BE SET TO -1 TO INDICATE THAT THIS     BCTIME........5200
C                  SPECIFICATION WAS MADE IN SUBROUTINE BCTIME.          BCTIME........5300
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........5400
C     UBC(IUP) = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE AT         BCTIME........5500
C                IU(TH) SPECIFIED CONCENTRATION OR TEMPERATURE NODE      BCTIME........5600
C                (WHERE IUP=IU+NPBC)                                     BCTIME........5700
C     IUBC(IUP) = ACTUAL NODE NUMBER OF IU(TH) SPECIFIED CONCENTRATION   BCTIME........5800
C                 OR TEMPERATURE NODE (WHERE IUP=IU+NPBC)                BCTIME........5900
C                 {WHEN NODE NUMBER I=IUBC(IU) IS NEGATIVE (I<0),        BCTIME........6000
C                 A VALUE MUST BE SPECIFIED FOR UBC.}                    BCTIME........6100
C     IBCUBC(IUP) = INDICATOR OF WHERE THIS CONCENTRATION OR TEMPERATURE BCTIME........6200
C                  SPECIFICATION WAS MADE. MUST BE SET TO -1 TO INDICATE BCTIME........6300
C                  THAT THIS SPECIFICATION WAS MADE IN SUBROUTINE        BCTIME........6400
C                  BCTIME.                                               BCTIME........6500
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........6600
C     IQSOP(IQP) = NODE NUMBER OF IQP(TH) FLUID SOURCE NODE.             BCTIME........6700
C                  {WHEN NODE NUMBER I=IQSOP(IQP) IS NEGATIVE (I<0),     BCTIME........6800
C                  VALUES MUST BE SPECIFIED FOR QIN AND UIN.}            BCTIME........6900
C     QIN(-I) = SPECIFIED FLUID SOURCE VALUE AT NODE (-I)                BCTIME........7000
C     UIN(-I) = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE OF ANY      BCTIME........7100
C               INFLOW OCCURRING AT FLUID SOURCE NODE (-I)               BCTIME........7200
C     IBCSOP(IQP) = INDICATOR OF WHERE THIS FLUID SOURCE SPECIFICATION   BCTIME........7300
C                   WAS MADE. MUST BE SET TO -1 TO INDICATE THAT THIS    BCTIME........7400
C                   SPECIFICATION WAS MADE IN SUBROUTINE BCTIME.         BCTIME........7500
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........7600
C     IQSOU(IQU) = NODE NUMBER OF IQU(TH) ENERGY OR                      BCTIME........7700
C                  SOLUTE MASS SOURCE NODE                               BCTIME........7800
C                  {WHEN NODE NUMBER I=IQSOU(IQU) IS NEGATIVE (I<0),     BCTIME........7900
C                  A VALUE MUST BE SPECIFIED FOR QUIN.}                  BCTIME........8000
C     QUIN(-I) = SPECIFIED ENERGY OR SOLUTE MASS SOURCE VALUE            BCTIME........8100
C                AT NODE (-I)                                            BCTIME........8200
C     IBCSOU(IQU) = INDICATOR OF WHERE THIS ENERGY OR SOLUTE MASS        BCTIME........8300
C                   SOURCE SPECIFICATION WAS MADE. MUST BE SET TO -1     BCTIME........8400
C                   TO INDICATE THAT THIS SPECIFICATION WAS MADE IN      BCTIME........8500
C                   SUBROUTINE BCTIME.                                   BCTIME........8600
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........8700
C                                                                        BCTIME........8800
C.....ADDITIONAL USEFUL VARIABLES                                        BCTIME........8900
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........9000
C     "FUNITS" ARE UNIT NUMBERS FOR INPUT AND OUTPUT FILES               BCTIME........9100
C         AS ASSIGNED IN THE INPUT FILE "SUTRA.FIL"                      BCTIME........9200
C                                                                        BCTIME........9300
C     X(I), Y(I), AND Z(I) ARE THE X-, Y-, AND Z-COORDINATES OF NODE I   BCTIME........9400
C     (FOR 2-D PROBLEMS, Z(I) IS THE MESH THICKNESS AT NODE I)           BCTIME........9500
C                                                                        BCTIME........9600
C     GRAVX, GRAVY AND GRAVZ ARE THE X-, Y-, AND Z-COMPONENTS OF THE     BCTIME........9700
C     GRAVITY VECTOR (FOR 2-D PROBLEMS, GRAVZ = 0)                       BCTIME........9800
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  BCTIME........9900
C                                                                        BCTIME.......10000
C                                                                        BCTIME.......10100
C.....NSOPI IS ACTUAL NUMBER OF FLUID SOURCE NODES                       BCTIME.......10200
      NSOPI=NSOP-1                                                       BCTIME.......10300
C.....NSOUI IS ACTUAL NUMBER OF ENERGY OR SOLUTE MASS SOURCE NODES       BCTIME.......10400
      NSOUI=NSOU-1                                                       BCTIME.......10500
C                                                                        BCTIME.......10600
C                                                                        BCTIME.......10700
C                                                                        BCTIME.......10800
C                                                                        BCTIME.......10900
C                                                                        BCTIME.......11000
C                                                                        BCTIME.......11100
C     ONE LINE TO AVOID TIDE AFFECTING EVAPORATION ONLY CASES
      TIDE=-1.D10
      IF(IPBCT) 50,240,240                                               BCTIME.......11200
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......11300
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......11400
C.....SECTION (1):  SET TIME-DEPENDENT SPECIFIED PRESSURES OR            BCTIME.......11500
C     CONCENTRATIONS (TEMPERATURES) OF INFLOWS AT SPECIFIED              BCTIME.......11600
C     PRESSURE NODES                                                     BCTIME.......11700
C                                                                        BCTIME.......11800
   50 CONTINUE                                                           BCTIME.......11900
C TA  --
C#  DATASET 13B: TIDE FLUCTUATION IN USUBS
C   TASP   -- SPRING TIDAL AMPLITUDE(M); 
C   TANE   -- NEAP TIDAL AMPLITUDE (M);
C   TPSP   -- TIDAL PERIOD OF SPRING TIDE(S); 
C   TPNP   -- TIDAL PERIOD OF NEAP TIDE(S); 
C   TM    -- MEAN TIDAL LEVEL[M]; 
C   RHOST -- THE DENSITY FOR TIDE WATER; 
C   SC    -- SALINITY OF THE SEAWATER; 
C   ITT   -- ITERRATION CRITERIA FOR BCTIME SHOULD BE LARGE OR EQUAL THAN 2
C         [TA]           [TP]           [TM]        [RHOST]        [SC]    [ITT]
C            0.0       4.32D+4          0.D+0        1025.0D+0    3.57D-02      2
C      WRITE(*,*) QET,UET,PET,UVM,NGT,ITE
C      WRITE(*,*) TASP,TANE,TPSP,TPNP,TM,RHOST,SC
C      TIDE=4.6+1.5*SIN(2*3.1415926*TSC/12/60/60)+0.5*SIN(2*3.1415926*IT*60/12.42/60/60)
C      TIDE=TM+TASP*SIN(2.D0*PI*IT*60/12/60/60)+TANP*SIN(2*PI*IT*60/12.42/60/60)
      TIDE=TM+TASP*SIN(2.D0*PI*TSEC/TPSP)+TANE*SIN(2.D0*PI*TSEC/TPNP)
C      TIDE=4.2D0+1.0D0*SIN(TSEC*3.1415926D0/360.D0/60.D0)       ! Chengji 2015-03-31
      IF (IT.EQ.1) THEN  
          OPEN(21,FILE='TIDE.DAT',STATUS='UNKNOWN')   
          WRITE(21,98)
   98     FORMAT('  IT',4X,'TIME(DAY)',3X,'TIDAL LEVEL (M)')
      ENDIF
      WRITE(21,99) IT, TSEC/3600./24.,TIDE
   99 FORMAT(I15,(1PE10.2,2X),500(1PE10.3,1X))
      
      DO 200 IP=1,NPBC                                                   BCTIME.......12000
      I=IPBC(IP)                                                         BCTIME.......12100
      IF(I) 100,200,200                                                  BCTIME.......12200
  100 CONTINUE                                                           BCTIME.......12300
C     NOTE: A FLOW AND TRANSPORT SOLUTION MUST OCCUR FOR ANY             BCTIME.......12400
C           TIME STEP IN WHICH PBC( ) CHANGES.                           BCTIME.......12500
C     PBC(IP) =  ((          ))                                          BCTIME.......12600
C     UBC(IP) =  ((          ))                                          BCTIME.......12700
C******************************************************************************************
      PBC(IP)=9.8D0*(1000.D0+SC*713.D0)*(TIDE-Y(IABS(I)))
  
      IF(Y(IABS(I)).LE.TIDE) THEN
          UBC(IP)=SC
      ELSE
          UBC(IP)=0.D0   
      ENDIF
      
      IF(PM1(IABS(I)).GT.0.AND.Y(IABS(I)).GT.TIDE) THEN
          PBC(IP)=0.D0
          CJGNUP(IP)=GNUP
          IF(Y(IABS(I))+PM1(IABS(I))/10245.GT.SEEPY) THEN
             SEEPY=Y(IABS(I))+PM1(IABS(I))/10245
             SEEPX=X(IABS(I))
          ENDIF
      ELSEIF(PM1(IABS(I)).GT.0.AND.Y(IABS(I)).LE.TIDE) THEN
          CJGNUP(IP)=GNUP
      ELSEIF(PM1(IABS(I)).LT.0.AND.Y(IABS(I)).GT.TIDE) THEN
          CJGNUP(IP)=0.D0
      ELSEIF(PM1(IABS(I)).LT.0.AND.Y(IABS(I)).LE.TIDE) THEN
          CJGNUP(IP)=GNUP
      ENDIF
C******************************************************************************************	  
C.....IBCPBC(IP) MUST BE SET TO -1 TO INDICATE THAT PBC(IP)              BCTIME.......12800
C        AND/OR UBC(IP) HAVE BEEN SET BY SUBROUTINE BCTIME.              BCTIME.......12900
      IBCPBC(IP) = -1                                                    BCTIME.......13000
  200 CONTINUE                                                           BCTIME.......13100
      WRITE(20094,'(I10,3E16.7)') IT,TIDE,SEEPX,SEEPY
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......13200
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......13300
C                                                                        BCTIME.......13400
C                                                                        BCTIME.......13500
C                                                                        BCTIME.......13600
C                                                                        BCTIME.......13700
C                                                                        BCTIME.......13800
C                                                                        BCTIME.......13900
  240 IF(IUBCT) 250,440,440                                              BCTIME.......14000
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......14100
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......14200
C.....SECTION (2):  SET TIME-DEPENDENT SPECIFIED                         BCTIME.......14300
C     CONCENTRATIONS (TEMPERATURES)                                      BCTIME.......14400
C                                                                        BCTIME.......14500
  250 CONTINUE                                                           BCTIME.......14600
      DO 400 IU=1,NUBC                                                   BCTIME.......14700
      IUP=IU+NPBC                                                        BCTIME.......14800
      I=IUBC(IUP)                                                        BCTIME.......14900
      IF(I) 300,400,400                                                  BCTIME.......15000
  300 CONTINUE                                                           BCTIME.......15100
C     NOTE: A TRANSPORT SOLUTION MUST OCCUR FOR ANY TIME STEP IN WHICH   BCTIME.......15200
C           UBC( ) CHANGES.  IN ADDITION, IF FLUID PROPERTIES ARE        BCTIME.......15300
C           SENSITIVE TO 'U', THEN A FLOW SOLUTION MUST OCCUR AS WELL.   BCTIME.......15400
C     UBC(IUP) =   ((          ))                                        BCTIME.......15500
      UBC(IUP) = SC

      IF(Y(IABS(I)).GT.TIDE) THEN
          CJGNUU(IUP)=0.D0
      ELSEIF(Y(IABS(I)).LE.TIDE) THEN
          CJGNUU(IUP)=GNUU
      ENDIF
C.....IBCUBC(IUP) MUST BE SET TO -1 TO INDICATE THAT UBC(IUP)            BCTIME.......15600
C        HAS BEEN SET BY SUBROUTINE BCTIME.                              BCTIME.......15700
      IBCUBC(IUP) = -1                                                   BCTIME.......15800
  400 CONTINUE                                                           BCTIME.......15900
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......16000
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......16100
C                                                                        BCTIME.......16200
C                                                                        BCTIME.......16300
C                                                                        BCTIME.......16400
C                                                                        BCTIME.......16500
C                                                                        BCTIME.......16600
C                                                                        BCTIME.......16700
  440 IF(IQSOPT) 450,640,640                                             BCTIME.......16800
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......16900
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......17000
C.....SECTION (3):  SET TIME-DEPENDENT FLUID SOURCES/SINKS,              BCTIME.......17100
C      OR CONCENTRATIONS (TEMPERATURES) OF SOURCE FLUID                  BCTIME.......17200
C                                                                        BCTIME.......17300
  450 CONTINUE                                                           BCTIME.......17400
      DO 600 IQP=1,NSOPI                                                 BCTIME.......17500
      I=IQSOP(IQP)                                                       BCTIME.......17600
      IF(I) 500,600,600                                                  BCTIME.......17700
  500 CONTINUE                                                           BCTIME.......17800
C      IF (IT.EQ.7194.AND.IABS(I).EQ.145) THEN
C        aabss=1.0D0
C      ENDIF
      IF(PM1(IABS(I)).LT.PET.AND.Y(IABS(I)).GE.TIDE) THEN
      CALL EVAPORATION (AET,PM1(IABS(I)),UM1(IABS(I)),RCIT(IABS(I))
     2,POR(IABS(I)),SW(IABS(I))
     3,NREG(IABS(I)),YY(IQP),SM(IABS(I))/SAREA(IQP))
C      0.004 M/DAY /3600/24  DAY/S *2M *1M * 1000 KG/M3    =[KG/S]
C     QET (M/S) * 2 (M) *1 (M) * 1000 KG/M3 = [KG/S] 
!            QIN(-I)=-QET*2.D0*1.D0*1.D3 
C     AET ENDS UP WITH A NEGATIVE VALUE!!!
            QIN(-I)=AET*SAREA(IQP)*1.D3 
            UIN(-I)=UET
      ELSE
            QIN(-I)=0.D0
            UIN(-I)=0.D0
      END IF
C     just wish to make one more commit
C     NOTE: A FLOW AND TRANSPORT SOLUTION MUST OCCUR FOR ANY             BCTIME.......17900
C           TIME STEP IN WHICH QIN( ) CHANGES.                           BCTIME.......18000
C     QIN(-I) =   ((           ))                                        BCTIME.......18100
C     NOTE: A TRANSPORT SOLUTION MUST OCCUR FOR ANY                      BCTIME.......18200
C           TIME STEP IN WHICH UIN( ) CHANGES.                           BCTIME.......18300
C     UIN(-I) =   ((           ))                                        BCTIME.......18400
C.....IBCSOP(IQP) MUST BE SET TO -1 TO INDICATE THAT QIN(-I)             BCTIME.......18500
C        AND/OR UIN(-I) HAVE BEEN SET BY SUBROUTINE BCTIME.              BCTIME.......18600
      IBCSOP(IQP) = -1                                                   BCTIME.......18700
  600 CONTINUE                                                           BCTIME.......18800
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......18900
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......19000
C                                                                        BCTIME.......19100
C                                                                        BCTIME.......19200
C                                                                        BCTIME.......19300
C                                                                        BCTIME.......19400
C                                                                        BCTIME.......19500
C                                                                        BCTIME.......19600
  640 IF(IQSOUT) 650,840,840                                             BCTIME.......19700
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......19800
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......19900
C.....SECTION (4):  SET TIME-DEPENDENT SOURCES/SINKS                     BCTIME.......20000
C     OF SOLUTE MASS OR ENERGY                                           BCTIME.......20100
C                                                                        BCTIME.......20200
  650 CONTINUE                                                           BCTIME.......20300
      DO 800 IQU=1,NSOUI                                                 BCTIME.......20400
      I=IQSOU(IQU)                                                       BCTIME.......20500
      IF(I) 700,800,800                                                  BCTIME.......20600
  700 CONTINUE                                                           BCTIME.......20700
C     NOTE: A TRANSPORT SOLUTION MUST OCCUR FOR ANY                      BCTIME.......20800
C           TIME STEP IN WHICH QUIN( ) CHANGES.                          BCTIME.......20900
C     QUIN(-I) =   ((           ))                                       BCTIME.......21000
C.....IBCSOU(IQU) MUST BE SET TO -1 TO INDICATE THAT QUIN(-I)            BCTIME.......21100
C        HAS BEEN SET BY SUBROUTINE BCTIME.                              BCTIME.......21200
      IBCSOU(IQU) = -1                                                   BCTIME.......21300
  800 CONTINUE                                                           BCTIME.......21400
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......21500
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  BCTIME.......21600
C                                                                        BCTIME.......21700
C                                                                        BCTIME.......21800
C                                                                        BCTIME.......21900
C                                                                        BCTIME.......22000
C                                                                        BCTIME.......22100
C                                                                        BCTIME.......22200
  840 CONTINUE                                                           BCTIME.......22300
C                                                                        BCTIME.......22400
      RETURN                                                             BCTIME.......22500
      END                                                                BCTIME.......22600
CC     SUBROUTINE        U  N  S  A  T              SUTRA VERSION 2.2     UNSAT..........100
CC                                                                        UNSAT..........200
CC *** PURPOSE :                                                          UNSAT..........300
CC ***  USER-PROGRAMMED SUBROUTINE GIVING:                                UNSAT..........400
CC ***  (1)  SATURATION AS A FUNCTION OF PRESSURE ( SW(PRES) )            UNSAT..........500
CC ***  (2)  DERIVATIVE OF SATURATION WITH RESPECT TO PRESSURE            UNSAT..........600
CC ***       AS A FUNCTION OF EITHER PRESSURE OR SATURATION               UNSAT..........700
CC ***       ( DSWDP(PRES), OR DSWDP(SW) )                                UNSAT..........800
CC ***  (3)  RELATIVE PERMEABILITY AS A FUNCTION OF EITHER                UNSAT..........900
CC ***       PRESSURE OR SATURATION ( REL(PRES) OR RELK(SW) )             UNSAT.........1000
CC ***                                                                    UNSAT.........1100
CC ***  CODE BETWEEN DASHED LINES MUST BE REPLACED TO GIVE THE            UNSAT.........1200
CC ***  PARTICULAR UNSATURATED RELATIONSHIPS DESIRED.                     UNSAT.........1300
CC ***                                                                    UNSAT.........1400
CC ***  DIFFERENT FUNCTIONS MAY BE GIVEN FOR EACH REGION OF THE MESH.     UNSAT.........1500
CC ***  REGIONS ARE SPECIFIED BY BOTH NODE NUMBER AND ELEMENT NUMBER      UNSAT.........1600
CC ***  IN INPUT DATA FILE FOR UNIT K1 (INP).                             UNSAT.........1700
CC                                                                        UNSAT.........1800
C      SUBROUTINE UNSAT(SW,DSWDP,RELK,PRES,KREG)                          UNSAT.........1900
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                UNSAT.........2000
C      DIMENSION KTYPE(2)                                                 UNSAT.........2100
C      COMMON /CONTRL/ GNUP,GNUU,UP,DTMULT,DTMAX,ME,ISSFLO,ISSTRA,ITCYC,  UNSAT.........2200
C     1   NPCYC,NUCYC,NPRINT,NBCFPR,NBCSPR,NBCPPR,NBCUPR,IREAD,           UNSAT.........2300
C     2   ISTORE,NOUMAT,IUNSAT,KTYPE                                      UNSAT.........2400
CC                                                                        UNSAT.........2500
CC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- UNSAT.........2600
CC     E X A M P L E   C O D I N G   FOR                                  UNSAT.........2700
CC     MESH WITH TWO REGIONS OF UNSATURATED PROPERTIES USING              UNSAT.........2800
CC     THREE PARAMETER-UNSATURATED FLOW RELATIONSHIPS OF                  UNSAT.........2900
CC     VAN GENUCHTEN(1980)                                                UNSAT.........3000
CC        RESIDUAL SATURATION, SWRES, GIVEN IN UNITS {L**0}               UNSAT.........3100
CC        PARAMETER, AA, GIVEN IN INVERSE PRESSURE UNITS {m*(s**2)/kg}    UNSAT.........3200
CC        PARAMETER, VN, GIVEN IN UNITS {L**0}                            UNSAT.........3300
CC                                                                        UNSAT.........3400
C      REAL SWRES,AA,VN,SWRM1,AAPVN,VNF,AAPVNN,DNUM,DNOM,SWSTAR           UNSAT.........3500
C      REAL SWRES1,SWRES2,AA1,AA2,VN1,VN2                                 UNSAT.........3600
CC                                                                        UNSAT.........3700
CC     DATA FOR REGION 1:                                                 UNSAT.........3800
C      DATA   SWRES1/0.30E0/,   AA1/5.0E-5/,   VN1/2.0E0/                 UNSAT.........3900
C      SAVE SWRES1, AA1, VN1                                              UNSAT.........4000
CC     DATA FOR REGION 2:                                                 UNSAT.........4100
C      DATA   SWRES2/0.30E0/,   AA2/5.0E-5/,   VN2/2.0E0/                 UNSAT.........4200
C      SAVE SWRES2, AA2, VN2                                              UNSAT.........4300
CC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- UNSAT.........4400
CC                                                                        UNSAT.........4500
CC *** BECAUSE THIS ROUTINE IS CALLED OFTEN FOR UNSATURATED FLOW RUNS,    UNSAT.........4600
CC *** EXECUTION TIME MAY BE SAVED BY CAREFUL CODING OF DESIRED           UNSAT.........4700
CC *** RELATIONSHIPS USING ONLY INTEGER AND SINGLE PRECISION VARIABLES!   UNSAT.........4800
CC *** RESULTS OF THE CALCULATIONS MUST THEN BE PLACED INTO DOUBLE        UNSAT.........4900
CC *** PRECISION VARIABLES SW, DSWDP AND RELK BEFORE LEAVING              UNSAT.........5000
CC *** THIS SUBROUTINE.                                                   UNSAT.........5100
CC                                                                        UNSAT.........5200
CC                                                                        UNSAT.........5300
CC*********************************************************************** UNSAT.........5400
CC*********************************************************************** UNSAT.........5500
CC                                                                        UNSAT.........5600
CC     SET PARAMETERS FOR CURRENT REGION, KREG                            UNSAT.........5700
C      GOTO(10,20),KREG                                                   UNSAT.........5800
C   10 SWRES=SWRES1                                                       UNSAT.........5900
C      AA=AA1                                                             UNSAT.........6000
C      VN=VN1                                                             UNSAT.........6100
C      GOTO 100                                                           UNSAT.........6200
C   20 SWRES=SWRES2                                                       UNSAT.........6300
C      AA=AA2                                                             UNSAT.........6400
C      VN=VN2                                                             UNSAT.........6500
C  100 CONTINUE                                                           UNSAT.........6600
CC                                                                        UNSAT.........6700
CC                                                                        UNSAT.........6800
CC*********************************************************************** UNSAT.........6900
CC*********************************************************************** UNSAT.........7000
CC.....SECTION (1):                                                       UNSAT.........7100
CC     SW VS. PRES   (VALUE CALCULATED ON EACH CALL TO UNSAT)             UNSAT.........7200
CC     CODING MUST GIVE A VALUE TO SATURATION, SW.                        UNSAT.........7300
CC                                                                        UNSAT.........7400
CC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  UNSAT.........7500
CC     THREE PARAMETER MODEL OF VAN GENUCHTEN(1980)                       UNSAT.........7600
C      SWRM1=1.E0-SWRES                                                   UNSAT.........7700
C      AAPVN=1.E0+(AA*(-PRES))**VN                                        UNSAT.........7800
C      VNF=(VN-1.E0)/VN                                                   UNSAT.........7900
C      AAPVNN=AAPVN**VNF                                                  UNSAT.........8000
C      S W   =   DBLE (SWRES+SWRM1/AAPVNN)                                UNSAT.........8100
CC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  UNSAT.........8200
CC*********************************************************************** UNSAT.........8300
CC*********************************************************************** UNSAT.........8400
CC                                                                        UNSAT.........8500
CC                                                                        UNSAT.........8600
CC                                                                        UNSAT.........8700
CC                                                                        UNSAT.........8800
CC                                                                        UNSAT.........8900
CC                                                                        UNSAT.........9000
C      IF(IUNSAT-2) 600,1200,1800                                         UNSAT.........9100
CC*********************************************************************** UNSAT.........9200
CC*********************************************************************** UNSAT.........9300
CC.....SECTION (2):                                                       UNSAT.........9400
CC     DSWDP VS. PRES, OR DSWDP VS. SW   (CALCULATED ONLY WHEN IUNSAT=1)  UNSAT.........9500
CC     CODING MUST GIVE A VALUE TO DERIVATIVE OF SATURATION WITH          UNSAT.........9600
CC     RESPECT TO PRESSURE, DSWDP.                                        UNSAT.........9700
CC                                                                        UNSAT.........9800
C  600 CONTINUE                                                           UNSAT.........9900
CC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  UNSAT........10000
C      DNUM=AA*(VN-1.E0)*SWRM1*(AA*(-PRES))**(VN-1.E0)                    UNSAT........10100
C      DNOM=AAPVN*AAPVNN                                                  UNSAT........10200
C      D S W D P   =   DBLE (DNUM/DNOM)                                   UNSAT........10300
CC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  UNSAT........10400
C      GOTO 1800                                                          UNSAT........10500
CC*********************************************************************** UNSAT........10600
CC*********************************************************************** UNSAT........10700
CC                                                                        UNSAT........10800
CC                                                                        UNSAT........10900
CC                                                                        UNSAT........11000
CC                                                                        UNSAT........11100
CC                                                                        UNSAT........11200
CC                                                                        UNSAT........11300
CC*********************************************************************** UNSAT........11400
CC*********************************************************************** UNSAT........11500
CC.....SECTION (3):                                                       UNSAT........11600
CC     RELK VS. P, OR RELK VS. SW   (CALCULATED ONLY WHEN IUNSAT=2)       UNSAT........11700
CC     CODING MUST GIVE A VALUE TO RELATIVE PERMEABILITY, RELK.           UNSAT........11800
CC                                                                        UNSAT........11900
C 1200 CONTINUE                                                           UNSAT........12000
CC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  UNSAT........12100
CC     GENERAL RELATIVE PERMEABILITY MODEL FROM VAN GENUCHTEN(1980)       UNSAT........12200
C      SWSTAR=(SW-SWRES)/SWRM1                                            UNSAT........12300
C      R E L K   =   DBLE (SQRT(SWSTAR)*                                  UNSAT........12400
C     1                   (1.E0-(1.E0-SWSTAR**(1.E0/VNF))**(VNF))**2)     UNSAT........12500
CC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  UNSAT........12600
CC                                                                        UNSAT........12700
CC*********************************************************************** UNSAT........12800
CC*********************************************************************** UNSAT........12900
CC                                                                        UNSAT........13000
CC                                                                        UNSAT........13100
CC                                                                        UNSAT........13200
CC                                                                        UNSAT........13300
CC                                                                        UNSAT........13400
CC                                                                        UNSAT........13500
C 1800 RETURN                                                             UNSAT........13600
CC                                                                        UNSAT........13700
C      END                                                                UNSAT........13800
