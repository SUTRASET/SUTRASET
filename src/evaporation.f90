! *** ASSIGNING EVAPORATION RATE THROUGH RELEVANT FUNCTIONS
! *** WHEN MET=1 THEN USING PENMAN (1948) EQUATION,
! *** THIS IS FURTHER MODIFIED BY KONUKCU (2007)
! *** WHEN MET=2 THEN USING PDV (1957) EQUATION.
! *** WHEN MET=3 CALCULATES ET BY A POTENTIAL ET
! *** SEE NOTEBOOK2 PAGE47 FOR FURTHER EXPLANATION 
SUBROUTINE EVAPORATION (AET,PC,CC,RHO,POR,SW,KREG,YY,DS,TSEC)
  USE M_PARAMS
  USE M_ET
  USE M_AEROR
  USE M_SURFR
  USE M_CTRL
  USE M_SALTR
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  COMMON /PARAMS/ COMPFL,COMPMA,DRWDU,CW,CS,RHOS,SIGMAW,SIGMAS,&
     RHOW0,URHOW0,VISC0,PRODF1,PRODS1,PRODF0,PRODS0,CHI1,CHI2 
!     DS   -- THICKNESS/DEPTH OF SALT MASS [M]
!     AET  -- ACTUAL EVAPORATION RATE [MM/DAY]
!     PC   -- PRESSURE AT THE SPECIFIED NODE [KG/M/S2]
!     CC   -- CONCENTRATION AT THE SPECIFIED NODE  [KG/KG]
!     RHO  -- DENSITY OF THE SPECIFIED NODE [KG/M2]
!     TMA  -- AIR TEMPERATURE [CELSIUS]
!     TMI  -- SOIL SURFACE TEMPERATURE [CELSIUS]
!     TM   -- MEAN DAILY AIR TEMPERATURE [CELSIUS]
!     DLT  -- SLOPE OF SATURATION VAPOR PRESSURE CURVEE/T[KPA/CELSIUS]
!     VLH  -- LATENT HEAT OF VAPORISATION [J/KG] 
!     ALF  -- ALBEDO [I]
!     RS   -- SHORT WAVE INCOMING RADIATION  [MJ/M2/DAY]
!     RH   -- RELATIVE HUMIDITY 0<RH<1
!     ED   -- SATURATED VAPOUR PRESSURE AT DEW POINT [KPA]
!     EMI  -- EMMISIVITY [I]
!     RN   -- NET RADIANT ENERGY 
!     GAM  -- PSYCHROMETRIC CONSTANT [KPA/CELSIUS]
!     EO(NOT ZERO)   -- SATURATED VAPOUR PRESSURE AT MEAN DAILY AIR
!                        TEMP. (TM) [PA]
!     FU   -- WIND FUNCTION [MJ/M2/KPA/DAY]
!     AP,BP-- PARAMETER IN WIND FUNCTION
!     U2   -- DAILY MEAN WIND SPEED AT 2M ABOVE GROUND [KM/DAY]
!     RC   -- UNIVERSAL GAS CONSTANT [J/MOL/K]
!     WMW   -- MOLECULAR WEIGHT OF WATER [KG/MOL]
!     TS   -- SOIL TEMPERATURE [CELSIUS]
!     DP   -- VAN'T HOFF DISSOCIATION FACTOR
!     STM  -- MOLECULAR WEIGHT OF NACL [KG/MOL]
!     FIO(NOT ZERO)  -- OSMOTIC POTENTIALS (M)
!     HO(NOT ZERO)   -- OSMOTIC PARAMETER
!     RHOW0   -- INITIAL WATER DENSITY (KG/M3)
!      TSK   -- TEMPERATURE AT THAT NODE
!      TMAK   -- TEMPERATURE AT THAT NODE  
! ... PEMAN METHOD ...      
      IF (MET.EQ.1) THEN ! USING PENMAN EVAPORATION EQUATION
      ELSEIF (MET.EQ.2) THEN                         
!       AIR TEMPERATURE TMAK=(TMA-273.15D0)

        TMAK  = 273.15D0+TMA
        TSK   = 273.15D0+TMI
        TERM1 = WMW/RHOW0/RC
!       RELATIVE HUMIDITY INDUCED BY OSMOTICAL POTENTIAL 
        HO = EXP(-WMW*2.D0*CHI(CC)*CC/STM)
!       RELATIVE HUMIDITY INDUCED BY MATRIC POTENTIAL
        HM = EXP(WMW*PC/RHOW0/RC/TSK)
!------- DFR, PWFR AND PCFR IS USED FOR OUTPUT IN .BCOF------
!        DFR=RSC                                       
!        PWFR=HO
!        PCFR=HM
!----------END PRINT------------------------
        SURF = SURFRSIS(PC,POR,SW,MSR,KREG,TSK,YY)
        RSC  = SALTRSIS(MSC,UVM,CC,DS)
        AET  = -TERM1*(PSAT(1,TSK)*HO*HM/TSK-PSAT(1,TMAK)*RH/TMAK)/  &
            (RAVT+SURF+RSC)
        IF (PRDDN.NE.0) AET=AET*DMAX1(DSIN(2.D0*PI/TSEC/PRDDN+OFSTDN),0.D0)
! WRITE (48,"( 7(1PE15.7))") SURF,SW,PC,TSK,HO,HM,AET
! SOME COMMENTS ON FORMAT IN FORTRAN 90 
! 7(1PE15.7) HERE 7 MEANS OUTPUT 7 TIMES
!                1P MEANS THE FIRST VALUE IS NON-ZERO 
!                   E.G., 1.3D1 NOT 0.13D2
!              EW.D EXPONENTIAL, TOTAL WIDTH, AND DECIMAL WIDTH
!                   IT IS ALSO ACCEPTABLE TO USE NUMBER TO REFERENCE THE
!                      FORMAT STYLE  
!1 FORMAT( 7(1PE20.7))
      ELSEIF (MET.EQ.3) THEN
!     CALCULATS ACTUAL EVAPORATION BASED ON A POTENTIAL ET      
!      SURF  =  SURFRSIS(PC,POR,SW,MSR,KREG,TSK,YY)
!      RSC   =  SALTRSIS(MSC,UVM,CC,DS)
!      AET   =  -QET*RAVT/ (RAVT+SURF+RSC)
!!!  ------------END OUTPUT---------------
!      ELSEIF (MET.EQ.4) THEN
      AET  = -QET
      ENDIF !MET
      RETURN
      END SUBROUTINE EVAPORATION
