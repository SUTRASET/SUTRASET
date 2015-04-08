C     MODULE            A  L  L  A  R  R           SUTRA VERSION 2.2     ALLARR.........100
C                                                                        ALLARR.........200
C *** PURPOSE :                                                          ALLARR.........300
C ***  TO DECLARE THE MAIN ALLOCATABLE ARRAYS.                           ALLARR.........400
C                                                                        ALLARR.........500
      MODULE ALLARR                                                      ALLARR.........600
      IMPLICIT NONE                                                      ALLARR.........700
      LOGICAL ALLO1, ALLO2, ALLO3                                        ALLARR.........800
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::                   ALLARR.........900
     1   PMAT, UMAT                                                      ALLARR........1000
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::                     ALLARR........1100
     1   PITER, UITER, PM1, DPDTITR, UM1, UM2, PVEL, SL, SR, X, Y, Z,    ALLARR........1200
     2   VOL, POR, CS1, CS2, CS3, SW, DSWDP, RHO, SOP, QIN, UIN, QUIN,   ALLARR........1300
     3   QINITR, RCIT, RCITM1, GNUP1, GNUU1                              ALLARR........1400
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::                     ALLARR........1500
     1   PVEC, UVEC                                                      ALLARR........1600
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::                     ALLARR........1700
     1   ALMAX, ALMIN, ATMAX, ATMIN, VMAG, VANG1,                        ALLARR........1800
     2   PERMXX, PERMXY, PERMYX, PERMYY, PANGL1                          ALLARR........1900
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::                     ALLARR........2000
     1   ALMID, ATMID, VANG2, PERMXZ, PERMYZ,                            ALLARR........2100
     2   PERMZX, PERMZY, PERMZZ, PANGL2, PANGL3                          ALLARR........2200
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::                     ALLARR........2300
     1   PBC, UBC, QPLITR                                                ALLARR........2400
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::                   ALLARR........2500
     1   GXSI, GETA, GZET                                                ALLARR........2600
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::                     ALLARR........2700
     1   FWK ,B                                                          ALLARR........2800
      INTEGER, DIMENSION(:), ALLOCATABLE ::                              ALLARR........2900
     1   IN, IQSOP, IQSOU, IPBC, IUBC,                                   ALLARR........3000
     2   NREG, LREG, IWK, IA, JA                                         ALLARR........3100
      INTEGER(1), DIMENSION(:), ALLOCATABLE ::                           ALLARR........3200
     1   IBCPBC, IBCUBC, IBCSOP, IBCSOU                                  ALLARR........3300
      INTEGER, DIMENSION(:), ALLOCATABLE ::                              ALLARR........3400
     1   IIDPBC,IIDUBC,IIDSOP,IIDSOU                                     ALLARR........3500
      CHARACTER*40, DIMENSION(:), ALLOCATABLE ::                         ALLARR........3600
     1   CIDBCS                                                          ALLARR........3700
      LOGICAL, DIMENSION(:), ALLOCATABLE ::                              ALLARR........3800
     1   BCSFL, BCSTR                                                    ALLARR........3900
      TYPE OBSDAT                                                        ALLARR........4000
         CHARACTER*40 :: NAME                                            ALLARR........4100
         CHARACTER*10 :: SCHED                                           ALLARR........4200
         CHARACTER*3 :: FRMT                                             ALLARR........4300
         INTEGER :: L                                                    ALLARR........4400
         DOUBLE PRECISION :: X, Y, Z                                     ALLARR........4500
         DOUBLE PRECISION :: XSI, ETA, ZET                               ALLARR........4600
      END TYPE OBSDAT                                                    ALLARR........4700
      TYPE (OBSDAT), DIMENSION (:), ALLOCATABLE :: OBSPTS                ALLARR........4800
C.....ARRAY SCHDLS IS DECLARED IN MODULE SCHDEF.                         ALLARR........4900
C                                                                        ALLARR........5000
      END MODULE ALLARR                                                  ALLARR........5100
C                                                                        ALLARR........5200
C     MODULE            L  L  D  E  F              SUTRA VERSION 2.2     LLDEF..........100
C                                                                        LLDEF..........200
C *** PURPOSE :                                                          LLDEF..........300
C ***  TO DEFINE THE DERIVED TYPE "LLD" (LINKED LIST OF                  LLDEF..........400
C ***  DOUBLE-PRECISION NUMBERS).                                        LLDEF..........500
C                                                                        LLDEF..........600
      MODULE LLDEF                                                       LLDEF..........700
      IMPLICIT NONE                                                      LLDEF..........800
C                                                                        LLDEF..........900
C.....DEFINE DERIVED TYPE LLD (DOUBLE-PRECISION LINKED LIST) WITH        LLDEF.........1000
C        THREE COMPONENTS: DVALU1, DVALU2 (DOUBLE-PRECISION NUMBERS),    LLDEF.........1100
C        AND NENT (POINTER TO NEXT ENTRY).                               LLDEF.........1200
      TYPE LLD                                                           LLDEF.........1300
         DOUBLE PRECISION :: DVALU1, DVALU2                              LLDEF.........1400
         TYPE (LLD), POINTER :: NENT                                     LLDEF.........1500
      END TYPE LLD                                                       LLDEF.........1600
C                                                                        LLDEF.........1700
      END MODULE LLDEF                                                   LLDEF.........1800
C                                                                        LLDEF.........1900
C     MODULE            E  X  P  I  N  T           SUTRA VERSION 2.2     EXPINT.........100
C                                                                        EXPINT.........200
C *** PURPOSE :                                                          EXPINT.........300
C ***  TO PROVIDE EXPLICIT INTERFACES FOR PROCEDURES THAT NEED THEM.     EXPINT.........400
C                                                                        EXPINT.........500
      MODULE EXPINT                                                      EXPINT.........600
      IMPLICIT NONE                                                      EXPINT.........700
C                                                                        EXPINT.........800
C.....EXPLICIT INTERFACE FOR SUBROUTINE LLD2AR                           EXPINT.........900
      INTERFACE                                                          EXPINT........1000
         SUBROUTINE LLD2AR(LSTLEN, DLIST, DARR1, DARR2)                  EXPINT........1100
            USE LLDEF                                                    EXPINT........1200
            INTEGER LSTLEN                                               EXPINT........1300
            TYPE (LLD), POINTER :: DLIST                                 EXPINT........1400
            DOUBLE PRECISION DARR1(*), DARR2(*)                          EXPINT........1500
         END SUBROUTINE LLD2AR                                           EXPINT........1600
      END INTERFACE                                                      EXPINT........1700
C                                                                        EXPINT........1800
C.....EXPLICIT INTERFACE FOR SUBROUTINE LLDINS                           EXPINT........1900
      INTERFACE                                                          EXPINT........2000
         SUBROUTINE LLDINS(LSTLEN, DLIST, DNUM1, DNUM2, DLAST)           EXPINT........2100
            USE LLDEF                                                    EXPINT........2200
            INTEGER LSTLEN                                               EXPINT........2300
            TYPE (LLD), POINTER :: DEN, DENPV, DENNW, DLIST, DLAST       EXPINT........2400
            DOUBLE PRECISION DNUM1, DNUM2                                EXPINT........2500
         END SUBROUTINE LLDINS                                           EXPINT........2600
      END INTERFACE                                                      EXPINT........2700
C                                                                        EXPINT........2800
C.....EXPLICIT INTERFACE FOR FUNCTION PUSWF                              EXPINT........2900
      INTERFACE                                                          EXPINT........3000
         FUNCTION PUSWF(L,XLOC,YLOC,ZLOC,SFRAC,PM1,UM1,PVEC,UVEC,        EXPINT........3100
     1      IN,LREG)                                                     EXPINT........3200
            DOUBLE PRECISION PUSWF(3),XLOC,YLOC,ZLOC,SFRAC               EXPINT........3300
            DOUBLE PRECISION PM1(NN),UM1(NN),PVEC(NN),UVEC(NN)           EXPINT........3400
            INTEGER L,IN(NIN),LREG(NE)                                   EXPINT........3500
            INTEGER NN,NE,NIN,NBI,NCBI,NB,NBHALF,NPBC,NUBC,              EXPINT........3600
     1         NSOP,NSOU,NBCN                                            EXPINT........3700
            COMMON /DIMS/ NN,NE,NIN,NBI,NCBI,NB,NBHALF,NPBC,NUBC,        EXPINT........3800
     1         NSOP,NSOU,NBCN,NCIDB                                      EXPINT........3900
         END FUNCTION                                                    EXPINT........4000
      END INTERFACE                                                      EXPINT........4100
C                                                                        EXPINT........4200
C.....EXPLICIT INTERFACE FOR FUNCTION DP3STR                             EXPINT........4300
      INTERFACE                                                          EXPINT........4400
         FUNCTION DP3STR(DPA)                                            EXPINT........4500
         DOUBLE PRECISION DPA(3)                                         EXPINT........4600
         CHARACTER DP3STR*45                                             EXPINT........4700
         END FUNCTION                                                    EXPINT........4800
      END INTERFACE                                                      EXPINT........4900
C                                                                        EXPINT........5000
C.....EXPLICIT INTERFACE FOR SUBROUTINE READIF                           EXPINT........5100
      INTERFACE                                                          EXPINT........5200
         SUBROUTINE READIF(KUU, NFB, INTFIL, ERRIN, CHERIN)              EXPINT........5300
            PARAMETER (KINMIN=10)                                        EXPINT........5400
            CHARACTER INTFIL*1000                                        EXPINT........5500
            CHARACTER*80 ERRCOD,ERRIN,CHERR(10)                          EXPINT........5600
            CHARACTER*80, DIMENSION(10), OPTIONAL :: CHERIN              EXPINT........5700
            CHARACTER*80 UNAME,FNAME                                     EXPINT........5800
            CHARACTER ERRF*3, FINS*80                                    EXPINT........5900
            LOGICAL IS                                                   EXPINT........6000
            DIMENSION INERR(10),RLERR(10)                                EXPINT........6100
            DIMENSION IUNIT(0:13)                                        EXPINT........6200
            DIMENSION FNAME(0:13)                                        EXPINT........6300
            COMMON /FNAMES/ UNAME,FNAME                                  EXPINT........6400
            COMMON /FUNIB/ NFBCS                                         EXPINT........6500
            COMMON /FUNITA/ IUNIT                                        EXPINT........6600
            COMMON /FUNITS/ K00,K0,K1,K2,K3,K4,K5,K6,K7,K8,K9,           EXPINT........6700
     1         K10,K11,K12,K13                                           EXPINT........6800
            COMMON /OBS/ NOBSN,NTOBS,NOBCYC,NOBLIN,NFLOMX                EXPINT........6900
         END SUBROUTINE READIF                                           EXPINT........7000
      END INTERFACE                                                      EXPINT........7100
C                                                                        EXPINT........7200
      END MODULE EXPINT                                                  EXPINT........7300
C                                                                        EXPINT........7400
C     MODULE            P  T  R  D  E  F           SUTRA VERSION 2.2     PTRDEF.........100
C                                                                        PTRDEF.........200
C *** PURPOSE :                                                          PTRDEF.........300
C ***  TO DEFINE POINTERS AND ARRAYS NEEDED TO CONSTRUCT THE             PTRDEF.........400
C ***  IA AND JA ARRAYS.                                                 PTRDEF.........500
C                                                                        PTRDEF.........600
      MODULE PTRDEF                                                      PTRDEF.........700
      IMPLICIT NONE                                                      PTRDEF.........800
C.....DEFINE DERIVED TYPE LNKLST (LINKED LIST) WITH TWO COMPONENTS:      PTRDEF.........900
C        NODNUM (NODE NUMBER) AND NENT (POINTER TO NEXT ENTRY).          PTRDEF........1000
      TYPE LNKLST                                                        PTRDEF........1100
         INTEGER :: NODNUM                                               PTRDEF........1200
         TYPE (LNKLST), POINTER :: NENT                                  PTRDEF........1300
      END TYPE LNKLST                                                    PTRDEF........1400
C.....DECLARE DENT, DENTPV, DENTPI, AND DENTNW AS GENERAL-PURPOSE        PTRDEF........1500
C        POINTERS OF TYPE LNKLST.                                        PTRDEF........1600
      TYPE (LNKLST), POINTER :: DENT, DENTPV, DENTPI, DENTNW             PTRDEF........1700
C.....DEFINE DERIVED TYPE IPOINT WITH ONE COMPONENT: A POINTER, PL,      PTRDEF........1800
C        OF TYPE LNKLST.                                                 PTRDEF........1900
      TYPE IPOINT                                                        PTRDEF........2000
         TYPE (LNKLST), POINTER :: PL                                    PTRDEF........2100
      END TYPE IPOINT                                                    PTRDEF........2200
C.....DECLARE HLIST, AN ARRAY OF POINTERS THAT WILL POINT TO THE HEAD    PTRDEF........2300
C        OF THE LINKED LIST OF NEIGHBORS FOR EACH NODE.                  PTRDEF........2400
      TYPE (IPOINT), ALLOCATABLE :: HLIST(:)                             PTRDEF........2500
C.....DECLARE ARRAY LLIST.                                               PTRDEF........2600
      INTEGER, DIMENSION(:), ALLOCATABLE :: LLIST                        PTRDEF........2700
C                                                                        PTRDEF........2800
      END MODULE PTRDEF                                                  PTRDEF........2900
C                                                                        PTRDEF........3000
C     MODULE            S  C  H  D  E  F           SUTRA VERSION 2.2     SCHDEF.........100
C                                                                        SCHDEF.........200
C *** PURPOSE :                                                          SCHDEF.........300
C ***  TO DEFINE POINTERS AND ARRAYS ASSOCIATED WITH SCHEDULES.          SCHDEF.........400
C                                                                        SCHDEF.........500
      MODULE SCHDEF                                                      SCHDEF.........600
      USE LLDEF                                                          SCHDEF.........700
      IMPLICIT NONE                                                      SCHDEF.........800
C                                                                        SCHDEF.........900
C.....DEFINE DERIVED TYPE SCHLST WITH FIVE COMPONENTS:                   SCHDEF........1000
C        NAME (SCHEDULE NAME)                                            SCHDEF........1100
C        LLEN (LIST LENGTH)                                              SCHDEF........1200
C        SLAST (POINTER TO TAIL OF LIST)                                 SCHDEF........1300
C        SLIST (POINTER TO HEAD OF LIST)                                 SCHDEF........1400
      TYPE SCHLST                                                        SCHDEF........1500
         CHARACTER*10 NAME                                               SCHDEF........1600
         INTEGER :: LLEN                                                 SCHDEF........1700
         TYPE (LLD), POINTER :: SLAST                                    SCHDEF........1800
         TYPE (LLD), POINTER :: SLIST                                    SCHDEF........1900
      END TYPE SCHLST                                                    SCHDEF........2000
C                                                                        SCHDEF........2100
C.....DECLARE SCHDLS, AN ARRAY OF INFORMATION ABOUT EACH SCHEDULE.       SCHDEF........2200
      TYPE (SCHLST), ALLOCATABLE :: SCHDLS(:)                            SCHDEF........2300
C                                                                        SCHDEF........2400
C.....DEFINE DERIVED TYPE OBSFM                                          SCHDEF........2500
      TYPE OBSFM                                                         SCHDEF........2600
         INTEGER :: ISCHED                                               SCHDEF........2700
         CHARACTER*3 :: FRMT                                             SCHDEF........2800
      END TYPE OBSFM                                                     SCHDEF........2900
C                                                                        SCHDEF........3000
C.....DEFINE DERIVED TYPE BCSFM                                          SCHDEF........3100
      TYPE BCSFM                                                         SCHDEF........3200
         INTEGER :: ISCHED                                               SCHDEF........3300
      END TYPE BCSFM                                                     SCHDEF........3400
C                                                                        SCHDEF........3500
C.....DECLARE MORE ARRAYS                                                SCHDEF........3600
      TYPE (OBSFM), ALLOCATABLE :: OFP(:)                                SCHDEF........3700
      TYPE (BCSFM), ALLOCATABLE :: BFP(:)                                SCHDEF........3800
      INTEGER, ALLOCATABLE :: IUNIO(:)                                   SCHDEF........3900
      CHARACTER*80, ALLOCATABLE :: FNAMO(:)                              SCHDEF........4000
      LOGICAL, ALLOCATABLE :: ONCK78(:)                                  SCHDEF........4100
C                                                                        SCHDEF........4200
      END MODULE SCHDEF                                                  SCHDEF........4300
C                                                                        SCHDEF........4400
C     MODULE            B  C  S  D  E  F           SUTRA VERSION 2.2     BCSDEF.........100
C                                                                        BCSDEF.........200
C *** PURPOSE :                                                          BCSDEF.........300
C ***  TO DEFINE ARRAYS ASSOCIATED WITH BCS FILES.                       BCSDEF.........400
C                                                                        BCSDEF.........500
      MODULE BCSDEF                                                      BCSDEF.........600
      USE LLDEF                                                          BCSDEF.........700
      IMPLICIT NONE                                                      BCSDEF.........800
C                                                                        BCSDEF.........900
C.....DECLARE ARRAYS                                                     BCSDEF........1000
      INTEGER, ALLOCATABLE :: IUNIB(:)                                   BCSDEF........1100
      CHARACTER*80, ALLOCATABLE :: FNAMB(:)                              BCSDEF........1200
      INTEGER, ALLOCATABLE :: LCNT(:)                                    BCSDEF........1300
      TYPE (LLD), ALLOCATABLE :: DENBCS(:)                               BCSDEF........1400
C                                                                        BCSDEF........1500
      END MODULE BCSDEF                                                  BCSDEF........1600
C                                                                        BCSDEF........1700
C     MODULE            F  I  N  D  E  F           SUTRA VERSION 2.2     FINDEF.........100
C                                                                        FINDEF.........200
C *** PURPOSE :                                                          FINDEF.........300
C ***  TO DEFINE ARRAYS ASSOCIATED FILE INSERTION.                       FINDEF.........400
C                                                                        FINDEF.........500
      MODULE FINDEF                                                      FINDEF.........600
      IMPLICIT NONE                                                      FINDEF.........700
C                                                                        FINDEF.........800
C.....DECLARE ARRAYS                                                     FINDEF.........900
      INTEGER, ALLOCATABLE :: NKS(:), KLIST(:,:)                         FINDEF........1000
      CHARACTER*80, ALLOCATABLE :: FNAIN(:,:)                            FINDEF........1100
C                                                                        FINDEF........1200
      END MODULE FINDEF                                                  FINDEF........1300
