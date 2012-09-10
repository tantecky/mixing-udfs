#include "cfx5ext.h"
#define DEBUG

dllexport(bubble_source)
      SUBROUTINE BUBBLE_SOURCE( 
     &  NLOC, NRET, NARG, RET, ARGS, CRESLT,CZ,DZ,IZ,LZ,RZ )
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      REAL COMPUTE_SOURCE
C-----Arguments
      INTEGER NLOC,NARG,NRET
      CHARACTER CRESLT*(*)
      REAL ARGS(NLOC,NARG), RET(NLOC,NRET)
      
c     ICLASS = ARGS(1,1)
c     RALFA = ARGS(1:NLOC,2)

c     F1 = ARGS(1:NLOC,3)
c     F2 = ARGS(1:NLOC,4)
c     F3 = ARGS(1:NLOC,5)
c     F4 = ARGS(1:NLOC,6)
c     F5 = ARGS(1:NLOC,7)
c     F6 = ARGS(1:NLOC,8)
c     F7 = ARGS(1:NLOC,9)
c     F8 = ARGS(1:NLOC,10)
c     F9 = ARGS(1:NLOC,11)
c     F10 = ARGS(1:NLOC,12)
c     F11 = ARGS(1:NLOC,13)
c     F12 = ARGS(1:NLOC,14)
      
      INTEGER IZ(*)
      CHARACTER CZ(*)*(1)
      DOUBLE PRECISION DZ(*)
      LOGICAL LZ(*)
      REAL RZ(*)
C-----Locale variables
      INTEGER ICLASS
      INTEGER ILOC
      REAL RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      
      RF(:, 1) = ARGS(:,3)
      RF(:, 2) = ARGS(:,4)
      RF(:, 3) = ARGS(:,5)
      RF(:, 4) = ARGS(:,6)
      RF(:, 5) = ARGS(:,7)
      RF(:, 6) = ARGS(:,8)
      RF(:, 7) = ARGS(:,9)
      RF(:, 8) = ARGS(:,10)
      RF(:, 9) = ARGS(:,11)
      RF(:, 10) = ARGS(:,12)
      RF(:, 11) = ARGS(:,13)
      RF(:, 12) = ARGS(:,14)
      
C-----Code
      ICLASS = INT(ARGS(1,1))
      
#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        CALL ERRMSG('Wrong ICLASS')
        CALL ERRMSG(ICLASS)
        STOP
      ENDIF
#endif

#ifdef DEBUG
      DO ILOC = 1, NLOC
        IF(ARGS(ILOC,2) .GT. 1.0D0 .OR. ARGS(ILOC,2)  .LT. 0.0D0) THEN
          CALL ERRMSG('Wrong air.Volume Fraction')
          CALL ERRMSG(ARGS(ILOC,2))
          STOP
        ENDIF
      END DO
#endif

      DO ILOC = 1, NLOC
        RET(ILOC,NRET) = COMPUTE_SOURCE(NLOC, ILOC, ICLASS, 
     *  ARGS(1,2), RF)
      END DO
      
      END
C=======================================================================
      REAL FUNCTION COMPUTE_SOURCE(NLOC, ILOC, ICLASS, RALFA, RF)
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      REAL RHO_G
      PARAMETER (RHO_G = 9.81)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, NUMBER_OF_CLASSES)
C-----Code
      COMPUTE_SOURCE = RHO_G*BUBBLE_CLASSES_VOL(ICLASS)*
     * (-DBRI(NLOC, ILOC, ICLASS, RALFA, RF))

      END 
C=======================================================================
      REAL FUNCTION NI(NLOC, ILOC, ICLASS, RALFA, RF)
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, NUMBER_OF_CLASSES)
C-----Code
      NI = RALFA(ILOC)*RF(ILOC,ICLASS)/BUBBLE_CLASSES_VOL(ICLASS)
      END
C=======================================================================
      REAL FUNCTION DBRI(NLOC, ILOC, ICLASS, RALFA, RF)
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      REAL NI
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, NUMBER_OF_CLASSES)
C-----Code
      IF(ICLASS .EQ. 1) THEN
        DBRI = 0.D0
      ELSE
        DBRI = NI(NLOC, ILOC, ICLASS, RALFA, RF)
      ENDIF

      END 
C=======================================================================
      BLOCKDATA
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12) 
C-----Locale variables
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
C     diameter of bubble classes
      DATA BUBBLE_CLASSES_DIA /0.5E-3, 1.0E-3, 2.0E-3, 3.0E-3, 4.0E-3
     * , 5.0E-3, 6.0E-3, 7.0E-3, 8.0E-3, 10.0E-3, 12.0E-3, 16.0E-3/
     
      DATA BUBBLE_CLASSES_VOL 
     */ 6.54498469497874D-11
     *, 5.23598775598299D-10
     *, 4.18879020478639D-09
     *, 1.41371669411541D-08
     *, 3.35103216382911D-08
     *, 6.54498469497874D-08
     *, 1.13097335529233D-07
     *, 1.79594380030217D-07
     *, 2.68082573106329D-07
     *, 5.23598775598299D-07
     *, 9.04778684233860D-07
     *, 2.14466058485063D-06 /

C-----Common blocks
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
      
      END
