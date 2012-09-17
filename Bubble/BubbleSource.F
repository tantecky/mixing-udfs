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
        IF(ARGS(ILOC,2) .GT. 1.0E0 .OR. ARGS(ILOC,2)  .LT. 0.0E0) THEN
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
      
      CRESLT = 'GOOD'
      
      END
C=======================================================================
      REAL FUNCTION COMPUTE_SOURCE(NLOC, ILOC, ICLASS, RALFA, RF)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      REAL RHO_G
      PARAMETER (RHO_G = 1.125E0)
C-----Called functions
      REAL BBRI
      REAL BAGI
      REAL DBRI
      REAL DAGI
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, 1:NUMBER_OF_CLASSES)
C-----Code

#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        CALL ERRMSG('Wrong COMPUTE_SOURCE - ICLASS')
        CALL ERRMSG(ICLASS)
        STOP
      ENDIF
#endif

      COMPUTE_SOURCE = RHO_G*BUBBLE_CLASSES_VOL(ICLASS)*
     * (
     *  BBRI(NLOC, ILOC, ICLASS, RALFA, RF)
     *  +BAGI(NLOC, ILOC, ICLASS, RALFA, RF)
     *  -DBRI(NLOC, ILOC, ICLASS, RALFA, RF)
     *  -DAGI(NLOC, ILOC, ICLASS, RALFA, RF)
     * )

      END
C=======================================================================
      REAL FUNCTION N(NLOC, ILOC, ICLASS, RALFA, RF)
      IMPLICIT NONE
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
      REAL RF(1:NLOC, 1:NUMBER_OF_CLASSES)
C-----Code

#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        CALL ERRMSG('Wrong N - ICLASS')
        CALL ERRMSG(ICLASS)
        STOP
      ENDIF
#endif

#ifdef DEBUG
      IF(RALFA(ILOC) .GT. 1.0E0 .OR. RALFA(ILOC)  .LT. 0.0E0) THEN
        CALL ERRMSG('Wrong air.Volume Fraction - N')
        CALL ERRMSG(RALFA(ILOC))
        STOP
        ENDIF
#endif

      N = RALFA(ILOC)*RF(ILOC,ICLASS)/BUBBLE_CLASSES_VOL(ICLASS)
      END

C=======================================================================
      REAL FUNCTION GK15(FCE, A, B, ICLASS, J)
      IMPLICIT NONE 
C-----Called functions
      REAL FCE
      EXTERNAL FCE 
C-----Arguments
      REAL A
      REAL B
      INTEGER ICLASS
      INTEGER J
C-----Locale variables
      REAL TRANS1
      REAL TRANS2
      REAL INTEGRAL
      INTEGER I
      REAL NODES(0:6)
      REAL WEIGHTS(0:7)
      
      DATA NODES
     */0.991455371120813E0,
     * 0.949107912342759E0,
     * 0.864864423359769E0,
     * 0.741531185599394E0,
     * 0.586087235467691E0,
     * 0.405845151377397E0,
     * 0.207784955007898E0
     */
     
      DATA WEIGHTS
     */ 0.022935322010529E0,
     * 0.063092092629979E0,
     * 0.104790010322250E0,
     * 0.140653259715525E0,
     * 0.169004726639267E0,
     * 0.190350578064785E0,
     * 0.204432940075298E0,
     * 0.209482141084728E0
     */
C-----Code
      INTEGRAL = 0.0E0
      TRANS1 = (B-A)/2.0E0
      TRANS2 = (A+B)/2.0E0
      
      DO I = 0, 6
        INTEGRAL = INTEGRAL +
     *  WEIGHTS(I)*((FCE(ICLASS, TRANS1*NODES(I) + TRANS2, J)) +
     *  (FCE(ICLASS, TRANS1*(-NODES(I)) + TRANS2, J)))
      ENDDO
  
  
      GK15 = TRANS1*(INTEGRAL +  WEIGHTS(7)*FCE(ICLASS, TRANS2, J))
      
      END
C=======================================================================
      REAL FUNCTION GAMMA_IJ(ICLASS, J)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Called functions
      REAL XI_BETA
      EXTERNAL XI_BETA
      REAL XI_MINUS_ONE_BETA
      EXTERNAL XI_MINUS_ONE_BETA
      REAL GK15
C-----Arguments
      INTEGER ICLASS
      INTEGER J
      
      IF(ICLASS .EQ. 1) THEN
        GAMMA_IJ = 
     *     GK15(XI_BETA, BUBBLE_CLASSES_VOL(ICLASS), 
     *     BUBBLE_CLASSES_VOL(ICLASS+1), ICLASS, J)     
        RETURN
      ENDIF
      
      IF(ICLASS .EQ. J) THEN
        GAMMA_IJ = 
     *     GK15(XI_MINUS_ONE_BETA, BUBBLE_CLASSES_VOL(ICLASS-1), 
     *     BUBBLE_CLASSES_VOL(ICLASS), ICLASS, J)
        RETURN
      ENDIF
      
      GAMMA_IJ =
     *GK15(XI_MINUS_ONE_BETA, BUBBLE_CLASSES_VOL(ICLASS-1), 
     *     BUBBLE_CLASSES_VOL(ICLASS), ICLASS, J)
     *     +
     *GK15(XI_BETA, BUBBLE_CLASSES_VOL(ICLASS), 
     *     BUBBLE_CLASSES_VOL(ICLASS+1), ICLASS, J)              
      
      END
C=======================================================================
      REAL FUNCTION BBRI(NLOC, ILOC, ICLASS, RALFA, RF)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      REAL N
      REAL GAMMA_IJ
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, NUMBER_OF_CLASSES)
C-----Locale variables
      INTEGER J
C-----Code
      
      BBRI = 0.0E0
      
      DO J = ICLASS, NUMBER_OF_CLASSES
        BBRI = BBRI + GAMMA_IJ(ICLASS, J)*N(NLOC, ILOC, J, RALFA, RF)
      ENDDO
      
#ifdef DEBUG
      IF(ISNAN(BBRI) .EQV. .TRUE.) THEN
        CALL ERRMSG('BBRI ISNAN')
        STOP
      ENDIF
#endif
      
      END
C=======================================================================
      REAL FUNCTION BAGI(NLOC, ILOC, ICLASS, RALFA, RF)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Called functions
      REAL N
      REAL XI
      REAL XI_MINUS_ONE
      INTEGER KRONECKER_D
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, NUMBER_OF_CLASSES)
C-----Locale variables
      INTEGER J
      INTEGER K
      REAL V
      REAL BRACKET_PRODUCT
      REAL AJK
      DATA AJK /1.0E-6/
C-----Code
#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        CALL ERRMSG('Wrong BAGI - ICLASS')
        CALL ERRMSG(ICLASS)
        STOP
      ENDIF
#endif

      BAGI = 0.0E0

      DO J = 1, (NUMBER_OF_CLASSES - 1)
        DO K = J, (NUMBER_OF_CLASSES - 1)
        
          V = BUBBLE_CLASSES_VOL(J)+BUBBLE_CLASSES_VOL(K)
          BRACKET_PRODUCT = 0.0E0
          
          IF(ICLASS .NE. 1) THEN
            IF(BUBBLE_CLASSES_VOL(ICLASS - 1) .LT. V .AND.
     *        V .LT. BUBBLE_CLASSES_VOL(ICLASS)) THEN
             BRACKET_PRODUCT = XI_MINUS_ONE(ICLASS, V)
            ENDIF
          ENDIF
          
          IF(ICLASS .NE. NUMBER_OF_CLASSES) THEN
            IF(BUBBLE_CLASSES_VOL(ICLASS) .LT. V .AND.
     *        V .LT. BUBBLE_CLASSES_VOL(ICLASS + 1)) THEN
              BRACKET_PRODUCT = XI(ICLASS, V)
            ENDIF
          ENDIF
          
        BAGI = BAGI + 
     *  BRACKET_PRODUCT*(1.0E0 - 0.5E0*KRONECKER_D(J,K))
     *  *N(NLOC, ILOC, J, RALFA, RF)*N(NLOC, ILOC, K, RALFA, RF)*AJK
          
        END DO
      END DO 

#ifdef DEBUG
      IF(ISNAN(BAGI) .EQV. .TRUE.) THEN
        CALL ERRMSG('BAGI ISNAN')
        STOP
      ENDIF
#endif
      END 
C=======================================================================
      REAL FUNCTION XI_BETA(I, V, J)
      IMPLICIT NONE
C-----Called functions
      REAL BETA
      REAL XI
C-----Arguments
      INTEGER J
      INTEGER I
      REAL V
      
      XI_BETA = XI(I,V)*BETA(J, V)
      END
C=======================================================================
      REAL FUNCTION XI_MINUS_ONE_BETA(I, V, J)
      IMPLICIT NONE
C-----Called functions
      REAL BETA
      REAL XI_MINUS_ONE
C-----Arguments
      INTEGER J
      INTEGER I
      REAL V
      
      XI_MINUS_ONE_BETA = XI_MINUS_ONE(I,V)*BETA(J, V)
      END
C=======================================================================
      REAL FUNCTION BETA(J, V)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER J
      REAL V
      

      IF(J .GT. NUMBER_OF_CLASSES .OR. J .LT. 1) THEN
        WRITE(*,*) ('Wrong XI - J')
        WRITE(*,*) (J)
        STOP
      ENDIF

      BETA = 60.0E0/BUBBLE_CLASSES_VOL(J) 
     * *(V/BUBBLE_CLASSES_VOL(J))**2
     * *(1.0E0 - V/BUBBLE_CLASSES_VOL(J))**2
      
      END
C=======================================================================
      REAL FUNCTION XI(I, V)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER I
      REAL V
      
#ifdef DEBUG
      IF(I .GT. NUMBER_OF_CLASSES .OR. I .LT. 1) THEN
        CALL ERRMSG('Wrong XI - I')
        CALL ERRMSG(I)
        STOP
      ENDIF
#endif
      XI = (BUBBLE_CLASSES_VOL(I+1) - V)
     */(BUBBLE_CLASSES_VOL(I+1) - BUBBLE_CLASSES_VOL(I))
      
      END
C=======================================================================
      REAL FUNCTION XI_MINUS_ONE(I, V)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER I
      REAL V
      
#ifdef DEBUG
      IF(I .GT. NUMBER_OF_CLASSES .OR. I .LT. 1) THEN
        CALL ERRMSG('Wrong XI_MINUS_ONE - I')
        CALL ERRMSG(I)
        STOP
      ENDIF
#endif
      XI_MINUS_ONE = (V - BUBBLE_CLASSES_VOL(I-1))
     */(BUBBLE_CLASSES_VOL(I) - BUBBLE_CLASSES_VOL(I-1))
      
      END
C=======================================================================
      INTEGER FUNCTION KRONECKER_D(J, K)
      IMPLICIT NONE
C-----Arguments
      INTEGER J
      INTEGER K
      
      IF(J .EQ. K) THEN
        KRONECKER_D = 1
      ELSE
        KRONECKER_D = 0
      ENDIF
      
      END
C=======================================================================
      REAL FUNCTION DBRI(NLOC, ILOC, ICLASS, RALFA, RF)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      REAL N
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, 1:NUMBER_OF_CLASSES)
C-----Code

#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        CALL ERRMSG('Wrong DBRI - ICLASS')
        CALL ERRMSG(ICLASS)
        STOP
      ENDIF
#endif

      IF(ICLASS .EQ. 1) THEN
        DBRI = 0.E0
      ELSE
        DBRI = N(NLOC, ILOC, ICLASS, RALFA, RF)
      ENDIF

#ifdef DEBUG
      IF(ISNAN(DBRI) .EQV. .TRUE.) THEN
        CALL ERRMSG('DBRI ISNAN')
        STOP
      ENDIF
#endif
      END 
C=======================================================================
      REAL FUNCTION DAGI(NLOC, ILOC, ICLASS, RALFA, RF)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      REAL N
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, 1:NUMBER_OF_CLASSES)
C-----Locale variables
      INTEGER J
      REAL AJK
      DATA AJK /1.0E-6/
C-----Code
      DAGI = 0.0E0

#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        CALL ERRMSG('Wrong DAGI - ICLASS')
        CALL ERRMSG(ICLASS)
        STOP
      ENDIF
#endif

      IF(ICLASS .EQ. NUMBER_OF_CLASSES) THEN
        RETURN
      ELSE
        DO J = 1, (NUMBER_OF_CLASSES - 1)
          DAGI = DAGI  
     *    + N(NLOC, ILOC, ICLASS, RALFA, RF)
     *      *N(NLOC, ILOC, J, RALFA, RF)*AJK
        END DO
      ENDIF

#ifdef DEBUG
      IF(ISNAN(DAGI) .EQV. .TRUE.) THEN
        CALL ERRMSG('DAGI ISNAN')
        STOP
      ENDIF
#endif
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
     */ 6.54498469497874E-11
     *, 5.23598775598299E-10
     *, 4.18879020478639E-09
     *, 1.41371669411541E-08
     *, 3.35103216382911E-08
     *, 6.54498469497874E-08
     *, 1.13097335529233E-07
     *, 1.79594380030217E-07
     *, 2.68082573106329E-07
     *, 5.23598775598299E-07
     *, 9.04778684233860E-07
     *, 2.14466058485063E-06 /

C-----Common blocks
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
      
      END
