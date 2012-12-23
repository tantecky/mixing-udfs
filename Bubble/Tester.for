#define REAL DOUBLE PRECISION
/*model constants*/
#define SIGMA_SURF_TENS (0.0728E0)
#define RHO_GAS (1.185E0)
#define RHO_LIQUID (997.0E0)
#define MU_LIQUID (0.0008899E0)
      
      PROGRAM TESTER
      IMPLICIT NONE
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      INTEGER NLOC
      INTEGER ICLASS
      PARAMETER (NLOC = 1)
      INTEGER NARG
      PARAMETER (NARG = 15)
      INTEGER NRET
      PARAMETER (NRET = 1)
      
      REAL G_DBRI
      REAL G_DAGI
      REAL G_BBRI
      REAL G_BAGI
      
      COMMON /C_DBRI/ G_DBRI
      COMMON /C_DAGI/ G_DAGI
      COMMON /C_BBRI/ G_BBRI
      COMMON /C_BAGI/ G_BAGI
      
      REAL ARGS(NLOC,NARG), RET(NLOC,NRET)
      REAL TOTSUM
      
      TOTSUM = 0.0E0
      
      DO ICLASS = 1, NUMBER_OF_CLASSES
        ARGS(1,1) = ICLASS
        ARGS(1:NLOC,2) = 0.1E0
        
        ARGS(1:NLOC,3) = 0.1E0
        ARGS(1:NLOC,4) = 0.1E0
        ARGS(1:NLOC,5) = 0.1E0
        ARGS(1:NLOC,6) = 0.1E0
        ARGS(1:NLOC,7) = 0.1E0
        ARGS(1:NLOC,8) = 0.1E0
        ARGS(1:NLOC,9) = 0.1E0
        ARGS(1:NLOC,10) = 0.1E0
        ARGS(1:NLOC,11) = 0.05E0
        ARGS(1:NLOC,12) = 0.05E0
        ARGS(1:NLOC,13) = 0.05E0
        ARGS(1:NLOC,14) = 0.05E0
        ARGS(1:NLOC,15) = 0.001E0
        
        CALL BUBBLE_SOURCE(NLOC, NRET, NARG, RET, ARGS)
        TOTSUM = TOTSUM + RET(NLOC,NRET)
      ENDDO
      
      WRITE(*,*) '==========================='
      WRITE(*,'(A, E20.10)') 'SUM_SOURCE: ', TOTSUM
      WRITE(*,'(A, E20.10)') 'BAGI-DAGI: ', (G_BAGI-G_DAGI)
      WRITE(*,'(A, E20.10)') 'BBRI-DBRI ', (G_BBRI-G_DBRI)
      
      END

      SUBROUTINE BUBBLE_SOURCE( 
     &  NLOC, NRET, NARG, RET, ARGS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      REAL COMPUTE_SOURCE
C-----Arguments
      INTEGER NLOC,NARG,NRET
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

      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong ICLASS')
        WRITE(*,*) (ICLASS)
        STOP
      ENDIF


      DO ILOC = 1, NLOC
        IF(ARGS(ILOC,2) .GT. 1.0E0 .OR. ARGS(ILOC,2)  .LT. 0.0E0) THEN
          WRITE(*,*) ('Wrong air.Volume Fraction')
          WRITE(*,*) (ARGS(ILOC,2))
          STOP
        ENDIF
      END DO


      DO ILOC = 1, NLOC
        RET(ILOC,NRET) = COMPUTE_SOURCE(NLOC, ILOC, ICLASS, 
     *  ARGS(1,2), RF, ARGS(1,15))   
      END DO 
     
      END
C=======================================================================
      REAL FUNCTION 
     *COMPUTE_SOURCE(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      REAL RHO_G
      PARAMETER (RHO_G = RHO_GAS)
C-----Called functions
      REAL BBRI
      REAL BAGI
      REAL DBRI
      REAL DAGI
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
      REAL G_DBRI
      REAL G_DAGI
      REAL G_BBRI
      REAL G_BAGI
      COMMON /C_DBRI/ G_DBRI
      COMMON /C_DAGI/ G_DAGI
      COMMON /C_BBRI/ G_BBRI
      COMMON /C_BAGI/ G_BAGI
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      REAL EPS(NLOC)
C-----Code


      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong COMPUTE_SOURCE - ICLASS')
        WRITE(*,*) (ICLASS)
        STOP
      ENDIF


      COMPUTE_SOURCE = RHO_G*BUBBLE_CLASSES_VOL(ICLASS)*
     * (
     *  BBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
     *  +BAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
     *  -DBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
     *  -DAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
     * )
      
      G_BBRI = G_BBRI + BUBBLE_CLASSES_VOL(ICLASS)*RHO_G*
     *BBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)    
               
      G_BAGI = G_BAGI + BUBBLE_CLASSES_VOL(ICLASS)*RHO_G*
     *BAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
               
      G_DBRI = G_DBRI + BUBBLE_CLASSES_VOL(ICLASS)*RHO_G*
     *DBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
               
      G_DAGI = G_DAGI + BUBBLE_CLASSES_VOL(ICLASS)*RHO_G*
     *DAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS) 
      
      WRITE(*,*) '---------------------------'
      WRITE(*,'(A, I0)') 'CLASS: ', ICLASS
      WRITE(*,'(A, F5.2)') 'VOLFRAC_G: ', RALFA(ILOC)
      WRITE(*,'(A, I0, A, F5.2)') 'F', ICLASS, ':', RF(ILOC,ICLASS)
      WRITE(*,
     * '(A, E10.5)') 'BBRI: ', BBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      WRITE(*,
     * '(A, E10.5)') 'BAGI: ', BAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'DBRI: ',-DBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'DAGI: ',-DAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'BAGI-DAGI: ',
     * BAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS) 
     *-DAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'BBRI-DBRI: ',
     * BBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS) 
     *-DBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      WRITE(*,'(A, E20.10)') 'SOURCE: ', COMPUTE_SOURCE
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


      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong N - ICLASS')
        WRITE(*,*) (ICLASS)
        STOP
      ENDIF



      IF(RALFA(ILOC) .GT. 1.0E0 .OR. RALFA(ILOC)  .LT. 0.0E0) THEN
        WRITE(*,*) ('Wrong air.Volume Fraction - N')
        WRITE(*,*) (RALFA(ILOC))
        STOP
      ENDIF


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
      REAL FUNCTION GAMMA_IJ(ICLASS, J, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
      REAL G_EPS
      COMMON /C_EPS/ G_EPS
C-----Called functions
      REAL XI_BETA
      EXTERNAL XI_BETA
      REAL XI_MINUS_ONE_BETA
      EXTERNAL XI_MINUS_ONE_BETA
      REAL GK15
C-----Arguments
      INTEGER ICLASS
      INTEGER J
      REAL EPS
C-----Code
      G_EPS = EPS

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
      REAL FUNCTION BBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Called functions
      REAL N
      REAL GAMMA_IJ
      REAL G_I
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, NUMBER_OF_CLASSES)
      REAL EPS(NLOC)
C-----Locale variables
      INTEGER J
C-----Code
      
      BBRI = 0.0E0
      
      DO J = ICLASS, NUMBER_OF_CLASSES
        BBRI = BBRI + G_I(J, EPS(ILOC))*
     *  GAMMA_IJ(ICLASS, J, EPS(ILOC))*N(NLOC, ILOC, J, RALFA, RF)
      ENDDO

      IF(ISNAN(BBRI) .EQV. .TRUE.) THEN
        WRITE(*,*) ('BBRI ISNAN')
        STOP
      ENDIF
      
      
      END
C=======================================================================
      REAL FUNCTION BAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
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
      REAL A_IJ
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, NUMBER_OF_CLASSES)
      REAL EPS(NLOC)
C-----Locale variables
      INTEGER J
      INTEGER K
      REAL V
      REAL BRACKET_PRODUCT
C-----Code

      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong BAGI - ICLASS')
        WRITE(*,*) (ICLASS)
        STOP
      ENDIF


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
     *  *N(NLOC, ILOC, J, RALFA, RF)*N(NLOC, ILOC, K, RALFA, RF)*
     *   A_IJ(J, K, EPS(ILOC))
          
        END DO
      END DO 


      IF(ISNAN(BAGI) .EQV. .TRUE.) THEN
        WRITE(*,*) ('BAGI ISNAN')
        STOP
      ENDIF
      
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

      IF(I .GT. NUMBER_OF_CLASSES .OR. I .LT. 1) THEN
        WRITE(*,*) ('Wrong XI - I')
        WRITE(*,*) (I)
        STOP
      ENDIF

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
      

      IF(I .GT. NUMBER_OF_CLASSES .OR. I .LT. 1) THEN
        WRITE(*,*) ('Wrong XI_MINUS_ONE - I')
        WRITE(*,*) (I)
        STOP
      ENDIF

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
      REAL FUNCTION A_IJ(I, J, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      REAL SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      REAL RHO_L
      PARAMETER (RHO_L = RHO_LIQUID)
      REAL H0
      PARAMETER (H0 = 1.0E-4)
      REAL HF
      PARAMETER (HF = 1.0E-8)
      REAL COAL_FACTOR
      PARAMETER (COAL_FACTOR = 1.0E0)
      REAL PI
      PARAMETER (PI = 3.1415926535897931E0)
C-----Common blocks
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments   
      INTEGER I
      INTEGER J
      REAL EPS
      REAL FREQ
      REAL EFF
      REAL R_IJ
      
      R_IJ = BUBBLE_CLASSES_DIA(I)*BUBBLE_CLASSES_DIA(J)
     *      / (BUBBLE_CLASSES_DIA(I)+BUBBLE_CLASSES_DIA(J)) 
      FREQ = SQRT(2.E0)/4.E0 * PI   
     *       *(BUBBLE_CLASSES_DIA(I)+BUBBLE_CLASSES_DIA(J))**2.E0
     *       * EPS**(1.E0/3.E0) * (BUBBLE_CLASSES_DIA(I)**(2.E0/3.E0) 
     *       + BUBBLE_CLASSES_DIA(J)**(2.E0/3.E0))**(1.E0/2.E0)
      EFF = EXP(-SQRT(RHO_L) * R_IJ**(5.E0/6.E0) * EPS**(1.E0/3.E0) 
     *   * LOG(H0/HF) /(4.E0*SQRT(SIGMA)))    
      
      A_IJ = COAL_FACTOR * FREQ * EFF
      
      END 
C=======================================================================
      REAL FUNCTION G_I(I, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      REAL SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      REAL RHO_L
      PARAMETER (RHO_L = RHO_LIQUID)
      REAL RHO_G
      PARAMETER (RHO_G = RHO_GAS)
      REAL MU_L
      PARAMETER (MU_L = MU_LIQUID)
      REAL BREAKUP_FACTOR
      PARAMETER (BREAKUP_FACTOR = 1.0E0)
      REAL P_ERF
      PARAMETER (P_ERF = 0.3275911E0)
C-----Common blocks
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments   
      INTEGER I
      REAL EPS
      REAL ERF_ARG
      REAL T_ERF
      REAL ERF

      IF(I .EQ. 1) THEN
        G_I = 0.E0
      ELSE
        ERF_ARG = SQRT(0.04E0 * SIGMA/(RHO_L * EPS**(2.E0/3.E0)
     *          * BUBBLE_CLASSES_DIA(I)**(5.E0/3.E0)) + 0.01E0 * MU_L
     *          / (SQRT(RHO_L*RHO_G) * EPS**(1.E0/3.E0)
     *          * BUBBLE_CLASSES_DIA(I)**(4.E0/3.E0)))
        T_ERF = 1.E0 / (1.E0 + P_ERF*ERF_ARG)
        ERF = 1.E0 - (0.254829592E0*T_ERF - 0.284496736E0*T_ERF**2.E0
     *      + 1.421413741E0*T_ERF**3.E0 - 1.453152027*T_ERF**4.E0
     *      + 1.061405429E0*T_ERF**5.E0) * EXP(-ERF_ARG**2.E0)
        
        G_I = BREAKUP_FACTOR * EPS**(1.E0/3.E0) * (1.0E0-ERF)
      ENDIF
      
      END      
C=======================================================================
      REAL FUNCTION DBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      REAL N
      REAL G_I
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      REAL EPS(NLOC)
C-----Code


      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong DBRI - ICLASS')
        WRITE(*,*) (ICLASS)
        STOP
      ENDIF



      DBRI = N(NLOC, ILOC, ICLASS, RALFA, RF)*
     *G_I(ICLASS, EPS(ILOC))


      IF(ISNAN(DBRI) .EQV. .TRUE.) THEN
        WRITE(*,*) ('DBRI ISNAN')
        STOP
      ENDIF
      

      END 
C=======================================================================
      REAL FUNCTION DAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      REAL N
      REAL A_IJ
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      REAL RALFA(NLOC)
      REAL RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      REAL EPS(NLOC)
C-----Locale variables
      INTEGER J
C-----Code
      DAGI = 0.0E0


      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong DAGI - ICLASS')
        WRITE(*,*) (ICLASS)
        STOP
      ENDIF


      IF(ICLASS .EQ. NUMBER_OF_CLASSES) THEN
        RETURN
      ELSE
        DO J = 1, (NUMBER_OF_CLASSES - 1)
          DAGI = DAGI  
     *    + N(NLOC, ILOC, ICLASS, RALFA, RF)
     *      *N(NLOC, ILOC, J, RALFA, RF)*A_IJ(ICLASS, J, EPS(ILOC))
        END DO
      ENDIF


      IF(ISNAN(DAGI) .EQV. .TRUE.) THEN
        WRITE(*,*) ('DAGI ISNAN')
        STOP
      ENDIF
      
      END 
C=======================================================================
      LOGICAL FUNCTION ISFINITE(X)
      REAL X
      REAL TMP
      
      TMP = X * 10.0E0
            
      IF(X .EQ.TMP) THEN
       ISFINITE = .FALSE.
      ELSE
        ISFINITE = .TRUE.
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
      REAL G_EPS
      REAL G_DBRI
      REAL G_DAGI
      REAL G_BBRI
      REAL G_BAGI
      
      DATA G_DBRI /0.0E0/
      DATA G_DAGI /0.0E0/
      DATA G_BBRI /0.0E0/
      DATA G_BAGI /0.0E0/
      
C     diameter of bubble classes
      DATA BUBBLE_CLASSES_DIA /0.1E-3, 1.0E-3, 2.0E-3, 3.0E-3, 4.0E-3
     * , 5.0E-3, 6.0E-3, 7.0E-3, 8.0E-3, 10.0E-3, 12.0E-3, 16.0E-3/
     
      DATA BUBBLE_CLASSES_VOL 
     */ 5.23598775598299e-13
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
      COMMON /C_EPS/ G_EPS
      
      COMMON /C_DBRI/ G_DBRI
      COMMON /C_DAGI/ G_DAGI
      COMMON /C_BBRI/ G_BBRI
      COMMON /C_BAGI/ G_BAGI
      
      END 
      

