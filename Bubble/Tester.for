#define REAL DOUBLE PRECISION

c#define MODEL_ALOPEA
c#define MODEL_LEHR
#define MODEL_MARTINEZ_BAZAN

#define DEBUG

C-------- constants
#define SIGMA_SURF_TENS (0.0728E0)
#define RHO_GAS (1.185E0)
#define RHO_LIQUID (997.0E0)
#define MU_LIQUID (0.0008899E0)
#define PI_CONST (3.1415926535897931E0)
#define BREAKUP_F (0.008E0)
c#define BREAKUP_F (1.0E0)
C-------- constants
      
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
        ARGS(1:NLOC,15) = 1.0E0
        
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
           
#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong ICLASS')
        WRITE(*,*) (ICLASS)
        CALL ABORT()
      ENDIF


      DO ILOC = 1, NLOC
        IF(ARGS(ILOC,2) .GT. 1.0E0 .OR. ARGS(ILOC,2)  .LT. 0.0E0) THEN
          WRITE(*,*) ('Wrong air.Volume Fraction')
          WRITE(*,*) (ARGS(ILOC,2))
          CALL ABORT()
        ENDIF
      END DO
#endif

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

#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong COMPUTE_SOURCE - ICLASS')
        WRITE(*,*) (ICLASS)
        CALL ABORT()
      ENDIF
#endif

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
     * '(A, E20.10)') 'BBRI: ', BBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'BAGI: ', BAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
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

#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong N - ICLASS')
        WRITE(*,*) (ICLASS)
        CALL ABORT()
      ENDIF

      IF(RALFA(ILOC) .GT. 1.0E0 .OR. RALFA(ILOC)  .LT. 0.0E0) THEN
        WRITE(*,*) ('Wrong air.Volume Fraction - N')
        WRITE(*,*) (RALFA(ILOC))
        CALL ABORT()
      ENDIF
#endif

      N = RALFA(ILOC)*RF(ILOC,ICLASS)/BUBBLE_CLASSES_VOL(ICLASS)
      END

C=======================================================================
      REAL FUNCTION GK15(FCE, A, B, ICLASS, J, BRANCH)
      IMPLICIT NONE 
C-----Called functions
      REAL FCE
      EXTERNAL FCE 
C-----Arguments
      REAL A
      REAL B
      INTEGER ICLASS
      INTEGER J
      INTEGER BRANCH
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
      
#ifdef DEBUG      
      IF(A .GE. B) THEN
        WRITE(*,*) ('GK15: A >= B')
        WRITE(*,*) 'A=',A
        WRITE(*,*) 'B=',B
        WRITE(*,*) 'TRANS1=',TRANS1
        WRITE(*,*) 'TRANS2=',TRANS2
        CALL ABORT()
      ENDIF
#endif
      
      DO I = 0, 6
        INTEGRAL = INTEGRAL +
     *  WEIGHTS(I)*((FCE(ICLASS, TRANS1*NODES(I) + TRANS2, J, BRANCH)) +
     *  (FCE(ICLASS, TRANS1*(-NODES(I)) + TRANS2, J, BRANCH)))
      ENDDO
  
  
      GK15 = 
     * TRANS1*(INTEGRAL +  WEIGHTS(7)*FCE(ICLASS, TRANS2, J, BRANCH))
      
      END
C=======================================================================
      REAL FUNCTION GK61(FCE, A, B, ICLASS, J, BRANCH)
      IMPLICIT NONE 
C-----Called functions
      REAL FCE
      EXTERNAL FCE 
C-----Arguments
      REAL A
      REAL B
      INTEGER ICLASS
      INTEGER J
      INTEGER BRANCH
C-----Locale variables
      REAL TRANS1
      REAL TRANS2
      REAL INTEGRAL
      INTEGER I
      REAL NODES(0:29)
      REAL WEIGHTS(0:30)
      
      DATA NODES
     * /0.9994844100504906375713259E0,
     * 0.9968934840746495402716301E0,
     * 0.9916309968704045948586284E0,
     * 0.9836681232797472099700326E0,
     * 0.9731163225011262683746939E0,
     * 0.9600218649683075122168710E0,
     * 0.9443744447485599794158313E0,
     * 0.9262000474292743258793243E0,
     * 0.9055733076999077985465226E0,
     * 0.8825605357920526815431165E0,
     * 0.8572052335460610989586585E0,
     * 0.8295657623827683974428981E0,
     * 0.7997278358218390830136689E0,
     * 0.7677774321048261949179773E0,
     * 0.7337900624532268047261711E0,
     * 0.6978504947933157969322924E0,
     * 0.6600610641266269613700537E0,
     * 0.6205261829892428611404776E0,
     * 0.5793452358263616917560249E0,
     * 0.5366241481420198992641698E0,
     * 0.4924804678617785749936931E0,
     * 0.4470337695380891767806099E0,
     * 0.4004012548303943925354762E0,
     * 0.3527047255308781134710372E0,
     * 0.3040732022736250773726771E0,
     * 0.2546369261678898464398051E0,
     * 0.2045251166823098914389577E0,
     * 0.1538699136085835469637947E0,
     * 0.1028069379667370301470968E0,
     * 0.0514718425553176958330252E0/
     
      DATA WEIGHTS
     * /0.0013890136986770076245516E0,
     * 0.0038904611270998840512672E0,
     * 0.0066307039159312921733198E0,
     * 0.0092732796595177634284411E0,
     * 0.0118230152534963417422329E0,
     * 0.0143697295070458048124514E0,
     * 0.0169208891890532726275723E0,
     * 0.0194141411939423811734090E0,
     * 0.0218280358216091922971675E0,
     * 0.0241911620780806013656864E0,
     * 0.0265099548823331016106017E0,
     * 0.0287540487650412928439788E0,
     * 0.0309072575623877624728843E0,
     * 0.0329814470574837260318142E0,
     * 0.0349793380280600241374997E0,
     * 0.0368823646518212292239111E0,
     * 0.0386789456247275929503487E0,
     * 0.0403745389515359591119953E0,
     * 0.0419698102151642461471475E0,
     * 0.0434525397013560693168317E0,
     * 0.0448148001331626631923556E0,
     * 0.0460592382710069881162717E0,
     * 0.0471855465692991539452615E0,
     * 0.0481858617570871291407795E0,
     * 0.0490554345550297788875282E0,
     * 0.0497956834270742063578116E0,
     * 0.0504059214027823468408931E0,
     * 0.0508817958987496064922975E0,
     * 0.0512215478492587721706563E0,
     * 0.0514261285374590259338629E0,
     * 0.0514947294294515675583404E0/
C-----Code
      INTEGRAL = 0.0E0
      TRANS1 = (B-A)/2.0E0
      TRANS2 = (A+B)/2.0E0
      
#ifdef DEBUG      
      IF(A .GE. B) THEN
        WRITE(*,*) ('GK61: A >= B')
        WRITE(*,*) 'A=',A
        WRITE(*,*) 'B=',B
        WRITE(*,*) 'TRANS1=',TRANS1
        WRITE(*,*) 'TRANS2=',TRANS2
        CALL ABORT()
      ENDIF
#endif
      
      DO I = 0, 29
        INTEGRAL = INTEGRAL +
     *  WEIGHTS(I)*((FCE(ICLASS, TRANS1*NODES(I) + TRANS2, J, BRANCH)) +
     *  (FCE(ICLASS, TRANS1*(-NODES(I)) + TRANS2, J, BRANCH)))
      ENDDO
  
  
      GK61 = 
     * TRANS1*(INTEGRAL +  WEIGHTS(30)*FCE(ICLASS, TRANS2, J, BRANCH))
      
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
C-----Called functions
      REAL XI_BETA
      EXTERNAL XI_BETA
      REAL XI_MINUS_ONE_BETA
      EXTERNAL XI_MINUS_ONE_BETA
C-----Arguments
      INTEGER ICLASS
      INTEGER J
      REAL EPS
      
#ifdef MODEL_ALOPEA
C-----Called functions
      REAL GK15
C-----Code
      IF(ICLASS .EQ. 1 .AND. ICLASS .EQ. J) THEN
        GAMMA_IJ = 0.0E0
        RETURN
      ENDIF
               
      IF(ICLASS .EQ. 1) THEN
        GAMMA_IJ = 
     *     GK15(XI_BETA, BUBBLE_CLASSES_VOL(ICLASS), 
     *     BUBBLE_CLASSES_VOL(ICLASS+1), ICLASS, J, 0)     
        RETURN
      ENDIF
      
      IF(ICLASS .EQ. J) THEN
        GAMMA_IJ = 
     *     GK15(XI_MINUS_ONE_BETA, BUBBLE_CLASSES_VOL(ICLASS-1), 
     *     BUBBLE_CLASSES_VOL(ICLASS), ICLASS, J, 0)
        RETURN
      ENDIF
      
      GAMMA_IJ =
     *GK15(XI_MINUS_ONE_BETA, BUBBLE_CLASSES_VOL(ICLASS-1), 
     *     BUBBLE_CLASSES_VOL(ICLASS), ICLASS, J, 0)
     *     +
     *GK15(XI_BETA, BUBBLE_CLASSES_VOL(ICLASS), 
     *     BUBBLE_CLASSES_VOL(ICLASS+1), ICLASS, J, 0)              

#elif defined MODEL_LEHR
C-----Common blocks
      REAL G_EPS
      COMMON /C_EPS/ G_EPS
C-----Locale variables
      REAL V0
      REAL V0HALF
      REAL VA
      REAL VB
C-----Called functions
      REAL GK61
C-----Code
      IF(ICLASS .EQ. 1 .AND. ICLASS .EQ. J) THEN
        GAMMA_IJ = 0.0E0
        RETURN
      ENDIF

      G_EPS = EPS
      V0 = BUBBLE_CLASSES_VOL(J)
      V0HALF = (V0/2.E0)
                 
      IF(ICLASS .EQ. 1) THEN
        VA = BUBBLE_CLASSES_VOL(ICLASS)
        VB = BUBBLE_CLASSES_VOL(ICLASS+1)
                
        IF(VA .LT. V0HALF .AND. VB .GT. V0HALF) THEN
          GAMMA_IJ = 
     *    GK61(XI_BETA, VA, V0HALF, ICLASS, J, 1)
          
          GAMMA_IJ = GAMMA_IJ +
     *    GK61(XI_BETA, V0HALF, VB, ICLASS, J, 2)
        ELSEIF(VA .GT. V0HALF) THEN
          GAMMA_IJ = 
     *    GK61(XI_BETA, VA, VB, ICLASS, J, 2) 
        ELSEIF(VB .LT. V0HALF) THEN
          GAMMA_IJ = 
     *    GK61(XI_BETA, VA, VB, ICLASS, J, 1) 
        ELSE
          WRITE(*,*) ('No solution - BETA')
          WRITE(*,*) (J)
          CALL ABORT()
        ENDIF
            
        RETURN
      ENDIF
      
      IF(ICLASS .EQ. J) THEN
        VA = BUBBLE_CLASSES_VOL(ICLASS-1)
        VB = BUBBLE_CLASSES_VOL(ICLASS)
        
        IF(VA .LT. V0HALF .AND. VB .GT. V0HALF) THEN
          GAMMA_IJ = 
     *    GK61(XI_MINUS_ONE_BETA, VA, V0HALF, ICLASS, J, 1)
          
          GAMMA_IJ = GAMMA_IJ +
     *    GK61(XI_MINUS_ONE_BETA, V0HALF, VB, ICLASS, J, 2)
        ELSEIF(VA .GT. V0HALF) THEN
          GAMMA_IJ = 
     *    GK61(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 2) 
        ELSEIF(VB .LT. V0HALF) THEN
          GAMMA_IJ = 
     *    GK61(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 1) 
        ELSE
          WRITE(*,*) ('No solution - BETA')
          WRITE(*,*) (J)
          CALL ABORT()
        ENDIF
            
        RETURN
      ENDIF
      
      VA = BUBBLE_CLASSES_VOL(ICLASS)
      VB = BUBBLE_CLASSES_VOL(ICLASS+1)
      
      IF(VA .LT. V0HALF .AND. VB .GT. V0HALF) THEN
         GAMMA_IJ = 
     *   GK61(XI_BETA, VA, V0HALF, ICLASS, J, 1)
          
         GAMMA_IJ = GAMMA_IJ +
     *   GK61(XI_BETA, V0HALF, VB, ICLASS, J, 2)
      ELSEIF(VA .GT. V0HALF) THEN
         GAMMA_IJ = 
     *   GK61(XI_BETA, VA, VB, ICLASS, J, 2) 
      ELSEIF(VB .LT. V0HALF) THEN
        GAMMA_IJ = 
     *   GK61(XI_BETA, VA, VB, ICLASS, J, 1) 
      ELSE
        WRITE(*,*) ('No solution - BETA')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
      
      VA = BUBBLE_CLASSES_VOL(ICLASS-1)
      VB = BUBBLE_CLASSES_VOL(ICLASS)
        
      IF(VA .LT. V0HALF .AND. VB .GT. V0HALF) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  GK61(XI_MINUS_ONE_BETA, VA, V0HALF, ICLASS, J, 1)
          
        GAMMA_IJ = GAMMA_IJ +
     *  GK61(XI_MINUS_ONE_BETA, V0HALF, VB, ICLASS, J, 2)
      ELSEIF(VA .GT. V0HALF) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  GK61(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 2) 
      ELSEIF(VB .LT. V0HALF) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  GK61(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 1) 
      ELSE
        WRITE(*,*) ('No solution - BETA')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF

#elif defined MODEL_MARTINEZ_BAZAN
C-----Symbolic constants
      REAL SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      REAL RHO_L
      PARAMETER (RHO_L = RHO_LIQUID)
      REAL BETA_PAR
      PARAMETER (BETA_PAR = 8.2E0)
      REAL PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      REAL G_LAMBDA
      REAL G_VMIN
      REAL G_VMAX
      COMMON /C_MB_PARS/ G_LAMBDA, G_VMIN, G_VMAX
      
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Called functions
      REAL GK15
C-----Locale variables
      REAL D0
      REAL DC
      REAL DMIN
      REAL DMAX
      REAL VA
      REAL VB
C-----Code
      IF(ICLASS .EQ. 1 .AND. ICLASS .EQ. J) THEN
        GAMMA_IJ = 0.0E0
        RETURN
      ENDIF
      
      D0 = BUBBLE_CLASSES_DIA(J)
      
      DC = 
     *(12.0E0*SIGMA/(BETA_PAR*RHO_L))**(3.E0/5.E0) / EPS**(2.E0/5.E0)
      
C CHECK CHECK CHECK!!!
      IF(DC .GT. D0) THEN
        GAMMA_IJ = 0.0E0
        RETURN
      ENDIF
     
      DMIN = (12.0E0*SIGMA/(BETA_PAR*RHO_L
     *     * D0))**(3.E0/2.E0) / EPS
      DMAX = D0 * (1.0E0 - (DMIN
     *     / D0)**3.E0)**(1.E0/3.E0)
      G_VMIN = PI * DMIN**3.E0 / 6.E0
      G_VMAX = PI * DMAX**3.E0 / 6.E0
      G_LAMBDA = DC / D0
      
      IF(ICLASS .EQ. 1) THEN
        VA = BUBBLE_CLASSES_VOL(ICLASS)
        VB = BUBBLE_CLASSES_VOL(ICLASS+1)
                
        IF(VB .LT. G_VMIN .OR. VA .GT. G_VMAX) THEN
          GAMMA_IJ = 0.0E0
C CHECK CHECK CHECK!!!
        ELSEIF(VA .LE. G_VMIN .AND. VB .GE. G_VMAX) THEN
          GAMMA_IJ = 1.0E0
        ELSEIF(G_VMIN .LT. VA .AND. G_VMAX .GT. VB) THEN
          GAMMA_IJ = 
     *    GK15(XI_BETA, VA, VB, ICLASS, J, 0)
        ELSEIF(G_VMIN .GT. VA .AND. VB .GT. G_VMIN) THEN
          GAMMA_IJ = 
     *    GK15(XI_BETA, G_VMIN, VB, ICLASS, J, 0)
        ELSEIF(VB .GT. G_VMAX .AND. G_VMAX .GT. VA) THEN
          GAMMA_IJ = 
     *    GK15(XI_BETA, VA, G_VMAX, ICLASS, J, 0)
        ELSE
          WRITE(*,*) ('No solution - BETA')
          WRITE(*,*) (J)
          CALL ABORT()
        ENDIF
            
        RETURN
      ENDIF
      
      IF(ICLASS .EQ. J) THEN
        VA = BUBBLE_CLASSES_VOL(ICLASS-1)
        VB = BUBBLE_CLASSES_VOL(ICLASS)
        
        IF(VB .LT. G_VMIN .OR. VA .GT. G_VMAX) THEN
          GAMMA_IJ = 0.0E0
C CHECK CHECK CHECK!!!
        ELSEIF(VA .LE. G_VMIN .AND. VB .GE. G_VMAX) THEN
          GAMMA_IJ = 1.0E0
        ELSEIF(G_VMIN .LT. VA .AND. G_VMAX .GT. VB) THEN
          GAMMA_IJ = 
     *    GK15(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 0)
        ELSEIF(G_VMIN .GT. VA .AND. VB .GT. G_VMIN) THEN
          GAMMA_IJ = 
     *    GK15(XI_MINUS_ONE_BETA, G_VMIN, VB, ICLASS, J, 0)
        ELSEIF(VB .GT. G_VMAX .AND. G_VMAX .GT. VA) THEN
          GAMMA_IJ = 
     *    GK15(XI_MINUS_ONE_BETA, VA, G_VMAX, ICLASS, J, 0)
        ELSE
          WRITE(*,*) ('No solution - BETA')
          WRITE(*,*) (J)
          CALL ABORT()
        ENDIF
            
        RETURN
      ENDIF
      
      VA = BUBBLE_CLASSES_VOL(ICLASS)
      VB = BUBBLE_CLASSES_VOL(ICLASS+1)
                
      IF(VB .LT. G_VMIN .OR. VA .GT. G_VMAX) THEN
        GAMMA_IJ = 0.0E0
C CHECK CHECK CHECK!!!
      ELSEIF(VA .LE. G_VMIN .AND. VB .GE. G_VMAX) THEN
        GAMMA_IJ = 1.0E0
      ELSEIF(G_VMIN .LT. VA .AND. G_VMAX .GT. VB) THEN
        GAMMA_IJ = 
     *  GK15(XI_BETA, VA, VB, ICLASS, J, 0)
      ELSEIF(G_VMIN .GT. VA .AND. VB .GT. G_VMIN) THEN
        GAMMA_IJ = 
     *  GK15(XI_BETA, G_VMIN, VB, ICLASS, J, 0)
      ELSEIF(VB .GT. G_VMAX .AND. G_VMAX .GT. VA) THEN
        GAMMA_IJ = 
     *  GK15(XI_BETA, VA, G_VMAX, ICLASS, J, 0)
      ELSE
        WRITE(*,*) ('No solution - BETA')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
      
      VA = BUBBLE_CLASSES_VOL(ICLASS-1)
      VB = BUBBLE_CLASSES_VOL(ICLASS)
        
      IF(VB .LT. G_VMIN .OR. VA .GT. G_VMAX) THEN
c        GAMMA_IJ = GAMMA_IJ + 0.0E0
        RETURN
C CHECK CHECK CHECK!!!
      ELSEIF(VA .LE. G_VMIN .AND. VB .GE. G_VMAX) THEN
        GAMMA_IJ = GAMMA_IJ + 1.0E0
      ELSEIF(G_VMIN .LT. VA .AND. G_VMAX .GT. VB) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  GK15(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 0)
      ELSEIF(G_VMIN .GT. VA .AND. VB .GT. G_VMIN) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  GK15(XI_MINUS_ONE_BETA, G_VMIN, VB, ICLASS, J, 0)
      ELSEIF(VB .GT. G_VMAX .AND. G_VMAX .GT. VA) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  GK15(XI_MINUS_ONE_BETA, VA, G_VMAX, ICLASS, J, 0)
      ELSE
        WRITE(*,*) ('No solution - BETA')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
    
#else
#error "Unknown model specified"
#endif 
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
     
#ifdef DEBUG
      CALL CHECK_FINITE(BBRI, __LINE__)
#endif
      
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
#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong BAGI - ICLASS')
        WRITE(*,*) (ICLASS)
        CALL ABORT()
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
     *  *N(NLOC, ILOC, J, RALFA, RF)*N(NLOC, ILOC, K, RALFA, RF)*
     *   A_IJ(J, K, EPS(ILOC))
          
        END DO
      END DO 
  
#ifdef DEBUG   
      CALL CHECK_FINITE(BAGI, __LINE__)
#endif 
      END 
C=======================================================================
#ifdef MODEL_MARTINEZ_BAZAN
      REAL FUNCTION BETA_DENOMINATOR(I, V, J, BRANCH)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      REAL PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      REAL G_LAMBDA
      REAL G_VMIN
      REAL G_VMAX
      COMMON /C_MB_PARS/ G_LAMBDA, G_VMIN, G_VMAX
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments
      INTEGER J
      INTEGER I
      INTEGER BRANCH
      REAL V
C-----Locale variables
      REAL D0
C-----Code 
#ifdef DEBUG
      IF(I .NE. 0 .OR.
     *J .GT. NUMBER_OF_CLASSES .OR. J .LT. 1
     *) THEN
        WRITE(*,*) ('Wrong I .OR. J')
        WRITE(*,*) 'I=',I,'J=',J
        CALL ABORT()
      ENDIF
      
      IF(BRANCH .NE. 0) THEN
        WRITE(*,*) ('Wrong BRANCH')
        WRITE(*,*) (BRANCH)
        CALL ABORT()
      ENDIF
#endif

      D0 = BUBBLE_CLASSES_DIA(J)
           
      BETA_DENOMINATOR = (((6.E0*V/PI)**(1.E0/3.E0)
     *     / D0)**(2.E0/3.E0) - G_LAMBDA**(5.E0/3.E0))
     *     * ((1 - ((6.E0*V/PI)**(1.E0/3.E0)
     *     / D0)**3.E0)**(2.E0/9.E0)
     *     - G_LAMBDA**(5.E0/3.E0))
      
      END
#endif      
C=======================================================================
      REAL FUNCTION XI_BETA(I, V, J, BRANCH)
      IMPLICIT NONE
C-----Called functions
      REAL BETA
      REAL XI
C-----Arguments
      INTEGER J
      INTEGER I
      INTEGER BRANCH
      REAL V
      
      XI_BETA = XI(I,V)*BETA(J, V, BRANCH)
      END
C=======================================================================
      REAL FUNCTION XI_MINUS_ONE_BETA(I, V, J, BRANCH)
      IMPLICIT NONE
C-----Called functions
      REAL BETA
      REAL XI_MINUS_ONE
C-----Arguments
      INTEGER J
      INTEGER I
      INTEGER BRANCH
      REAL V
      
      XI_MINUS_ONE_BETA = XI_MINUS_ONE(I,V)*BETA(J, V, BRANCH)
      END
C=======================================================================
#ifdef MODEL_ALOPEA
      REAL FUNCTION BETA(J, V, BRANCH)
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
      INTEGER BRANCH

#ifdef DEBUG
      IF(J .GT. NUMBER_OF_CLASSES .OR. J .LT. 1) THEN
        WRITE(*,*) ('Wrong XI - J')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
      
      IF(BRANCH .NE. 0) THEN
        WRITE(*,*) ('Wrong BRANCH')
        WRITE(*,*) (BRANCH)
        CALL ABORT()
      ENDIF
#endif

      BETA = 60.0E0/BUBBLE_CLASSES_VOL(J) 
     * *(V/BUBBLE_CLASSES_VOL(J))**2
     * *(1.0E0 - V/BUBBLE_CLASSES_VOL(J))**2
     
#ifdef DEBUG    
      CALL CHECK_FINITE(BETA, __LINE__)
#endif
      
      END
C=======================================================================
#elif defined MODEL_LEHR
      REAL FUNCTION BETA(J, V, BRANCH)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      REAL SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      REAL RHO_L
      PARAMETER (RHO_L = (RHO_LIQUID))
c      REAL P_ERF
c      PARAMETER (P_ERF = 0.3275911E0)
      REAL PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
      REAL G_EPS
      COMMON /C_EPS/ G_EPS
C-----Arguments
      INTEGER J
      REAL V
      INTEGER BRANCH
C-----Locale variables
      REAL ERF_ARG
c      REAL T_ERF
      REAL ERF
      REAL V0
      REAL D0

#ifdef DEBUG                  
      IF(J .GT. NUMBER_OF_CLASSES .OR. J .LT. 1) THEN
        WRITE(*,*) ('Wrong XI - J')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
      
      IF(BRANCH .NE. 1 .AND. BRANCH .NE. 2) THEN
        WRITE(*,*) ('Wrong BRANCH')
        WRITE(*,*) (BRANCH)
        CALL ABORT()
      ENDIF
#endif

      V0 = BUBBLE_CLASSES_VOL(J) 
      D0 = BUBBLE_CLASSES_DIA(J)
      
      ERF_ARG = 3.E0/2.E0 * LOG(2.E0**(1.E0/15.E0) * D0
     *        * RHO_L**(3.E0/5.E0) * G_EPS**(2.E0/5.E0)
     *        / SIGMA**(3.E0/5.E0))
c      T_ERF = 1.E0 / (1.E0 + P_ERF*ERF_ARG)
c      ERF = 1.E0 - (0.254829592E0*T_ERF - 0.284496736E0*T_ERF**2.E0
c     *    + 1.421413741E0*T_ERF**3.E0 - 1.453152027*T_ERF**4.E0
c     *    + 1.061405429E0*T_ERF**5.E0) * EXP(-ERF_ARG**2.E0)

      ERF = DERF(ERF_ARG)

      IF(BRANCH .EQ. 1) THEN
      BETA = 
     *   1.E0/(SQRT(PI)*V) * EXP(-9.E0/4.E0 * (LOG(2.0E0**(2.E0/5.E0)
     *   * RHO_L**(3.E0/5.E0) * (6.E0*V/PI)**(1.E0/3.E0)
     *   * G_EPS**(2.E0/5.E0) / SIGMA**(3.E0/5.E0)))**2.E0) / (1.E0+ERF) 
      ELSEIF(BRANCH .EQ. 2) THEN
        BETA = 1.E0 / (SQRT(PI) * (V0 - V))
     *   * EXP(-9.E0/4.E0 * (LOG(2.0E0**(2.E0/5.E0) * RHO_L**(3.E0/5.E0)
     *   * (6.E0*(V0 - V)/PI)**(1.E0/3.E0)
     *   * G_EPS**(2.E0/5.E0) / SIGMA**(3.E0/5.E0)))**2.E0) / (1.E0+ERF)   
      ELSE
        WRITE(*,*) ('Wrong BRANCH')
        WRITE(*,*) (BRANCH)
        CALL ABORT()
      ENDIF
      
#ifdef DEBUG    
      CALL CHECK_FINITE(BETA, __LINE__)
#endif
       
      END
C=======================================================================
#elif defined MODEL_MARTINEZ_BAZAN
      REAL FUNCTION BETA(J, V, BRANCH)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      REAL PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      REAL G_LAMBDA
      REAL G_VMIN
      REAL G_VMAX
      COMMON /C_MB_PARS/ G_LAMBDA, G_VMIN, G_VMAX
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments
      INTEGER J
      REAL V
      INTEGER BRANCH
C-----Called functions
      REAL BETA_DENOMINATOR
      EXTERNAL BETA_DENOMINATOR
      REAL GK15
C-----Locale variables
      REAL D0
      REAL NOMINATOR
      REAL DENOMINATOR
C-----Code
#ifdef DEBUG                  
      IF(J .GT. NUMBER_OF_CLASSES .OR. J .LT. 1) THEN
        WRITE(*,*) ('Wrong XI - J')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
      
      IF(BRANCH .NE. 0) THEN
        WRITE(*,*) ('Wrong BRANCH')
        WRITE(*,*) (BRANCH)
        CALL ABORT()
      ENDIF
#endif

      D0 = BUBBLE_CLASSES_DIA(J)
      
      NOMINATOR = (((6.E0*V/PI)**(1.E0/3.E0)
     *     / D0)**(2.E0/3.E0) - G_LAMBDA**(5.E0/3.E0))
     *     * ((1 - ((6.E0*V/PI)**(1.E0/3.E0)
     *     / D0)**3.E0)**(2.E0/9.E0)
     *     - G_LAMBDA**(5.E0/3.E0))
     
      DENOMINATOR = GK15(BETA_DENOMINATOR, G_VMIN, G_VMAX, 0, J, 0)
          
      BETA = NOMINATOR / DENOMINATOR
      
#ifdef DEBUG    
      CALL CHECK_FINITE(BETA, __LINE__)
#endif
      END      
#else
#error "Unknown model specified"
#endif 
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
        WRITE(*,*) ('Wrong XI - I')
        WRITE(*,*) (I)
        CALL ABORT()
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
        WRITE(*,*) ('Wrong XI_MINUS_ONE - I')
        WRITE(*,*) (I)
        CALL ABORT()
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
      PARAMETER (PI = PI_CONST)
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

      REAL BREAKUP_FACTOR
      PARAMETER (BREAKUP_FACTOR = BREAKUP_F)
C-----Common blocks
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments   
      INTEGER I
      REAL EPS

#ifdef MODEL_ALOPEA
      REAL P_ERF
      PARAMETER (P_ERF = 0.3275911E0)
      REAL RHO_G
      PARAMETER (RHO_G = RHO_GAS)
      REAL MU_L
      PARAMETER (MU_L = MU_LIQUID)
      REAL ERF_ARG
      REAL T_ERF
      REAL ERF
C-----Code 
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
      
#elif defined MODEL_LEHR
C-----Code 
      IF(I .EQ. 1) THEN
        G_I = 0.0E0
      ELSE
        G_I = BREAKUP_FACTOR * 0.5E0 *BUBBLE_CLASSES_DIA(I)**(5.E0/3.E0)
     *      * EPS**(19.E0/15.E0) * RHO_L**(7.E0/5.E0)
     *      / SIGMA**(7.E0/5.E0) * EXP(-SQRT(2.0E0) * SIGMA**(9.E0/5.E0)
     *      /(BUBBLE_CLASSES_DIA(I)**3.E0 * RHO_L**(9.E0/5.E0)
     *      * EPS**(6.E0/5.E0)))
      ENDIF
      
#elif defined MODEL_MARTINEZ_BAZAN

C-----Symbolic constants
      REAL K_G
      PARAMETER (K_G = 0.25E0)
      REAL BETA_PAR
      PARAMETER (BETA_PAR = 8.2E0)
C-----Locale variables  
      REAL DISRUPT
      REAL CONFINE
C-----Code 
      DISRUPT = BETA_PAR * (EPS * BUBBLE_CLASSES_DIA(I))**(2.E0/3.E0)
      CONFINE = 12.E0 * SIGMA / (RHO_L * BUBBLE_CLASSES_DIA(I))
      
      IF(I .EQ. 1 .OR. CONFINE .GT. DISRUPT) THEN
        G_I = 0.0E0
      ELSE
        G_I = BREAKUP_FACTOR * K_G * SQRT(DISRUPT - CONFINE)
     *    / BUBBLE_CLASSES_DIA(I)
      ENDIF
      
#else
#error "Unknown model specified"
#endif      
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

#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong DBRI - ICLASS')
        WRITE(*,*) (ICLASS)
        CALL ABORT()
      ENDIF
#endif

      DBRI = N(NLOC, ILOC, ICLASS, RALFA, RF)*
     *G_I(ICLASS, EPS(ILOC))

#ifdef DEBUG    
      CALL CHECK_FINITE(DBRI, __LINE__)
#endif
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

#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong DAGI - ICLASS')
        WRITE(*,*) (ICLASS)
        CALL ABORT()
      ENDIF
#endif

      IF(ICLASS .EQ. NUMBER_OF_CLASSES) THEN
        RETURN
      ELSE
        DO J = 1, (NUMBER_OF_CLASSES - 1)
          DAGI = DAGI  
     *    + N(NLOC, ILOC, ICLASS, RALFA, RF)
     *      *N(NLOC, ILOC, J, RALFA, RF)*A_IJ(ICLASS, J, EPS(ILOC))
        END DO
      ENDIF

#ifdef DEBUG
      CALL CHECK_FINITE(DAGI, __LINE__)
#endif
      END 
C=======================================================================
#ifdef DEBUG
      SUBROUTINE CHECK_FINITE(X, LINE)
      IMPLICIT NONE
C-----Arguments
      REAL X
      INTEGER LINE


      IF(ISNAN(X) .OR. ABS(X) .GE. HUGE(X)) THEN
        WRITE(*,*) 'Variable is NOT a finite number:',X
        WRITE(*,*) 'Throw by line:', LINE
        CALL ABORT()
      ENDIF
           
      END
#endif
C=======================================================================
      BLOCKDATA
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12) 
C-----Locale variables
      REAL BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      REAL G_DBRI
      REAL G_DAGI
      REAL G_BBRI
      REAL G_BAGI
#ifdef MODEL_LEHR
      REAL G_EPS
#elif defined MODEL_MARTINEZ_BAZAN
      REAL G_LAMBDA
      REAL G_VMIN
      REAL G_VMAX
#endif
      
      DATA G_DBRI /0.0E0/
      DATA G_DAGI /0.0E0/
      DATA G_BBRI /0.0E0/
      DATA G_BAGI /0.0E0/
      
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
#ifdef MODEL_LEHR
      COMMON /C_EPS/ G_EPS
#elif defined MODEL_MARTINEZ_BAZAN
      COMMON /C_MB_PARS/ G_LAMBDA, G_VMIN, G_VMAX
#endif
      
      COMMON /C_DBRI/ G_DBRI
      COMMON /C_DAGI/ G_DAGI
      COMMON /C_BBRI/ G_BBRI
      COMMON /C_BAGI/ G_BAGI
            
      END