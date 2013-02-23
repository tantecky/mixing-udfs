#include "cfx5ext.h"

c#define DEBUG
C-------- used model
c#define MODEL_ALOPEA
#define MODEL_LEHR
c#define MODEL_MARTINEZ_BAZAN
C-------- used model

C-------- constants
#define SIGMA_SURF_TENS (0.0728D0)
#define RHO_GAS (1.185D0)
#define RHO_LIQUID (997.0D0)
#define MU_LIQUID (0.0008899D0)
#define PI_CONST (3.1415926535897931D0)
c#define BREAKUP_F (0.008D0)
#define BREAKUP_F (1.0D0)
C-------- constants

C--------CFX double precision solver check
#ifdef cfd_version
#ifndef DOUBLE_PRECISION
#error "The code has to be compiled as double precison (-double)"
#endif
#endif
C--------CFX double precision solver check

dllexport(bubble_source)
      SUBROUTINE BUBBLE_SOURCE( 
     &  NLOC, NRET, NARG, RET, ARGS, CRESLT,CZ,DZ,IZ,LZ,RZ )
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      DOUBLE PRECISION COMPUTE_SOURCE
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
c      DOUBLE PRECISION RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      
c      RF(:, 1) = ARGS(:,3)
c      RF(:, 2) = ARGS(:,4)
c      RF(:, 3) = ARGS(:,5)
c      RF(:, 4) = ARGS(:,6)
c      RF(:, 5) = ARGS(:,7)
c      RF(:, 6) = ARGS(:,8)
c      RF(:, 7) = ARGS(:,9)
c      RF(:, 8) = ARGS(:,10)
c      RF(:, 9) = ARGS(:,11)
c      RF(:, 10) = ARGS(:,12)
c      RF(:, 11) = ARGS(:,13)
c      RF(:, 12) = ARGS(:,14)
      
C-----Code
      ICLASS = INT(ARGS(1,1))
           
#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong ICLASS')
        WRITE(*,*) (ICLASS)
        CALL ABORT()
      ENDIF


      DO ILOC = 1, NLOC
        IF(ARGS(ILOC,2) .GT. 1.0D0 .OR. ARGS(ILOC,2)  .LT. 0.0D0) THEN
          WRITE(*,*) ('Wrong air.Volume Fraction')
          WRITE(*,*) (ARGS(ILOC,2))
          CALL ABORT()
        ENDIF
      END DO
#endif

      DO ILOC = 1, NLOC
        RET(ILOC,NRET) = COMPUTE_SOURCE(NLOC, ILOC, ICLASS, 
     *  ARGS(1,2), ARGS(:,3:14), ARGS(1,15))   
      END DO 

      CRESLT = 'GOOD'
     
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION 
     *COMPUTE_SOURCE(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      DOUBLE PRECISION RHO_G
      PARAMETER (RHO_G = RHO_GAS)
C-----Called functions
      DOUBLE PRECISION BBRI
      DOUBLE PRECISION BAGI
      DOUBLE PRECISION DBRI
      DOUBLE PRECISION DAGI
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION RALFA(NLOC)
      DOUBLE PRECISION RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      DOUBLE PRECISION EPS(NLOC)
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

      END
C=======================================================================
      DOUBLE PRECISION FUNCTION N(NLOC, ILOC, ICLASS, RALFA, RF)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION RALFA(NLOC)
      DOUBLE PRECISION RF(1:NLOC, 1:NUMBER_OF_CLASSES)
C-----Code

#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong N - ICLASS')
        WRITE(*,*) (ICLASS)
        CALL ABORT()
      ENDIF

      IF(RALFA(ILOC) .GT. 1.0D0 .OR. RALFA(ILOC)  .LT. 0.0D0) THEN
        WRITE(*,*) ('Wrong air.Volume Fraction - N')
        WRITE(*,*) (RALFA(ILOC))
        CALL ABORT()
      ENDIF
#endif

      N = RALFA(ILOC)*RF(ILOC,ICLASS)/BUBBLE_CLASSES_VOL(ICLASS)
      END
C=======================================================================
#ifdef MODEL_LEHR
      DOUBLE PRECISION FUNCTION INT_LEHR
     *                          (FCE, A, B, ICLASS, J, BRANCH, EPS)
      IMPLICIT NONE 
C-----Called functions
      DOUBLE PRECISION FCE
      EXTERNAL FCE 
      DOUBLE PRECISION GK15
      DOUBLE PRECISION GK61
C-----Arguments
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      INTEGER ICLASS
      INTEGER J
      INTEGER BRANCH
      DOUBLE PRECISION EPS
C-----Code
      IF(EPS .LT. 1.0D0) THEN
        INT_LEHR = GK15(FCE, A, B, ICLASS, J, BRANCH)
      ELSE
        INT_LEHR = GK61(FCE, A, B, ICLASS, J, BRANCH)
      ENDIF
      
      END
#endif 
C=======================================================================
      DOUBLE PRECISION FUNCTION G7(FCE, A, B, ICLASS, J, BRANCH)
      IMPLICIT NONE 
C-----Called functions
      DOUBLE PRECISION FCE
      EXTERNAL FCE 
C-----Arguments
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      INTEGER ICLASS
      INTEGER J
      INTEGER BRANCH
C-----Locale variables
      DOUBLE PRECISION TRANS1
      DOUBLE PRECISION TRANS2
      DOUBLE PRECISION INTEGRAL
      INTEGER I
      DOUBLE PRECISION NODES(0:2)
      DOUBLE PRECISION WEIGHTS(0:3)
      
      DATA NODES
     */0.949107912342759D0,
     * 0.741531185599394D0,
     * 0.405845151377397D0
     */
     
      DATA WEIGHTS
     */0.129484966168870D0,
     * 0.279705391489277D0,
     * 0.381830050505119D0,
     * 0.417959183673469D0
     */
C-----Code
      INTEGRAL = 0.0D0
      TRANS1 = (B-A)/2.0D0
      TRANS2 = (A+B)/2.0D0
      
#ifdef DEBUG      
      IF(A .GE. B) THEN
        WRITE(*,*) ('G7: A >= B')
        WRITE(*,*) 'A=',A
        WRITE(*,*) 'B=',B
        WRITE(*,*) 'TRANS1=',TRANS1
        WRITE(*,*) 'TRANS2=',TRANS2
        CALL ABORT()
      ENDIF
#endif
      
      DO I = 0, 2
        INTEGRAL = INTEGRAL +
     *  WEIGHTS(I)*((FCE(ICLASS, TRANS1*NODES(I) + TRANS2, J, BRANCH)) +
     *  (FCE(ICLASS, TRANS1*(-NODES(I)) + TRANS2, J, BRANCH)))
      ENDDO
  
  
      G7 = 
     * TRANS1*(INTEGRAL +  WEIGHTS(3)*FCE(ICLASS, TRANS2, J, BRANCH))
      
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION GK15(FCE, A, B, ICLASS, J, BRANCH)
      IMPLICIT NONE 
C-----Called functions
      DOUBLE PRECISION FCE
      EXTERNAL FCE 
C-----Arguments
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      INTEGER ICLASS
      INTEGER J
      INTEGER BRANCH
C-----Locale variables
      DOUBLE PRECISION TRANS1
      DOUBLE PRECISION TRANS2
      DOUBLE PRECISION INTEGRAL
      INTEGER I
      DOUBLE PRECISION NODES(0:6)
      DOUBLE PRECISION WEIGHTS(0:7)
      
      DATA NODES
     */0.991455371120813D0,
     * 0.949107912342759D0,
     * 0.864864423359769D0,
     * 0.741531185599394D0,
     * 0.586087235467691D0,
     * 0.405845151377397D0,
     * 0.207784955007898D0
     */
     
      DATA WEIGHTS
     */0.022935322010529D0,
     * 0.063092092629979D0,
     * 0.104790010322250D0,
     * 0.140653259715525D0,
     * 0.169004726639267D0,
     * 0.190350578064785D0,
     * 0.204432940075298D0,
     * 0.209482141084728D0
     */
C-----Code
      INTEGRAL = 0.0D0
      TRANS1 = (B-A)/2.0D0
      TRANS2 = (A+B)/2.0D0
      
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
      DOUBLE PRECISION FUNCTION GK31(FCE, A, B, ICLASS, J, BRANCH)
      IMPLICIT NONE 
C-----Called functions
      DOUBLE PRECISION FCE
      EXTERNAL FCE 
C-----Arguments
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      INTEGER ICLASS
      INTEGER J
      INTEGER BRANCH
C-----Locale variables
      DOUBLE PRECISION TRANS1
      DOUBLE PRECISION TRANS2
      DOUBLE PRECISION INTEGRAL
      INTEGER I
      DOUBLE PRECISION NODES(0:14)
      DOUBLE PRECISION WEIGHTS(0:15)
      
      DATA NODES
     */0.998002298693397D0,
     * 0.987992518020485D0,
     * 0.967739075679139D0,
     * 0.937273392400705D0,
     * 0.897264532344081D0,
     * 0.848206583410427D0,
     * 0.790418501442465D0,
     * 0.724417731360170D0,
     * 0.650996741297416D0,
     * 0.570972172608538D0,
     * 0.485081863640239D0,
     * 0.394151347077563D0,
     * 0.299180007153168D0,
     * 0.201194093997434D0,
     * 0.101142066918717D0               
     */
     
      DATA WEIGHTS
     */0.005377479872923D0,
     * 0.015007947329316D0,
     * 0.025460847326715D0,
     * 0.035346360791375D0,
     * 0.044589751324764D0,
     * 0.053481524690928D0,
     * 0.062009567800670D0,
     * 0.069854121318728D0,
     * 0.076849680757720D0,
     * 0.083080502823133D0,
     * 0.088564443056211D0,
     * 0.093126598170825D0,
     * 0.096642726983623D0,
     * 0.099173598721791D0,
     * 0.100769845523875D0,
     * 0.101330007014791D0                                        
     */
C-----Code
      INTEGRAL = 0.0D0
      TRANS1 = (B-A)/2.0D0
      TRANS2 = (A+B)/2.0D0
      
#ifdef DEBUG      
      IF(A .GE. B) THEN
        WRITE(*,*) ('GK31: A >= B')
        WRITE(*,*) 'A=',A
        WRITE(*,*) 'B=',B
        WRITE(*,*) 'TRANS1=',TRANS1
        WRITE(*,*) 'TRANS2=',TRANS2
        CALL ABORT()
      ENDIF
#endif
      
      DO I = 0, 14
        INTEGRAL = INTEGRAL +
     *  WEIGHTS(I)*((FCE(ICLASS, TRANS1*NODES(I) + TRANS2, J, BRANCH)) +
     *  (FCE(ICLASS, TRANS1*(-NODES(I)) + TRANS2, J, BRANCH)))
      ENDDO
  
  
      GK31 = 
     * TRANS1*(INTEGRAL +  WEIGHTS(15)*FCE(ICLASS, TRANS2, J, BRANCH))
      
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION GK61(FCE, A, B, ICLASS, J, BRANCH)
      IMPLICIT NONE 
C-----Called functions
      DOUBLE PRECISION FCE
      EXTERNAL FCE 
C-----Arguments
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      INTEGER ICLASS
      INTEGER J
      INTEGER BRANCH
C-----Locale variables
      DOUBLE PRECISION TRANS1
      DOUBLE PRECISION TRANS2
      DOUBLE PRECISION INTEGRAL
      INTEGER I
      DOUBLE PRECISION NODES(0:29)
      DOUBLE PRECISION WEIGHTS(0:30)
      
      DATA NODES
     */0.9994844100504906375713259D0,
     * 0.9968934840746495402716301D0,
     * 0.9916309968704045948586284D0,
     * 0.9836681232797472099700326D0,
     * 0.9731163225011262683746939D0,
     * 0.9600218649683075122168710D0,
     * 0.9443744447485599794158313D0,
     * 0.9262000474292743258793243D0,
     * 0.9055733076999077985465226D0,
     * 0.8825605357920526815431165D0,
     * 0.8572052335460610989586585D0,
     * 0.8295657623827683974428981D0,
     * 0.7997278358218390830136689D0,
     * 0.7677774321048261949179773D0,
     * 0.7337900624532268047261711D0,
     * 0.6978504947933157969322924D0,
     * 0.6600610641266269613700537D0,
     * 0.6205261829892428611404776D0,
     * 0.5793452358263616917560249D0,
     * 0.5366241481420198992641698D0,
     * 0.4924804678617785749936931D0,
     * 0.4470337695380891767806099D0,
     * 0.4004012548303943925354762D0,
     * 0.3527047255308781134710372D0,
     * 0.3040732022736250773726771D0,
     * 0.2546369261678898464398051D0,
     * 0.2045251166823098914389577D0,
     * 0.1538699136085835469637947D0,
     * 0.1028069379667370301470968D0,
     * 0.0514718425553176958330252D0/
     
      DATA WEIGHTS
     */0.0013890136986770076245516D0,
     * 0.0038904611270998840512672D0,
     * 0.0066307039159312921733198D0,
     * 0.0092732796595177634284411D0,
     * 0.0118230152534963417422329D0,
     * 0.0143697295070458048124514D0,
     * 0.0169208891890532726275723D0,
     * 0.0194141411939423811734090D0,
     * 0.0218280358216091922971675D0,
     * 0.0241911620780806013656864D0,
     * 0.0265099548823331016106017D0,
     * 0.0287540487650412928439788D0,
     * 0.0309072575623877624728843D0,
     * 0.0329814470574837260318142D0,
     * 0.0349793380280600241374997D0,
     * 0.0368823646518212292239111D0,
     * 0.0386789456247275929503487D0,
     * 0.0403745389515359591119953D0,
     * 0.0419698102151642461471475D0,
     * 0.0434525397013560693168317D0,
     * 0.0448148001331626631923556D0,
     * 0.0460592382710069881162717D0,
     * 0.0471855465692991539452615D0,
     * 0.0481858617570871291407795D0,
     * 0.0490554345550297788875282D0,
     * 0.0497956834270742063578116D0,
     * 0.0504059214027823468408931D0,
     * 0.0508817958987496064922975D0,
     * 0.0512215478492587721706563D0,
     * 0.0514261285374590259338629D0,
     * 0.0514947294294515675583404D0/
C-----Code
      INTEGRAL = 0.0D0
      TRANS1 = (B-A)/2.0D0
      TRANS2 = (A+B)/2.0D0
      
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
      DOUBLE PRECISION FUNCTION GAMMA_IJ(ICLASS, J, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Called functions
      DOUBLE PRECISION XI_BETA
      EXTERNAL XI_BETA
      DOUBLE PRECISION XI_MINUS_ONE_BETA
      EXTERNAL XI_MINUS_ONE_BETA
C-----Arguments
      INTEGER ICLASS
      INTEGER J
      DOUBLE PRECISION EPS
      
#ifdef MODEL_ALOPEA
C-----Called functions
      DOUBLE PRECISION G7
C-----Code
      IF(ICLASS .EQ. 1 .AND. ICLASS .EQ. J) THEN
        GAMMA_IJ = 0.0D0
        RETURN
      ENDIF
               
      IF(ICLASS .EQ. 1) THEN
        GAMMA_IJ = 
     *     G7(XI_BETA, BUBBLE_CLASSES_VOL(ICLASS), 
     *     BUBBLE_CLASSES_VOL(ICLASS+1), ICLASS, J, 0)     
        RETURN
      ENDIF
      
      IF(ICLASS .EQ. J) THEN
        GAMMA_IJ = 
     *     G7(XI_MINUS_ONE_BETA, BUBBLE_CLASSES_VOL(ICLASS-1), 
     *     BUBBLE_CLASSES_VOL(ICLASS), ICLASS, J, 0)
        RETURN
      ENDIF
      
      GAMMA_IJ =
     *G7(XI_MINUS_ONE_BETA, BUBBLE_CLASSES_VOL(ICLASS-1), 
     *     BUBBLE_CLASSES_VOL(ICLASS), ICLASS, J, 0)
     *     +
     *G7(XI_BETA, BUBBLE_CLASSES_VOL(ICLASS), 
     *     BUBBLE_CLASSES_VOL(ICLASS+1), ICLASS, J, 0)              

#elif defined MODEL_LEHR
C-----Common blocks
      DOUBLE PRECISION G_EPS
      COMMON /C_EPS/ G_EPS
C-----Locale variables
      DOUBLE PRECISION V0
      DOUBLE PRECISION V0HALF
      DOUBLE PRECISION VA
      DOUBLE PRECISION VB
C-----Called functions
      DOUBLE PRECISION INT_LEHR
C-----Code
      IF(ICLASS .EQ. 1 .AND. ICLASS .EQ. J) THEN
        GAMMA_IJ = 0.0D0
        RETURN
      ENDIF

      G_EPS = EPS
      V0 = BUBBLE_CLASSES_VOL(J)
      V0HALF = (V0/2.D0)
                 
      IF(ICLASS .EQ. 1) THEN
        VA = BUBBLE_CLASSES_VOL(ICLASS)
        VB = BUBBLE_CLASSES_VOL(ICLASS+1)
                
        IF(VA .LT. V0HALF .AND. VB .GT. V0HALF) THEN
          GAMMA_IJ = 
     *    INT_LEHR(XI_BETA, VA, V0HALF, ICLASS, J, 1, EPS)
          
          GAMMA_IJ = GAMMA_IJ +
     *    INT_LEHR(XI_BETA, V0HALF, VB, ICLASS, J, 2, EPS)
        ELSEIF(VA .GT. V0HALF) THEN
          GAMMA_IJ = 
     *    INT_LEHR(XI_BETA, VA, VB, ICLASS, J, 2, EPS) 
        ELSEIF(VB .LT. V0HALF) THEN
          GAMMA_IJ = 
     *    INT_LEHR(XI_BETA, VA, VB, ICLASS, J, 1, EPS) 
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
     *    INT_LEHR(XI_MINUS_ONE_BETA, VA, V0HALF, ICLASS, J, 1, EPS)
          
          GAMMA_IJ = GAMMA_IJ +
     *    INT_LEHR(XI_MINUS_ONE_BETA, V0HALF, VB, ICLASS, J, 2, EPS)
        ELSEIF(VA .GT. V0HALF) THEN
          GAMMA_IJ = 
     *    INT_LEHR(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 2, EPS) 
        ELSEIF(VB .LT. V0HALF) THEN
          GAMMA_IJ = 
     *    INT_LEHR(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 1, EPS) 
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
     *   INT_LEHR(XI_BETA, VA, V0HALF, ICLASS, J, 1, EPS)
          
         GAMMA_IJ = GAMMA_IJ +
     *   INT_LEHR(XI_BETA, V0HALF, VB, ICLASS, J, 2, EPS)
      ELSEIF(VA .GT. V0HALF) THEN
         GAMMA_IJ = 
     *   INT_LEHR(XI_BETA, VA, VB, ICLASS, J, 2, EPS) 
      ELSEIF(VB .LT. V0HALF) THEN
        GAMMA_IJ = 
     *   INT_LEHR(XI_BETA, VA, VB, ICLASS, J, 1, EPS) 
      ELSE
        WRITE(*,*) ('No solution - BETA')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
      
      VA = BUBBLE_CLASSES_VOL(ICLASS-1)
      VB = BUBBLE_CLASSES_VOL(ICLASS)
        
      IF(VA .LT. V0HALF .AND. VB .GT. V0HALF) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  INT_LEHR(XI_MINUS_ONE_BETA, VA, V0HALF, ICLASS, J, 1, EPS)
          
        GAMMA_IJ = GAMMA_IJ +
     *  INT_LEHR(XI_MINUS_ONE_BETA, V0HALF, VB, ICLASS, J, 2, EPS)
      ELSEIF(VA .GT. V0HALF) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  INT_LEHR(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 2, EPS) 
      ELSEIF(VB .LT. V0HALF) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  INT_LEHR(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 1, EPS) 
      ELSE
        WRITE(*,*) ('No solution - BETA')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF

#elif defined MODEL_MARTINEZ_BAZAN
C-----Symbolic constants
      DOUBLE PRECISION SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      DOUBLE PRECISION RHO_L
      PARAMETER (RHO_L = RHO_LIQUID)
      DOUBLE PRECISION BETA_PAR
      PARAMETER (BETA_PAR = 8.2D0)
      DOUBLE PRECISION PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      DOUBLE PRECISION G_LAMBDA
      DOUBLE PRECISION G_VMIN
      DOUBLE PRECISION G_VMAX
      COMMON /C_MB_PARS/ G_LAMBDA, G_VMIN, G_VMAX
      
      DOUBLE PRECISION BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Called functions
      DOUBLE PRECISION G7
C-----Locale variables
      DOUBLE PRECISION D0
      DOUBLE PRECISION DC
      DOUBLE PRECISION DMIN
      DOUBLE PRECISION DMAX
      DOUBLE PRECISION VA
      DOUBLE PRECISION VB
      DOUBLE PRECISION HELPVARIABLE
C-----Code
      IF(ICLASS .EQ. 1 .AND. ICLASS .EQ. J) THEN
        GAMMA_IJ = 0.0D0
        RETURN
      ENDIF
      
      D0 = BUBBLE_CLASSES_DIA(J)
      
      DC = 
     *(12.0D0*SIGMA/(BETA_PAR*RHO_L))**(3.D0/5.D0) / EPS**(2.D0/5.D0)
      
      IF(DC .GT. D0) THEN
        GAMMA_IJ = 0.0D0
        RETURN
      ENDIF
     
      DMIN = (12.0D0*SIGMA/(BETA_PAR*RHO_L
     *     * D0))**(3.D0/2.D0) / EPS
      DMAX = D0 * (1.0D0 - (DMIN
     *     / D0)**3.D0)**(1.D0/3.D0)
      G_VMIN = PI * DMIN**3.D0 / 6.D0
      G_VMAX = PI * DMAX**3.D0 / 6.D0
      G_LAMBDA = DC / D0
      
C-----A treacherous feature of the M-B model is that for some combinations of low values of EPS (cca 0.2 - 0.5, in our case) and bubble diameter, VMIN > VMAX may occur.
C-----Simultaneously, the integration BETA curve is located below x-axis. As a result, the conditions contained in if-tree are confused resulting mostly in GAMMA_IJ = 0.
C-----To remedy this behavior, VMIN and VMAX must be swapped when VMIN > VMAX. The result of the integration is then negative because of the BETA curve below x-axis. This
C-----is corrected via ABS(int(expression)) of all integrals.         
            
      IF(G_VMIN .GT. 0.D0 .AND. G_VMAX .GT. 0.D0 .AND.
     *  G_VMIN .GT. G_VMAX) THEN
        HELPVARIABLE = G_VMIN
        G_VMIN = G_VMAX
        G_VMAX = HELPVARIABLE
      ELSE
        G_VMIN = G_VMIN
        G_VMAX = G_VMAX
      ENDIF
            
      IF(ICLASS .EQ. 1) THEN
        VA = BUBBLE_CLASSES_VOL(ICLASS)
        VB = BUBBLE_CLASSES_VOL(ICLASS+1)
                
        IF(VB .LT. G_VMIN .OR. VA .GT. G_VMAX) THEN
          GAMMA_IJ = 0.0D0
        ELSEIF(VA .LE. G_VMIN .AND. VB .GE. G_VMAX) THEN
          GAMMA_IJ =
     *    DABS(G7(XI_BETA, G_VMIN, G_VMAX, ICLASS, J, 0))
        ELSEIF(VA .LT. G_VMIN .AND. VB .LT. G_VMAX) THEN
          GAMMA_IJ = 
     *    DABS(G7(XI_BETA, G_VMIN, VB, ICLASS, J, 0))
        ELSEIF(VA .GT. G_VMIN .AND. VB .GT. G_VMAX) THEN
          GAMMA_IJ = 
     *    DABS(G7(XI_BETA, VA, G_VMAX, ICLASS, J, 0))
        ELSEIF(VA .GT. G_VMIN .AND. VB .LT. G_VMAX) THEN
          GAMMA_IJ = 
     *    DABS(G7(XI_BETA, VA, VB, ICLASS, J, 0))
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
          GAMMA_IJ = 0.0D0
        ELSEIF(VA .LE. G_VMIN .AND. VB .GE. G_VMAX) THEN
          GAMMA_IJ =
     *    DABS(G7(XI_MINUS_ONE_BETA, G_VMIN, G_VMAX, ICLASS, J, 0))
        ELSEIF(VA .LT. G_VMIN .AND. VB .LT. G_VMAX) THEN
          GAMMA_IJ = 
     *    DABS(G7(XI_MINUS_ONE_BETA, G_VMIN, VB, ICLASS, J, 0))
        ELSEIF(VA .GT. G_VMIN .AND. VB .GT. G_VMAX) THEN
          GAMMA_IJ = 
     *    DABS(G7(XI_MINUS_ONE_BETA, VA, G_VMAX, ICLASS, J, 0))
        ELSEIF(VA .GT. G_VMIN .AND. VB .LT. G_VMAX) THEN
          GAMMA_IJ = 
     *    DABS(G7(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 0))
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
        GAMMA_IJ = 0.0D0
      ELSEIF(VA .LE. G_VMIN .AND. VB .GE. G_VMAX) THEN
        GAMMA_IJ =
     *  DABS(G7(XI_BETA, G_VMIN, G_VMAX, ICLASS, J, 0))
      ELSEIF(VA .LT. G_VMIN .AND. VB .LT. G_VMAX) THEN
        GAMMA_IJ = 
     *  DABS(G7(XI_BETA, G_VMIN, VB, ICLASS, J, 0))
      ELSEIF(VA .GT. G_VMIN .AND. VB .GT. G_VMAX) THEN
        GAMMA_IJ = 
     *  DABS(G7(XI_BETA, VA, G_VMAX, ICLASS, J, 0))
      ELSEIF(VA .GT. G_VMIN .AND. VB .LT. G_VMAX) THEN
        GAMMA_IJ = 
     *  DABS(G7(XI_BETA, VA, VB, ICLASS, J, 0))
      ELSE          
        WRITE(*,*) ('No solution - BETA')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
      
      VA = BUBBLE_CLASSES_VOL(ICLASS-1)
      VB = BUBBLE_CLASSES_VOL(ICLASS)
        
      IF(VB .LT. G_VMIN .OR. VA .GT. G_VMAX) THEN
c        GAMMA_IJ = GAMMA_IJ + 0.0D0
        RETURN
      ELSEIF(VA .LE. G_VMIN .AND. VB .GE. G_VMAX) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  DABS(G7(XI_MINUS_ONE_BETA, G_VMIN, G_VMAX, ICLASS, J, 0))
      ELSEIF(VA .LT. G_VMIN .AND. VB .LT. G_VMAX) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  DABS(G7(XI_MINUS_ONE_BETA, G_VMIN, VB, ICLASS, J, 0))
      ELSEIF(VA .GT. G_VMIN .AND. VB .GT. G_VMAX) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  DABS(G7(XI_MINUS_ONE_BETA, VA, G_VMAX, ICLASS, J, 0))
      ELSEIF(VA .GT. G_VMIN .AND. VB .LT. G_VMAX) THEN
        GAMMA_IJ = GAMMA_IJ +
     *  DABS(G7(XI_MINUS_ONE_BETA, VA, VB, ICLASS, J, 0))
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
      DOUBLE PRECISION FUNCTION BBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Called functions
      DOUBLE PRECISION N
      DOUBLE PRECISION GAMMA_IJ
      DOUBLE PRECISION G_I
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION RALFA(NLOC)
      DOUBLE PRECISION RF(1:NLOC, NUMBER_OF_CLASSES)
      DOUBLE PRECISION EPS(NLOC)
C-----Locale variables
      INTEGER J
C-----Code
      
      BBRI = 0.0D0
      
      DO J = ICLASS, NUMBER_OF_CLASSES
        BBRI = BBRI + G_I(J, EPS(ILOC))*
     *  GAMMA_IJ(ICLASS, J, EPS(ILOC))*N(NLOC, ILOC, J, RALFA, RF)
      ENDDO
     
#ifdef DEBUG
      CALL CHECK_FINITE(BBRI, __LINE__)
#endif
      
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION BAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Called functions
      DOUBLE PRECISION N
      DOUBLE PRECISION XI
      DOUBLE PRECISION XI_MINUS_ONE
      INTEGER KRONECKER_D
      DOUBLE PRECISION A_IJ
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION RALFA(NLOC)
      DOUBLE PRECISION RF(1:NLOC, NUMBER_OF_CLASSES)
      DOUBLE PRECISION EPS(NLOC)
C-----Locale variables
      INTEGER J
      INTEGER K
      DOUBLE PRECISION V
      DOUBLE PRECISION BRACKET_PRODUCT
C-----Code
#ifdef DEBUG
      IF(ICLASS .GT. NUMBER_OF_CLASSES .OR. ICLASS .LT. 1) THEN
        WRITE(*,*) ('Wrong BAGI - ICLASS')
        WRITE(*,*) (ICLASS)
        CALL ABORT()
      ENDIF
#endif

      BAGI = 0.0D0

      DO J = 1, (NUMBER_OF_CLASSES - 1)
        DO K = J, (NUMBER_OF_CLASSES - 1)
        
          V = BUBBLE_CLASSES_VOL(J)+BUBBLE_CLASSES_VOL(K)
          BRACKET_PRODUCT = 0.0D0
          
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
     *  BRACKET_PRODUCT*(1.0D0 - 0.5D0*KRONECKER_D(J,K))
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
      DOUBLE PRECISION FUNCTION BETA_DENOMINATOR(I, V, J, BRANCH)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      DOUBLE PRECISION PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      DOUBLE PRECISION G_LAMBDA
      DOUBLE PRECISION G_VMIN
      DOUBLE PRECISION G_VMAX
      COMMON /C_MB_PARS/ G_LAMBDA, G_VMIN, G_VMAX
      DOUBLE PRECISION BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments
      INTEGER J
      INTEGER I
      INTEGER BRANCH
      DOUBLE PRECISION V
C-----Locale variables
      DOUBLE PRECISION D0
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
           
      BETA_DENOMINATOR = (((6.D0*V/PI)**(1.D0/3.D0)
     *     / D0)**(2.D0/3.D0) - G_LAMBDA**(5.D0/3.D0))
     *     * ((1 - ((6.D0*V/PI)**(1.D0/3.D0)
     *     / D0)**3.D0)**(2.D0/9.D0)
     *     - G_LAMBDA**(5.D0/3.D0))
      
      END
#endif      
C=======================================================================
      DOUBLE PRECISION FUNCTION XI_BETA(I, V, J, BRANCH)
      IMPLICIT NONE
C-----Called functions
      DOUBLE PRECISION BETA
      DOUBLE PRECISION XI
C-----Arguments
      INTEGER J
      INTEGER I
      INTEGER BRANCH
      DOUBLE PRECISION V
      
      XI_BETA = XI(I,V)*BETA(J, V, BRANCH)
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION XI_MINUS_ONE_BETA(I, V, J, BRANCH)
      IMPLICIT NONE
C-----Called functions
      DOUBLE PRECISION BETA
      DOUBLE PRECISION XI_MINUS_ONE
C-----Arguments
      INTEGER J
      INTEGER I
      INTEGER BRANCH
      DOUBLE PRECISION V
      
      XI_MINUS_ONE_BETA = XI_MINUS_ONE(I,V)*BETA(J, V, BRANCH)
      END
C=======================================================================
#ifdef MODEL_ALOPEA
      DOUBLE PRECISION FUNCTION BETA(J, V, BRANCH)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER J
      DOUBLE PRECISION V
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

      BETA = 60.0D0/BUBBLE_CLASSES_VOL(J) 
     * *(V/BUBBLE_CLASSES_VOL(J))**2
     * *(1.0D0 - V/BUBBLE_CLASSES_VOL(J))**2
     
#ifdef DEBUG    
      CALL CHECK_FINITE(BETA, __LINE__)
#endif
      
      END
C=======================================================================
#elif defined MODEL_LEHR
      DOUBLE PRECISION FUNCTION BETA(J, V, BRANCH)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      DOUBLE PRECISION SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      DOUBLE PRECISION RHO_L
      PARAMETER (RHO_L = (RHO_LIQUID))
      DOUBLE PRECISION P_ERF
      PARAMETER (P_ERF = 0.3275911D0)
      DOUBLE PRECISION PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
      DOUBLE PRECISION BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
      DOUBLE PRECISION G_EPS
      COMMON /C_EPS/ G_EPS
C-----Arguments
      INTEGER J
      DOUBLE PRECISION V
      INTEGER BRANCH
C-----Locale variables
      DOUBLE PRECISION ERF_ARG
      DOUBLE PRECISION T_ERF
      DOUBLE PRECISION ERF
      DOUBLE PRECISION V0
      DOUBLE PRECISION D0

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
      
      ERF_ARG = 3.D0/2.D0 * DLOG(2.D0**(1.D0/15.D0) * D0
     *        * RHO_L**(3.D0/5.D0) * G_EPS**(2.D0/5.D0)
     *        / SIGMA**(3.D0/5.D0))
      T_ERF = 1.D0 / (1.D0 + P_ERF*ERF_ARG)
      ERF = 1.D0 - (0.254829592D0*T_ERF - 0.284496736D0*T_ERF**2.D0
     *    + 1.421413741D0*T_ERF**3.D0 - 1.453152027*T_ERF**4.D0
     *    + 1.061405429D0*T_ERF**5.D0) * DEXP(-ERF_ARG**2.D0)

      IF(BRANCH .EQ. 1) THEN
      BETA = 
     *   1.D0/(DSQRT(PI)*V) * DEXP(-9.D0/4.D0 * (DLOG(2.0D0**(2.D0/5.D0)
     *   * RHO_L**(3.D0/5.D0) * (6.D0*V/PI)**(1.D0/3.D0)
     *   * G_EPS**(2.D0/5.D0) 
     *   / SIGMA**(3.D0/5.D0)))**2.D0) / (1.D0+ERF) 
      ELSEIF(BRANCH .EQ. 2) THEN
        BETA = 1.D0 / (DSQRT(PI) * (V0 - V))
     *   * DEXP(-9.D0/4.D0 * (DLOG(2.0D0**(2.D0/5.D0)*RHO_L**(3.D0/5.D0)
     *   * (6.D0*(V0 - V)/PI)**(1.D0/3.D0)
     *   * G_EPS**(2.D0/5.D0) 
     *   / SIGMA**(3.D0/5.D0)))**2.D0) / (1.D0+ERF)   
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
      DOUBLE PRECISION FUNCTION BETA(J, V, BRANCH)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      DOUBLE PRECISION PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      DOUBLE PRECISION G_LAMBDA
      DOUBLE PRECISION G_VMIN
      DOUBLE PRECISION G_VMAX
      COMMON /C_MB_PARS/ G_LAMBDA, G_VMIN, G_VMAX
      DOUBLE PRECISION BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments
      INTEGER J
      DOUBLE PRECISION V
      INTEGER BRANCH
C-----Called functions
      DOUBLE PRECISION BETA_DENOMINATOR
      EXTERNAL BETA_DENOMINATOR
      DOUBLE PRECISION G7
C-----Locale variables
      DOUBLE PRECISION D0
      DOUBLE PRECISION NOMINATOR
      DOUBLE PRECISION DENOMINATOR
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
      
      NOMINATOR = (((6.D0*V/PI)**(1.D0/3.D0)
     *     / D0)**(2.D0/3.D0) - G_LAMBDA**(5.D0/3.D0))
     *     * ((1 - ((6.D0*V/PI)**(1.D0/3.D0)
     *     / D0)**3.D0)**(2.D0/9.D0)
     *     - G_LAMBDA**(5.D0/3.D0))
     
      DENOMINATOR = G7(BETA_DENOMINATOR, G_VMIN, G_VMAX, 0, J, 0)
          
      BETA = 2.D0 * NOMINATOR / DENOMINATOR
      
#ifdef DEBUG    
      CALL CHECK_FINITE(BETA, __LINE__)
#endif
      END      
#else
#error "Unknown model specified"
#endif 
C=======================================================================
      DOUBLE PRECISION FUNCTION XI(I, V)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER I
      DOUBLE PRECISION V
      
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
      DOUBLE PRECISION FUNCTION XI_MINUS_ONE(I, V)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER I
      DOUBLE PRECISION V
      
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
      DOUBLE PRECISION FUNCTION A_IJ(I, J, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      DOUBLE PRECISION SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      DOUBLE PRECISION RHO_L
      PARAMETER (RHO_L = RHO_LIQUID)
      DOUBLE PRECISION H0
      PARAMETER (H0 = 1.0D-4)
      DOUBLE PRECISION HF
      PARAMETER (HF = 1.0D-8)
      DOUBLE PRECISION COAL_FACTOR
      PARAMETER (COAL_FACTOR = 1.0D0)
      DOUBLE PRECISION PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments   
      INTEGER I
      INTEGER J
      DOUBLE PRECISION EPS
      DOUBLE PRECISION FREQ
      DOUBLE PRECISION EFF
      DOUBLE PRECISION R_IJ
      
      R_IJ = BUBBLE_CLASSES_DIA(I)*BUBBLE_CLASSES_DIA(J)
     *      / (BUBBLE_CLASSES_DIA(I)+BUBBLE_CLASSES_DIA(J)) 
      FREQ = DSQRT(2.D0)/4.D0 * PI   
     *       *(BUBBLE_CLASSES_DIA(I)+BUBBLE_CLASSES_DIA(J))**2.D0
     *       * EPS**(1.D0/3.D0) * (BUBBLE_CLASSES_DIA(I)**(2.D0/3.D0) 
     *       + BUBBLE_CLASSES_DIA(J)**(2.D0/3.D0))**(1.D0/2.D0)
      EFF = DEXP(-DSQRT(RHO_L) * R_IJ**(5.D0/6.D0) * EPS**(1.D0/3.D0) 
     *   * DLOG(H0/HF) /(4.D0*DSQRT(SIGMA)))    
      
      A_IJ = COAL_FACTOR * FREQ * EFF
      
      END 
C=======================================================================
      DOUBLE PRECISION FUNCTION G_I(I, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      DOUBLE PRECISION SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      DOUBLE PRECISION RHO_L
      PARAMETER (RHO_L = RHO_LIQUID)

      DOUBLE PRECISION BREAKUP_FACTOR
      PARAMETER (BREAKUP_FACTOR = BREAKUP_F)
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments   
      INTEGER I
      DOUBLE PRECISION EPS

#ifdef MODEL_ALOPEA
      DOUBLE PRECISION P_ERF
      PARAMETER (P_ERF = 0.3275911D0)
      DOUBLE PRECISION RHO_G
      PARAMETER (RHO_G = RHO_GAS)
      DOUBLE PRECISION MU_L
      PARAMETER (MU_L = MU_LIQUID)
      DOUBLE PRECISION ERF_ARG
      DOUBLE PRECISION T_ERF
      DOUBLE PRECISION ERF
C-----Code 
      IF(I .EQ. 1) THEN
        G_I = 0.D0
      ELSE
        ERF_ARG = DSQRT(0.04D0 * SIGMA/(RHO_L * EPS**(2.D0/3.D0)
     *          * BUBBLE_CLASSES_DIA(I)**(5.D0/3.D0)) + 0.01D0 * MU_L
     *          / (DSQRT(RHO_L*RHO_G) * EPS**(1.D0/3.D0)
     *          * BUBBLE_CLASSES_DIA(I)**(4.D0/3.D0)))
        T_ERF = 1.D0 / (1.D0 + P_ERF*ERF_ARG)
        ERF = 1.D0 - (0.254829592D0*T_ERF - 0.284496736D0*T_ERF**2.D0
     *      + 1.421413741D0*T_ERF**3.D0 - 1.453152027*T_ERF**4.D0
     *      + 1.061405429D0*T_ERF**5.D0) * DEXP(-ERF_ARG**2.D0)
        
        G_I = BREAKUP_FACTOR * 2.52 * EPS**(1.D0/3.D0) * (1.0D0-ERF)
      ENDIF
      
#elif defined MODEL_LEHR
C-----Code 
      IF(I .EQ. 1) THEN
        G_I = 0.0D0
      ELSE
        G_I = BREAKUP_FACTOR * 0.5D0 *BUBBLE_CLASSES_DIA(I)**(5.D0/3.D0)
     *      * EPS**(19.D0/15.D0) * RHO_L**(7.D0/5.D0)
     *      / SIGMA**(7.D0/5.D0) * DEXP(-DSQRT(2.0D0)*SIGMA**(9.D0/5.D0)
     *      /(BUBBLE_CLASSES_DIA(I)**3.D0 * RHO_L**(9.D0/5.D0)
     *      * EPS**(6.D0/5.D0)))
      ENDIF
      
#elif defined MODEL_MARTINEZ_BAZAN

C-----Symbolic constants
      DOUBLE PRECISION K_G
      PARAMETER (K_G = 0.25D0)
      DOUBLE PRECISION BETA_PAR
      PARAMETER (BETA_PAR = 8.2D0)
C-----Locale variables  
      DOUBLE PRECISION DISRUPT
      DOUBLE PRECISION CONFINE
C-----Code 
      DISRUPT = BETA_PAR * (EPS * BUBBLE_CLASSES_DIA(I))**(2.D0/3.D0)
      CONFINE = 12.D0 * SIGMA / (RHO_L * BUBBLE_CLASSES_DIA(I))
      
      IF(I .EQ. 1 .OR. CONFINE .GT. DISRUPT) THEN
        G_I = 0.0D0
      ELSE
        G_I = BREAKUP_FACTOR * K_G * DSQRT(DISRUPT - CONFINE)
     *    / BUBBLE_CLASSES_DIA(I)
      ENDIF
      
#else
#error "Unknown model specified"
#endif      
      END      
C=======================================================================
      DOUBLE PRECISION FUNCTION DBRI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      DOUBLE PRECISION N
      DOUBLE PRECISION G_I
C-----Common blocks
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_VOL/ BUBBLE_CLASSES_VOL
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION RALFA(NLOC)
      DOUBLE PRECISION RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      DOUBLE PRECISION EPS(NLOC)
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
      DOUBLE PRECISION FUNCTION DAGI(NLOC, ILOC, ICLASS, RALFA, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      DOUBLE PRECISION N
      DOUBLE PRECISION A_IJ
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION RALFA(NLOC)
      DOUBLE PRECISION RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      DOUBLE PRECISION EPS(NLOC)
C-----Locale variables
      INTEGER J
C-----Code
      DAGI = 0.0D0

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
      DOUBLE PRECISION X
      INTEGER LINE


      IF(ISNAN(X) .OR. DABS(X) .GE. HUGE(X)) THEN
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
      DOUBLE PRECISION BUBBLE_CLASSES_VOL(1:NUMBER_OF_CLASSES)
      DOUBLE PRECISION BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
#ifdef MODEL_LEHR
      DOUBLE PRECISION G_EPS
#elif defined MODEL_MARTINEZ_BAZAN
      DOUBLE PRECISION G_LAMBDA
      DOUBLE PRECISION G_VMIN
      DOUBLE PRECISION G_VMAX
#endif    
      
C     diameter of bubble classes
      DATA BUBBLE_CLASSES_DIA /0.5D-3, 1.0D-3, 2.0D-3, 3.0D-3, 4.0D-3
     * , 5.0D-3, 6.0D-3, 7.0D-3, 8.0D-3, 10.0D-3, 12.0D-3, 16.0D-3/
     
      DATA BUBBLE_CLASSES_VOL 
     */ 6.54498469497872D-11
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
#ifdef MODEL_LEHR
      COMMON /C_EPS/ G_EPS
#elif defined MODEL_MARTINEZ_BAZAN
      COMMON /C_MB_PARS/ G_LAMBDA, G_VMIN, G_VMAX
#endif
      END