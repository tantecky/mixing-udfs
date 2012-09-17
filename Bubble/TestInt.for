      PROGRAM TESTER
      IMPLICIT NONE
      
      REAL GK15
      REAL TESTINT
      EXTERNAL TESTINT
      
      
      WRITE(*,*) GK15(TESTINT, 0.003E0, 0.004E0, 1, 1)
      END
      
      REAL FUNCTION TESTINT(I, V)
      IMPLICIT NONE
C-----Symbolic constants
C-----Arguments
      INTEGER I
      REAL Vi,Vj, Vi1
      REAL V
      
      Vi = 0.004E0
      Vi1= 0.002E0
      Vj =  0.008E0
      
      TESTINT = 
     *(V-Vi1)/(Vi-Vi1) * 60E0/Vj * (V/Vj)**2E0 * (1.E0-V/Vj)**2E0
   
   
      
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
