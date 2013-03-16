#define DEBUG

C-------- used model
c#define MODEL_ALOPEA
c#define MODEL_LEHR
#define MODEL_MARTINEZ_BAZAN
C-------- used model

C-------- constants
#define SIGMA_SURF_TENS (0.0728D0)
#define RHO_GAS (1.185D0)
#define RHO_LIQUID (997.0D0)
#define MU_LIQUID (0.0008899D0)
#define PI_CONST (3.1415926535897931D0)
#define BREAKUP_F (1.0D0)
C-------- constants

C--------CFX double precision solver check
#ifdef cfd_version
#ifndef DOUBLE_PRECISION
#error "The code has to be compiled as double precison (-double)"
#endif
#endif
C--------CFX double precision solver check
      
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
#ifdef DEBUG
      DOUBLE PRECISION G_DBRI
      DOUBLE PRECISION G_DAGI
      DOUBLE PRECISION G_BBRI
      DOUBLE PRECISION G_BAGI
      
      COMMON /C_DBRI/ G_DBRI
      COMMON /C_DAGI/ G_DAGI
      COMMON /C_BBRI/ G_BBRI
      COMMON /C_BAGI/ G_BAGI
#endif
      
      DOUBLE PRECISION ARGS(NLOC,NARG), RET(NLOC,NRET)
      DOUBLE PRECISION TOTSUM
      
      TOTSUM = 0.0D0
      
      DO ICLASS = 1, NUMBER_OF_CLASSES
        ARGS(1,1) = ICLASS
        ARGS(1:NLOC,2) = 0.1D0
        
        ARGS(1:NLOC,3) = 0.1D0
        ARGS(1:NLOC,4) = 0.1D0
        ARGS(1:NLOC,5) = 0.1D0
        ARGS(1:NLOC,6) = 0.1D0
        ARGS(1:NLOC,7) = 0.1D0
        ARGS(1:NLOC,8) = 0.1D0
        ARGS(1:NLOC,9) = 0.1D0
        ARGS(1:NLOC,10) = 0.1D0
        ARGS(1:NLOC,11) = 0.05D0
        ARGS(1:NLOC,12) = 0.05D0
        ARGS(1:NLOC,13) = 0.05D0
        ARGS(1:NLOC,14) = 0.05D0
        ARGS(1:NLOC,15) = 1.0D0
        
        CALL BUBBLE_SOURCE(NLOC, NRET, NARG, RET, ARGS)
        TOTSUM = TOTSUM + RET(NLOC,NRET)
      ENDDO
#ifdef DEBUG
      WRITE(*,*) '==========================='
      WRITE(*,'(A, E20.10)') 'SUM_SOURCE: ', TOTSUM
      WRITE(*,'(A, E20.10)') 'BAGI-DAGI: ', (G_BAGI-G_DAGI)
      WRITE(*,'(A, E20.10)') 'BBRI-DBRI ', (G_BBRI-G_DBRI)
#endif   
      END

      SUBROUTINE BUBBLE_SOURCE( 
     &  NLOC, NRET, NARG, RET, ARGS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      DOUBLE PRECISION COMPUTE_SOURCE
C-----Arguments
      INTEGER NLOC,NARG,NRET
      DOUBLE PRECISION ARGS(NLOC,NARG), RET(NLOC,NRET)
      
c     ICLASS = ARGS(1,1)
c     ALPHA_G = ARGS(1:NLOC,2)

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
     
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION 
     *COMPUTE_SOURCE(NLOC, ILOC, ICLASS, ALPHA_G, RF, EPS)
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
      DOUBLE PRECISION VB(1:NUMBER_OF_CLASSES)
      COMMON /C_VB/ VB
      DOUBLE PRECISION G_DBRI
      DOUBLE PRECISION G_DAGI
      DOUBLE PRECISION G_BBRI
      DOUBLE PRECISION G_BAGI
      COMMON /C_DBRI/ G_DBRI
      COMMON /C_DAGI/ G_DAGI
      COMMON /C_BBRI/ G_BBRI
      COMMON /C_BAGI/ G_BAGI
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION ALPHA_G(NLOC)
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

      COMPUTE_SOURCE = 
     *  RHO_G * ALPHA_G(ILOC)**2.D0 *
     * (
     *  BAGI(NLOC, ILOC, ICLASS, RF, EPS)
     * -DAGI(NLOC, ILOC, ICLASS, RF, EPS)
     * )
     * +                    
     *  RHO_G * ALPHA_G(ILOC)* 
     * (
     *  BBRI(NLOC, ILOC, ICLASS, RF, EPS)    
     * -DBRI(NLOC, ILOC, ICLASS, RF, EPS)
     * )
     
#ifdef DEBUG
      G_BBRI = G_BBRI + RHO_G * ALPHA_G(ILOC) *
     *BBRI(NLOC, ILOC, ICLASS, RF, EPS)    
               
      G_BAGI = G_BAGI + RHO_G * ALPHA_G(ILOC)**2.D0 *
     *BAGI(NLOC, ILOC, ICLASS, RF, EPS)
               
      G_DBRI = G_DBRI + RHO_G * ALPHA_G(ILOC) *
     *DBRI(NLOC, ILOC, ICLASS, RF, EPS)
               
      G_DAGI = G_DAGI + RHO_G * ALPHA_G(ILOC)**2.D0 *
     *DAGI(NLOC, ILOC, ICLASS, RF, EPS) 
      
      WRITE(*,*) '---------------------------'
      WRITE(*,'(A, I0)') 'CLASS: ', ICLASS
      WRITE(*,'(A, F5.2)') 'VOLFRAC_G: ', ALPHA_G(ILOC)
      WRITE(*,'(A, I0, A, F5.2)') 'F', ICLASS, ':', RF(ILOC,ICLASS)
      WRITE(*,
     * '(A, E20.10)') 'BBRI: ', BBRI(NLOC, ILOC, ICLASS, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'BAGI: ', BAGI(NLOC, ILOC, ICLASS, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'DBRI: ',-DBRI(NLOC, ILOC, ICLASS, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'DAGI: ',-DAGI(NLOC, ILOC, ICLASS, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'BAGI-DAGI: ',
     * BAGI(NLOC, ILOC, ICLASS, RF, EPS) 
     *-DAGI(NLOC, ILOC, ICLASS, RF, EPS)
      WRITE(*,
     * '(A, E20.10)') 'BBRI-DBRI: ',
     * BBRI(NLOC, ILOC, ICLASS, RF, EPS) 
     *-DBRI(NLOC, ILOC, ICLASS, RF, EPS)
      WRITE(*,'(A, E20.10)') 'SOURCE: ', COMPUTE_SOURCE
#endif
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION G7(FCE, A, B, ICLASS, J)
      IMPLICIT NONE 
C-----Called functions
      DOUBLE PRECISION FCE
      EXTERNAL FCE 
C-----Arguments
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      INTEGER ICLASS
      INTEGER J
C-----Locale variables
      DOUBLE PRECISION TRANS1
      DOUBLE PRECISION TRANS2
      DOUBLE PRECISION INTEGRAL
      INTEGER L
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
      
      DO L = 0, 2
        INTEGRAL = INTEGRAL +
     *  WEIGHTS(L)*((FCE(ICLASS, TRANS1*NODES(L) + TRANS2, J)) +
     *  (FCE(ICLASS, TRANS1*(-NODES(L)) + TRANS2, J)))
      ENDDO

  
      G7 = 
     * TRANS1*(INTEGRAL +  WEIGHTS(3)*FCE(ICLASS, TRANS2, J))
      
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION GK15(FCE, A, B, ICLASS, J)
      IMPLICIT NONE 
C-----Called functions
      DOUBLE PRECISION FCE
      EXTERNAL FCE 
C-----Arguments
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      INTEGER ICLASS
      INTEGER J
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
     *  WEIGHTS(I)*((FCE(ICLASS, TRANS1*NODES(I) + TRANS2, J)) +
     *  (FCE(ICLASS, TRANS1*(-NODES(I)) + TRANS2, J)))
      ENDDO
  
  
      GK15 = 
     * TRANS1*(INTEGRAL +  WEIGHTS(7)*FCE(ICLASS, TRANS2, J))
      
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION G_IJ(I, J, EPS)
      IMPLICIT NONE
C-----Called functions
      DOUBLE PRECISION G_I
      DOUBLE PRECISION BETA
C-----Arguments
      INTEGER I
      INTEGER J
      DOUBLE PRECISION EPS
      
      G_IJ = G_I(I, EPS) * BETA(I, J, EPS)
      
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION BBRI(NLOC, ILOC, ICLASS, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION VB(1:NUMBER_OF_CLASSES)
      COMMON /C_VB/ VB
C-----Called functions
      DOUBLE PRECISION G_IJ
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION RF(1:NLOC, NUMBER_OF_CLASSES)
      DOUBLE PRECISION EPS(NLOC)
C-----Locale variables
      INTEGER J
C-----Code
      
      BBRI = 0.D0
      
      DO J = ICLASS + 1, NUMBER_OF_CLASSES
        BBRI = BBRI + G_IJ(J, ICLASS, EPS(ILOC)) * RF(ILOC, J)
      ENDDO
     
#ifdef DEBUG
      CALL CHECK_FINITE(BBRI, __LINE__)
#endif
      
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION BAGI(NLOC, ILOC, ICLASS, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION VB(1:NUMBER_OF_CLASSES)
      COMMON /C_VB/ VB
C-----Called functions
      DOUBLE PRECISION A_IJ
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
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

      DO J = 1, ICLASS
        DO K = 1, ICLASS
        
          V = VB(J) + VB(K)
          BRACKET_PRODUCT = 0.D0
          
          IF(ICLASS .NE. 1) THEN
            IF(VB(ICLASS - 1) .LT. V .AND.
     *        V .LT. VB(ICLASS)) THEN
             BRACKET_PRODUCT = (V - VB(ICLASS-1))
     *       /(VB(ICLASS) - VB(ICLASS-1))
            ENDIF
          ENDIF

c-----The following condition ensures that bubbles of the biggest class N do not coalesce. If one or two did, a bubble of size bigger than VB(N) would have to be formed;
c-----the following statement prevents this event. In order not to form a bubble bigger than VB(N) from smaller bubbles, it is necessary to pay attention to the restriction
c-----that 2*VB(N-1) < VB(N). This limitation is possible to remove simply by splitting a particle bigger than V(N) into a number of particles of size V(N) - BRACKET_PRODUCT = V/V(N). 

          IF(ICLASS .NE. NUMBER_OF_CLASSES) THEN
            IF(VB(ICLASS) .LT. V .AND.
     *        V .LT. VB(ICLASS + 1)) THEN
             BRACKET_PRODUCT = (VB(ICLASS + 1) - V)
     *       /(VB(ICLASS+1) - VB(ICLASS))
            ENDIF
          ENDIF
          
        BAGI = BAGI + BRACKET_PRODUCT * 0.5D0 
     *       * (VB(J) + VB(K)) / (VB(J) * VB(K))
     *       * RF(ILOC, J) * RF(ILOC, K) * A_IJ(J, K, EPS(ILOC))
          
        END DO
      END DO 
  
#ifdef DEBUG   
      CALL CHECK_FINITE(BAGI, __LINE__)
#endif 
      END 
C=======================================================================
#ifdef MODEL_ALOPEA
      DOUBLE PRECISION FUNCTION BETA(J, I, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION VB(1:NUMBER_OF_CLASSES)
      COMMON /C_VB/ VB
C-----Arguments
      INTEGER I
      INTEGER J
      DOUBLE PRECISION EPS

#ifdef DEBUG
      IF(J .GT. NUMBER_OF_CLASSES .OR. J .LT. 1) THEN
        WRITE(*,*) ('Wrong XI - J')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
#endif

      BETA = 60.D0 * (VB(I)/VB(J))**2.D0 * (1.D0 - VB(I)/VB(J))**2.D0
     
#ifdef DEBUG    
      CALL CHECK_FINITE(BETA, __LINE__)
#endif
      
      END
C=======================================================================
#elif defined MODEL_LEHR
      DOUBLE PRECISION FUNCTION BETA(J, I, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      DOUBLE PRECISION SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      DOUBLE PRECISION RHO_L
      PARAMETER (RHO_L = RHO_LIQUID)
      DOUBLE PRECISION P_ERF
      PARAMETER (P_ERF = 0.3275911D0)
      DOUBLE PRECISION PI
      PARAMETER (PI = PI_CONST)
C-----Common blocks
      DOUBLE PRECISION VB(1:NUMBER_OF_CLASSES)
      COMMON /C_VB/ VB
      DOUBLE PRECISION DB(1:NUMBER_OF_CLASSES)
      COMMON /C_DB/ DB
C-----Arguments
      INTEGER I
      INTEGER J
      DOUBLE PRECISION EPS
C-----Locale variables
      DOUBLE PRECISION ERF_ARG
      DOUBLE PRECISION T_ERF
      DOUBLE PRECISION ERF

#ifdef DEBUG                  
      IF(J .GT. NUMBER_OF_CLASSES .OR. J .LT. 1) THEN
        WRITE(*,*) ('Wrong XI - J')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
#endif
      
      ERF_ARG = 3.D0/2.D0 * DLOG(2.D0**(1.D0/15.D0) * DB(J)
     *     * RHO_L**(3.D0/5.D0) * EPS**(2.D0/5.D0) / SIGMA**(3.D0/5.D0))
            
      IF(ERF_ARG .GT. 0.D0) THEN         
      T_ERF = 1.D0 / (1.D0 + P_ERF*ERF_ARG)
      ERF = 1.D0 - (0.254829592D0*T_ERF - 0.284496736D0*T_ERF**2.D0
     *    + 1.421413741D0*T_ERF**3.D0 - 1.453152027*T_ERF**4.D0
     *    + 1.061405429D0*T_ERF**5.D0) * DEXP(-ERF_ARG**2.D0)
      ELSE
      T_ERF = 1.D0 / (1.D0 + P_ERF*(-ERF_ARG))
      ERF = -(1.D0 - (0.254829592D0*T_ERF - 0.284496736D0*T_ERF**2.D0
     *    + 1.421413741D0*T_ERF**3.D0 - 1.453152027*T_ERF**4.D0
     *    + 1.061405429D0*T_ERF**5.D0) * DEXP(-ERF_ARG**2.D0))
      ENDIF

      IF(VB(I) .GT. 0.D0 .AND. VB(I) .LT. VB(J)/2.D0) THEN
      BETA = 1.D0 / (DSQRT(PI) * VB(I)/VB(J))
     *   * DEXP(-9.D0/4.D0 * (DLOG(2.0D0**(2.D0/5.D0)
     *   * RHO_L**(3.D0/5.D0) * DB(I) * EPS**(2.D0/5.D0) 
     *   / SIGMA**(3.D0/5.D0)))**2.D0) / (1.D0 + ERF + 1.D-10) 
      ELSEIF(VB(I) .GT. VB(J)/2.D0 .AND. VB(I) .LT. VB(J)) THEN
        BETA = 1.D0 / (DSQRT(PI) * (VB(J) - VB(I)) / VB(J))
     *   * DEXP(-9.D0/4.D0 * (DLOG(2.D0**(2.D0/5.D0) *RHO_L**(3.D0/5.D0)
     *   * (6.D0*(VB(J) - VB(I))/PI)**(1.D0/3.D0) * EPS**(2.D0/5.D0) 
     *   / SIGMA**(3.D0/5.D0)))**2.D0) / (1.D0 + ERF + 1.D-10)   
      ELSE
        WRITE(*,*) ('Error - Lehr BETA')
        WRITE(*,*) (BETA)
        CALL ABORT()
      ENDIF
      
#ifdef DEBUG    
      CALL CHECK_FINITE(BETA, __LINE__)
#endif
       
      END
C=======================================================================
#elif defined MODEL_MARTINEZ_BAZAN
      DOUBLE PRECISION FUNCTION BETA(J, I, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      DOUBLE PRECISION SIGMA
      PARAMETER (SIGMA = SIGMA_SURF_TENS)
      DOUBLE PRECISION RHO_L
      PARAMETER (RHO_L = RHO_LIQUID)
      DOUBLE PRECISION BETA_PAR
      PARAMETER (BETA_PAR = 8.2D0)
C-----Common blocks
      DOUBLE PRECISION DMIN
      DOUBLE PRECISION DMAX
      DOUBLE PRECISION DC
      DOUBLE PRECISION G_LAMBDA
      DOUBLE PRECISION G_DSTARMIN
      DOUBLE PRECISION G_DSTARMAX
      COMMON /C_MB_PARS/ G_LAMBDA, G_DSTARMIN, G_DSTARMAX
      DOUBLE PRECISION DB(1:NUMBER_OF_CLASSES)
      COMMON /C_DB/ DB
      DOUBLE PRECISION VB(1:NUMBER_OF_CLASSES)
      COMMON /C_VB/ VB
C-----Arguments
      INTEGER I
      INTEGER J
      DOUBLE PRECISION EPS
C-----Called functions
      DOUBLE PRECISION BETA_DENOMINATOR
      EXTERNAL BETA_DENOMINATOR
      DOUBLE PRECISION G7
C-----Locale variables
      DOUBLE PRECISION NOMINATOR
      DOUBLE PRECISION DENOMINATOR
C-----Code
#ifdef DEBUG                  
      IF(J .GT. NUMBER_OF_CLASSES .OR. J .LT. 1) THEN
        WRITE(*,*) ('Wrong XI - J')
        WRITE(*,*) (J)
        CALL ABORT()
      ENDIF
#endif

      DC = (12.D0*SIGMA/(BETA_PAR*RHO_L))**(3.D0/5.D0) /EPS**(2.D0/5.D0)
      DMIN = (12.D0*SIGMA / (BETA_PAR*RHO_L*DB(J)))**(3.D0/2.D0) / EPS
      DMAX = DB(J) * (1.0D0 - (DMIN / DB(J))**3.D0)**(1.D0/3.D0)
      G_DSTARMIN = DMIN / DB(J)
      G_DSTARMAX = DMAX / DB(J)
      G_LAMBDA = DC / DB(J)
      
      IF(DB(I) .LT. DMIN .OR. DB(I) .GT. DMAX) THEN
      BETA = 0.D0
      ELSE      
      NOMINATOR = ((DB(I) / DB(J))**(2.D0/3.D0) - G_LAMBDA**(5.D0/3.D0))
     *  *((1.D0-(DB(I)/DB(J))**3.D0)**(2.D0/9.D0)-G_LAMBDA**(5.D0/3.D0))
     
      DENOMINATOR = G7(BETA_DENOMINATOR, G_DSTARMIN, G_DSTARMAX, 0, J)

C-----A treacherous feature of the M-B model is that for some combinations of low values of EPS (cca 0.2 - 0.5, in our case) and bubble diameter, DMIN > DMAX may occur.
C-----Simultaneously, the BETA curve is then located below the x-axis and the BETA values are therefore negative. Fortunately, in MUSIG method (in contrast with the CM),
C-----this behavior does not matter because if DMIN > DMAX and the curve is under the x-axis, the result of the integration is positive. 
          
      BETA = 2.D0 * NOMINATOR / DENOMINATOR
      ENDIF
      
#ifdef DEBUG    
      CALL CHECK_FINITE(BETA, __LINE__)
#endif
      END      
#else
#error "Unknown model specified"
#endif 
C=======================================================================
#ifdef MODEL_MARTINEZ_BAZAN
      DOUBLE PRECISION FUNCTION BETA_DENOMINATOR(I, DSTAR, J)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      DOUBLE PRECISION G_LAMBDA
      DOUBLE PRECISION G_DSTARMIN
      DOUBLE PRECISION G_DSTARMAX
      COMMON /C_MB_PARS/ G_LAMBDA, G_DSTARMIN, G_DSTARMAX
      DOUBLE PRECISION DB(1:NUMBER_OF_CLASSES)
      COMMON /C_DB/ DB
C-----Arguments
      INTEGER I
      INTEGER J
      DOUBLE PRECISION DSTAR
C-----Code 
#ifdef DEBUG
      IF(I .NE. 0 .OR.
     *J .GT. NUMBER_OF_CLASSES .OR. J .LT. 1
     *) THEN
        WRITE(*,*) ('Wrong I .OR. J')
        WRITE(*,*) 'I=',I,'J=',J
        CALL ABORT()
      ENDIF
#endif

      BETA_DENOMINATOR = (DSTAR**(2.D0/3.D0) - G_LAMBDA**(5.D0/3.D0))
     *   *((1.D0 - DSTAR**3.D0)**(2.D0/9.D0) - G_LAMBDA**(5.D0/3.D0))
      
      END
#endif      
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
      DOUBLE PRECISION DB(1:NUMBER_OF_CLASSES)
      COMMON /C_DB/ DB
C-----Arguments   
      INTEGER I
      INTEGER J
      DOUBLE PRECISION EPS
      DOUBLE PRECISION FREQ
      DOUBLE PRECISION EFF
      DOUBLE PRECISION R_IJ
      
      R_IJ = DB(I) * DB(J) / (DB(I) + DB(J)) 
      FREQ = DSQRT(2.D0) / 4.D0 * PI   
     *       * (DB(I) + DB(J))**2.D0 * EPS**(1.D0/3.D0)
     *       * (DB(I)**(2.D0/3.D0) + DB(J)**(2.D0/3.D0))**(1.D0/2.D0)
      EFF = DEXP(-DSQRT(RHO_L) * R_IJ**(5.D0/6.D0) * EPS**(1.D0/3.D0) 
     *      * DLOG(H0/HF) /(4.D0 * DSQRT(SIGMA)))    
      
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
      DOUBLE PRECISION DB(1:NUMBER_OF_CLASSES)
      COMMON /C_DB/ DB
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
        ERF_ARG = DSQRT(0.04D0 * SIGMA / (RHO_L * EPS**(2.D0/3.D0)
     *   * DB(I)**(5.D0/3.D0)) + 0.01D0 * MU_L
     *   / (DSQRT(RHO_L*RHO_G) * EPS**(1.D0/3.D0) * DB(I)**(4.D0/3.D0)))
        T_ERF = 1.D0 / (1.D0 + P_ERF*ERF_ARG)
        ERF = 1.D0 - (0.254829592D0*T_ERF - 0.284496736D0*T_ERF**2.D0
     *      + 1.421413741D0*T_ERF**3.D0 - 1.453152027*T_ERF**4.D0
     *      + 1.061405429D0*T_ERF**5.D0) * DEXP(-ERF_ARG**2.D0)
        
        G_I = BREAKUP_FACTOR * EPS**(1.D0/3.D0) * (1.D0 - ERF)
      ENDIF
      
#elif defined MODEL_LEHR
C-----Code 
      IF(I .EQ. 1) THEN
        G_I = 0.D0
      ELSE
        G_I = BREAKUP_FACTOR * 0.5D0 * DB(I)**(5.D0/3.D0)
     *      * EPS**(19.D0/15.D0) * RHO_L**(7.D0/5.D0)
     *      / SIGMA**(7.D0/5.D0) * DEXP(-DSQRT(2.D0) *SIGMA**(9.D0/5.D0)
     *      /(DB(I)**3.D0 * RHO_L**(9.D0/5.D0) * EPS**(6.D0/5.D0)))
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
      DISRUPT = BETA_PAR * (EPS * DB(I))**(2.D0/3.D0)
      CONFINE = 12.D0 * SIGMA / (RHO_L * DB(I))
      
      IF(I .EQ. 1 .OR. CONFINE .GT. DISRUPT) THEN
        G_I = 0.D0
      ELSE
        G_I = BREAKUP_FACTOR * K_G * DSQRT(DISRUPT - CONFINE)
     *      / DB(I)
      ENDIF
      
#else
#error "Unknown model specified"
#endif      
      END      
C=======================================================================
      DOUBLE PRECISION FUNCTION DBRI(NLOC, ILOC, ICLASS, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      DOUBLE PRECISION G_IJ
C-----Common blocks
      DOUBLE PRECISION VB(1:NUMBER_OF_CLASSES)
      COMMON /C_VB/ VB
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      DOUBLE PRECISION EPS(NLOC)
C-----Locale variables
      INTEGER J
C-----Code
      
      DBRI = 0.0D0
      
      DO J = 1, ICLASS - 1
        DBRI = DBRI + G_IJ(ICLASS, J, EPS(ILOC)) * RF(ILOC,ICLASS)   
      ENDDO
     
#ifdef DEBUG
      CALL CHECK_FINITE(DBRI, __LINE__)
#endif
      
      END
C=======================================================================
      DOUBLE PRECISION FUNCTION DAGI(NLOC, ILOC, ICLASS, RF, EPS)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Called functions
      DOUBLE PRECISION A_IJ
C-----Common blocks
      DOUBLE PRECISION VB(1:NUMBER_OF_CLASSES)
      COMMON /C_VB/ VB
C-----Arguments
      INTEGER NLOC
      INTEGER ILOC
      INTEGER ICLASS
      DOUBLE PRECISION RF(1:NLOC, 1:NUMBER_OF_CLASSES)
      DOUBLE PRECISION EPS(NLOC)
C-----Locale variables
      INTEGER J
C-----Code
      DAGI = 0.D0

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
        DO J = 1, NUMBER_OF_CLASSES - 1
          DAGI = DAGI + 1.D0 / VB(J)
     *     * RF(ILOC, ICLASS) * RF(ILOC, J) * A_IJ(ICLASS, J, EPS(ILOC))
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
        WRITE(*,*) 'Thrown by line:', LINE
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
      DOUBLE PRECISION VB(1:NUMBER_OF_CLASSES)
      DOUBLE PRECISION DB(1:NUMBER_OF_CLASSES)
      DOUBLE PRECISION G_DBRI
      DOUBLE PRECISION G_DAGI
      DOUBLE PRECISION G_BBRI
      DOUBLE PRECISION G_BAGI
#ifdef MODEL_MARTINEZ_BAZAN
      DOUBLE PRECISION G_LAMBDA
      DOUBLE PRECISION G_DSTARMIN
      DOUBLE PRECISION G_DSTARMAX
#endif
#ifdef DEBUG  
      DATA G_DBRI /0.0D0/
      DATA G_DAGI /0.0D0/
      DATA G_BBRI /0.0D0/
      DATA G_BAGI /0.0D0/
#endif
      
C-----Diameters of bubble classes
      DATA DB /0.5D-3, 1.0D-3, 2.0D-3, 3.0D-3, 4.0D-3
     * , 5.0D-3, 6.0D-3, 7.0D-3, 8.0D-3, 10.0D-3, 12.0D-3, 16.0D-3/
     
      DATA VB 
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
      COMMON /C_DB/ DB
      COMMON /C_VB/ VB
#ifdef MODEL_MARTINEZ_BAZAN
      COMMON /C_MB_PARS/ G_LAMBDA, G_DSTARMIN, G_DSTARMAX
#endif
#ifdef DEBUG  
      COMMON /C_DBRI/ G_DBRI
      COMMON /C_DAGI/ G_DAGI
      COMMON /C_BBRI/ G_BBRI
      COMMON /C_BAGI/ G_BAGI
#endif           
      END