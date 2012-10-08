C=======================================================================
      EPS = ARGS(1,15)
      SIGMA = ARGS(1,16)
      RHO_L = ARGS(1,17)
      RHO_G = ARGS(1,18)
      MU_L = ARGS(1,19)
C=======================================================================
      REAL FUNCTION A_IJ(I, J, EPS, SIGMA, RHO_L)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments   
      INTEGER I
      INTEGER J
      REAL EPS
      REAL SIGMA
      REAL RHO_L
      REAL FREQ
      REAL EFF
      REAL RIJ
      REAL H0
      PARAMETER (H0 = 1.0E-4)
      REAL HF
      PARAMETER (HF = 1.0E-8)
      REAL COAL_FACTOR
      PARAMETER (COAL_FACTOR = 1.0E0)
      REAL PI
      PARAMETER (PI = 3.1415926535897931E0)
      
      RIJ = BUBBLE_CLASSES_DIA(I)*BUBBLE_CLASSES_DIA(J)
     *      / (BUBBLE_CLASSES_DIA(I)+BUBBLE_CLASSES_DIA(J)) 
      FREQ = SQRT(2)/4 * PI   
     *       *(BUBBLE_CLASSES_DIA(I)+BUBBLE_CLASSES_DIA(J))**2
     *       * EPS**(1/3) * (BUBBLE_CLASSES_DIA(I)**(2/3) 
     *       + BUBBLE_CLASSES_DIA(J)**(2/3))**(1/2)
      EFF = EXP(-SQRT(RHO_L) * R_IJ**(5/6) * EPS**(1/3) * LOG(H0/HF)
     *      /(4*SQRT(SIGMA)))
      
      A_IJ = COAL_FACTOR * FREQ * EFF
      
      END
C=======================================================================
      REAL FUNCTION G_I(I, EPS, SIGMA, RHO_L, RHO_G, MU_L)
      IMPLICIT NONE
C-----Symbolic constants
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Common blocks
      REAL BUBBLE_CLASSES_DIA(1:NUMBER_OF_CLASSES)
      COMMON /C_BUBBLE_CLASSES_DIA/ BUBBLE_CLASSES_DIA
C-----Arguments   
      INTEGER I
      REAL EPS
      REAL SIGMA
      REAL RHO_L
      REAL RHO_G
      REAL MU_L
      REAL ERF_ARG
      REAL T_ERF
      REAL ERF
      REAL COAL_FACTOR
      PARAMETER (BREAKUP_FACTOR = 1.0E0)
      REAL P_ERF
      PARAMETER (P_ERF = 0.3275911E0)

      ERF_ARG = SQRT(0.04 * SIGMA/(RHO_L * EPS**(2/3)
    *          * BUBBLE_CLASSES_DIA(I)**(5/3)) + 0.01 * MU_L
    *          / (SQRT(RHO_L*RHO_G) * EPS**(1/3)
    *          * BUBBLE_CLASSES_DIA(I)**(4/3)))
      T_ERF = 1 / (1 + P_ERF*ERF_ARG)
      ERF = 1 - (0.254829592*T_ERF - 0.284496736*T_ERF**2
     *      + 1.421413741*T_ERF**3 - 1.453152027*T_ERF**4
     *      + 1.061405429*T_ERF**5) * EXP(-ERF_ARG**2)
      
      G_I = BREAKUP_FACTOR * EPS**(1/3) * (1-ERF)
      
      END
C=======================================================================