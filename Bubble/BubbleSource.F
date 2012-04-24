#include "cfx5ext.h"
dllexport(bubble_source)
      SUBROUTINE BUBBLE_SOURCE( 
     &  NLOC, NRET, NARG, RET, ARGS, CRESLT, CZ,DZ,IZ,LZ,RZ )
CC
CD User routine: template for user CEL function
CC
CC --------------------
CC        Input
CC --------------------
CC
CC  NLOC   - size of current locale
CC  NRET   - number of components in result
CC  NARG   - number of arguments in call
CC  ARGS() - (NLOC,NARG) argument values
CC
CC --------------------
CC      Modified
CC --------------------
CC
CC  Stacks possibly.
CC
CC --------------------
CC        Output
CC --------------------
CC
CC  RET()  - (NLOC,NRET) return values
CC  CRESLT - 'GOOD' for success
CC
CC --------------------
CC       Details
CC --------------------
CC
CC  Print 'Hello World!', and return value of argument.
CC
CC
CC======================================================================
C
C ------------------------------
C     Preprocessor includes
C ------------------------------
C
C
C ------------------------------
C        Global Parameters
C ------------------------------
C
C
C ------------------------------
C        Argument list
C ------------------------------
C
      INTEGER NLOC,NARG,NRET
C
      CHARACTER CRESLT*(*)
C
      REAL ARGS(NLOC,NARG), RET(NLOC,NRET)

C
      INTEGER IZ(*)
      CHARACTER CZ(*)*(1)
      DOUBLE PRECISION DZ(*)
      LOGICAL LZ(*)
      REAL RZ(*)
C
C ------------------------------
C        External routines
C ------------------------------
C
C
C ------------------------------
C        Local Parameters
C ------------------------------
C
C
C ------------------------------
C        Local Variables
C ------------------------------
C
      INTEGER ILOC
C
C ------------------------------
C        Stack pointers
C ------------------------------
C
C=======================================================================
C
C ---------------------------
C    Executable Statements
C ---------------------------
C
C
       CALL SET_A_0( RET, NLOC*NRET )


C       CALL USER_PRINT_INTR('NLOC', NLOC)
C      CALL USER_PRINT_INTR('NARG', NARG)
C      CALL USER_PRINT_INTR('NRET', NRET)
C       CALL USER_PRINT_INTR('ARGS', ARGS(NLOC,1))
       
       DO ILOC=1,NLOC
         CALL USER_PRINT_INTR('ARGS', ARGS(ILOC,1))
       END DO


C Set success flag.
      CRESLT = 'GOOD'
C
C=======================================================================
      END


