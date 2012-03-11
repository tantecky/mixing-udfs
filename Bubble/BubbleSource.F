#include "cfx5ext.h"
dllexport(BUBBLE_SOURCE)
      SUBROUTINE BUBBLE_SOURCE (
     & NLOC,NRET,NARG,RET,ARGS,CRESLT,CZ,DZ,IZ,LZ,RZ)
C------------------------------
C        Argument list
      INTEGER NLOC, NRET, NARG
      CHARACTER CRESLT*(*)
      REAL RET(1:NLOC,1:NRET), ARGS(1:NLOC,1:NARG)

      INTEGER IZ(*)
      CHARACTER CZ(*)*(1)
      DOUBLE PRECISION DZ(*)
      LOGICAL LZ(*)
      REAL RZ(*)
C ------------------------------   

      CALL USER_PRINT_INTR('NLOC', NLOC)
      CALL USER_PRINT_INTR('NARG', NARG)
   
      
      END
