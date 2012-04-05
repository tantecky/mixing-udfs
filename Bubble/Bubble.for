#include "cfx5ext.h"
#define __stack_point__ INTEGER*8
dllexport(bubble_source)
      SUBROUTINE BUBBLE_SOURCE( 
     &  NLOC, NRET, NARG, RET, ARGS, CRESLT, CZ,DZ,IZ,LZ,RZ )
      IMPLICIT NONE

      INTEGER NLOC,NARG,NRET
      

C
      CHARACTER CRESLT*(*)
C
      REAL ARGS(NLOC,NARG), RET(NLOC,NRET)
C
      __stack_point__ pVar
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
      REAL ComputeSource
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
      INTEGER bClass
      bClass = INT(ARGS(1,1))
      
      RET(1,1) = ComputeSource(bClass, CRESLT, pVar,CZ,DZ,IZ,LZ,RZ)
C
C Set success flag.
      CRESLT = 'GOOD'
C
C=======================================================================
      END






cdllexport(bubble_source)
c      SUBROUTINE BUBBLE_SOURCE (
c     & NLOC,NRET,NARG,RET,ARGS,CRESLT,CZ,DZ,IZ,LZ,RZ)
c      REAL ComputeSource
cC------------------------------
cC        Argument list
cC      NLOC: Number of locations in space 
cC            over which the calculations have to be performed.
c      INTEGER NLOC
cC      NRET: Number of return variables. This is always 1 in CFX,
cC            but is included to allow future extensions.
c      INTEGER NRET
cC      NARG: Number of arguments passed to the function.
c      INTEGER NARG
      
c      CHARACTER CRESLT*(*)
cC      RET(1:NLOC,1:NRET): Return variables (at each point in space).
c      REAL RET(1:NLOC,1:NRET), ARGS(1:NLOC,1:NARG)

c      __stack_point__ pVar
c      INTEGER IZ(*)
c      CHARACTER CZ(*)*(1)
c      DOUBLE PRECISION DZ(*)
c      LOGICAL LZ(*)
c      REAL RZ(*)
cC ------------------------------  
c      INTEGER bClass
c      bClass = INT(ARGS(1,1))
      
c      RET(1,1) = ComputeSource(bClass, CRESLT, pVar,CZ,DZ,IZ,LZ,RZ)
      
      

cC      CALL USER_PRINT_INTR('NLOC', NLOC)
cC      CALL USER_PRINT_INTR('NARG', NARG)
   
      
c      END
           
      REAL FUNCTION ComputeSource(bClass, CRESLT, pVar, 
     *CZ,DZ,IZ,LZ,RZ)
      IMPLICIT NONE
      
      
      
      REAL ComputeBubbleVol
      REAL Dbr
      
      __stack_point__ pVar
      INTEGER IZ(*)
      CHARACTER CZ(*)*(1)
      DOUBLE PRECISION DZ(*)
      LOGICAL LZ(*)
      REAL RZ(*)
      
      INTEGER bClass
      REAL density_g
      REAL volFrac
      REAL fi
      CHARACTER CRESLT*(*)
      
      CALL USER_GETVAR('air.Volume Fraction', CRESLT,
     * pVar,CZ,DZ,IZ,LZ,RZ)
     
      IF (CRESLT .NE. 'GOOD' ) THEN
        CALL ERRMSG('Error air.Volume Fraction')
        CALL ERRMSG(CRESLT)
        STOP
      ENDIF
      
      volFrac = RZ(pVar)
      
      IF(volFrac .GT. 1.01 .OR. volFrac .LT. -0.01) THEN
        CALL ERRMSG('Error air.Volume wrong num')
        CALL ERRMSG(volFrac)
        STOP
      ENDIF
      
      CALL USER_GETVAR('air.Density', CRESLT,
     * pVar,CZ,DZ,IZ,LZ,RZ)
     
      IF (CRESLT .NE. 'GOOD' ) THEN
        CALL ERRMSG('Error in air.Density')
        CALL ERRMSG(CRESLT)
        STOP
      ENDIF

      density_g = RZ(pVar)
      
      IF(density_g .GT. 1.5 .OR. density_g .LT. 1.0) THEN
        CALL ERRMSG('Error air.Density wrong num')
        CALL ERRMSG(density_g)
        STOP
      ENDIF
      
      CALL USER_GETVAR('air.f1', CRESLT,
     * pVar,CZ,DZ,IZ,LZ,RZ)
     
      IF (CRESLT .NE. 'GOOD' ) THEN
        CALL ERRMSG('Error in fi')
        CALL ERRMSG(CRESLT)
        STOP
      ENDIF
      
      fi = RZ(pVar)
      
      ComputeSource = density_g*ComputeBubbleVol(bClass)*(
     *-Dbr(bClass, volFrac, fi)
     *)
      
      END
      
      REAL FUNCTION Dbr(bClass, volFrac, fi)
      IMPLICIT NONE
      
      REAL ComputeBubbleVol
      
      
      INTEGER bClass
      REAL volFrac
      REAL ni
      REAL fi
      
          
      ni = volFrac*fi/ComputeBubbleVol(bClass) 
      
      Dbr=1.0E0*ni
      
      END
      
      
      REAL FUNCTION ComputeBubbleVol(bClass)
      IMPLICIT NONE
      
      REAL BubbleClasses(1:12)
      COMMON /C_BubbleClasses/ BubbleClasses
      
      INTEGER bClass
      
      ComputeBubbleVol = 3.14159265*BubbleClasses(bClass)**3.0E0/6.0E0
      
      END
      
      BLOCKDATA
      IMPLICIT NONE
      
      REAL BubbleClasses(1:12)
      
      DATA BubbleClasses /0.5E-3, 1.0E-3, 2.0E-3, 3.0E-3, 4.0E-3
     * , 5.0E-3, 6.0E-3, 7.0E-3, 8.0E-3, 10.0E-3, 12.0E-3, 16.0E-3/
     
      COMMON /C_BubbleClasses/ BubbleClasses
      
      END
