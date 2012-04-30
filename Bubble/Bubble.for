#include "cfx5ext.h"
#define DEBUG
dllexport(bubble_source)
      SUBROUTINE BUBBLE_SOURCE( 
     &  NLOC, NRET, NARG, RET, ARGS, CRESLT,CZ,DZ,IZ,LZ,RZ )
      IMPLICIT NONE
#include "stack_point.h"
C-----Volane funkce
      REAL ComputeSource
C-----Parametry podprogramu
      INTEGER NLOC,NARG,NRET
      CHARACTER CRESLT*(*)
      REAL ARGS(NLOC,NARG), RET(NLOC,NRET)
      INTEGER IZ(*)
      CHARACTER CZ(*)*(1)
      DOUBLE PRECISION DZ(*)
      LOGICAL LZ(*)
      REAL RZ(*)
C-----Lokalni promene
      __stack_point__ pVar
      INTEGER bClass
      INTEGER ILOC
      bClass = INT(ARGS(1,1))
C-----Kod
C vynuluje vse pro jistotu
      CALL SET_A_0( RET, NLOC*NRET )
      
C      RET(1,1) = ComputeSource(bClass, CRESLT, pVar,CZ,DZ,IZ,LZ,RZ)
      CRESLT = 'GOOD'

      END
C=======================================================================
      REAL FUNCTION ComputeSource(bClass, CRESLT, pVar, 
     *CZ,DZ,IZ,LZ,RZ)
      IMPLICIT NONE
C-----Symbolicke konstanty
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
C-----Volane funkce
      REAL ComputeBubbleVol
      REAL Dbr
C-----Common bloky
      CHARACTER NameOfClasses(1:NUMBER_OF_CLASSES)*10
      COMMON /C_NameOfClasses/ NameOfClasses
C-----Parametry podprogramu
      __stack_point__ pVar
      INTEGER IZ(*)
      CHARACTER CZ(*)*(1)
      DOUBLE PRECISION DZ(*)
      LOGICAL LZ(*)
      REAL RZ(*)
      CHARACTER CRESLT*(*)
C-----Lokalni promene
      INTEGER bClass
      REAL density_g
      REAL volFrac
      REAL fi
 
      CALL USER_GETVAR('air.Volume Fraction', CRESLT,
     * pVar,CZ,DZ,IZ,LZ,RZ)
     
      IF (CRESLT .NE. 'GOOD' ) THEN
        CALL ERRMSG('Error air.Volume Fraction CRESLT: '//CRESLT)
        CALL ERRMSG(CRESLT)
        STOP
      ENDIF
      
      volFrac = RZ(pVar)

#ifdef DEBUG
      IF(volFrac .GT. 1.01 .OR. volFrac .LT. -0.01) THEN
        CALL ERRMSG('Error air.Volume wrong num')
        CALL ERRMSG(volFrac)
        STOP
      ENDIF
#endif
      
      CALL USER_GETVAR('air.Density', CRESLT,
     * pVar,CZ,DZ,IZ,LZ,RZ)
     
      IF (CRESLT .NE. 'GOOD' ) THEN
        CALL ERRMSG('Error in air.Density CRESLT: '//CRESLT)
        STOP
      ENDIF

      density_g = RZ(pVar)
      
#ifdef DEBUG
      IF(density_g .GT. 1.5 .OR. density_g .LT. 1.0) THEN
        CALL ERRMSG('Error air.Density wrong num')
        CALL ERRMSG(density_g)
        STOP
      ENDIF
#endif
      
      CALL USER_GETVAR
     *(NameOfClasses(bClass)(1:LEN_TRIM(NameOfClasses(bClass)))
     *, CRESLT,
     * pVar,CZ,DZ,IZ,LZ,RZ)
     
      IF (CRESLT .NE. 'GOOD' ) THEN
        CALL ERRMSG('Error in fi CRESLT: '//CRESLT)
        CALL ERRMSG(CRESLT)
        STOP
      ENDIF
      
      fi = RZ(pVar)
      
      ComputeSource = density_g*ComputeBubbleVol(bClass)*(
     *-Dbr(bClass, volFrac, fi)
     *)
      
      END
C=======================================================================
      REAL FUNCTION Dbr(bClass, volFrac, fi)
      IMPLICIT NONE
C-----Volane funkce
      REAL ComputeBubbleVol
C-----Parametry podprogramu
      INTEGER bClass
      REAL volFrac
      REAL fi
C-----Lokalni promene
      REAL ni

      ni = volFrac*fi/ComputeBubbleVol(bClass) 
      
      Dbr = 1.0E0*ni
      
      END
C=======================================================================
      REAL FUNCTION ComputeBubbleVol(bClass)
      IMPLICIT NONE
C-----Symbolicke konstanty
      REAL PI
      PARAMETER (PI = 3.1415926535897931D0)
C-----Common bloky
      REAL BubbleClasses(1:12)
      COMMON /C_BubbleClasses/ BubbleClasses
C-----Parametry podprogramu
      INTEGER bClass
      
      ComputeBubbleVol = PI*BubbleClasses(bClass)**3.0E0/6.0E0
      
      END
C=======================================================================
      BLOCKDATA
      IMPLICIT NONE
C-----Symbolicke konstanty
      INTEGER NUMBER_OF_CLASSES
      PARAMETER (NUMBER_OF_CLASSES = 12)
      
      REAL BubbleClasses(1:NUMBER_OF_CLASSES)
      INTEGER NumberOfClasses
      CHARACTER NameOfClasses(1:NUMBER_OF_CLASSES)*10
      
      DATA NumberOfClasses /NUMBER_OF_CLASSES/
      DATA BubbleClasses /0.5E-3, 1.0E-3, 2.0E-3, 3.0E-3, 4.0E-3
     * , 5.0E-3, 6.0E-3, 7.0E-3, 8.0E-3, 10.0E-3, 12.0E-3, 16.0E-3/
     
      DATA NameOfClasses(1) /'air.f1'/
      DATA NameOfClasses(2) /'air.f2'/
      DATA NameOfClasses(3) /'air.f3'/
      DATA NameOfClasses(4) /'air.f4'/
      DATA NameOfClasses(5) /'air.f5'/
      DATA NameOfClasses(6) /'air.f6'/
      DATA NameOfClasses(7) /'air.f7'/
      DATA NameOfClasses(8) /'air.f8'/
      DATA NameOfClasses(9) /'air.f9'/
      DATA NameOfClasses(10) /'air.f10'/
      DATA NameOfClasses(11) /'air.f11'/
      DATA NameOfClasses(12) /'air.f12'/
    
      COMMON /C_BubbleClasses/ BubbleClasses
      COMMON /C_NumberOfClasses/ NumberOfClasses
      COMMON /C_NameOfClasses/ NameOfClasses
      
      END
