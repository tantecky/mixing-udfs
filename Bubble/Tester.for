      PROGRAM TESTER
      IMPLICIT NONE
      
      character NameOfClasses(1:12)*50
      integer i
      
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
      
      DO i=1,12
        WRITE(*,*) NameOfClasses(i)(1:LEN_TRIM(NameOfClasses(i)))
      ENDDO
c      WRITE(*,*) matname(1)(1:LEN_TRIM(matname(1))), matname(2) 
      
      END
