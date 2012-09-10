      PROGRAM TESTER
      IMPLICIT NONE
      
c      character NameOfClasses(1:12)*50
c      integer i
      
c      DATA NameOfClasses(1) /'air.f1'/
c      DATA NameOfClasses(2) /'air.f2'/
c      DATA NameOfClasses(3) /'air.f3'/
c      DATA NameOfClasses(4) /'air.f4'/
c      DATA NameOfClasses(5) /'air.f5'/
c      DATA NameOfClasses(6) /'air.f6'/
c      DATA NameOfClasses(7) /'air.f7'/
c      DATA NameOfClasses(8) /'air.f8'/
c      DATA NameOfClasses(9) /'air.f9'/
c      DATA NameOfClasses(10) /'air.f10'/
c      DATA NameOfClasses(11) /'air.f11'/
c      DATA NameOfClasses(12) /'air.f12'/
      
c      DO i=1,12
c        WRITE(*,*) NameOfClasses(i)(1:LEN_TRIM(NameOfClasses(i)))
c      ENDDO
c      WRITE(*,*) matname(1)(1:LEN_TRIM(matname(1))), matname(2) 


      INTEGER test(1:5)
      INTEGER test2(1:5)
      
      DATA test /1,2,3,4,5/
      DATA test2 /10,20,30,40,50/
      
      
      CALL LULZ(test,test2)

      
      END
      
      SUBROUTINE LULZ(A,B)
      IMPLICIT NONE
      
      INTEGER A(5)
      INTEGER B(5)
      
      INTEGER C(2,5)
      C(1,:) = A
      C(2,:) = B
       
       
      WRITE(*,*) C(2,1)
      
      END
