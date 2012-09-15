MODULE prmat

 IMPLICIT NONE

 CONTAINS

  SUBROUTINE printmatr(A,Z,S)
  
   REAL*8 A(:,:)
   INTEGER Z,S,I,P

   P = 1

  DO WHILE (P.NE.S)
   IF ((S-P).GE.6) THEN
    DO I = 1,Z
     WRITE(*,'(7E17.5E2)') A(I,P:P+6)
    END DO
    P = P + 6
   ELSE IF ((S-P).EQ.5) THEN
    DO I = 1,Z
     WRITE(*,'(6E17.5E2)') A(I,P:P+5)
    END DO
    P = P + 5 
   ELSE IF ((S-P).EQ.4) THEN
    DO I = 1,Z
     WRITE(*,'(5E217.5E2)') A(I,P:P+4)
    END DO
    P = P + 4
   ELSE IF ((S-P).EQ.3) THEN
    DO I = 1,Z
     WRITE(*,'(4E17.5E2)') A(I,P:P+3)
    END DO
    P = P + 3
   ELSE IF ((S-P).EQ.2) THEN
    DO I = 1,Z
     WRITE(*,'(3E17.5E2)') A(I,P:P+2)
    END DO
    P = P + 2
   ELSE IF ((S-P).EQ.1) THEN
    DO I = 1,Z
     WRITE(*,'(2E17.5E2)') A(I,P:P+1)
    END DO
    P = P + 1
   END IF
  END DO

  END SUBROUTINE printmatr

  SUBROUTINE printveci(V,D)

   INTEGER V(:)
   INTEGER D,I

   DO I = 1,D
    WRITE(*,*) V(I)
   END DO
  
  END SUBROUTINE printveci

END MODULE
   
