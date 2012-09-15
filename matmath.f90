MODULE matmath
 
 USE prmat
 IMPLICIT NONE
  
 CONTAINS

  SUBROUTINE matmultr(A,ZA,SA,B,ZB,SB,C)

   REAL*8 A(:,:)
   REAL*8 B(:,:)
   INTEGER ZA,SA,ZB,SB,I,J,K
   REAL*8, ALLOCATABLE :: C(:,:)

   ALLOCATE( C(ZA,SB) )
   I = 0
   J = 0
   K = 0

   DO I=1,ZA
    DO J = 1,SB
     DO K = 1,SA
      C(I,J) = C(I,J) + A(I,K)*B(K,J)
     END DO
    END DO
   END DO

 END SUBROUTINE matmultr

 SUBROUTINE matcomp(A,ZA,SA,B,ZB,SB,C,D)
 
  REAL*8 A(:,:), B(:,:)
  INTEGER ZA,SA,ZB,SB,D,I,J,K,G
  REAL*8, ALLOCATABLE :: C(:,:)

  D = SA + SB
  I = 0
  J = 0
  K = 0
  
  !!! Getting number of different columns !!!
  DO I=1,SB
   DO J=1,SA
    IF((B(1,I).EQ.A(1,J)).AND.(B(2,I).EQ.A(2,J)).AND.(B(3,I).EQ.A(3,J)).AND.(B(4,I).EQ.A(4,J))) THEN
     D = D - 1
    END IF
   END DO
  END DO
  I = 0
  J = 0
  ALLOCATE( C(4,D) ) 
  
  !!! Storing A in C !!!
  DO I=1,SA
   DO J = 1,4
    C(J,I) = A(J,I)
   END DO
  END DO
  I = 0
  J = 0
  IF(D.EQ.SA) THEN
   GO TO 100
  END IF

  !!! Getting other positions !!!
  K = SA
  DO I=1,SB
   G = 0
   DO J=1,SA
    IF((B(1,I).EQ.A(1,J)).AND.(B(2,I).EQ.A(2,J)).AND.(B(3,I).EQ.A(3,J)).AND.(B(4,I).EQ.A(4,J))) THEN
     G = G + 1
    ELSE
     EXIT
    END IF
   END DO
  IF(G.EQ.0) THEN
   K = K + 1
   C(1,K) = B(1,I)
   C(2,K) = B(2,I)
   C(3,K) = B(3,I)
   C(4,K) = B(4,I)
  END IF
  END DO
100 END SUBROUTINE matcomp

END MODULE
