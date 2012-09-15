MODULE symmetry

 USE prmat
 USE pointgroups
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! contains following subroutines:                                      !
!  -pointgroup(): for SYSDIM=0; reads pointgroup and atom positions;   !
!  -readatomposca(NUMINAT,ELEM,POSMAT): reads number of inequivalent   ! 
!   atoms [NUMINAT], elements [ELEM()] and position of atoms           !
!   [POSMAT(3,NUMINAT)] as cartesian coordinates in Angstrom           !
!  -readposmatc(NUMINAT,ELEM,POSMAT): reads elements [ELEM(NUMINAT)]   !
!   and positions [POSMAT(3,NUMINAT)] of inequivalent atoms in         !
!   Cartesian inputs.                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IMPLICIT NONE
 
 CONTAINS

  SUBROUTINE pointgroup(NUMINAT,POSMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calls readatomposca()                                              !
  !       [pointgroup](NUMINAT,ELEM,POSMAT,NUMAT)                      !
  !                                                                    !
  ! variables: -PG: point group number                                 ! 
  !            -NUMINAT: number of inequivalent atoms                  !
  !            -NUMAT: number of all atoms                             !
  !            -CORDTYPE: information about type of inp. coord.        !
  !            -ELEM(:): vector containing the element numbers         !
  !            -POSMAT(3,:): matrix containing the atom positions      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER PG, NUMINAT, NUMAT
   CHARACTER*8 CORDTYPE
   REAL*8, ALLOCATABLE :: POSMAT(:,:)

   READ(*,*) PG
   READ(*,*) NUMINAT
   READ(*,*) CORDTYPE
   WRITE(*,*) "Point group:",PG
   WRITE(*,*) "Number of inequivalent atoms:",NUMINAT

   IF(CORDTYPE.EQ."CARTA") THEN
    WRITE(*,*) "Atom positions given in Cartesian coordinates [Ang]."
    CALL readatomposca(NUMINAT,POSMAT)
   END IF

   IF(PG.EQ.1) THEN
    WRITE(*,*) "Point group C1."
    GO TO 100
   ELSE IF(PG.EQ.2) THEN
    WRITE(*,*) "Point group Ci."
    CALL ci(NUMINAT,POSMAT,NUMAT)
   ELSE IF(PG.EQ.3) THEN
    WRITE(*,*) "Point group C2(x)."
    CALL c2x(NUMINAT,POSMAT,NUMAT)
   ELSE IF(PG.EQ.15) THEN
    WRITE(*,*) "Point group C2v(z)."
    CALL c2vz(NUMINAT,POSMAT,NUMAT)
   END IF
   
100  END SUBROUTINE pointgroup

  SUBROUTINE readatomposca(NUMINAT,POSMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calls readposmatc()                                                !
  !                                                                    !
  ! variables: -NUMINAT: number of inequivalent atoms                  !
  !            -ELEM(:): vector containing the element numbers         !
  !            -POSMAT(3,:): matrix containing the atom positions      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER NUMINAT
   REAL*8, ALLOCATABLE :: POSMAT(:,:)

   ALLOCATE( POSMAT(5,NUMINAT) )
   
   POSMAT(1:5,1:NUMINAT) = 0
   
   CALL readposmatc(NUMINAT,POSMAT)
   CALL printmatr(POSMAT,5,NUMINAT)
   
  END SUBROUTINE readatomposca

  SUBROUTINE readposmatc(NUMINAT,POSMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! reads cartesian geometry inputs                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER NUMINAT, I, E
   REAL*8 POSMAT(:,:), X, Y, Z

   DO I=1,NUMINAT
    READ(*,*) E, X, Y, Z
    POSMAT(1,I) = I
    POSMAT(2,I) = E
    POSMAT(3,I) = X
    POSMAT(4,I) = Y
    POSMAT(5,I) = Z
   END DO
  
  END SUBROUTINE readposmatc
  
END MODULE
