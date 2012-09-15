MODULE readtheinput

 USE symmetry
 USE basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reads the inputfile in the following order:                          !
!   - title: first line of input file                                  !
!   - geometry: first block of input file; ended with END              !
!                                                                      ! 
! contains following modules:                                          !
!   - initialinput(): reads/writes title and calls readgeominput()     !
!   - readsgeominput(): reads dimensionality (SYSDIM)                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE initialinput()

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! variables: - TITLE: comment line at beginning of inputfile          !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  CHARACTER*8 TITLE
  INTEGER NUMINAT
  REAL*8,ALLOCATABLE :: POSMAT(:,:)
  REAL*8,ALLOCATABLE :: BASMAT(:,:,:)

  READ(*,*) TITLE
  WRITE(*,*) TITLE

  CALL readgeominput(NUMINAT,POSMAT)
  CALL readbasisinput(NUMINAT,BASMAT)

 END SUBROUTINE initialinput

 SUBROUTINE readgeominput(NUMINAT,POSMAT)
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! variables: -SYSDIM: dimensionality of system (0)                    !
 !                                                                     !
 ! calls different subroutines for specific dimensionality:            !
 !  - SYSDIM = 0: calls pointgroup() for molecules                     !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  INTEGER SYSDIM,NUMINAT
  REAL*8,ALLOCATABLE :: POSMAT(:,:)

  READ(*,*) SYSDIM
  IF(SYSDIM.EQ.0) THEN
   WRITE(*,*) "Dimensionality of the System: 0. Calculating a MOLECULE."
   CALL pointgroup(NUMINAT,POSMAT)
  END IF

 END SUBROUTINE readgeominput

END MODULE readtheinput
