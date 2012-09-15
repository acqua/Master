MODULE pointgroups
 
 USE matmath
 USE prmat
 
 IMPLICIT NONE

 CONTAINS
  SUBROUTINE ci(NUMINAT,POSMAT,NUMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calls: -matmultr                                                   !
  !        -matcomp                                                    !
  !                                                                    !
  ! variables: -NUMINAT (h)                                            !
  !            -ELEM(NUMINAT (h)                                       !
  !            -POSMAT(3:NUMINAT) (h)                                  !
  !            -NUMAT: number of all atoms                             !
  !            -CIPOS: atoms by CI matrix                              !
  !                                                                    !
  ! transformation matrices:                                           !
  !  -CI: inversion on center (0,0,0)                                  !
  !            ( -1,  0,  0)                                           !            
  !       CI = (  0, -1,  0)                                           !
  !            (  0,  0, -1)                                           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER NUMINAT, NUMAT,SMAT1
   REAL*8, ALLOCATABLE :: POSMAT(:,:)
   REAL*8 CIM(5,5)
   REAL*8, ALLOCATABLE :: CIPOS(:,:)
   REAL*8, ALLOCATABLE :: MAT1(:,:)

  !!! initialization !!!

   CIM(1:5,1:5) = 0

  !!! construction of transformation matrices !!!

   CIM(1,1) = 1
   CIM(2,2) =  1
   CIM(3,3) = -1
   CIM(4,4) = -1
   CIM(5,5) = -1

  !!! C2 rotation aroung z axis, mirror plane xz, mirror plane yz !!!

   CALL matmultr(CIM,5,5,POSMAT,5,NUMINAT,CIPOS)

  !!! comparing pos matrices !!!

   CALL matcomp(POSMAT,5,NUMINAT,CIPOS,5,NUMINAT,MAT1,SMAT1)
   DEALLOCATE( CIPOS,POSMAT )
   NUMAT = SMAT1
  !!! storing POSMAT !!!

   ALLOCATE( POSMAT(5,NUMAT) )
   POSMAT = MAT1
   DEALLOCATE( MAT1 )
   CALL printmatr(POSMAT,5,NUMAT)

  END SUBROUTINE ci
  
  SUBROUTINE c2x(NUMINAT,POSMAT,NUMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calls: -matmultr                                                   !
  !        -matcomp                                                    !
  !                                                                    !
  ! variables: -NUMINAT (h)                                            !
  !            -ELEM(NUMINAT (h)                                       !
  !            -POSMAT(3:NUMINAT) (h)                                  !
  !            -NUMAT: number of all atoms                             !
  !            -CIPOS: atoms by CI matrix                              !
  !                                                                    !
  ! transformation matrices:                                           !
  !  -CI: inversion on center (0,0,0)                                  !
  !             (  1,  0,  0)                                          !            
  !       C2X = (  0, -1,  0)                                          !
  !             (  0,  0, -1)                                          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER NUMINAT, NUMAT,SMAT1
   REAL*8, ALLOCATABLE :: POSMAT(:,:)
   REAL*8 C2(5,5)
   REAL*8, ALLOCATABLE :: C2XPOS(:,:)
   REAL*8, ALLOCATABLE :: MAT1(:,:)

  !!! initialization !!!

   C2(1:5,1:5) = 0
   
  !!! construction of transformation matrices !!!

   C2(1,1) =  1
   C2(2,2) =  1
   C2(3,3) =  1
   C2(4,4) = -1
   C2(5,5) = -1

  !!! C2 rotation aroung z axis, mirror plane xz, mirror plane yz !!!

   CALL matmultr(C2,5,5,POSMAT,5,NUMINAT,C2XPOS)

  !!! comparing pos matrices !!!

   CALL matcomp(POSMAT,5,NUMINAT,C2XPOS,5,NUMINAT,MAT1,SMAT1)
   DEALLOCATE( C2XPOS,POSMAT )
   NUMAT = SMAT1

  !!! storing POSMAT !!!

   ALLOCATE( POSMAT(5,NUMAT) )
   POSMAT = MAT1
   DEALLOCATE( MAT1 )
   CALL printmatr(POSMAT,5,NUMAT)

  END SUBROUTINE c2x

  SUBROUTINE c2vz(NUMINAT,POSMAT,NUMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calls: -matmultr                                                   !
  !        -matcomp                                                    !
  !                                                                    !
  ! variables: -NUMINAT (h)                                            !
  !            -ELEM(NUMINAT) (h)                                      !
  !            -POSMAT(3:NUMINAT) (h)                                  !
  !            -NUMAT: number of all atoms                             !
  !            -C2ZPOS: atoms by C2Z matrix                            !
  !                                                                    !
  ! transformation matrices:                                           !
  !  -C2Z: rotation about 180° around z axis                           !
  !              ( -1,  0,  0)                                         !
  !        C2Z = (  0, -1,  0)                                         !
  !              (  0,  0,  1)                                         !
  !  -SXZ: mirror plane along xz                                       !
  !              (  1,  0,  0)                                         !
  !        SXZ = (  0, -1,  0)                                         !
  !              (  0,  0,  1)                                         !
  !  -SYZ: mirror plane along yz                                       !
  !              ( -1,  0,  0)                                         !
  !        SYZ = (  0,  1,  0)                                         !
  !              (  0,  0,  1)                                         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   INTEGER NUMINAT, NUMAT,SMAT1,SMAT2,SMAT3
   REAL*8, ALLOCATABLE :: POSMAT(:,:)
   REAL*8 C2Z(5,5), SXZ(5,5), SYZ(5,5)
   REAL*8, ALLOCATABLE :: C2ZPOS(:,:)
   REAL*8, ALLOCATABLE :: SXZPOS(:,:)
   REAL*8, ALLOCATABLE :: SYZPOS(:,:)
   REAL*8, ALLOCATABLE :: MAT1(:,:)
   REAL*8, ALLOCATABLE :: MAT2(:,:)
   REAL*8, ALLOCATABLE :: MAT3(:,:)

  !!! initialization !!!
  
   C2Z(1:5,1:5) = 0
   SXZ(1:5,1:5) = 0
   SYZ(1:5,1:5) = 0

  !!! construction of transformation matrices !!!
   
   C2Z(1,1) =  1
   C2Z(2,2) =  1
   C2Z(3,3) = -1
   C2Z(4,4) = -1
   C2Z(5,5) =  1
   SXZ(1,1) =  1
   SXZ(2,2) =  1
   SXZ(3,3) =  1
   SXZ(4,4) = -1
   SXZ(5,5) =  1
   SYZ(1,1) =  1
   SYZ(2,2) =  1
   SYZ(3,3) = -1
   SYZ(4,4) =  1
   SYZ(5,5) =  1

  !!! C2 rotation aroung z axis, mirror plane xz, mirror plane yz !!!
  
   CALL matmultr(C2Z,5,5,POSMAT,5,NUMINAT,C2ZPOS)
   CALL matmultr(SXZ,5,5,POSMAT,5,NUMINAT,SXZPOS)
   CALL matmultr(SYZ,5,5,POSMAT,5,NUMINAT,SYZPOS)

  !!! comparing pos matrices !!!

   CALL matcomp(POSMAT,5,NUMINAT,C2ZPOS,5,NUMINAT,MAT1,SMAT1)
   DEALLOCATE( C2ZPOS )
   CALL matcomp(MAT1,5,SMAT1,SXZPOS,5,NUMINAT,MAT2,SMAT2)
   DEALLOCATE( MAT1,SXZPOS )
   CALL matcomp(MAT2,5,SMAT2,SYZPOS,5,NUMINAT,MAT3,NUMAT)
   DEALLOCATE( MAT2,SYZPOS,POSMAT )

  !!! storing POSMAT !!!

   ALLOCATE( POSMAT(5,NUMAT) )
   POSMAT = MAT3
   DEALLOCATE( MAT3 )
   CALL printmatr(POSMAT,5,NUMAT)

  END SUBROUTINE c2vz

END MODULE
