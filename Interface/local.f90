!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This driver should demonstrate how these libraries might be joined using
!! the local matrix multiplication level of granularity.
!! You can run it with the following command:
!!    ./bin/localdriv Interface/mata.mtx Interface/matb.mtx Interface/matc.mtx \
!!                    Interface/blocks.inp
!! The point of modification is the subroutine MultiplyPairList.
PROGRAM LocalDriver
  USE DataTypesModule, ONLY : NTREAL
  USE DMatrixModule, ONLY : Matrix_ldr, ConstructEmptyMatrix, DestructMatrix
  USE SMatrixModule, ONLY : Matrix_lsr, ConstructMatrixFromFile, DestructMatrix, PrintMatrix
  USE SMatrixAlgebraModule, ONLY : IncrementMatrix, MatrixNorm
  IMPLICIT NONE
  !! A data type which describes the pair of matrices to multiply
  TYPE MultPair_t
     TYPE(Matrix_ldr), POINTER :: amat
     TYPE(Matrix_ldr), POINTER :: bmat
     TYPE(Matrix_ldr), POINTER :: cmat
     INTEGER :: M, N, K
  END TYPE MultPair_t
  !! The file names.
  CHARACTER(len=80) :: mata_file, matb_file, matc_file, blocking_file
  !! The input and output matrices.
  TYPE(Matrix_lsr) :: mata, matb, matc, matc_computed
  !! The matrix once we split it by blocks.
  TYPE(Matrix_ldr), DIMENSION(:,:), ALLOCATABLE, TARGET :: splita
  TYPE(Matrix_ldr), DIMENSION(:,:), ALLOCATABLE, TARGET :: splitb
  TYPE(Matrix_ldr), DIMENSION(:,:), ALLOCATABLE, TARGET :: splitc
  !! A list of block sizes
  INTEGER, DIMENSION(:), ALLOCATABLE :: block_list
  !! The list of multiplication pairs to send to DBCSR
  TYPE(MultPair_t), DIMENSION(:), ALLOCATABLE :: pair_list
  !! Temporary variables
  REAL(NTREAL) :: normval
  INTEGER :: num_blocks
  INTEGER :: II, JJ, KK, CC

  !! Read in input matrices from file.
  CALL get_command_argument(1, mata_file)
  CALL ConstructMatrixFromFile(mata, mata_file)
  CALL get_command_argument(2, matb_file)
  CALL ConstructMatrixFromFile(matb, matb_file)

  !! Read the check matrix from file
  CALL get_command_argument(3, matc_file)
  CALL ConstructMatrixFromFile(matc, matc_file)

  !! Read a second file which has the blocking parameters.
  CALL get_command_argument(4, blocking_file)
  CALL read_blocks(block_list, blocking_file)
  num_blocks = SIZE(block_list)

  !! Split the matrix into blocks using those parameters.
  CALL CSR_To_Dense2d(mata, block_list, splita)
  CALL CSR_To_Dense2d(matb, block_list, splitb)
  !! Allocate C as an empty matrix
  ALLOCATE(splitc(num_blocks, num_blocks))
  DO II = 1, num_blocks
     DO JJ = 1, num_blocks
        CALL ConstructEmptyMatrix(splitc(II,JJ), splita(II,JJ)%rows, &
             & splita(II,JJ)%columns)
        splitc(II,JJ)%DATA = 0
     END DO
  END DO

  !! Convert to a list of blocks to multiply
  ALLOCATE(pair_list(num_blocks*num_blocks*num_blocks))
  CC = 1 !! because one based indexing is a pain.
  DO II = 1, num_blocks
     DO JJ = 1, num_blocks
        DO KK = 1, num_blocks
           pair_list(CC)%M = splita(II,KK)%rows
           pair_list(CC)%N = splitb(KK,JJ)%columns
           pair_list(CC)%K = splita(II,KK)%columns
           pair_list(CC)%amat => splita(II,KK)
           pair_list(CC)%bmat => splitb(KK,JJ)
           pair_list(CC)%cmat => splitc(II,JJ)
           CC = CC + 1
        END DO
     END DO
  END DO

  !! Perform the multiplication
  CALL MultiplyPairList(pair_list)

  !! Merge back
  CALL Dense2d_To_CSR(splitc, matc_computed)

  ! !! Check the result.
  CALL IncrementMatrix(matc_computed, matc, alpha_in=-1.0_NTREAL)
  normval = MatrixNorm(matc)
  WRITE(*,*) "Error:", normval
  IF (normval .GT. 1e-10) THEN
     STOP 1
  END IF

  !! Cleanup
  DO II = 1, num_blocks
     DO JJ = 1, num_blocks
        CALL DestructMatrix(splita(II,JJ))
        CALL DestructMatrix(splitb(II,JJ))
        CALL DestructMatrix(splitc(II,JJ))
     END DO
  END DO
  CALL DestructMatrix(mata)
  CALL DestructMatrix(matb)
  CALL DestructMatrix(matc)
  CALL DestructMatrix(matc_computed)
  DEALLOCATE(splita)
  DEALLOCATE(splitb)
  DEALLOCATE(splitc)
  DEALLOCATE(pair_list)
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read in the blocking information from a file.
  !! The file is as follows. The first line is the number of blocks. Then
  !! a blank space. After that, all the block sizes are separated by line.
  SUBROUTINE read_blocks(block_list, blocking_file)
    !> The list of block offsets to allocate and read in.
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: block_list
    !> The name of the file with the blocking information.
    CHARACTER(len=80), INTENT(IN) :: blocking_file
    INTEGER :: num_blocks
    INTEGER :: II
    INTEGER, PARAMETER :: file_handler = 16

    OPEN(file_handler, file=blocking_file, status='old')
    READ(file_handler,*) num_blocks
    ALLOCATE(block_list(num_blocks))

    DO II = 1, num_blocks
       READ(file_handler,*) block_list(II)
    END DO
    CLOSE(file_handler)
  END SUBROUTINE read_blocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a matrix from NTPoly CSR to a 2D array of dense matrices.
  SUBROUTINE CSR_To_Dense2d(mat, block_list, mat_array)
    USE DMatrixModule, ONLY : ConstructMatrixDFromS
    USE SMatrixModule, ONLY : SplitMatrix, DestructMatrix
    !> The matrix to convert.
    TYPE(Matrix_lsr), INTENT(IN) :: mat
    !> A list of block sizes
    INTEGER, DIMENSION(:), INTENT(IN) :: block_list
    !> The 2d matrix split by the block list
    TYPE(Matrix_ldr), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: mat_array
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: split_mat_sparse

    !! Split to a sparse matrix
    num_blocks = SIZE(block_list)
    ALLOCATE(split_mat_sparse(num_blocks, num_blocks))
    CALL SplitMatrix(mat, num_blocks, num_blocks, split_mat_sparse, &
         & block_list, block_list)

    !! Convert to an array of dense matrices
    ALLOCATE(mat_array(num_blocks, num_blocks))
    DO II = 1, num_blocks
       DO JJ = 1, num_blocks
          CALL ConstructMatrixDFromS(split_mat_sparse(II,JJ), &
               & mat_array(II,JJ))
       END DO
    END DO

    !! Cleanup
    DO II = 1, num_blocks
       DO JJ = 1, num_blocks
          CALL DestructMatrix(split_mat_sparse(II,JJ))
       END DO
    END DO
    DEALLOCATE(split_mat_sparse)
  END SUBROUTINE CSR_To_Dense2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply all of the matrix pairs in a list.
  !! This is the part you work on ;)
  !! I've provided a reference implementation using a call to BLAS to check.
  SUBROUTINE MultiplyPairList(pair_list)
    !> The pair list to work on.
    TYPE(MultPair_t), DIMENSION(:), INTENT(INOUT) :: pair_list
    INTEGER :: M, N, K
    INTEGER :: LDA, LDB, LDC
    INTEGER :: II
    DOUBLE PRECISION, PARAMETER :: ALPHA = 1.0
    DOUBLE PRECISION, PARAMETER :: BETA = 1.0
    CHARACTER, PARAMETER :: TRANSA = 'N'
    CHARACTER, PARAMETER :: TRANSB = 'N'

    !! As a demonstration, here is how it would be done with a loop over calls
    !! to blas dgemm.
    DO II = 1, SIZE(pair_list)
       M = pair_list(II)%M
       N = pair_list(II)%N
       K = pair_list(II)%K
       LDA = M
       LDB = K
       LDC = M

       CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, pair_list(II)%amat%data, &
            & LDA, pair_list(II)%bmat%data, LDB, BETA, &
            & pair_list(II)%cmat%data, LDC)
    END DO
  END SUBROUTINE MultiplyPairList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a matrix from a 2d array of dense matirces to NTPoly CSR.
  SUBROUTINE Dense2d_To_CSR(mat_array, mat)
    USE DMatrixModule, ONLY : ConstructMatrixSFromD
    USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr
    USE SMatrixModule, ONLY : ComposeMatrix, DestructMatrix
    !> The 2d matrix split by the block list
    TYPE(Matrix_ldr), DIMENSION(:,:), INTENT(IN) :: mat_array
    !> The matrix to convert.
    TYPE(Matrix_lsr), INTENT(OUT) :: mat
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: split_mat_sparse
    TYPE(MatrixMemoryPool_lr), DIMENSION(:,:), ALLOCATABLE :: pool

    !! Convert to a 2D array of sparse matrices.
    num_blocks = SIZE(mat_array, DIM=1)
    ALLOCATE(split_mat_sparse(num_blocks, num_blocks))
    DO II = 1, num_blocks
       DO JJ = 1, num_blocks
          CALL ConstructMatrixSFromD(mat_array(II,JJ), split_mat_sparse(II,JJ))
       END DO
    END DO

    !! Merge
    CALL ComposeMatrix(split_mat_sparse, num_blocks, num_blocks, mat)

    !! Cleanup
    DO II = 1, num_blocks
       DO JJ = 1, num_blocks
          CALL DestructMatrix(split_mat_sparse(II,JJ))
       END DO
    END DO
    DEALLOCATE(split_mat_sparse)
  END SUBROUTINE Dense2d_To_CSR
END PROGRAM LocalDriver
