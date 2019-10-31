!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This driver should demonstrate how these libraries might be joined using
!! the local matrix multiplication level of granularity.
!! You can run it with the following command:
!!    ./bin/localdriv Interface/input.mtx Interface/output.mtx \
!!                    Interface/blocks.inp
PROGRAM LocalDriver
  USE DMatrixModule, ONLY : Matrix_ldr
  USE SMatrixModule, ONLY : Matrix_lsr, ConstructMatrixFromFile, &
       & SplitMatrix, DestructMatrix
  IMPLICIT NONE
  CHARACTER(len=80) :: input_file, blocking_file
  TYPE(Matrix_lsr) :: mat, mat2
  INTEGER, DIMENSION(:), ALLOCATABLE :: block_list

  !! Read in the matrix from file.
  CALL get_command_argument(1, input_file)
  CALL ConstructMatrixFromFile(mat, input_file)

  !! Read the check matrix from file
  CALL get_command_argument(2, input_file)
  CALL ConstructMatrixFromFile(mat2, input_file)

  !! Read a second file which has the blocking parameters.
  CALL get_command_argument(3, blocking_file)
  CALL read_blocks(block_list, blocking_file)

  !! Split the matrix into blocks using those parameters.
  ! ALLOCATE(split_mat(num_blocks, num_blocks))
  ! CALL SplitMatrix(mat, num_blocks, num_blocks, split_mat, &
  !      & block_list, block_list)

  !! Convert to a list of dense matrices
  ! ALLOCATE()

  !! Convert to a list of blocks to multiply

  !! Perform the multiplication

  !! Merge back

  !! Check the result.

CONTAINS
  !> Read in the blocking information from a file.
  SUBROUTINE read_blocks(block_list, blocking_file)
    !> The list of block offsets to allocate and read in.
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: block_list
    CHARACTER(len=80), INTENT(IN) :: blocking_file
    INTEGER :: num_blocks
    INTEGER :: II

  END SUBROUTINE read_blocks
END PROGRAM LocalDriver
