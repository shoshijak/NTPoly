!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This driver should demonstrate how these libraries might be joined using
!! the local matrix multiplication level of granularity.
PROGRAM LocalDriver
  USE SMatrixModule, ONLY : Matrix_lsr, ConstructMatrixFromFile, &
       & SplitMatrix, DestructMatrix
  IMPLICIT NONE
  CHARACTER(len=80) :: input_file, blocking_file
  TYPE(Matrix_lsr) :: mat
  INTEGER, DIMENSION(:), ALLOCATABLE :: block_list

  !! Read in the matrix from file.
  CALL get_command_argument(1, input_file)
  CALL ConstructMatrixFromFile(mat, input_file)

  !! Read a second file which has the blocking parameters.
  CALL get_command_argument(1, blocking_file)

  !! Split the matrix into blocks using those parameters.
  CALL

  !! Cleanup
  CALL DestructMatrix(mat)
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
