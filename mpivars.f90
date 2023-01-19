  !     .. header file of MPI variables
module mpivars

    use mpi

    implicit none

    !     .. variables

    integer :: rank              ! current node number
    integer :: size              ! total number of nodes 
    integer  :: numproc           ! size  total number of nodes
    integer :: ierr
    integer :: stat(MPI_STATUS_SIZE)
    integer :: source
    integer :: dest
    integer :: tag
    parameter(tag = 0)
    integer ::  err

integer  flag_solver


end module mpivars
