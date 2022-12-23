program game_of_life
    use mpi_f08

    implicit none
    
    ! Variables
    integer :: height, width
    integer :: max_gen, gen
    integer :: n_ranks, my_rank, root, north_rank, south_rank, east_rank, west_rank
    integer :: n_rows, n_cols, row, col, ib, ie, jb, je
    logical, dimension(:, :), pointer :: old_world, new_world, tmp_world
    type(MPI_Comm) :: comm
    type(MPI_Datatype) :: a_row, a_col
    type(MPI_Status) :: status

    ! Initializing MPI, getting number of processes and IDs
    call MPI_Init()
    call MPI_Comm_rank( comm, my_rank)
    call MPI_Comm_size( comm, n_ranks)

    ! Setting root rank to 0 and checking ranks to processes
    root = 0
    if (my_rank == root) read *, height, width, max_gen, n_rows, n_cols

    ! Broadcasting n_rows, n_cols, height, width & max_gen
    call MPI_Bcast(n_rows, 1, MPI_INTEGER, root, comm)
    call MPI_Bcast(n_cols, 1, MPI_INTEGER, root, comm)
    call MPI_Bcast(height, 1, MPI_INTEGER, root, comm)
    call MPI_Bcast(width, 1, MPI_INTEGER, root, comm)
    call MPI_Bcast(max_gen, 1, MPI_INTEGER, root, comm)

    !Error message when number of processes != number of subgrids
    if ((n_rows * n_cols) /= n_ranks) then
        print "(a)", "Incorrect number of processes"
        call MPI_Abort( comm, MPI_ERR_TOPOLOGY )
    end if

    ! Getting coordinates for current rank
    call get_coords(my_rank, n_rows, n_cols, row, col)

    ! Setting rank neighbours
    north_rank = get_rank( row - 1, col,     n_rows, n_cols )
    south_rank = get_rank( row + 1, col,     n_rows, n_cols )
    west_rank  = get_rank( row,     col - 1, n_rows, n_cols )
    east_rank  = get_rank( row,     col + 1, n_rows, n_cols )

    ! Partitioning grid in subgrids
    call partition( row, n_rows, height, ib, ie )
    call partition( col, n_cols,  width, jb, je )

    ! Allocating world map variables
    allocate(old_world(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(new_world(ib - 1:ie + 1, jb - 1:je + 1))

    ! Creating MPI types
    call MPI_Type_contiguous(ie - ib + 3, MPI_LOGICAL, a_col)
    call MPI_Type_commit(a_col)
    block
        integer(kind=MPI_ADDRESS_KIND) :: lb, real_extent
        type(MPI_Datatype) :: a_tmp_row

        call MPI_Type_vector(je - jb + 1, 1, ie - ib + 3, MPI_LOGICAL, a_tmp_row)
        call MPI_Type_get_extent( MPI_LOGICAL, lb, real_extent )
        call MPI_Type_create_resized( a_tmp_row, lb, real_extent, a_row )
        call MPI_Type_commit( a_row )
    end block

    ! Reading input map and distributing it
    call read_map(old_world, height, width)

    ! Updating ghost layer borders
    call update_borders(old_world)

    ! Actual Game of Life
    do gen = 1, max_gen
        !printing Generation number and map before each iteration (when in root)
        if (my_rank == root) print "(a, i0)", "Generation ", gen
        call print_map(old_world, height, width)
        !calculating next map
        call next_gen(old_world, new_world)
        !updating borders
        call update_borders(new_world)
        !if (my_rank == root) call wait_cls(100)
        !aborting if world is still
        call MPI_Bcast(world_is_still( old_world, new_world ), 1, MPI_LOGICAL, root, comm)
        if (world_is_still( old_world, new_world )) exit
        !swapping new map to old map
        tmp_world => old_world;  old_world => new_world;  new_world => tmp_world 
    end do

    ! Deallocating pointers
    if (associated(old_world)) deallocate(old_world)
    if (associated(new_world)) deallocate(new_world)

    !Finallizing MPI
    call MPI_Type_free(a_row)   !two data types
    call MPI_Type_free(a_col)
    call MPI_Finalize()         !actual environment

contains

    function world_is_still( old_map, new_map ) result(still)
        logical, dimension(:, :), pointer, intent(in) :: old_map, new_map
        logical :: still
        logical, dimension(n_ranks)                   :: all_is_still

        still = all( old_map .eqv. new_map )
        call MPI_Gather(still, 1, MPI_LOGICAL, all_is_still,   1, MPI_LOGICAL, root, comm)
        if (my_rank == root) then
            still = all ( all_is_still .eqv. .True.)
        else
            still = .false.
        endif
    end function world_is_still

    subroutine get_coords( rank, n_rows, n_cols, row, col )
        integer, intent(in)    :: rank, n_rows, n_cols
        integer, intent(inout) :: row, col

        row = modulo(rank, n_rows)
        col = (rank - row) / n_rows
        if (0 <= col .and. col < n_cols) then
            return
        else
            print "(a, 2(i0, a))", "get_coords: rank ", rank, &
                " is outside the column range [0, ", n_cols, ")."
            call MPI_Abort( comm, MPI_ERR_TOPOLOGY )
        end if
    end subroutine get_coords

    
    integer function get_rank(row, col, n_rows, n_cols)
        integer, intent(in) :: row, col, n_rows, n_cols
        integer :: row_, col_

        row_ = row
        col_ = col

        if (row_ < 0) then
            row_ = n_rows - abs(row_)
        else if (row_ >= n_rows) then
            row_ = row_ - n_rows
        end if

        if (col_ < 0) then
            col_ = n_cols - abs(col_)
        else if (col_ >= n_cols) then
            col_ = col_ - n_cols
        end if

        get_rank = row_ + col_ * n_rows
        return
    end function get_rank

    subroutine update_borders( map )
        logical, dimension(:, :), pointer, intent(inout) :: map
    
        ! synchronize rows with north and south
        call MPI_Sendrecv( map(ib,jb), 1, a_row, north_rank, 1, &
                           map(ie+1,jb), 1, a_row, south_rank, 1, comm, status )
        call MPI_Sendrecv( map(ie,jb), 1, a_row, south_rank, 2, &
                           map(ib-1,jb), 1, a_row, north_rank, 2, comm, status )
    
        ! synchronize columns with east and west
        call MPI_Sendrecv( map(ib-1,jb), 1, a_col, west_rank, 3, &
                           map(ib-1,je+1), 1, a_col, east_rank, 3, comm, status )
        call MPI_Sendrecv( map(ib-1,je), 1, a_col, east_rank, 4, &
                           map(ib-1,jb-1), 1, a_col, west_rank, 4, comm, status )
    end subroutine update_borders
    
    subroutine read_map( map, h, w )
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer, intent(in) :: h, w
        character(len=:), allocatable :: line       ! variable to hold a line
        logical,          allocatable :: temp(:)    ! temporary array to hold a portion of the line
        integer :: i, j, rb, re, cb, ce             ! indices for looping
        integer :: current_row, current_col, dst    ! variables for holding partition info
        
        ! root reads map
        if (my_rank == root) then
            allocate(character(len=w) :: line)
            do current_row = 0, n_rows - 1
                call partition(current_row, n_rows, h, rb, re)
                ! loop over the rows
                do i = rb, re
                    ! read in a line
                    read *, line
                    ! loop over the columns
                    do current_col = 0, n_cols - 1
                        call partition(current_col, n_cols, w, cb, ce)
                        dst = get_rank(current_row, current_col, n_rows, n_cols)
                        allocate(temp(ce - cb + 1))
                        ! loop over the map partition
                        do j = cb, ce 
                            ! convert map characters to logical values
                            select case (line(j:j))
                            case ('X')
                                temp(j - cb + 1) = .True.
                            case ('.')
                                temp(j - cb + 1) = .False.
                            case default
                                stop "read_map: wrong input character `" // line(j:j) // "`"
                            end select
                        end do
                        ! if partition is from rank, copy temp to map
                        if (dst == root) then
                            map(i, cb : ce)  = temp
                        ! if partition is not from rank, send temp to rank
                        else
                            call MPI_Send( temp, ce - cb + 1, MPI_LOGICAL, dst, 0,  comm )
                        end if
                        ! deallocate temp array
                        if (allocated( temp )) deallocate(temp)
                    end do
                end do
            end do
            ! deallocate line
            if (allocated( line )) deallocate(line)
        else
            ! loop over rows
            do i = ib, ie
                ! receive row from root
                call MPI_Recv(map(i,jb), 1, a_row, root, 0, comm, status)
            end do
        end if
    end subroutine read_map

    subroutine print_map(map, h, w)
        logical, dimension(:, :), pointer, intent(in) :: map
        integer, intent(in)                           :: h, w
        character(len=:), allocatable :: line
        logical,          allocatable :: temp(:)
        integer :: i, j, rb, re, cb, ce

        if (my_rank == root) then
            block
                integer :: current_row
                integer :: current_col
                integer :: dst
                allocate(character(len=w) :: line)

                do current_row = 0, n_rows - 1
                    call partition(current_row, n_rows, h, rb, re)
                    do i = rb, re
                        do current_col = 0, n_cols - 1
                            call partition(current_col, n_cols, w, cb, ce)
                            allocate(temp(ce - cb + 1))
                            dst = get_rank(current_row, current_col, n_rows, n_cols)
                            if (dst == root) then
                                do j = cb, ce 
                                    line(j:j) = merge ( 'X', '.', map(i,j))
                                end do
                            else 
                                call MPI_Recv(temp, ce-cb+1, MPI_LOGICAL, dst, 0, comm, status)
                                do j = cb, ce 
                                    line(j:j) = merge ( 'X', '.', temp(j-cb+1))
                                end do
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do
                        print "(a)", line
                    end do
                end do
                print *
                if (allocated( line )) deallocate(line)
            end block
       else
           do i = ib, ie
               call MPI_Send( map(i, jb), 1, a_row, root, 0, comm)
           end do
       end if
        

    end subroutine print_map
    
    subroutine next_gen(old_map, new_map)
        logical, dimension(:, :), pointer, intent(inout) :: old_map, new_map
        integer :: i, j
        integer :: c ! the number of live neighbors

        do j = jb,je
            do i = ib,ie
                c = count(old_map(i - 1:i + 1, j - 1:j + 1))
                if (old_map(i, j)) then ! cell is live
                    new_map(i, j) = merge(.true., .false., 3 <= c .and. c <= 4)
                else ! cell is dead
                    new_map(i, j) = merge(.true., .false., c == 3)
                end if
            end do
        end do
        return
    end subroutine next_gen

    ! Divides 2D grid in the closest to equal subgrids
    subroutine partition( id, n_ids, size, b, e )
        integer, intent(in)    :: id, n_ids, size
        integer, intent(inout) :: b, e
        integer :: remainder, quotient

        remainder = modulo( size, n_ids )
        quotient  = (size - remainder) / n_ids
        b = 1 + quotient * (id    ) + min( remainder, id     )
        e =     quotient * (id + 1) + min( remainder, id + 1 )
    end subroutine partition

end program game_of_life
