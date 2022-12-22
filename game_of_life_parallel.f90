program game_of_life
    use mpi_f08
    implicit none
    
    ! Variables
    integer :: height, width
    integer :: max_gen, gen
    integer :: n_ranks, my_rank, root, north_rank, south_rank, east_rank, west_rank
    integer :: n_rows, n_cols, row, col, ib, ie, jb, je
    logical, dimension(:, :), pointer :: old_world, new_world, tmp_world
    logical :: ghost_flag, root_world_still_flag
    integer(kind=MPI_ADDRESS_KIND) :: lb, real_extent
    type(MPI_Datatype) :: a_row, a_col, a_tmp_row
    type(MPI_Status) :: status

    ! Initializing MPI, getting number of processes and IDs
    call MPI_Init()
    call MPI_Comm_rank( MPI_COMM_WORLD, my_rank)
    call MPI_Comm_size( MPI_COMM_WORLD, n_ranks)

    ! Setting root rank to 0 and checking ranks to processes
    root = 0
    if (my_rank == root) then
        read *, height, width, max_gen, n_rows, n_cols

        !error message when number of processes != number of subgrids
        if ((n_rows * n_cols) /= n_ranks) then
            print "(a)", "Incorrect number of processes"
            call MPI_Abort( MPI_COMM_WORLD, MPI_ERR_TOPOLOGY )
        end if
    end if

    ! Broadcasting n_rows, n_cols, height, width & max_gen
    call MPI_Bcast(n_rows, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(n_cols, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(height, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(width, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(max_gen, 1, MPI_INTEGER, root, MPI_COMM_WORLD)

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

    ! Defining MPI types
    call MPI_Type_contiguous(ie - ib + 3, MPI_LOGICAL, a_col)
    call MPI_Type_commit(a_col)
    call MPI_Type_vector(je - jb + 1, 1, ie - ib + 3, MPI_LOGICAL, a_tmp_row)
    call MPI_Type_get_extent( MPI_LOGICAL, lb, real_extent )
    call MPI_Type_create_resized( a_tmp_row, lb, real_extent, a_row )
    call MPI_Type_commit( a_row )

    ! Reading input map and broadcasting it
    call read_map(old_world, height, width)

    ! Updating ghost layer borders
    call update_borders(old_world)

    ! Actual Game of Life
    do gen = 1, max_gen
        !printing Generation number and map before each iteration
        if (my_rank == root) print "(a, i0)", "Generation ", gen
        call print_map(old_world, height, width)
        !calculating next map
        call next_gen(old_world, new_world)
        !updating borders
        call update_borders(new_world)
        !if (my_rank == root) call wait_cls(100)
        !aborting if world is still
        if (world_is_still( old_world, new_world )) exit
        !swapping new map to old map
        tmp_world => old_world;  old_world => new_world;  new_world => tmp_world 
    end do

    ! Cleaning memory from worlds
    if (associated( old_world )) deallocate(old_world)
    if (associated( new_world )) deallocate(new_world)

    !Finallizing MPI
    call MPI_Type_free(a_row)   !two data types
    call MPI_Type_free(a_col)
    call MPI_Finalize()         !actual environment

contains

    logical function world_is_still(old_map, new_map)
        logical, dimension(:, :), pointer, intent(in) :: old_map, new_map
        logical :: still

        still = all(old_map .eqv. new_map)
        call MPI_Allreduce(still, world_is_still, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
        return
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
            call MPI_Abort( MPI_COMM_WORLD, MPI_ERR_TOPOLOGY )
        end if
    end subroutine get_coords

    
    integer function get_rank(row, col, n_rows, n_cols)
        integer, intent(in) :: row, col, n_rows, n_cols

        if (row < 0) then
            row = n_rows - abs(row)
        else if (row >= n_rows) then
            row = row - n_rows
        end if

        if (col < 0) then
            col = n_cols - abs(col)
        else if (col >= n_cols) then
            col = col - n_cols
        end if

        get_rank = row + col * n_rows
        return
    end function get_rank

    subroutine update_borders(map)
        logical, dimension(:, :), pointer, intent(inout) :: map
    
        ! rows
        ! send top border to south rank and receive bottom border from north rank
        call MPI_Sendrecv(map(ie + 1, jb), 1, a_row, south_rank, 1, &
                         map(ib, jb), 1, a_row, north_rank, 1, &
                         MPI_COMM_WORLD, status)
        ! send bottom border to north rank and receive top border from south rank
        call MPI_Sendrecv(map(ib - 1, jb), 1, a_row, north_rank, 2, &
                         map(ie, jb), 1, a_row, south_rank, 2, &
                         MPI_COMM_WORLD, status)
    
        ! columns
        ! send left border to east rank and receive right border from west rank
        call MPI_Sendrecv(map(ib - 1, je + 1), 1, a_col, east_rank, 3, &
                         map(ib - 1, jb), 1, a_col, west_rank, 3, &
                         MPI_COMM_WORLD, status)
        ! send right border to west rank and receive left border from east rank
        call MPI_Sendrecv(map(ib - 1, jb - 1), 1, a_col, west_rank, 3, &
                         map(ib - 1, je), 1, a_col, east_rank, 3, &
                         MPI_COMM_WORLD, status)
    end subroutine update_borders

subroutine read_map(map, h, w)
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer, intent(in) :: h, w
        integer :: i, j 

        if (my_rank == root) then 

            block
                logical, dimension(:), allocatable :: t_row
                character(len=:), allocatable :: line
                integer :: dst
                allocate(t_row(0:w+1))
                allocate(character(len=w) :: line)

                ! We use a row and a column loop to study all processes
                do row = 0, n_rows - 1 
                    ! The rows indexes of each rank in the array are obtained
                    call partition(row, n_rows, h, rb, re)
                    do i = rb, re
                        ! The grid of the world is read line by line from the 
                        ! input file
                        read *, line(:)

                        do j = 1, w
                            ! The lines are transformed into logical values 
                            select case (line(j:j)) 
                            case ('X')
                                t_row(j:j) = .true.
                            case ('.')
                                t_row(j:j) = .false.
                            case default
                                stop "read_map: wrong input character `" // &
                                &line(j:j) // "`"
                            end select
                        end do !j

                        do col = 0, n_cols - 1
                            ! The column indexes of each rank in the array are 
                            ! obtained
                            call partition(col, n_cols, w, cb, ce)
                            ! The rank associated to each row and column
                            ! combination are determined
                            dst = get_rank(row, col, n_rows, n_cols)
                            ! If the lines belong to root, they are added 
                            ! directly in the map
                            if (dst == root) then
                                map(i, cb:ce) = t_row(cb:ce)
                            else !dst /= root
                                ! Part of the line, which consists on logical 
                                ! elements, is distributed to the appropriate 
                                ! process
                                call MPI_Send(t_row(cb-1), ce - cb + 3, &
                                &MPI_LOGICAL, dst, 0, MPI_COMM_WORLD)
                            end if
                        end do ! col
                    end do !i
                end do !row
                if (allocated(t_row)) deallocate(t_row)
                if (allocated(line)) deallocate(line)

            end block
            
        else  ! my_rank /= root

            ! The proper lines that are associated to each process are received
            rb = ib 
            re = ie
            
            do i = rb, re
                call MPI_Recv(map(i, jb-1), 1, a_row, root, 0, MPI_COMM_WORLD, status)
            end do
        end if 

        return 
        
    end subroutine read_map

    subroutine print_map(map,  h, w)
        integer, intent(in) :: h, w
        integer :: i, j
        logical, dimension(:, :), pointer, intent(inout) :: map

        if (my_rank == root) then

            print "(a, i0)", "Generation ", gen

            block
                logical, dimension(:), allocatable :: t_row
                character(len=:), allocatable :: line_character
                logical, dimension(:), allocatable :: line_logical

                integer :: src
                allocate(character(len=w) :: line_character)
                allocate(line_logical(1:w))
                allocate(t_row(0:w+1))
                
                ! This section is very similar to the read_map one, though in 
                ! this case, we have to determine the rank from which each line
                ! is being sent
                do row = 0, n_rows - 1 
                    call partition(row, n_rows, h, rb, re)

                    do i = rb, re 
                        do col = 0, n_cols - 1
                            call partition(col, n_cols, w, cb, ce)
                            ! The source is determined
                            src = get_rank(row, col, n_rows, n_cols) 
                            if (src == root) then
                                 line_logical(cb:ce) = map(i,cb:ce)
                            else
                                ! The lines are received from each process
                                call MPI_Recv(t_row(cb-1), ce - cb + 3, &
                                &MPI_LOGICAL, src, 0, MPI_COMM_WORLD,status)
                                line_logical(cb:ce) = t_row(cb:ce)
                            end if 
                        end do ! col

                       ! The grid is printed line by line once the True o False 
                       ! values are converted to characters
                        do j=1,w 
                            select case(line_logical(j))
                            case(.true.)
                                line_character(j:j) = 'X'
                            case(.false.)
                                line_character(j:j) = '.'
                        end select
                        end do !j

                        print '(a)', line_character(:)
                    end do !i
                end do !row

                print *

                if (allocated(line_character)) deallocate(line_character)
                if (allocated(line_logical)) deallocate(line_logical)
                if (allocated(t_row)) deallocate(t_row)

            end block

        else ! my_rank/=root

            rb = ib
            re = ie
            do i = rb, re
                ! Each process sends its sub-grid line by line to the root 
                ! process
                call MPI_Send(map(i, jb-1), 1, a_row, root, 0, MPI_COMM_WORLD)
            end do
        
        end if 
        return 
        
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

    ! Wait specified number of ms and then clear the terminal screen.
    subroutine wait_cls( ms )
        integer, intent(in) :: ms
        integer :: tick, tack
        real :: rate

        call system_clock( count=tick, count_rate=rate )
        do
            call system_clock( count=tack )
            if (real( tack - tick ) / rate >= ms * 1e-3) exit
        end do
        ! Clear the terminal screen using console escape code ^[2J.
        print "(2a)", achar( 27 ), '[2J'
    end subroutine wait_cls

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
