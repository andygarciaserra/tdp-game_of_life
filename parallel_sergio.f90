program game_of_life
    use mpi_f08
    implicit none
    integer :: height, width
    integer :: max_gen, gen
    integer :: n_rows, n_cols, row, col
    integer :: my_rank, n_ranks, root, north_rank, south_rank, east_rank, west_rank
    integer :: ib, ie, jb, je
    logical, dimension(:, :), pointer :: old_world, new_world, tmp_world
    logical :: ghost_flag, root_world_still_flag
    type(MPI_Datatype) :: a_row, a_col
    type(MPI_Status)   :: status

    call MPI_Init()
    call MPI_Comm_rank( MPI_COMM_WORLD, my_rank)
    call MPI_Comm_size( MPI_COMM_WORLD, n_ranks)
    
    root = 0
    if (my_rank == root) then
        read *, height, width, max_gen, n_rows, n_cols

        if ((n_rows * n_cols) /= n_ranks) then
            print "(a)", "Incorrect number of processes"
            call MPI_Abort( MPI_COMM_WORLD, MPI_ERR_TOPOLOGY )
        end if
    end if
    
    call MPI_Bcast(n_rows, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(n_cols, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(height, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(width, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(max_gen, 1, MPI_INTEGER, root, MPI_COMM_WORLD)

    call get_coords( my_rank, n_rows, n_cols, row, col)
    north_rank = get_rank(row - 1, col, n_rows, n_cols)
    south_rank = get_rank(row + 1, col, n_rows, n_cols)
    west_rank  = get_rank(row, col - 1, n_rows, n_cols)
    east_rank  = get_rank(row, col + 1, n_rows, n_cols)


    call partition( row, n_rows, height, ib, ie )
    call partition( col, n_cols, width, jb, je )

    allocate(old_world(ib - 1 : ie + 1, jb - 1 : je + 1))
    allocate(new_world(ib - 1 : ie + 1, jb - 1 : je + 1))

    ! Definitions of MPI types
    call MPI_Type_contiguous(ie - ib + 3, MPI_LOGICAL, a_col) ! column with ghost layers
    CALL MPI_Type_commit(a_col)
    block
        type(MPI_Datatype) :: a_tmp_row
        integer(kind=MPI_ADDRESS_KIND) :: lb, real_extent

        call MPI_Type_vector(je - jb + 1, 1, ie - ib + 3, MPI_LOGICAL, a_tmp_row)
        call MPI_Type_get_extent( MPI_LOGICAL, lb, real_extent )
        call MPI_Type_create_resized( a_tmp_row, lb, real_extent, a_row )
        call MPI_Type_commit( a_row )
    end block

    call read_map( old_world, height, width )
    ghost_flag = .True.
    call update_borders( old_world )

    do gen = 1, max_gen
        if (my_rank == root) print "(a, i0)", "Generation ", gen
        call print_map( old_world, height, width )
        call next_gen( old_world, new_world )
        call update_borders( new_world )
        root_world_still_flag = world_is_still( old_world, new_world )
        call MPI_Bcast(root_world_still_flag, 1, MPI_LOGICAL, root, MPI_COMM_WORLD)
        if (root_world_still_flag) exit
        !if (my_rank == root) call wait_cls( 100 )
        ! Swap maps
        tmp_world => old_world;  old_world => new_world;  new_world => tmp_world
    end do

    if (associated( old_world )) deallocate(old_world)
    if (associated( new_world )) deallocate(new_world)

    call MPI_Type_free( a_row )
    call MPI_Type_free( a_col )
    call MPI_Finalize()

contains

function world_is_still( old_map, new_map ) result(still)
    logical, dimension(:, :), pointer, intent(in)   :: old_map, new_map
    logical                                         :: still
    logical, dimension(n_ranks)                     :: all_is_still

    still = all( old_map .eqv. new_map )
    call MPI_Gather(still, 1, MPI_LOGICAL, all_is_still,   1, MPI_LOGICAL, root, MPI_COMM_WORLD)
    if (my_rank == root) then
        still = all ( all_is_still .eqv. .True.)
    else
        still = .false.
    endif
end function world_is_still

    subroutine update_borders( map )
        logical, dimension(:, :), pointer, intent(inout) :: map

        ! synchronize rows
        call MPI_Sendrecv( map(ib,jb),     1, a_row, north_rank, 1, &
                           map(ie + 1,jb), 1, a_row, south_rank, 1, MPI_COMM_WORLD, status)
        call MPI_Sendrecv( map(ie,jb),     1, a_row, south_rank, 2, &
                           map(ib-1, jb),   1, a_row, north_rank, 2, MPI_COMM_WORLD, status)

        ! syncrhonize cols (big cols)
        call MPI_Sendrecv( map(ib - 1, jb), 1, a_col, west_rank, 3, &
                           map(ib - 1, je + 1), 1, a_col, east_rank, 3, MPI_COMM_WORLD, status)
        call MPI_Sendrecv( map(ib - 1, je), 1, a_col, east_rank, 3, &
                           map(ib - 1, jb - 1), 1, a_col, west_rank, 3, MPI_COMM_WORLD, status)


    end subroutine update_borders

    subroutine read_map( map, h, w )
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer, intent(in) :: h, w
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
                        read *, line
                        do current_col = 0, n_cols - 1
                            call partition(current_col, n_cols, w, cb, ce)
                            dst = get_rank(current_row, current_col, n_rows, n_cols)
                            allocate(temp(ce - cb + 1))
                            do j = cb, ce 
                                select case (line(j:j))
                                case ('X')
                                    temp(j - cb + 1) = .True.
                                case ('.')
                                    temp(j - cb + 1) = .False.
                                case default
                                    stop "read_map: wrong input character `" // line(j:j) // "`"
                                end select
                            end do
                            if (dst == root) then
                                map(i, cb : ce)  = temp
                            else
                                call MPI_Send( temp, ce - cb + 1, MPI_LOGICAL, dst, 0,  MPI_COMM_WORLD )
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do
                    end do
                end do
                if (allocated( line )) deallocate(line)
            end block
        else
            do i = ib, ie
                call MPI_Recv(map(i,jb), 1, a_row, root, 0, MPI_COMM_WORLD, status)
            end do
        end if
    end subroutine read_map

    subroutine barrier_print_map (map, ghost_flag)
        logical, dimension(:, :), pointer, intent(in) :: map
        logical,                           intent(in) :: ghost_flag

        character(len=:), allocatable :: line
        integer :: rank
        integer :: i, j
        integer :: ghost_index
        
        ! ghost_index controls if ghost layers are printed or not   
        ghost_index = merge ( 1, 0, ghost_flag)

        do rank = 0, n_ranks
            if (rank == my_rank) then
                print *, "Process: ", my_rank
                allocate(character(len=je-jb+1 + (2*ghost_index)) :: line)
                do i = ib - ghost_index, ie + ghost_index
                    do j = jb - ghost_index, je + ghost_index
                        line(j-jb+1+ghost_index:j-jb+1+ghost_index) = merge ( 'X', '.', map(i,j))
                    end do
                    print "(a)", line
                end do
                print *
                if (allocated( line )) deallocate(line)
            end if
            call MPI_Barrier( MPI_COMM_WORLD ) 
        end do
    end subroutine

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
                                call MPI_Recv(temp, ce-cb+1, MPI_LOGICAL, dst, 0, MPI_COMM_WORLD, status)
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
               call MPI_Send( map(i, jb), 1, a_row, root, 0, MPI_COMM_WORLD)
           end do
       end if
        

    end subroutine

    subroutine next_gen( old_map, new_map )
        logical, dimension(:, :), pointer, intent(inout) :: old_map, new_map
        integer :: i, j
        integer :: c ! the number of live neighbors

        do j = jb, je
            do i = ib, ie
                c = count( old_map(i - 1:i + 1, j - 1:j + 1) )
                if (old_map(i, j)) then ! cell is live
                    new_map(i, j) = merge( .true., .false., 3 <= c .and. c <= 4 )
                else ! cell is dead
                    new_map(i, j) = merge( .true., .false., c == 3 )
                end if
            end do
        end do
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

    ! Parallel subroutines ----------------------------------------------------
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

    integer function get_rank( row, col, n_rows, n_cols )
        integer, intent(in) ::  row, col, n_rows, n_cols
        integer :: aux_row, aux_col
        aux_row = row
        aux_col = col

        if (      0 <= col .and. col < n_cols &
            .and. 0 <= row .and. row < n_rows) then
                get_rank = row + col * n_rows
        else ! case when we apply toroidal topology
            if (row < 0) then
                aux_row = n_rows - 1
            else if (row >= n_rows) then 
                aux_row = 0
            end if

            if (col < 0) then
                aux_col = n_cols - 1
            else if (col >= n_cols) then 
                aux_col = 0
            end if

            get_rank = aux_row + aux_col * n_rows
        end if
    end function get_rank

    subroutine partition (id, n_ids, size, b, e)
        integer, intent(in)    :: id, n_ids, size
        integer, intent(inout) :: b, e
        integer :: remainder, quotient

        remainder = modulo( size, n_ids )
        quotient  = (size - remainder) / n_ids
        b = 1 + quotient * (id    ) + min( remainder, id     )
        e =     quotient * (id + 1) + min( remainder, id + 1 )
    end subroutine partition

end program game_of_life
