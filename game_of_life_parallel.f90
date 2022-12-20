program game_of_life
    implicit none
    integer :: height, width
    integer :: max_gen, gen
    logical, dimension(:, :), pointer :: old_world, new_world, tmp_world

    read *, height, width, max_gen !, n_rows, n_cols
    allocate(old_world(0:height + 1, 0:width + 1))
    allocate(new_world(0:height + 1, 0:width + 1))
    call read_map( old_world, height, width )
    call update_borders( old_world, height, width )

    do gen = 1, max_gen
        print "(a, i0)", "Generation ", gen
        call print_map( old_world, height, width )
        call next_gen( old_world, new_world, height, width )
        call update_borders( new_world, height, width )
        call wait_cls( 100 )
        if (world_is_still( old_world, new_world )) exit
        ! Swap maps
        tmp_world => old_world;  old_world => new_world;  new_world => tmp_world
    end do

    if (associated( old_world )) deallocate(old_world)
    if (associated( new_world )) deallocate(new_world)

contains

    logical function world_is_still( old_map, new_map )
        logical, dimension(:, :), pointer, intent(in) :: old_map, new_map

        world_is_still = all( old_map .eqv. new_map )
    end function world_is_still

    subroutine update_borders( map, h, w )
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer, intent(in) :: h, w

        ! Inner rows
        map(0,     1:w) = map(h, 1:w)
        map(h + 1, 1:w) = map(1, 1:w)
        ! Full columns
        map(0:h + 1, 0    ) = map(0:h + 1, w)
        map(0:h + 1, w + 1) = map(0:h + 1, 1)
    end subroutine update_borders

    subroutine read_map( map, h, w )
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer, intent(in) :: h, w
        character(len=:), allocatable :: line
        integer :: i, j
        
        allocate(character(len=w) :: line)
        do i = 1, h
            read *, line
            do j = 1, w
                select case (line(j:j))
                case ('X')
                    map(i, j) = .true.
                case ('.')
                    map(i, j) = .false.
                case default
                    stop "read_map: wrong input character `" // line(j:j) // "`"
                end select
            end do
        end do
        if (allocated( line )) deallocate(line)
    end subroutine read_map

    subroutine print_map( map, h, w )
        logical, dimension(:, :), pointer, intent(in) :: map
        integer, intent(in) :: h, w
        character(len=:), allocatable :: line
        integer :: i, j

        allocate(character(len=w) :: line)
        do i = 1, h
            do j = 1, w
                line(j:j) = merge( 'X', '.', map(i, j) )
            end do
            print "(a)", line
        end do
        print *
        if (allocated( line )) deallocate(line)
    end subroutine print_map

    subroutine next_gen( old_map, new_map, h, w )
        logical, dimension(:, :), pointer, intent(inout) :: old_map, new_map
        integer, intent(in) :: h, w
        integer :: i, j
        integer :: c ! the number of live neighbors

        do j = 1, w
            do i = 1, h
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

end program game_of_life
