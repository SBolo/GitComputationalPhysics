PROGRAM ch1201

implicit none

    integer :: x, i, j
    integer :: y

    open(unit = 1, file = 'DATA');

    DO x = 1, 10
        DO y = 1, 10

        write(1, fmt = 100) x, y
        100 format(' ', I2, ' ', I2)

        END DO
    END DO

        close(unit = 1)

    do i = 1, 3
        do j = 1, 3
            print 200, i, j
            200 format('x = ',' ',i3,'y = ',' ',i3)
        end do
    end do

END PROGRAM ch1201