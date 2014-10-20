program funzione_if

implicit none

real :: r = 3., s
integer :: i

    if ( r == 0.74 ) then
    
        do i = 1, 10
            print 100, "Numero: ", i
        end do

    else

        print 200, "Hai sbagliato!"

    end if

    
    100 format(1X,A8,I2)
    200 format (1X,A14)
    stop

end program funzione_if