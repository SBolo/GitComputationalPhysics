program funzione_if

implicit none

real :: r = .74, s
integer :: i

    if ( r == 0.74 ) then
    
        do i = 1, 10, 2 !2 è lo step: può essere solo intero
                if(i==1.OR.i==3) cycle !se è verificato, salta l'istruzione e torna a inizio ciclo sommando STEP all'indice

                if(i==7) exit !esce dal ciclo quando i = 7
            print 100, "Numero: ", i
        end do

    else

        print 200, "Hai sbagliato!"

    end if

    
    100 format(1X,A8,I2)
    200 format (1X,A14)
    stop

end program funzione_if