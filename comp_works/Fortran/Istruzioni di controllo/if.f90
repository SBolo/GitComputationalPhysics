program funzione_if

implicit none

real :: r = 0.74, s

    write(*,100) "Inserisci un numero: "
    100 format (1X,A20)
    read(*,*)
    

    if(s == r) then
        write(*,10), "Hai indovinato il numero!"
        10 format(1X,A30)
    else
        write(*,10), "Non hai indovinato il numero!"
    end if

    stop

end program funzione_if