program funzione_if

implicit none

real :: r = .74, s
integer :: i


    print 200, "Inserisci un numero: "
    read(*,*) i
    print 200, "I cubi minori di 10.000 sono: "
    do
        i=i**3
        print 100, i
        if(i>10000) exit
    end do

    
    100 format(1X,I7)
    200 format (1X,A21)
    stop

end program funzione_if