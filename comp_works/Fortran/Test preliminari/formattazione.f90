MODULE mymod
IMPLICIT NONE
INTEGER, PARAMETER :: single=4 !singola precisione
INTEGER, PARAMETER :: double=8 !doppia precisione
END MODULE

program product

USE mymod
implicit none

integer :: i = 2, j = 7, k
real(kind=double) :: f = 79.34
complex :: z1 = (1.,2.), z2 = (4., 7.)

!evita la dichiarazione di oggetti come impliciti: regola di merda, scrivilo e basta

    k = i*j

    print*, "Modulo: ", SQRT(z1**2+z2**2)
    print*, "Esponenziale: ", exp(f)
    print*, "Seno iperbolico: ", sinh(3.)
    print*, "Logaritmo naturale: ", log(2.72)

    ! format => serve a formattare il testo. 
    ! 1 = numero che definisce la formattazione
    ! 1x = uno spazio all'inizio della riga
    ! I1 = seguito da un intero da un carattere

    print*, i==2 !restituisce T, true

    write(*,100) "form: ", f
    100 format(1X,A5,1X,F5.2) !F5.2 = numero di 5 caratteri (compreso il .) con 2 cifre decimali

    print*, i/=j

end program product