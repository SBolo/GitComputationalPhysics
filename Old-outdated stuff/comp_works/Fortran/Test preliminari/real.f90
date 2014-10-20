program product

implicit none

integer :: i = 2, j = 7, k
real :: f = 79.34

!evita la dichiarazione di oggetti come impliciti: regola di merda, scrivilo e basta

    k = i*j

    print*, "Risultato 2: ", k, f, k*f !altro modo per printare
    print*, "Elevamento a potenza ", k**2, 2./3. !elevamento a potenza come gnuplot + divisione col punto

end program product