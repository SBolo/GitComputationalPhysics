program product

implicit none

integer :: i = 2, j = 7, k
real :: f = 79.34
complex :: z1 = (1.,2.), z2 = (4., 7.)

!evita la dichiarazione di oggetti come impliciti: regola di merda, scrivilo e basta

    k = i*j

    print*, "Prodotto complesso: ", z1*z2
    print*, "Prodotto complesso e reale: ", z1*f
    print*, "Modulo: ", SQRT(z1**2+z2**2)

end program product