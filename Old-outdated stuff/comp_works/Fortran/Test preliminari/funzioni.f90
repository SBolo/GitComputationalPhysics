program product

implicit none

integer :: i = 2, j = 7, k
real :: f = 79.34
complex :: z1 = (1.,2.), z2 = (4., 7.)

!evita la dichiarazione di oggetti come impliciti: regola di merda, scrivilo e basta

    k = i*j

    print*, "Modulo: ", SQRT(z1**2+z2**2)
    print*, "Esponenziale: ", exp(f)
    print*, "Seno iperbolico: ", sinh(3.)
    print*, "Logaritmo naturale: ", log(2.72)
    write(*,*) i==j !restituisce F, ovvero false
    print*, i==2 !restituisce T, true

    print*, i/=j

end program product