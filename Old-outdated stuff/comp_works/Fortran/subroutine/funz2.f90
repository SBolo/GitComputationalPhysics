MODULE params
!funziona come una struttura, in cui però puoi anche definire variabili
!non male, perchè puoi assegnare direttamente un valore a x, y, z già nel modulo
!una volta dichiarate nel modulo, le variabili sono disponibili nel main
implicit none !solito, da ricordarsi
save !meglio metterlo per essere sicuri che si salvino

    real :: x, y, z
    real :: k

END MODULE params

PROGRAM main
USE params !per usare il modulo nel programma effettivo

IMPLICIT NONE


REAL :: mainsum

x = 3
y = 2
k = 2


mainsum = bisection(0.,2.)

    print*, mainsum


CONTAINS !da mettere quando specifico funzioni
! le funzioni vanno in fondo

    FUNCTION f(x)

        IMPLICIT NONE

        REAL :: f, x
        f = x**2 - 1

    END FUNCTION


    FUNCTION bisection(min, max)
       REAL :: bisection, sign, m, min, max, a, b

a = min
b = max

    do

            sign = f(min) * f(max)
        if(sign > 0) then
                    print*, "Non va bene"
                    EXIT

        else
                m = (min + max)/2.
                print*, "Siamo nell'if giusto", m
    

                if (f(m)*f(min) < 0) then
                        max = m
                        print*, "Siamo qui"
                else
                        min = m
                        print*, "Siamo la"
                end if

        end if

    bisection = m

    if( (max - min) < 0.001 ) EXIT
    end do

    END FUNCTION


END PROGRAM