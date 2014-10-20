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


mainsum = ciao(x,y,k)

    print*, mainsum


CONTAINS !da mettere quando specifico funzioni
! le funzioni vanno in fondo

    FUNCTION somma(a,b,c)
    !la funzione viene trattata come una variabile, quindi bisogna definirla come tale

        IMPLICIT NONE

        REAL :: somma, a, b, c
        somma = a+b+c

    END FUNCTION


    FUNCTION ciao(a,b,c)
        REAL :: ciao, a, b, c

    ciao = somma(a,b,c) + 3 ! posso usare funzioni semplicemente all'interno di altre
    ! anche se sono a tutti gli effetti delle variabili

    END FUNCTION


END PROGRAM