PROGRAM array

! Creiamo un array e maneggiamolo

IMPLICIT NONE !da mettere sempre per buona norma. Evita la dichiarazione implicita di variabili quando non vengono definite
              !a inizio programma

INTEGER :: dim, i, j
REAL, DIMENSION(-3:3) :: vec
REAL, DIMENSION(1:6, 1:6) :: vec2
REAL, DIMENSION(1:6, 1:6) :: vec3
REAL, DIMENSION(1:3) :: arg = (/3., 3.2, 3.3/)
!
    vec(-3) = 1
    vec2(1,1) = 5

    PRINT *, "Ciao!", arg(1), arg(2), arg(3)

    vec3 = vec2; ! copio un array in un altro, tanta roba!
! posso fare tutte le operazioni che desidero su array e matrici tranne che il prodotto di matrici, il * 
! è la moltiplicazione posto a posto

    PRINT *, 'Ho allocato un vettore con posizioni negative', vec(-3)
    PRINT *, 'Ho copiato vec2 in vec3', vec3(1,1)

    PRINT *,' Scrivi tre numeri separati da uno spazio.'
    READ *, vec(1) , vec(2), vec(3)

    PRINT *, ' Risultato'
    PRINT *, (vec(1) - vec(2))*vec(3)

    DO i = 0, 6, 2
        DO j = 0, 6, 2
            vec2(i,j) = i+j
        END DO
    END DO

    PRINT *, "Hai caricato l'array ed il sesto posto è occupato da", vec2(6,6)
!


END PROGRAM array