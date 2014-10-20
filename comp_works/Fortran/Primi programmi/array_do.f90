PROGRAM array

! Creiamo un array e maneggiamolo

IMPLICIT NONE !da mettere sempre per buona norma
INTEGER :: i
INTEGER, PARAMETER :: N = 10
    ! Se voglio fissare la dimensione di un array, serve che l'intero sia param
REAL, DIMENSION(0:N) :: vec !scelgo l'elemento iniziale dell'array
!
    DO i = 0, N
    ! funziona esattamente come il for, ma l'incremento posso ometterlo
    ! viene sottointeso che sia 1
        vec(i) = i + 2
        PRINT *, vec(i)
    ENDDO
!
END PROGRAM array