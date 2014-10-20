PROGRAM array

! Creiamo un array e maneggiamolo

IMPLICIT NONE !da mettere sempre per buona norma
INTEGER :: i, N
REAL, DIMENSION(:), ALLOCATABLE :: vec !scelgo l'elemento iniziale dell'array
! lo avverto che voglio tenere lo spazio per allocare successivamente un array
!

    PRINT *, ' Dimensione array'
    READ *, N
    ALLOCATE(vec(0:N)) ! alloca l'array ora che ne conosce la dimensione

    DO i = 0, N
    ! funziona esattamente come il for, ma l'incremento posso ometterlo
    ! viene sottointeso che sia 1
        vec(i) = i + 2
        PRINT *, vec(i)
    ENDDO
!
END PROGRAM array