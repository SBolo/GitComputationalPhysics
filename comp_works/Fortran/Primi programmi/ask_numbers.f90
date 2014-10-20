PROGRAM ask_name

! This program reads in and prints out a name !

IMPLICIT NONE !da mettere sempre per buona norma

REAL :: N1, N2, N3 !char da 20 caratteri che chiamo First_Name
REAL, PARAMETER :: C = 4. ! variabile a cui non ho accesso in sovrascrittura
    ! se tento di sovrascrivere parameter il compilatore mi restituisce errore
!
    PRINT *,' Scrivi tre numeri separati da uno spazio.'
    READ *, N1 , N2, N3

    N1 = C * (N1 + N2 + N3)
    N2 = (N1 + N2 + N3)/3.
    PRINT *,N1, N2
!
END PROGRAM ask_name