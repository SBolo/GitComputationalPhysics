PROGRAM ask_name

! This program reads in and prints out a name !

IMPLICIT NONE !da mettere sempre per buona norma

CHARACTER*20 :: First_Name !char da 20 caratteri che chiamo First_Name
!
    PRINT *,' Type in your first name.'
    PRINT *,' up to 20 characters'
    READ *,First_Name
    PRINT *,First_Name
!
END PROGRAM ask_name