Module Data

  IMPLICIT NONE

  Real :: Lx, Ly, D, dx, dy, dt , epsilon
  Integer :: Nx, Ny, R
  Character (len = 32) :: SystType
  

CONTAINS

  Subroutine ReadData()

    Character (len = 32) :: file_name

    print *, "Nom du fichier ?"
    read *, file_name

    ! Ouverture fichier data.txt
    open (unit=11,file=file_name,action="read",status="old")

      Read (unit=11,fmt=*) Lx
      Read (unit=11,fmt=*) Ly
      Read (unit=11,fmt=*) Nx
      Read (unit=11,fmt=*) Ny
      Read (unit=11,fmt=*) D
      Read (unit=11,fmt=*) SystType
      Read (unit=11,fmt=*) epsilon
      Read (unit=11,fmt=*) R

    close(unit=11)                   ! fermeture fu fichier

    dx = Lx/(Nx+1)
    dy = Ly/(Ny+1)

    dt = 0.01

  End Subroutine

END MODULE Data
